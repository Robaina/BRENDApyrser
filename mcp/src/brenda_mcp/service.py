"""
Service layer for the BRENDA MCP server.

Wraps the :class:`brendapyrser.BRENDA` API with:

* a lazily-loaded, process-wide cached database (the BRENDA JSON is large, so we
  parse it once on first use), and
* helpers that turn the rich Python objects into compact, JSON-serialisable
  summaries — statistics and histograms rather than thousands of raw numbers —
  so tool results stay small enough to be cheap for an LLM to read.

Everything here is transport-agnostic; ``server.py`` only adds the MCP tool
definitions on top.
"""

from __future__ import annotations

import math
import os
import statistics
import sys
import threading
from typing import Any, Iterable, Optional

from brendapyrser import BRENDA

# Environment variable pointing at the BRENDA JSON database. Accepts a plain
# ``.json`` file, a ``.json.gz``, or the ``.json.tar.gz`` archive BRENDA ships.
DB_PATH_ENV = "BRENDA_DATABASE_PATH"

# Units for each queryable quantity (mirrors brendapyrser.constants.units, with
# the temperature/pH conditions added).
UNITS: dict[str, str] = {
    "km": "mM",
    "kcat": "1/s",
    "ki": "mM",
    "kcat_km": "mM^-1 s^-1",
    "specific_activity": "umol/min/mg",
    "temperature_optimum": "degC",
    "temperature_range": "degC",
    "temperature_stability": "degC",
    "ph_optimum": "pH",
    "ph_range": "pH",
    "ph_stability": "pH",
}

# parameter name -> Reaction property holding an EnzymePropertyDict
_KINETIC_PROPERTY = {
    "km": "KMvalues",
    "kcat": "Kcatvalues",
    "ki": "KIvalues",
    "kcat_km": "KKMvalues",
}

_brenda: Optional[BRENDA] = None
_lock = threading.Lock()


def _log(msg: str) -> None:
    # NEVER write to stdout: it is the MCP stdio channel. Logs go to stderr.
    print(f"[brenda-mcp] {msg}", file=sys.stderr, flush=True)


def get_brenda(path: Optional[str] = None) -> BRENDA:
    """
    Return the shared :class:`BRENDA` instance, parsing the database on first
    use. ``path`` overrides the ``BRENDA_DATABASE_PATH`` environment variable
    (used by the test suite to point at the small fixture).
    """
    global _brenda
    if _brenda is not None and path is None:
        return _brenda

    db_path = path or os.environ.get(DB_PATH_ENV)
    if not db_path:
        raise RuntimeError(
            f"No BRENDA database configured. Set the {DB_PATH_ENV} environment "
            "variable to the path of the BRENDA JSON file "
            "(e.g. data/brenda_2026_1.json or the .json.tar.gz archive)."
        )
    if not os.path.exists(db_path):
        raise RuntimeError(f"BRENDA database not found at: {db_path}")

    with _lock:
        # Re-check inside the lock in case another thread loaded it meanwhile.
        if _brenda is not None and path is None:
            return _brenda
        _log(f"loading BRENDA database from {db_path} (this can take a while) ...")
        brenda = BRENDA(db_path)
        _log(
            f"loaded release {brenda.release!r}: "
            f"{len(brenda.reactions)} enzyme entries."
        )
        if path is None:
            _brenda = brenda
        return brenda


# --------------------------------------------------------------------------- #
# Statistics helpers                                                          #
# --------------------------------------------------------------------------- #
def _percentile(sorted_vals: list[float], pct: float) -> float:
    """Linear-interpolated percentile of an already-sorted list."""
    if not sorted_vals:
        return float("nan")
    if len(sorted_vals) == 1:
        return sorted_vals[0]
    rank = (pct / 100) * (len(sorted_vals) - 1)
    lo = math.floor(rank)
    hi = math.ceil(rank)
    if lo == hi:
        return sorted_vals[lo]
    frac = rank - lo
    return sorted_vals[lo] * (1 - frac) + sorted_vals[hi] * frac


def _histogram(sorted_vals: list[float], bins: int) -> list[dict[str, float]]:
    lo, hi = sorted_vals[0], sorted_vals[-1]
    width = (hi - lo) / bins
    edges = [lo + i * width for i in range(bins + 1)]
    counts = [0] * bins
    for v in sorted_vals:
        idx = int((v - lo) / width) if width else 0
        if idx >= bins:  # the maximum value lands in the last bin
            idx = bins - 1
        counts[idx] += 1
    return [
        {"min": round(edges[i], 4), "max": round(edges[i + 1], 4), "count": counts[i]}
        for i in range(bins)
    ]


def _coerce(
    values: Iterable[Any], *, max_value: Optional[float], drop_negative: bool
) -> list[float]:
    """Coerce to finite floats, dropping the parser's -999 sentinel (negatives)
    and anything above ``max_value`` (used to trim implausible outliers, exactly
    as the example notebook does with ``values < 1000``)."""
    out: list[float] = []
    for v in values:
        try:
            f = float(v)
        except (TypeError, ValueError):
            continue
        if not math.isfinite(f):
            continue
        if drop_negative and f < 0:  # -999 sentinel for unparseable values
            continue
        if max_value is not None and f > max_value:
            continue
        out.append(f)
    return out


def summarize_values(
    values: Iterable[Any],
    *,
    bins: int = 10,
    max_value: Optional[float] = None,
    drop_negative: bool = True,
) -> dict[str, Any]:
    """Reduce a sequence of numbers to a compact summary (count + descriptive
    stats + optional histogram). Returns ``{"count": 0}`` when nothing parses."""
    clean = sorted(_coerce(values, max_value=max_value, drop_negative=drop_negative))
    n = len(clean)
    if n == 0:
        return {"count": 0}
    stats: dict[str, Any] = {
        "count": n,
        "min": round(clean[0], 6),
        "max": round(clean[-1], 6),
        "mean": round(statistics.fmean(clean), 6),
        "median": round(statistics.median(clean), 6),
        "p25": round(_percentile(clean, 25), 6),
        "p75": round(_percentile(clean, 75), 6),
    }
    if bins and n >= 2 and clean[-1] > clean[0]:
        stats["histogram"] = _histogram(clean, bins)
    return stats


# --------------------------------------------------------------------------- #
# Matching helpers (forgiving compound / organism matching for LLM callers)    #
# --------------------------------------------------------------------------- #
def _organism_matches(query: str, species: Iterable[str]) -> bool:
    q = query.lower()
    return any(q in (s or "").lower() for s in species)


def _match_compound_keys(prop: dict, compound: Optional[str]) -> list[str]:
    """Resolve a user-supplied compound name to actual keys of an
    EnzymePropertyDict: exact, then case-insensitive, then substring."""
    if compound is None:
        return list(prop.keys())
    if compound in prop:
        return [compound]
    low = compound.lower()
    ci = [k for k in prop if k.lower() == low]
    if ci:
        return ci
    return [k for k in prop if low in k.lower()]


def _collect_property_values(
    prop: dict,
    *,
    compound: Optional[str],
    organism: Optional[str],
) -> tuple[list[float], dict[str, list[float]]]:
    """Walk an EnzymePropertyDict ({compound: [{value, species, ...}, ...]}) and
    return (all_values, {compound: [values]}) after applying compound/organism
    filters at the record level."""
    by_compound: dict[str, list[float]] = {}
    for key in _match_compound_keys(prop, compound):
        vals: list[float] = []
        for entry in prop.get(key, []):
            if organism and not _organism_matches(organism, entry.get("species", [])):
                continue
            vals.append(entry.get("value"))
        if vals:
            by_compound[key] = vals
    all_values = [v for vals in by_compound.values() for v in vals]
    return all_values, by_compound


# --------------------------------------------------------------------------- #
# Reaction-level summaries                                                     #
# --------------------------------------------------------------------------- #
def reaction_overview(r) -> dict[str, Any]:
    """Compact overview of one enzyme, including a map of which data fields are
    populated so the agent knows what is worth querying next."""
    temp = r.temperature
    ph = r.PH
    return {
        "ec_number": r.ec_number,
        "name": r.name,
        "systematic_name": r.systematic_name,
        "reaction": r.reaction_str,
        "reaction_type": r.reaction_type,
        "synonyms": r.synonyms[:15],
        "n_organisms": len(r.organisms),
        "data_available": {
            "km_compounds": len(r.KMvalues),
            "kcat_compounds": len(r.Kcatvalues),
            "ki_compounds": len(r.KIvalues),
            "kcat_km_compounds": len(r.KKMvalues),
            "specific_activity_records": len(r.specificActivities),
            "temperature_optimum_records": len(temp.get("optimum", [])),
            "ph_optimum_records": len(ph.get("optimum", [])),
            "cofactors": len(r.cofactors),
            "inhibitors": len(r.inhibitors),
            "activators": len(r.activators),
            "metals": len(r.metals),
            "natural_substrate_product_pairs": len(r.substratesAndProducts),
            "references": len(r.references),
        },
    }


def enzyme_kinetics(
    r,
    parameter: str,
    *,
    compound: Optional[str],
    organism: Optional[str],
    max_value: Optional[float],
    include_values: bool,
    values_limit: int = 200,
) -> dict[str, Any]:
    """KM / kcat / Ki / kcat/KM / specific-activity summary for one enzyme."""
    result: dict[str, Any] = {
        "ec_number": r.ec_number,
        "name": r.name,
        "parameter": parameter,
        "unit": UNITS.get(parameter, ""),
        "compound_filter": compound,
        "organism_filter": organism,
    }

    if parameter == "specific_activity":
        records = r.specificActivities
        vals = [
            rec.get("value")
            for rec in records
            if not organism or _organism_matches(organism, rec.get("species", []))
        ]
        result["stats"] = summarize_values(vals, max_value=max_value)
        if include_values:
            result["values"] = _coerce(vals, max_value=max_value, drop_negative=True)[
                :values_limit
            ]
        return result

    if parameter not in _KINETIC_PROPERTY:
        raise ValueError(
            f"Unknown kinetic parameter {parameter!r}. "
            f"Valid: {sorted(list(_KINETIC_PROPERTY) + ['specific_activity'])}"
        )

    prop = getattr(r, _KINETIC_PROPERTY[parameter])
    all_values, by_compound = _collect_property_values(
        prop, compound=compound, organism=organism
    )
    result["stats"] = summarize_values(all_values, max_value=max_value)
    # Per-compound breakdown (most-measured compounds first), compact.
    breakdown = []
    for cmp_name, vals in sorted(by_compound.items(), key=lambda kv: -len(kv[1])):
        s = summarize_values(vals, bins=0, max_value=max_value)
        if s["count"]:
            breakdown.append(
                {
                    "compound": cmp_name,
                    "count": s["count"],
                    "median": s["median"],
                    "min": s["min"],
                    "max": s["max"],
                }
            )
    result["by_compound"] = breakdown[:25]
    if include_values:
        result["values"] = _coerce(all_values, max_value=max_value, drop_negative=True)[
            :values_limit
        ]
    return result


def enzyme_conditions(
    r,
    prop_name: str,
    condition: str,
    *,
    organism: Optional[str],
) -> dict[str, Any]:
    """Temperature or pH summary (optimum / range / stability) for one enzyme."""
    source = r.temperature if prop_name == "temperature" else r.PH
    records = source.get(condition, [])
    param_key = f"{prop_name}_{condition}"
    result: dict[str, Any] = {
        "ec_number": r.ec_number,
        "name": r.name,
        "property": prop_name,
        "condition": condition,
        "unit": UNITS.get(param_key, ""),
        "organism_filter": organism,
    }
    kept = [
        rec
        for rec in records
        if not organism or _organism_matches(organism, rec.get("species", []))
    ]
    if condition == "range":
        # values are [low, high] pairs
        pairs = [
            rec.get("value")
            for rec in kept
            if isinstance(rec.get("value"), (list, tuple)) and len(rec["value"]) == 2
        ]
        lows = [p[0] for p in pairs]
        highs = [p[1] for p in pairs]
        result["count"] = len(pairs)
        result["low"] = summarize_values(lows, bins=0)
        result["high"] = summarize_values(highs, bins=0)
        result["ranges"] = [
            [round(float(lo), 3), round(float(hi), 3)]
            for lo, hi in pairs
            if lo >= 0 and hi >= 0
        ][:100]
    else:
        result["stats"] = summarize_values([rec.get("value") for rec in kept])
    return result


def enzyme_compounds(r, kind: str, *, limit: int = 100) -> dict[str, Any]:
    """Cofactors / inhibitors / activators / metals / substrate-product pairs /
    synonyms for one enzyme."""
    result: dict[str, Any] = {"ec_number": r.ec_number, "name": r.name, "kind": kind}
    actuators = {
        "cofactors": r.cofactors,
        "inhibitors": r.inhibitors,
        "activators": r.activators,
        "metals": r.metals,
    }
    if kind in actuators:
        d = actuators[kind]
        result["total"] = len(d)
        result["items"] = [
            {"compound": name, "organisms": (info.get("species") or [])[:8]}
            for name, info in list(d.items())[:limit]
        ]
    elif kind == "substrates_products":
        pairs = r.substratesAndProducts
        result["total"] = len(pairs)
        result["items"] = pairs[:limit]
    elif kind == "synonyms":
        syn = r.synonyms
        result["total"] = len(syn)
        result["items"] = syn[:limit]
    else:
        raise ValueError(
            "kind must be one of: cofactors, inhibitors, activators, metals, "
            "substrates_products, synonyms"
        )
    return result


def enzyme_organisms(r, *, limit: int = 100) -> dict[str, Any]:
    orgs = r.organisms
    return {
        "ec_number": r.ec_number,
        "name": r.name,
        "total": len(orgs),
        "organisms": orgs[:limit],
    }


def enzyme_references(r, *, limit: int = 25) -> dict[str, Any]:
    refs = r.references  # {id: citation_string}
    items = list(refs.values())
    return {
        "ec_number": r.ec_number,
        "name": r.name,
        "total": len(items),
        "references": items[:limit],
    }


# --------------------------------------------------------------------------- #
# Database-level searches and distributions                                   #
# --------------------------------------------------------------------------- #
def search_enzymes(query: str, *, limit: int = 25) -> dict[str, Any]:
    """Substring search over EC number, recommended name and synonyms."""
    q = query.lower().strip()
    hits = []
    for r in get_brenda().reactions:
        name = r.name or ""
        if (
            q in r.ec_number.lower()
            or q in name.lower()
            or any(q in (s or "").lower() for s in r.synonyms)
        ):
            hits.append({"ec_number": r.ec_number, "name": name})
            if len(hits) >= limit:
                break
    return {"query": query, "count": len(hits), "results": hits}


def find_enzymes_by_compound(
    compound: str, *, role: str = "any", limit: int = 25
) -> dict[str, Any]:
    rl = get_brenda().reactions
    if role == "substrate":
        matches = rl.filter_by_substrate(compound)
    elif role == "product":
        matches = rl.filter_by_product(compound)
    elif role == "any":
        matches = rl.filter_by_compound(compound)
    else:
        raise ValueError("role must be one of: substrate, product, any")
    return {
        "compound": compound,
        "role": role,
        "total": len(matches),
        "results": [
            {"ec_number": r.ec_number, "name": r.name} for r in matches[:limit]
        ],
    }


def find_enzymes_by_organism(organism: str, *, limit: int = 25) -> dict[str, Any]:
    matches = get_brenda().reactions.filter_by_organism(organism)
    return {
        "organism": organism,
        "total": len(matches),
        "results": [
            {"ec_number": r.ec_number, "name": r.name} for r in matches[:limit]
        ],
    }


def parameter_distribution(
    parameter: str,
    *,
    organism: Optional[str] = None,
    compound: Optional[str] = None,
    max_value: Optional[float] = None,
    bins: int = 10,
) -> dict[str, Any]:
    """Aggregate one kinetic/condition parameter across the whole database (or a
    genus/species subset), returning summary statistics + a histogram. This is
    the tool behind questions like "median KM in BRENDA" or "median optimal
    temperature in the genus Thermotoga"."""
    brenda = get_brenda()
    reactions = (
        brenda.reactions.filter_by_organism(organism) if organism else brenda.reactions
    )

    values: list[Any] = []
    n_scanned = 0
    for r in reactions:
        n_scanned += 1
        if parameter in _KINETIC_PROPERTY:
            prop = getattr(r, _KINETIC_PROPERTY[parameter])
            vals, _ = _collect_property_values(
                prop, compound=compound, organism=organism
            )
            values.extend(vals)
        elif parameter == "specific_activity":
            for rec in r.specificActivities:
                if not organism or _organism_matches(organism, rec.get("species", [])):
                    values.append(rec.get("value"))
        elif parameter in ("temperature_optimum", "ph_optimum"):
            prop_name, condition = parameter.split("_", 1)
            source = r.temperature if prop_name == "temperature" else r.PH
            for rec in source.get(condition, []):
                if not organism or _organism_matches(organism, rec.get("species", [])):
                    values.append(rec.get("value"))
        else:
            raise ValueError(
                "parameter must be one of: km, kcat, ki, kcat_km, "
                "specific_activity, temperature_optimum, ph_optimum"
            )

    return {
        "parameter": parameter,
        "unit": UNITS.get(parameter, ""),
        "organism_filter": organism,
        "compound_filter": compound,
        "max_value_filter": max_value,
        "enzymes_scanned": n_scanned,
        "stats": summarize_values(values, bins=bins, max_value=max_value),
    }


def database_info() -> dict[str, Any]:
    brenda = get_brenda()
    return {
        "release": brenda.release,
        "schema_version": brenda.schema_version,
        "n_enzymes": len(brenda.reactions),
        "copyright": " ".join(brenda.copyright.split()),
    }
