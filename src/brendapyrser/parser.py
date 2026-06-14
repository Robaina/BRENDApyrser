#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Object to parse and manipulate the BRENDA database.

Since the 2024.x releases, BRENDA is distributed as a single structured JSON
document (schema https://www.brenda-enzymes.org/schemas/docs/2.0.0) instead of
the legacy flat text file. This module parses that JSON while preserving the
historical public API (``BRENDA``, ``Reaction``, ``ReactionList`` and their
properties) and exposing the additional structure now available (UniProt/GenBank
accessions, structured citations with PMIDs, reversibility, etc.).
"""

from __future__ import annotations

import gzip
import io
import json
import re
import tarfile
from importlib import metadata
from pathlib import Path

import numpy as np
import pandas as pd

from .constants import fields, units

meta = metadata.metadata("brendapyrser")
__version__ = meta["Version"]
__author__ = meta["Author"]


_VALUE_UNIT_RE = re.compile(r"^(?P<num>.*?)\s*\{(?P<unit>.*)\}\s*$")
_REVERSIBILITY_RE = re.compile(r"\{(ir|r)\}")

# BRENDA still embeds legacy annotations inside reaction-string values:
#   {ir}/{r} reversibility, |...| pipe comments, and <1,2,3> numeric reference
# tags. These must be stripped before parsing substrates/products so they do not
# leak into compound names. The reference-tag pattern is numeric-only on purpose
# so it never touches chemical notation such as glycosidic bonds, e.g. "(1->4)".
_BRACE_RE = re.compile(r"\{.*?\}")
_PIPE_RE = re.compile(r"\|([^|]*)\|")
_REF_TAG_RE = re.compile(r"<[0-9][0-9,]*>")
_PROTEIN_TAG_RE = re.compile(r"#\d+(?:,\d+)*#")
_MULTISPACE_RE = re.compile(r"\s{2,}")


def _load_database(path_to_database) -> dict:
    """
    Load a BRENDA JSON database, transparently handling a plain ``.json`` file,
    a gzip-compressed ``.json.gz`` file, or a ``.json.tar.gz`` archive (the
    format in which BRENDA currently distributes the database).
    """
    path = Path(path_to_database)
    name = path.name.lower()

    if name.endswith((".tar.gz", ".tgz")):
        with tarfile.open(path, "r:gz") as tar:
            members = [
                m
                for m in tar.getmembers()
                if m.isfile() and m.name.lower().endswith(".json")
            ]
            if not members:
                raise ValueError(f"No .json member found in archive '{path}'")
            if len(members) > 1:
                raise ValueError(
                    f"Multiple .json members found in archive '{path}': "
                    f"{[m.name for m in members]}"
                )
            fobj = tar.extractfile(members[0])
            with io.TextIOWrapper(fobj, encoding="utf-8") as text:
                return json.load(text)

    if name.endswith(".gz"):
        with gzip.open(path, "rt", encoding="utf-8") as fh:
            return json.load(fh)

    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def _split_value_unit(value: str):
    """
    Split a BRENDA value such as ``"0.045 {NADH}"`` into ``("0.045", "NADH")``.
    Returns ``(number_part, brace_content)`` where ``brace_content`` is ``None``
    when the value carries no ``{...}`` annotation.
    """
    match = _VALUE_UNIT_RE.match(value or "")
    if match:
        return match.group("num").strip(), match.group("unit").strip()
    return (value or "").strip(), None


def _clean_reaction_string(value: str) -> str:
    """
    Strip legacy annotations (``|...|`` comments, ``{...}`` reversibility and
    ``<1,2>`` reference tags) from a reaction-string value, leaving only the
    reaction equation itself. Chemical notation such as ``(1->4)`` glycosidic
    bonds or polymer subscripts (``[...]n+m``) is preserved.
    """
    s = _PIPE_RE.sub(" ", value or "")
    s = _BRACE_RE.sub(" ", s)
    s = _REF_TAG_RE.sub(" ", s)
    return _MULTISPACE_RE.sub(" ", s).strip()


def _split_reaction_side(side: str) -> list:
    """
    Split one side of a reaction equation into compounds, splitting only on
    surrounding-space ``+`` so that ionic ``+`` (e.g. ``H+``, ``NAD+``) and
    polymer subscripts (``[...]n+m``) stay attached to the compound name.
    ``?`` and ``more`` placeholders are dropped.
    """
    return [
        token.strip()
        for token in side.split(" + ")
        if token.strip() not in ("", "?", "more")
    ]


class BRENDA:
    """
    Provides methods to parse the BRENDA database (https://www.brenda-enzymes.org/)
    """

    def __init__(self, path_to_database):
        database = _load_database(path_to_database)
        self.__release = database.get("release", "")
        self.__schema_version = database.get("version", "")
        data = database.get("data", {})
        # Every key is an EC number except the "spontaneous" pseudo-entry, which
        # is not an enzyme and is therefore excluded from the reaction list.
        self.__reactions = [
            Reaction(entry) for key, entry in data.items() if key != "spontaneous"
        ]
        self.__copyright = """Copyrighted by Dietmar Schomburg, Techn. University
        Braunschweig, GERMANY. Distributed under the License as stated
        at http:/www.brenda-enzymes.org"""
        self.__fields = fields
        self.__units = units

    def _repr_html_(self):
        """This method is executed automatically by Jupyter to print html!"""
        return """
        <table>
            <tr>
                <td><strong>Number of Enzymes</strong></td><td>{n_ec}</td>
            </tr><tr>
                <td><strong>BRENDA release</strong></td><td>{release}</td>
            </tr><tr>
                <td><strong>BRENDA copyright</strong></td><td>{cr}</td>
            </tr><tr>
                <td><strong>Brendapyrser version</strong></td><td>{parser}</td>
            </tr><tr>
                <td><strong>Author</strong></td><td>{author}</td>
            </tr>
        </table>
        """.format(
            n_ec=len(self.__reactions),
            release=self.__release,
            cr=self.__copyright,
            parser=__version__,
            author=__author__,
        )

    @property
    def fields(self):
        return self.__fields

    @property
    def units(self):
        return self.__units

    @property
    def reactions(self):
        return ReactionList(self.__reactions)

    @property
    def copyright(self):
        return self.__copyright

    @property
    def release(self):
        """BRENDA release the database was downloaded from, e.g. '2026.1'."""
        return self.__release

    @property
    def schema_version(self):
        """Version of the BRENDA JSON schema used by the loaded database."""
        return self.__schema_version

    def getOrganisms(self) -> list:
        """
        Get list of all represented species in BRENDA
        """
        species = set()
        for rxn in self.__reactions:
            species.update(rxn.organisms)
        species.discard("")
        return list({s for s in species if "no activity" not in s})

    def getKMcompounds(self) -> list:
        """
        Get list of all substrates in BRENDA with KM data
        """
        cpds = set()
        for rxn in self.__reactions:
            cpds.update(rxn.KMvalues.keys())
        cpds.discard("")
        return list(cpds)


class ReactionList(list):
    # Make ReactionList slicing return ReactionList object
    def __init__(self, seq=None):
        super(self.__class__, self).__init__(seq)

    def __getslice__(self, start, stop):
        return self.__class__(super(self.__class__, self).__getslice__(start, stop))

    def __getitem__(self, key):
        if isinstance(key, slice):
            return self.__class__(super(self.__class__, self).__getitem__(key))
        else:
            return super(self.__class__, self).__getitem__(key)

    def get_by_id(self, id: str):
        try:
            return [rxn for rxn in self if rxn.ec_number == id][0]
        except Exception:
            raise ValueError(f"Enzyme with EC {id} not found in database")

    def get_by_name(self, name: str):
        try:
            return [rxn for rxn in self if rxn.name.lower() == name.lower()][0]
        except Exception:
            raise ValueError(f"Enzyme {name} not found in database")

    def filter_by_substrate(self, substrate: str) -> list[Reaction]:
        """
        Filter reactions by a specific substrate
        """
        return [
            rxn
            for rxn in self
            if any(
                [substrate in mets["substrates"] for mets in rxn.substratesAndProducts]
            )
        ]

    def filter_by_product(self, product: str) -> list[Reaction]:
        """
        Filter reactions by a specific product
        """
        return [
            rxn
            for rxn in self
            if any([product in mets["products"] for mets in rxn.substratesAndProducts])
        ]

    def filter_by_compound(self, compound: str) -> list[Reaction]:
        """
        Filter reactions by a substrate or product
        """
        return [
            rxn
            for rxn in self
            if any(
                [
                    (compound in mets["substrates"] or compound in mets["products"])
                    for mets in rxn.substratesAndProducts
                ]
            )
        ]

    def filter_by_organism(self, species: str):
        def is_contained(p, S):
            return any([p in s.lower() for s in S])

        return self.__class__(
            [rxn for rxn in self if is_contained(species.lower(), rxn.organisms)]
        )


class EnzymeDict(dict):
    def filter_by_organism(self, species: str):
        filtered_dict = {}

        def is_contained(p, S):
            return any([p in s for s in S])

        for k in self.keys():
            filtered_values = [
                v for v in self[k] if is_contained(species, v["species"])
            ]
            if len(filtered_values) > 0:
                filtered_dict[k] = filtered_values
        return self.__class__(filtered_dict)

    def get_values(self):
        return [v["value"] for k in self.keys() for v in self[k]]


class EnzymePropertyDict(EnzymeDict):
    def filter_by_compound(self, compound: str):
        try:
            return self.__class__({compound: self[compound]})
        except Exception:
            return self.__class__({compound: []})


class EnzymeConditionDict(EnzymeDict):
    def filter_by_condition(self, condition: str):
        try:
            return self.__class__({condition: self[condition]})
        except Exception:
            raise KeyError(
                f'Invalid condition, valid conditions are: {", ".join(list(self.keys()))}'
            )


class Reaction:
    def __init__(self, entry: dict):
        self.__entry = entry
        self.__ec_number = entry.get("id", "")
        self.__name = entry.get("recommended_name", "")
        self.__systematic_name = entry.get("systematic_name", "")
        self.__proteins_raw = entry.get("protein", {})
        self.__references_raw = entry.get("reference", {})

    # ------------------------------------------------------------------ #
    # Internal helpers                                                    #
    # ------------------------------------------------------------------ #
    @staticmethod
    def __eval_range_value(v):
        """Evaluate a single numeric value, averaging ``a-b`` ranges."""
        try:
            if not re.search(r"\d-\d", v):
                return float(v)
            return float(np.mean([float(s) for s in v.split("-")]))
        except Exception:
            return -999

    @staticmethod
    def __eval_range_pair(v):
        """Evaluate a numeric range ``a-b`` into ``[a, b]``."""
        try:
            return [float(s) for s in v.split("-")]
        except Exception:
            return [-999, -999]

    @staticmethod
    def __format_citation(ref: dict) -> str:
        """Render a structured reference record as a citation string."""
        authors = "; ".join(ref.get("authors", []))
        citation = (
            f"{authors}: {ref.get('title', '')}. {ref.get('journal', '')} "
            f"({ref.get('year', '')}) {ref.get('vol', '')}, {ref.get('pages', '')}."
        ).strip()
        pmid = ref.get("pmid", "")
        if pmid:
            citation += f" {{Pubmed:{pmid}}}"
        return citation

    @staticmethod
    def __reversibility(value: str):
        """Return 'reversible'/'irreversible'/None from a {r}/{ir} annotation."""
        match = _REVERSIBILITY_RE.search(value or "")
        if not match:
            return None
        return "irreversible" if match.group(1) == "ir" else "reversible"

    def __organisms_for(self, protein_ids: list) -> list:
        """Map BRENDA protein ids to their organism names."""
        return list(
            {
                self.__proteins_raw[pid]["organism"]
                for pid in protein_ids
                if pid in self.__proteins_raw
            }
        )

    def __citations_for(self, ref_ids: list) -> list:
        """Map BRENDA reference ids to full citation strings."""
        return [
            self.__format_citation(self.__references_raw[rid])
            for rid in ref_ids
            if rid in self.__references_raw
        ]

    def __getDictOfEnzymeActuators(self, field: str) -> EnzymePropertyDict:
        res = {}
        for record in self.__entry.get(field, []):
            value = (record.get("value", "") or "").strip()
            if value != "more":
                res[value] = {
                    "species": self.__organisms_for(record.get("proteins", [])),
                    "meta": record.get("comment", ""),
                    "refs": self.__citations_for(record.get("references", [])),
                }
        return EnzymePropertyDict(res)

    def __getDictOfEnzymeProperties(self, field: str) -> EnzymePropertyDict:
        res = {}
        for record in self.__entry.get(field, []):
            number, substrate = _split_value_unit(record.get("value", ""))
            substrate = substrate if substrate is not None else ""
            if substrate == "more":
                continue
            res.setdefault(substrate, []).append(
                {
                    "value": self.__eval_range_value(number),
                    "species": self.__organisms_for(record.get("proteins", [])),
                    "meta": record.get("comment", ""),
                    "refs": self.__citations_for(record.get("references", [])),
                }
            )
        return EnzymePropertyDict(res)

    def __extractTempOrPHData(self, field: str, is_range: bool) -> list:
        values = []
        for record in self.__entry.get(field, []):
            raw = record.get("value", "")
            value = (
                self.__eval_range_pair(raw)
                if is_range
                else self.__eval_range_value(raw)
            )
            values.append(
                {
                    "value": value,
                    "species": self.__organisms_for(record.get("proteins", [])),
                    "meta": record.get("comment", ""),
                    "refs": record.get("references", []),
                }
            )
        return values

    def field(self, name: str) -> list:
        """
        Return the raw, enriched records for any BRENDA data field as a list of
        ``{value, organisms, comment, references}`` dicts, with protein and
        reference ids resolved to organism names and full citations. Useful for
        fields without a dedicated property (e.g. ``cloned``, ``application``).
        """
        enriched = []
        for record in self.__entry.get(name, []):
            if not isinstance(record, dict):
                continue
            enriched.append(
                {
                    "value": record.get("value", ""),
                    "organisms": self.__organisms_for(record.get("proteins", [])),
                    "comment": record.get("comment", ""),
                    "references": self.__citations_for(record.get("references", [])),
                }
            )
        return enriched

    # ------------------------------------------------------------------ #
    # Display helpers                                                     #
    # ------------------------------------------------------------------ #
    def printReactionSummary(self):
        data = {
            "EC number": self.__ec_number,
            "Name": self.__name,
            "Systematic name": self.__systematic_name,
            "Reaction type": self.reaction_type,
            "Reaction": self.reaction_str,
        }
        return pd.DataFrame.from_dict(data, orient="index", columns=[""])

    def _repr_html_(self):
        """This method is executed automatically by Jupyter to print html!"""
        return """
        <table>
            <tr>
                <td><strong>Enzyme identifier</strong></td><td>{ec}</td>
            </tr><tr>
                <td><strong>Name</strong></td><td>{name}</td>
            </tr><tr>
                <td><strong>Systematic name</strong></td><td>{sys_name}</td>
            </tr><tr>
                <td><strong>Reaction type</strong></td><td>{rxn_type}</td>
            </tr><tr>
                <td><strong>Reaction</strong></td><td>{rxn_str}</td>
            </tr>
        </table>
        """.format(
            ec=self.__ec_number,
            name=self.__name,
            sys_name=self.__systematic_name,
            rxn_type=self.reaction_type,
            rxn_str=self.reaction_str,
        )

    # ------------------------------------------------------------------ #
    # Public API (preserved)                                             #
    # ------------------------------------------------------------------ #
    @property
    def summary(self):
        return self.printReactionSummary()

    @property
    def ec_number(self):
        return self.__ec_number

    @property
    def name(self):
        return self.__name

    @property
    def systematic_name(self):
        return self.__systematic_name

    @property
    def reaction_str(self) -> str:
        reactions = self.__entry.get("reaction", [])
        return reactions[0].get("value", "") if reactions else ""

    @property
    def mechanism(self) -> list[str]:
        return [r.get("value", "") for r in self.__entry.get("reaction", [])]

    @property
    def reaction_type(self) -> list[str]:
        return [r.get("value", "") for r in self.__entry.get("reaction_type", [])]

    @property
    def cofactors(self):
        return self.__getDictOfEnzymeActuators("cofactor")

    @property
    def metals(self):
        return self.__getDictOfEnzymeActuators("metals_ions")

    @property
    def inhibitors(self):
        return self.__getDictOfEnzymeActuators("inhibitor")

    @property
    def activators(self):
        return self.__getDictOfEnzymeActuators("activating_compound")

    @property
    def KMvalues(self):
        return self.__getDictOfEnzymeProperties("km_value")

    @property
    def KIvalues(self):
        return self.__getDictOfEnzymeProperties("ki_value")

    @property
    def KKMvalues(self):
        return self.__getDictOfEnzymeProperties("kcat_km_value")

    @property
    def Kcatvalues(self):
        return self.__getDictOfEnzymeProperties("turnover_number")

    @property
    def specificActivities(self):
        res = []
        for record in self.__entry.get("specific_activity", []):
            number, _ = _split_value_unit(record.get("value", ""))
            res.append(
                {
                    "value": self.__eval_range_value(number),
                    "species": self.__organisms_for(record.get("proteins", [])),
                    "meta": record.get("comment", ""),
                    "refs": self.__citations_for(record.get("references", [])),
                    "specific_info": "",
                }
            )
        return res

    @property
    def substratesAndProducts(self) -> list:
        """
        Returns list of dicts with evaluated "natural" substrates and products
        of the enzyme across organisms.
        """
        substrates, products, res = [], [], []
        for record in self.__entry.get("natural_substrates_products", []):
            rxn = _clean_reaction_string(record.get("value", ""))
            if rxn.count("=") != 1:
                continue
            subs_side, prods_side = rxn.split("=")
            subs = sorted(_split_reaction_side(subs_side))
            prods = sorted(_split_reaction_side(prods_side))
            if subs and prods and subs not in substrates:
                substrates.append(subs)
                products.append(prods)
                res.append({"substrates": subs, "products": prods})
        return res

    @property
    def substrates_products(self) -> list:
        """
        All catalogued substrate/product pairs (BRENDA ``substrates_products``),
        enriched with reversibility, organisms and citations.
        """
        res = []
        for record in self.__entry.get("substrates_products", []):
            value = record.get("value", "")
            # |...| pipe comments are not pulled into the JSON `comment` field by
            # BRENDA, so fold them in here rather than discard them.
            pipe_notes = [
                _PROTEIN_TAG_RE.sub("", _REF_TAG_RE.sub("", note)).strip()
                for note in _PIPE_RE.findall(value)
            ]
            pipe_notes = [note for note in pipe_notes if note]
            comment = record.get("comment", "")
            if pipe_notes:
                extra = "; ".join(pipe_notes)
                comment = f"{comment}; {extra}" if comment else extra
            res.append(
                {
                    "value": _clean_reaction_string(value),
                    "reversibility": self.__reversibility(value),
                    "organisms": self.__organisms_for(record.get("proteins", [])),
                    "comment": comment,
                    "references": self.__citations_for(record.get("references", [])),
                }
            )
        return res

    @property
    def synonyms(self) -> list:
        return [r.get("value", "") for r in self.__entry.get("synonyms", [])]

    @property
    def molecular_weight(self) -> list:
        return self.field("molecular_weight")

    @property
    def subunits(self) -> list:
        return self.field("subunits")

    @property
    def pi_value(self) -> list:
        return self.field("pi_value")

    @property
    def ic50_value(self) -> list:
        return self.field("ic50_value")

    @property
    def source_tissue(self) -> list:
        return self.field("source_tissue")

    @property
    def localization(self) -> list:
        return self.field("localization")

    @property
    def temperature(self):
        return EnzymeConditionDict(
            {
                "optimum": self.__extractTempOrPHData("temperature_optimum", False),
                "range": self.__extractTempOrPHData("temperature_range", True),
                "stability": self.__extractTempOrPHData("temperature_stability", False),
            }
        )

    @property
    def PH(self):
        return EnzymeConditionDict(
            {
                "optimum": self.__extractTempOrPHData("ph_optimum", False),
                "range": self.__extractTempOrPHData("ph_range", True),
                "stability": self.__extractTempOrPHData("ph_stability", False),
            }
        )

    @property
    def proteins(self) -> dict:
        """
        Returns a dict listing all proteins for the given EC number, keyed by
        BRENDA protein id. Preserves the historical ``name``/``proteinID``/``refs``
        keys and adds ``organism``/``accessions``/``source``/``comment``.
        """
        result = {}
        for pid, p in self.__proteins_raw.items():
            accessions = p.get("accessions", [])
            result[pid] = {
                "name": p.get("organism", ""),
                "organism": p.get("organism", ""),
                "proteinID": accessions[0] if accessions else "",
                "accessions": accessions,
                "source": p.get("source", ""),
                "comment": p.get("comment", ""),
                "refs": p.get("references", []),
            }
        return result

    @property
    def organisms(self) -> list:
        """
        Returns a list containing all represented species in the database for this reaction
        """
        organisms = list({p.get("organism", "") for p in self.__proteins_raw.values()})
        organisms.sort()
        return organisms

    @property
    def references(self) -> dict:
        """
        Bibliography cited for this EC number as ``{id: citation_string}``.
        See :pyattr:`bibliography` for the structured records.
        """
        return {
            rid: self.__format_citation(ref)
            for rid, ref in self.__references_raw.items()
        }

    @property
    def bibliography(self) -> dict:
        """
        Structured bibliography as ``{id: {title, authors, journal, year, vol,
        pages, pmid}}``.
        """
        return self.__references_raw
