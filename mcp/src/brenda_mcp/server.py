"""
BRENDA MCP server.

Exposes the BRENDApyrser API as Model Context Protocol tools so that any
MCP-compatible agent (Claude Desktop, Claude Code, or a custom client driving
Claude / DeepSeek / etc.) can query the BRENDA enzyme database in natural
language.

Run it (stdio transport, the default):

    BRENDA_DATABASE_PATH=/path/to/brenda_2026_1.json brenda-mcp
    # or, without installing the package:
    PYTHONPATH=src BRENDA_DATABASE_PATH=... python -m brenda_mcp.server

The database is parsed lazily on the first tool call and cached for the life of
the process.
"""

from __future__ import annotations

from typing import Literal, Optional

from mcp.server.fastmcp import FastMCP

from . import service

INSTRUCTIONS = """\
Tools to query the BRENDA enzyme database (https://www.brenda-enzymes.org) via
the BRENDApyrser parser. Enzymes are identified by EC number (e.g. "2.7.1.40").

Typical workflow:
  * Don't know the EC number? Call `search_enzymes` with a name fragment.
  * `get_enzyme` gives an overview, including a `data_available` map showing
    which kinetic/condition data is populated — use it to decide what to query.
  * `get_enzyme_kinetics` / `get_enzyme_conditions` return summary statistics
    (count, median, quartiles, histogram), optionally filtered by compound and
    organism, rather than raw value dumps.
  * `compute_parameter_distribution` aggregates a parameter across the whole
    database or a genus/species subset (e.g. median optimal temperature in the
    genus Thermotoga vs. the whole database).

Kinetic values carry units: KM/Ki in mM, kcat (turnover number) in 1/s, kcat/KM
in mM^-1 s^-1, specific activity in umol/min/mg, temperature in degC.
"""

mcp = FastMCP("brenda", instructions=INSTRUCTIONS)

# Parameter vocabularies, surfaced to the model as JSON-schema enums.
KineticParam = Literal["km", "kcat", "ki", "kcat_km", "specific_activity"]
ConditionProperty = Literal["temperature", "ph"]
ConditionKind = Literal["optimum", "range", "stability"]
CompoundKind = Literal[
    "cofactors", "inhibitors", "activators", "metals", "substrates_products", "synonyms"
]
CompoundRole = Literal["substrate", "product", "any"]
DistributionParam = Literal[
    "km",
    "kcat",
    "ki",
    "kcat_km",
    "specific_activity",
    "temperature_optimum",
    "ph_optimum",
]


@mcp.tool()
def get_database_info() -> dict:
    """Return BRENDA release, JSON schema version, number of enzyme entries and
    the copyright notice. Useful as a first call to confirm the database is
    loaded and see which release the answers come from."""
    return service.database_info()


@mcp.tool()
def search_enzymes(query: str, limit: int = 25) -> dict:
    """Find enzymes by a substring of their EC number, recommended name, or any
    synonym. Use this when you have an enzyme name (e.g. "pyruvate kinase") but
    not its EC number.

    Args:
        query: Text to search for, case-insensitive (e.g. "alcohol dehydrogenase").
        limit: Maximum number of matches to return.
    """
    return service.search_enzymes(query, limit=limit)


@mcp.tool()
def get_enzyme(ec_number: str) -> dict:
    """Overview of a single enzyme: name, systematic name, catalysed reaction,
    reaction type, synonyms, number of source organisms, and a `data_available`
    map of how many records exist for each kinetic/condition field. Call this
    before the more specific tools to see what is worth querying.

    Args:
        ec_number: EC number, e.g. "2.7.1.40".
    """
    r = service.get_brenda().reactions.get_by_id(ec_number)
    return service.reaction_overview(r)


@mcp.tool()
def get_enzyme_kinetics(
    ec_number: str,
    parameter: KineticParam,
    compound: Optional[str] = None,
    organism: Optional[str] = None,
    max_value: Optional[float] = None,
    include_values: bool = False,
) -> dict:
    """Summary statistics for a kinetic parameter of one enzyme, with an
    optional per-compound breakdown. Returns count, min/max, mean, median and
    quartiles (plus a histogram) — not the raw values, unless `include_values`
    is set.

    Args:
        ec_number: EC number, e.g. "2.7.1.40".
        parameter: One of km, kcat, ki, kcat_km, specific_activity.
        compound: Restrict to a substrate/inhibitor (e.g. "phosphoenolpyruvate").
            Matched case-insensitively, falling back to substring match.
        organism: Restrict to records from organisms whose name contains this
            text (e.g. "Bos taurus", "Escherichia coli").
        max_value: Drop values above this threshold (to trim implausible
            outliers, e.g. max_value=1000 for KM in mM).
        include_values: If true, also return the (capped) list of raw values.
    """
    r = service.get_brenda().reactions.get_by_id(ec_number)
    return service.enzyme_kinetics(
        r,
        parameter,
        compound=compound,
        organism=organism,
        max_value=max_value,
        include_values=include_values,
    )


@mcp.tool()
def get_enzyme_conditions(
    ec_number: str,
    property: ConditionProperty,
    condition: ConditionKind = "optimum",
    organism: Optional[str] = None,
) -> dict:
    """Temperature or pH data for one enzyme. `condition` selects the optimum,
    the active range, or the stability range. Returns summary statistics
    (ranges are returned as [low, high] pairs plus stats on the bounds).

    Args:
        ec_number: EC number, e.g. "2.7.1.40".
        property: "temperature" or "ph".
        condition: "optimum", "range", or "stability".
        organism: Restrict to records from organisms whose name contains this text.
    """
    r = service.get_brenda().reactions.get_by_id(ec_number)
    return service.enzyme_conditions(r, property, condition, organism=organism)


@mcp.tool()
def get_enzyme_compounds(ec_number: str, kind: CompoundKind, limit: int = 100) -> dict:
    """List the cofactors, inhibitors, activators, metals/ions, natural
    substrate-product pairs, or synonyms recorded for one enzyme.

    Args:
        ec_number: EC number, e.g. "2.7.1.40".
        kind: cofactors, inhibitors, activators, metals, substrates_products, or synonyms.
        limit: Maximum number of items to return.
    """
    r = service.get_brenda().reactions.get_by_id(ec_number)
    return service.enzyme_compounds(r, kind, limit=limit)


@mcp.tool()
def get_enzyme_organisms(ec_number: str, limit: int = 100) -> dict:
    """List the source organisms in which this enzyme has been characterised.

    Args:
        ec_number: EC number, e.g. "2.7.1.40".
        limit: Maximum number of organism names to return.
    """
    r = service.get_brenda().reactions.get_by_id(ec_number)
    return service.enzyme_organisms(r, limit=limit)


@mcp.tool()
def get_enzyme_references(ec_number: str, limit: int = 25) -> dict:
    """Return the literature citations (with PubMed IDs where available) for one
    enzyme.

    Args:
        ec_number: EC number, e.g. "2.7.1.40".
        limit: Maximum number of citations to return.
    """
    r = service.get_brenda().reactions.get_by_id(ec_number)
    return service.enzyme_references(r, limit=limit)


@mcp.tool()
def find_enzymes_by_compound(
    compound: str, role: CompoundRole = "any", limit: int = 25
) -> dict:
    """Find enzymes that act on a given compound. `role` restricts to reactions
    where the compound is a substrate, a product, or either. Uses exact
    compound-name matching against BRENDA's natural substrate/product lists.

    Args:
        compound: Exact compound name, e.g. "phosphoenolpyruvate".
        role: "substrate", "product", or "any".
        limit: Maximum number of enzymes to return.
    """
    return service.find_enzymes_by_compound(compound, role=role, limit=limit)


@mcp.tool()
def find_enzymes_by_organism(organism: str, limit: int = 25) -> dict:
    """Find enzymes characterised in a given organism or taxon. Matching is a
    case-insensitive substring of the organism name, so a genus like
    "Thermotoga" returns enzymes from every species in that genus.

    Args:
        organism: Organism or taxon name (e.g. "Escherichia coli", "Thermotoga").
        limit: Maximum number of enzymes to return.
    """
    return service.find_enzymes_by_organism(organism, limit=limit)


@mcp.tool()
def compute_parameter_distribution(
    parameter: DistributionParam,
    organism: Optional[str] = None,
    compound: Optional[str] = None,
    max_value: Optional[float] = None,
    bins: int = 10,
) -> dict:
    """Aggregate one parameter across the whole database, or across a genus /
    species subset, returning summary statistics and a histogram. This answers
    population-level questions such as "what is the median KM in BRENDA?" or "do
    enzymes in the genus Thermotoga have higher optimal temperatures than the
    database as a whole?". Note: this scans many enzymes and may take several
    seconds on the full database.

    Args:
        parameter: km, kcat, ki, kcat_km, specific_activity, temperature_optimum, or ph_optimum.
        organism: Restrict to a taxon (case-insensitive substring, e.g. "Thermotoga").
        compound: Restrict kinetic parameters to a specific substrate/inhibitor.
        max_value: Drop values above this threshold (trims outliers; e.g. 1000 for KM).
        bins: Number of histogram bins.
    """
    return service.parameter_distribution(
        parameter,
        organism=organism,
        compound=compound,
        max_value=max_value,
        bins=bins,
    )


def main() -> None:
    """Console-script entry point: run the server over stdio."""
    import logging
    import os

    # Quiet the per-request INFO chatter from the MCP runtime unless debugging.
    # (Our own progress messages in service._log go straight to stderr.)
    if not os.environ.get("BRENDA_MCP_DEBUG"):
        logging.getLogger("mcp").setLevel(logging.WARNING)
    mcp.run()


if __name__ == "__main__":
    main()
