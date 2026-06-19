#!/usr/bin/env python
"""
No-API-key smoke test for the BRENDA MCP server.

Launches the server over stdio against the small committed test fixture
(``tests/fixtures/brenda_sample.json``, EC 1.1.1.304 + 6.6.99.99), then connects
as an MCP client and exercises every tool, asserting on known fixture values. No
LLM and no API keys are involved — this validates the MCP plumbing and the tool
logic end to end.

    python tests/smoke_test.py
"""

from __future__ import annotations

import asyncio
import json
import sys
from datetime import timedelta
from pathlib import Path

from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client

HERE = Path(__file__).resolve().parent
MCP_DIR = HERE.parent
SRC_DIR = MCP_DIR / "src"
REPO_ROOT = MCP_DIR.parent
FIXTURE = REPO_ROOT / "tests" / "fixtures" / "brenda_sample.json"

TIMEOUT = timedelta(seconds=120)


async def call(session: ClientSession, name: str, **arguments):
    result = await session.call_tool(name, arguments, read_timeout_seconds=TIMEOUT)
    assert not getattr(
        result, "isError", False
    ), f"{name} returned an error: {result.content}"
    text = "\n".join(
        getattr(b, "text", "") for b in result.content if getattr(b, "text", "")
    )
    return json.loads(text)


def check(label: str, condition: bool, detail: str = "") -> None:
    status = "PASS" if condition else "FAIL"
    print(f"  [{status}] {label}" + (f" — {detail}" if detail else ""))
    if not condition:
        raise AssertionError(f"{label}: {detail}")


async def main() -> None:
    if not FIXTURE.exists():
        sys.exit(f"Fixture not found: {FIXTURE}")

    server_params = StdioServerParameters(
        command=sys.executable,
        args=["-m", "brenda_mcp.server"],
        env={
            "PYTHONPATH": str(SRC_DIR),
            "BRENDA_DATABASE_PATH": str(FIXTURE),
            "PATH": __import__("os").environ.get("PATH", ""),
        },
    )

    async with stdio_client(server_params) as (read, write):
        async with ClientSession(read, write) as session:
            await session.initialize()
            tools = (await session.list_tools()).tools
            print(
                f"Server exposes {len(tools)} tools: {', '.join(t.name for t in tools)}\n"
            )
            check("expected tool count", len(tools) == 11, f"got {len(tools)}")

            info = await call(session, "get_database_info")
            check("database release", info["release"] == "2026.1", info["release"])
            check("enzyme count", info["n_enzymes"] == 2, str(info["n_enzymes"]))

            found = await call(session, "search_enzymes", query="reductase")
            ecs = {h["ec_number"] for h in found["results"]}
            check("search finds 1.1.1.304", "1.1.1.304" in ecs, str(ecs))

            ov = await call(session, "get_enzyme", ec_number="1.1.1.304")
            check("overview name", "reductase" in ov["name"].lower(), ov["name"])
            check(
                "overview reports km compounds",
                ov["data_available"]["km_compounds"] > 0,
                str(ov["data_available"]["km_compounds"]),
            )

            km = await call(
                session, "get_enzyme_kinetics", ec_number="1.1.1.304", parameter="km"
            )
            check(
                "km stats present",
                km["stats"]["count"] > 0,
                str(km["stats"].get("count")),
            )
            check("km unit is mM", km["unit"] == "mM", km["unit"])
            check(
                "km by_compound present",
                len(km["by_compound"]) > 0,
                f"{len(km['by_compound'])} compounds",
            )

            km_nadh = await call(
                session,
                "get_enzyme_kinetics",
                ec_number="1.1.1.304",
                parameter="km",
                compound="NADH",
                include_values=True,
            )
            check(
                "km filtered by NADH",
                km_nadh["stats"]["count"] > 0,
                str(km_nadh["stats"].get("count")),
            )
            check(
                "km values returned when requested",
                "values" in km_nadh,
                str(km_nadh.get("values", [])[:3]),
            )

            topt = await call(
                session,
                "get_enzyme_conditions",
                ec_number="1.1.1.304",
                property="temperature",
                condition="optimum",
            )
            check(
                "temperature optimum count == 3",
                topt["stats"]["count"] == 3,
                str(topt["stats"].get("count")),
            )
            check(
                "temperature optimum median == 30",
                topt["stats"]["median"] == 30,
                str(topt["stats"].get("median")),
            )

            cof = await call(
                session, "get_enzyme_compounds", ec_number="1.1.1.304", kind="cofactors"
            )
            cof_names = {c["compound"] for c in cof["items"]}
            check("cofactors include NADH", "NADH" in cof_names, str(cof_names))

            orgs = await call(session, "get_enzyme_organisms", ec_number="1.1.1.304")
            check("organisms listed", orgs["total"] > 0, str(orgs["total"]))

            refs = await call(session, "get_enzyme_references", ec_number="1.1.1.304")
            check("references listed", refs["total"] > 0, str(refs["total"]))

            by_org = await call(
                session, "find_enzymes_by_organism", organism="Staphylococcus"
            )
            by_org_ecs = {h["ec_number"] for h in by_org["results"]}
            check(
                "find_enzymes_by_organism(Staphylococcus) finds 1.1.1.304",
                "1.1.1.304" in by_org_ecs,
                str(by_org_ecs),
            )

            dist = await call(
                session,
                "compute_parameter_distribution",
                parameter="temperature_optimum",
            )
            check(
                "distribution over fixture has values",
                dist["stats"]["count"] >= 3,
                str(dist["stats"].get("count")),
            )
            check(
                "distribution scanned 2 enzymes",
                dist["enzymes_scanned"] == 2,
                str(dist["enzymes_scanned"]),
            )

    print("\nAll smoke-test checks passed.")


if __name__ == "__main__":
    asyncio.run(main())
