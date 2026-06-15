# BRENDA MCP server

A [Model Context Protocol](https://modelcontextprotocol.io) server that exposes
the [BRENDApyrser](../) API as tools, so any MCP-compatible agent — Claude
Desktop, Claude Code, or a custom client driving Claude / DeepSeek / etc. — can
query the BRENDA enzyme database in natural language.

It wraps `brendapyrser.BRENDA` and returns **compact summaries** (counts,
medians, quartiles, histograms) rather than raw value dumps, keeping tool
results cheap for an LLM to read.

```
mcp/
├── pyproject.toml            # installable package: `brenda-mcp`
├── requirements.txt          # server + example deps
├── .env.example              # copy to .env and fill in keys (gitignored)
├── src/brenda_mcp/
│   ├── server.py             # FastMCP server + tool definitions
│   └── service.py            # DB loader, stats/formatting helpers
├── examples/
│   └── enzyme_kinetics_agent.py   # LLM agent (Claude + DeepSeek) demo
└── tests/
    └── smoke_test.py         # no-API-key end-to-end check against the fixture
```

## Tools

| Tool | What it returns |
|------|-----------------|
| `get_database_info` | BRENDA release, schema version, number of enzymes |
| `search_enzymes` | Enzymes matching a name/synonym/EC substring |
| `get_enzyme` | Overview of one enzyme + a `data_available` map |
| `get_enzyme_kinetics` | KM / kcat / Ki / kcat·KM / specific-activity stats (by compound, by organism) |
| `get_enzyme_conditions` | Temperature / pH optimum, range, or stability |
| `get_enzyme_compounds` | Cofactors / inhibitors / activators / metals / substrate-product pairs / synonyms |
| `get_enzyme_organisms` | Source organisms for an enzyme |
| `get_enzyme_references` | Literature citations (with PubMed IDs) |
| `find_enzymes_by_compound` | Enzymes acting on a compound (as substrate/product/either) |
| `find_enzymes_by_organism` | Enzymes characterised in a taxon (e.g. a genus) |
| `compute_parameter_distribution` | Database-wide or genus-wide aggregate of a parameter |

## Prerequisites

1. **The BRENDA database (JSON).** Due to BRENDA's license it cannot be shipped;
   download `brenda_<release>.json.tar.gz` from
   <https://www.brenda-enzymes.org/download.php>. The server accepts the
   `.json`, `.json.gz`, or `.json.tar.gz` form directly.
2. **A Python ≥ 3.10 environment with `brendapyrser` installed.** On this
   machine that is the `brendapyrser-dev` conda env (Python 3.11):
   `/home/robaina/miniconda3/envs/brendapyrser-dev/bin/python`.

## Install

Install the server (and example) dependencies into that environment:

```bash
cd mcp
pip install -r requirements.txt          # mcp, brendapyrser, anthropic, openai, python-dotenv
# optional — register the `brenda-mcp` console script:
pip install -e .
```

## Configure

```bash
cp .env.example .env        # .env is gitignored
```

Edit `.env`:

```ini
BRENDA_DATABASE_PATH=/abs/path/to/brenda_2026_1.json
ANTHROPIC_API_KEY=sk-ant-...
DEEPSEEK_API_KEY=sk-...
```

## Verify (no API keys needed)

Runs the server against the small committed fixture and exercises every tool:

```bash
python tests/smoke_test.py
# → "All smoke-test checks passed."
```

## Run the example agent

An LLM decides which BRENDA tools to call to answer an enzyme-kinetics question
(the default question mirrors `docs/examples.ipynb`: pyruvate kinase's KM for
phosphoenolpyruvate, and Thermotoga's optimal temperature vs. the database).

```bash
python examples/enzyme_kinetics_agent.py --provider claude
python examples/enzyme_kinetics_agent.py --provider deepseek
python examples/enzyme_kinetics_agent.py --provider claude --question "Which enzymes use phosphoenolpyruvate as a substrate?"
```

The example launches the MCP server itself (as a stdio subprocess using the same
Python interpreter), connects as an MCP client, converts the MCP tool schemas to
each provider's tool format, and runs a manual tool-use loop. The first tool
call parses the full database (can take a while); subsequent calls are fast.

> **Models.** Claude defaults to `claude-opus-4-8` with adaptive thinking;
> DeepSeek defaults to `deepseek-v4-pro`. Override via `ANTHROPIC_MODEL` /
> `DEEPSEEK_MODEL` in `.env`.

## Register with an MCP client

### Claude Desktop / Claude Code

Add to your MCP config (`claude_desktop_config.json`, or Claude Code's
`.mcp.json`). Use the absolute interpreter path so `brendapyrser` + `mcp` are on
the path, and put `src/` on `PYTHONPATH`:

```json
{
  "mcpServers": {
    "brenda": {
      "command": "/home/robaina/miniconda3/envs/brendapyrser-dev/bin/python",
      "args": ["-m", "brenda_mcp.server"],
      "env": {
        "PYTHONPATH": "/home/robaina/Documents/Hapdera/Hapdera-Projects/BRENDApyrser/mcp/src",
        "BRENDA_DATABASE_PATH": "/home/robaina/Documents/Hapdera/Hapdera-Projects/BRENDApyrser/data/brenda_2026_1.json"
      }
    }
  }
}
```

If you ran `pip install -e .`, you can use the console script instead:
`"command": ".../bin/brenda-mcp"` with no `args`/`PYTHONPATH` (still set
`BRENDA_DATABASE_PATH`).

## Notes

- **Lazy loading.** The database is parsed on the first tool call and cached for
  the life of the process. Expect the first call to take a while on the full
  709 MB JSON; set `BRENDA_DATABASE_PATH` to the fixture for fast iteration.
- **stdout is sacred.** The stdio transport uses stdout for protocol traffic;
  all server logging goes to stderr. Set `BRENDA_MCP_DEBUG=1` to restore the
  MCP runtime's verbose per-request logging.
- **DeepSeek** is reached through its OpenAI-compatible endpoint, so the example
  uses the `openai` SDK pointed at `DEEPSEEK_BASE_URL`. The same MCP server and
  tool set serve both providers unchanged.
