#!/usr/bin/env python
"""
Example: an LLM agent that answers enzyme-kinetics questions through the BRENDA
MCP server.

It reproduces the kind of analyses in ``docs/examples.ipynb`` — looking up
pyruvate kinase's KM for phosphoenolpyruvate, and comparing the optimal
temperature of a thermophilic genus against the whole database — but instead of
hand-written pandas/matplotlib code, an LLM decides which BRENDA tools to call.

The same MCP toolbox is driven by two providers so you can compare them:

* Claude   — via the official `anthropic` SDK (`claude-opus-4-8` by default)
* DeepSeek — via the OpenAI-compatible `openai` SDK

Both connect to the *same* MCP server (``brenda_mcp.server``) launched as a
subprocess over stdio; only the model and the tool-schema adapter differ.

Usage:
    cp .env.example .env      # then fill in your API keys
    python examples/enzyme_kinetics_agent.py --provider claude
    python examples/enzyme_kinetics_agent.py --provider deepseek
    python examples/enzyme_kinetics_agent.py --provider claude --question "..."
"""

from __future__ import annotations

import argparse
import asyncio
import json
import os
import sys
from datetime import timedelta
from pathlib import Path

from dotenv import load_dotenv
from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client

HERE = Path(__file__).resolve().parent
MCP_DIR = HERE.parent
SRC_DIR = MCP_DIR / "src"
ENV_PATH = MCP_DIR / ".env"

load_dotenv(ENV_PATH)

# Per-tool-call timeout: the first call parses the (large) BRENDA database, and
# compute_parameter_distribution scans the whole database, so be generous.
TOOL_TIMEOUT = timedelta(seconds=600)
MAX_TURNS = 14  # guard against runaway tool loops

SYSTEM_PROMPT = (
    "You are an enzymology research assistant with access to the BRENDA enzyme "
    "database through tools. Answer the user's question using those tools rather "
    "than prior knowledge. Identify enzymes by EC number, report quantitative "
    "results with their units (KM/Ki in mM, kcat in 1/s, temperature in degC), "
    "and state the sample sizes (n) behind any median or mean. If a value is a "
    "database-wide or genus-wide aggregate, say so. Be concise."
)

DEFAULT_QUESTION = (
    "I'm studying pyruvate kinase (EC 2.7.1.40). "
    "(1) What is the typical KM for phosphoenolpyruvate across organisms, and "
    "how does it compare to the KM measured specifically in Bos taurus? "
    "(2) Separately: is it true that enzymes from the thermophilic genus "
    "Thermotoga have a higher optimal temperature than enzymes in BRENDA as a "
    "whole? Use the data to support your answer."
)


# --------------------------------------------------------------------------- #
# MCP toolbox: connect to the server and adapt its tools to each provider     #
# --------------------------------------------------------------------------- #
class MCPToolbox:
    def __init__(self, session: ClientSession):
        self.session = session
        self.tools = []

    async def load(self) -> None:
        self.tools = (await self.session.list_tools()).tools

    async def call(self, name: str, arguments: dict) -> str:
        result = await self.session.call_tool(
            name, arguments, read_timeout_seconds=TOOL_TIMEOUT
        )
        parts = [
            getattr(b, "text", "") for b in result.content if getattr(b, "text", "")
        ]
        text = "\n".join(parts).strip()
        if not text and getattr(result, "structuredContent", None):
            text = json.dumps(result.structuredContent)
        if getattr(result, "isError", False):
            return f"ERROR: {text or '(unknown tool error)'}"
        return text or "(no content)"

    def anthropic_tools(self) -> list[dict]:
        return [
            {
                "name": t.name,
                "description": t.description or "",
                "input_schema": t.inputSchema,
            }
            for t in self.tools
        ]

    def openai_tools(self) -> list[dict]:
        return [
            {
                "type": "function",
                "function": {
                    "name": t.name,
                    "description": t.description or "",
                    "parameters": t.inputSchema,
                },
            }
            for t in self.tools
        ]


def _log_tool_call(name: str, arguments: dict) -> None:
    print(f"  → tool: {name}({json.dumps(arguments, ensure_ascii=False)})")


# --------------------------------------------------------------------------- #
# Claude backend (official anthropic SDK, manual agentic loop)                #
# --------------------------------------------------------------------------- #
async def run_claude(toolbox: MCPToolbox, question: str) -> str:
    from anthropic import AsyncAnthropic

    model = os.environ.get("ANTHROPIC_MODEL", "claude-opus-4-8")
    effort = os.environ.get("ANTHROPIC_EFFORT", "medium")
    client = AsyncAnthropic()  # reads ANTHROPIC_API_KEY from the environment
    tools = toolbox.anthropic_tools()
    messages: list[dict] = [{"role": "user", "content": question}]

    for _ in range(MAX_TURNS):
        resp = await client.messages.create(
            model=model,
            max_tokens=16000,
            system=SYSTEM_PROMPT,
            tools=tools,
            # Adaptive thinking lets Claude decide how much to reason between
            # tool calls; effort trades thoroughness against token cost.
            thinking={"type": "adaptive", "display": "summarized"},
            output_config={"effort": effort},
            messages=messages,
        )

        for block in resp.content:
            if block.type == "thinking" and getattr(block, "thinking", ""):
                print(f"  [thinking] {block.thinking.strip()[:300]}")
            elif block.type == "text" and block.text.strip():
                print(f"  [claude] {block.text.strip()}")

        if resp.stop_reason != "tool_use":
            return "".join(b.text for b in resp.content if b.type == "text").strip()

        messages.append({"role": "assistant", "content": resp.content})
        tool_results = []
        for block in resp.content:
            if block.type == "tool_use":
                _log_tool_call(block.name, block.input)
                output = await toolbox.call(block.name, block.input)
                tool_results.append(
                    {"type": "tool_result", "tool_use_id": block.id, "content": output}
                )
        messages.append({"role": "user", "content": tool_results})

    return "(stopped: reached the maximum number of tool-use turns)"


# --------------------------------------------------------------------------- #
# DeepSeek backend (OpenAI-compatible API, manual function-calling loop)       #
# --------------------------------------------------------------------------- #
async def run_deepseek(toolbox: MCPToolbox, question: str) -> str:
    from openai import AsyncOpenAI

    model = os.environ.get("DEEPSEEK_MODEL", "deepseek-v4-pro")
    client = AsyncOpenAI(
        api_key=os.environ["DEEPSEEK_API_KEY"],
        base_url=os.environ.get("DEEPSEEK_BASE_URL", "https://api.deepseek.com"),
    )
    tools = toolbox.openai_tools()
    messages: list[dict] = [
        {"role": "system", "content": SYSTEM_PROMPT},
        {"role": "user", "content": question},
    ]

    for _ in range(MAX_TURNS):
        resp = await client.chat.completions.create(
            model=model, messages=messages, tools=tools
        )
        msg = resp.choices[0].message

        if msg.content and msg.content.strip():
            print(f"  [deepseek] {msg.content.strip()}")

        assistant: dict = {"role": "assistant", "content": msg.content or ""}
        if msg.tool_calls:
            assistant["tool_calls"] = [
                {
                    "id": tc.id,
                    "type": "function",
                    "function": {
                        "name": tc.function.name,
                        "arguments": tc.function.arguments,
                    },
                }
                for tc in msg.tool_calls
            ]
        messages.append(assistant)

        if not msg.tool_calls:
            return (msg.content or "").strip()

        for tc in msg.tool_calls:
            try:
                arguments = json.loads(tc.function.arguments or "{}")
            except json.JSONDecodeError:
                arguments = {}
            _log_tool_call(tc.function.name, arguments)
            output = await toolbox.call(tc.function.name, arguments)
            messages.append({"role": "tool", "tool_call_id": tc.id, "content": output})

    return "(stopped: reached the maximum number of tool-use turns)"


# --------------------------------------------------------------------------- #
# Driver                                                                       #
# --------------------------------------------------------------------------- #
async def main_async(provider: str, question: str) -> None:
    if not os.environ.get("BRENDA_DATABASE_PATH"):
        sys.exit(
            "BRENDA_DATABASE_PATH is not set. Copy .env.example to .env and point "
            "it at your BRENDA JSON file."
        )

    # Launch the MCP server as a subprocess, using the *same* interpreter (so it
    # shares this environment's `brendapyrser` + `mcp`), with src/ on PYTHONPATH.
    server_params = StdioServerParameters(
        command=sys.executable,
        args=["-m", "brenda_mcp.server"],
        env={**os.environ, "PYTHONPATH": str(SRC_DIR)},
    )

    async with stdio_client(server_params) as (read, write):
        async with ClientSession(read, write) as session:
            await session.initialize()
            toolbox = MCPToolbox(session)
            await toolbox.load()

            print(
                f"Connected to BRENDA MCP server: {len(toolbox.tools)} tools available"
            )
            print(f"Provider: {provider}\n")
            print(f"Question:\n  {question}\n")
            print("--- agent trace ---")

            try:
                if provider == "claude":
                    answer = await run_claude(toolbox, question)
                else:
                    answer = await run_deepseek(toolbox, question)
            except Exception as exc:  # noqa: BLE001 — surface a clean message, not a stack trace
                print(f"\nLLM call failed ({type(exc).__name__}): {exc}")
                return

            print("\n--- final answer ---")
            print(answer)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--provider",
        choices=["claude", "deepseek"],
        default="claude",
        help="Which LLM backend to use (default: claude).",
    )
    parser.add_argument(
        "--question",
        default=DEFAULT_QUESTION,
        help="Question to ask the agent (defaults to the pyruvate-kinase example).",
    )
    args = parser.parse_args()

    key_var = "ANTHROPIC_API_KEY" if args.provider == "claude" else "DEEPSEEK_API_KEY"
    if not os.environ.get(key_var):
        sys.exit(f"{key_var} is not set. Add it to {ENV_PATH} (see .env.example).")

    asyncio.run(main_async(args.provider, args.question))


if __name__ == "__main__":
    main()
