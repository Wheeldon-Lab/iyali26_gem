# GEM Curation MCP Server

A [FastMCP](https://github.com/jlowin/fastmcp) server that exposes
[COBRApy](https://cobrapy.readthedocs.io/) diagnostics for genome-scale metabolic
models (GEMs) as MCP tools, so Claude can inspect and reason about a model
directly. It is the **diagnose** half of an eventual *diagnose → repair → validate*
loop for autonomous model curation.

It defaults to the canonical **iYali26** model (*Yarrowia lipolytica*),
`model.xml` in the project root, but can load any SBML (`.xml`) or JSON
(`.json`) model.

> The model files are **not** part of this repository — they live in the GEM
> project directory. `load_model` resolves a relative path against that project
> root, so point it at your local `model.xml` (or any model you output
> from a curation run).

The loaded model is **cached in server state**, so every tool operates on the live
`cobra.Model` object rather than re-parsing the SBML on each call. Every tool
returns structured JSON (never printed text), and the "no model loaded yet" case
comes back as a clean error dict.

## Tools

| Tool | What it does |
|------|--------------|
| `load_model(path="model.xml")` | Load a model, cache it, return a summary (id, reaction/metabolite/gene counts, objective, compartments, and **annotation coverage** — how many metabolites carry a formula/charge). Relative paths resolve against the project root. |
| `list_mass_charge_imbalances(include_boundary=False, limit=None)` | The core diagnostic. Every reaction whose mass or charge doesn't balance, with the per-element and per-charge discrepancy. Boundary/exchange reactions are excluded by default (they are intentionally open). Also returns a **`diagnosability`** block (see below). |
| `get_reaction(reaction_id)` | Full detail on one reaction: stoichiometry, each metabolite's formula/charge/compartment, bounds, reversibility, GPR, genes, subsystem, and its own mass/charge balance. |
| `get_metabolite(metabolite_id)` | A metabolite's formula, charge, compartment, and every reaction it participates in (with coefficient and reactant/product role). |
| `run_fba()` | Flux balance analysis: objective value and the nonzero fluxes, sorted by descending magnitude. |

## Install

The project already has **COBRApy** (0.30.0) and the **GLPK** solver. Only
`fastmcp` needs installing into the same Python environment:

```bash
/opt/miniconda3/bin/pip install fastmcp
```

(Verified with `fastmcp` 3.4.2, `cobra` 0.30.0, Python 3.13.)

## Run

Directly, for a quick check (it speaks MCP over stdio and will wait for a client):

```bash
/opt/miniconda3/bin/python "/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem/mcp_server/gem_mcp_server.py"
```

In normal use you don't run it by hand — Claude Desktop launches it for you (below).

## Connect to Claude Desktop

Claude Desktop reads MCP servers from `claude_desktop_config.json`:

- **macOS:** `~/Library/Application Support/Claude/claude_desktop_config.json`
- **Windows:** `%APPDATA%\Claude\claude_desktop_config.json`

Add (or merge) this `gem-curation` entry. Both paths are absolute so Claude can
launch the server from anywhere:

```json
{
  "mcpServers": {
    "gem-curation": {
      "command": "/opt/miniconda3/bin/python",
      "args": [
        "/Users/david/Desktop/Lab/Ian wheeldon/code/iyali26_gem/mcp_server/gem_mcp_server.py"
      ]
    }
  }
}
```

Then **fully quit and reopen Claude Desktop**. The five tools appear under the
tools (🔨) menu. Ask Claude to `load_model` first, then e.g. "list the mass/charge
imbalances" or "show me reaction R5".

> Using the standalone `fastmcp` package means `command` must point at the Python
> interpreter that has `fastmcp` installed — here `/opt/miniconda3/bin/python`.
> If you install it elsewhere, update `command` to match.

## Diagnosability

`check_mass_balance()` can only balance a reaction when **every** metabolite in
it has a formula (for mass) and a charge (for charge). When a formula is missing
it returns an empty result — indistinguishable from "balanced". So a count of
`0` imbalances can mean *either* "clean" *or* "couldn't check anything".

To keep the agent honest, `load_model` returns `annotation_coverage` (how many
metabolites have a formula / charge), and `list_mass_charge_imbalances` returns
a `diagnosability` block:

```
"diagnosability": {
  "num_reactions_checked":   2102,
  "num_mass_checkable":      0,      # every reaction has a formula-less metabolite
  "num_mass_uncheckable":    2102,
  "num_charge_checkable":    2102,
  "num_charge_uncheckable":  0,
  "warning": "Mass-balance diagnosis is meaningless for this model: ..."
}
```

**Always read `diagnosability.warning` before trusting a low imbalance count.**

## Notes

- On load, COBRApy prints warnings to **stderr** for metabolites whose formulas
  contain `*` or parentheses (polymer / pseudo-formulas such as biomass
  constituents). These are harmless, do not reach the agent, and do not affect
  the tool results. For any reaction that can't be element-balanced because of
  such a formula, the imbalance is reported via a `balance_error` field instead
  of crashing.
- **`model.xml`** (the default) is the canonical, deterministic iYali26 baseline.
  The current baseline has 2,300 reactions, 1,865 metabolites, and 1,073 genes.
  Formula coverage is incomplete, so always inspect the returned diagnosability
  counts before interpreting a low imbalance count.

## Not yet implemented

By design, this first pass is diagnosis only. Repair (editing charges/formulas/
stoichiometry) and validation (Memote scoring) tools come next, once this runs
against the real model.
