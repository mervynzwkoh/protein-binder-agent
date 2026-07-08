from urllib import response
import os
import anthropic
from dotenv import load_dotenv
from tools import (
    resolve_target, get_structure, get_known_ligands,
    load_candidate_library, dock_candidates,
)
from report import generate_report

load_dotenv()
client = anthropic.Anthropic()

# Shared state between tool calls - avoids routing large data (sequences,
# candidate SMILES lists) through the LLM's own context window.
_state = {}


def _tool_resolve_target(name_or_id: str) -> dict:
    result = resolve_target(name_or_id)
    _state["target"] = result
    return {"accession": result["accession"], "gene_name": result["gene_name"]}


def _tool_get_structure() -> dict:
    target = _state["target"]
    result = get_structure(target["accession"], target["sequence"])
    _state["structure"] = result
    return result


def _tool_get_known_ligands(max_results: int = 50) -> dict:
    target = _state["target"]
    ligands = get_known_ligands(target["accession"], max_results)
    _state["known_ligands"] = ligands
    return {"count": len(ligands), "sample": ligands[:5]}


def _tool_load_candidate_library(max_candidates: int = 20) -> dict:
    known_ids = {l["chembl_id"] for l in _state.get("known_ligands", [])}
    candidates = load_candidate_library(exclude_chembl_ids=known_ids, max_candidates=max_candidates)
    _state["candidates"] = candidates
    return {"count": len(candidates), "sample": candidates[:5]}


def _tool_dock_candidates() -> list[dict]:
    structure = _state["structure"]
    candidates = _state["candidates"]
    results = dock_candidates(structure["path"], candidates)
    _state["docking_results"] = results
    return results

def _tool_dock_known_ligands() -> list[dict]:
    """Dock the known-ligand calibration set - used to check whether Vina
    actually separates real binders from the candidate pool for this target."""
    structure = _state["structure"]
    known_ligands = _state["known_ligands"]

    # Multiple assay measurements can point at the same molecule - dock each
    # unique compound once, not once per measurement.
    unique_ligands = {l["chembl_id"]: l for l in known_ligands}.values()

    results = dock_candidates(structure["path"], list(unique_ligands))
    _state["calibration_docking_results"] = results
    return results

def _tool_generate_report(narrative: str) -> str:
    """Save the current run's findings as a Markdown report with a 3D pose visualization.
    `narrative` should be your own interpretation/rationale in prose - it will be
    embedded in the report verbatim."""
    return generate_report(_state, narrative)

TOOL_FUNCTIONS = {
    "resolve_target": _tool_resolve_target,
    "get_structure": _tool_get_structure,
    "get_known_ligands": _tool_get_known_ligands,
    "load_candidate_library": _tool_load_candidate_library,
    "dock_candidates": _tool_dock_candidates,
    "dock_known_ligands": _tool_dock_known_ligands,
    "generate_report": _tool_generate_report,
}

TOOLS = [
    {
        "name": "resolve_target",
        "description": "Resolve a gene name or UniProt ID to a canonical target accession.",
        "input_schema": {
            "type": "object",
            "properties": {"name_or_id": {"type": "string"}},
            "required": ["name_or_id"],
        },
    },
    {
        "name": "get_structure",
        "description": "Get a 3D structure for the currently resolved target (experimental PDB, or ESMFold prediction if none exists). Call resolve_target first.",
        "input_schema": {"type": "object", "properties": {}},
    },
    {
        "name": "get_known_ligands",
        "description": "Get known bioactive ligands (IC50 data) for the current target from ChEMBL, as a calibration/validation set.",
        "input_schema": {
            "type": "object",
            "properties": {"max_results": {"type": "integer", "description": "Default 50"}},
        },
    },
    {
        "name": "load_candidate_library",
        "description": "Load a random sample of approved-drug candidates (Lipinski-filtered, excluding any already-known ligands for this target) to screen as repurposing candidates.",
        "input_schema": {
            "type": "object",
            "properties": {"max_candidates": {"type": "integer", "description": "Default 20"}},
        },
    },
    {
        "name": "dock_candidates",
        "description": "Dock the currently loaded candidates against the current target's structure using AutoDock Vina. Requires get_structure and load_candidate_library to have been called first.",
        "input_schema": {"type": "object", "properties": {}},
    },
    {
        "name": "dock_known_ligands",
        "description": "Dock the known-ligand calibration set against the current target's structure. Call this to check whether the docking setup separates real binders from the candidate pool - do this before presenting candidate rankings as meaningful.",
        "input_schema": {"type": "object", "properties": {}},
    },
    {
        "name": "generate_report",
        "description": "Save your findings as a persisted Markdown report with a 3D visualization of the top pose. Call this last, after dock_candidates, passing your own interpretation as the narrative argument.",
        "input_schema": {
            "type": "object",
            "properties": {"narrative": {"type": "string", "description": "Your interpretation and rationale, in prose"}},
            "required": ["narrative"],
        },
    },
]

SYSTEM_PROMPT = """You are a computational drug discovery triage assistant. Given a protein \
target, you autonomously: resolve it, retrieve its structure, pull known ligands as a \
calibration reference, load a candidate pool of approved drugs (repurposing candidates, \
already excluding known binders), and dock them.

Always call dock_known_ligands before presenting the candidate ranking as meaningful - \
compare the two score distributions directly and state whether known ligands score \
notably better than candidates, which is the actual evidence for whether this docking \
setup works for this target.

When reasoning about results: compare candidate docking scores to the known-ligand \
scores as a rough calibration check - if known ligands don't score notably better than \
random, treat that as a signal the docking setup may be unreliable for this target, and \
say so explicitly rather than presenting the ranking uncritically. State whether the \
structure used was experimental or ESMFold-predicted, since that materially affects \
confidence. Always end with a ranked shortlist and a one-sentence rationale per candidate, \
plus an explicit caveat that docking scores are a coarse triage signal, not proof of binding.

Once you've reasoned through the calibration check and ranking, always finish by calling generate_report \
with your full interpretation as the narrative argument - this is where your complete \
analysis belongs. Your final chat response after that should be brief: 2-3 sentences \
confirming what you found and pointing to the saved report, not a repeat of the full \
analysis."""


def run_agent(user_prompt: str) -> str:
    global _state
    _state = {}  # reset state for a fresh run
    messages = [{"role": "user", "content": user_prompt}]

    while True:
        response = client.messages.create(
            model="claude-sonnet-5",
            max_tokens=4096,
            system=SYSTEM_PROMPT,
            tools=TOOLS,
            messages=messages,
        )
        messages.append({"role": "assistant", "content": response.content})
        print(f"[agent] stop_reason={response.stop_reason}, output_tokens={response.usage.output_tokens}")

        if response.stop_reason != "tool_use":
            return "".join(block.text for block in response.content if block.type == "text")

        tool_results = []
        for block in response.content:
            if block.type == "tool_use":
                print(f"[agent] calling {block.name}({block.input})")
                try:
                    result = TOOL_FUNCTIONS[block.name](**block.input)
                    content = str(result)
                except Exception as e:
                    content = f"ERROR: {e}"
                tool_results.append({"type": "tool_result", "tool_use_id": block.id, "content": content})
        messages.append({"role": "user", "content": tool_results})


if __name__ == "__main__":
    print(run_agent("Find and rank the best repurposing candidates for EGFR."))