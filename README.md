# Target-to-Binder Triage Agent

An autonomous Claude agent that triages a druggable protein target end-to-end: resolves it, retrieves its structure (experimental or ESMFold-predicted), pulls known bioactive ligands as a calibration set, screens a randomly-sampled pool of approved drugs via AutoDock Vina, and produces a ranked, evidence-checked shortlist with an explicit self-assessment of its own reliability.

Built to demonstrate agent orchestration (Claude tool-use) applied to a real, methodologically-grounded computational biology pipeline.

## How it works

The agent doesn't follow a fixed script — it decides which tools to call and in what order, based on a system prompt that requires it to validate its own pipeline before presenting results:

1. **Resolve target** — gene name/UniProt ID → canonical accession + sequence (UniProt)
2. **Get structure** — experimental PDB structure if one exists, otherwise an ESMFold prediction
3. **Get known ligands** — real bioactive compounds for this target from ChEMBL (the calibration set)
4. **Load candidate library** — a randomly-sampled, Lipinski-filtered pool of approved drugs, explicitly excluding anything already known to bind this target
5. **Dock known ligands** — validates whether AutoDock Vina actually separates real binders from background for this specific target
6. **Dock candidates** — screens the repurposing candidates against the same receptor
7. **Generate report** — a saved Markdown report plus an interactive 3D visualization of the top pose

## Example output

See `examples/EGFR_report.md` for a full run against EGFR, including the calibration comparison and ranked shortlist.

## Quickstart

```powershell
git clone https://github.com/mervynzwkoh/protein-binder-agent/tree/main
cd protein-binder-agent
python -m venv venv
.\venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

Download AutoDock Vina's Windows executable from the [official releases page](https://github.com/ccsb-scripps/AutoDock-Vina/releases) and place it at `tools/vina.exe`.

## Usage

Run the agent on any protein target, from the project root, as a module:

```powershell
python -m src.agent EGFR
python -m src.agent DRD2
python -m src.agent P00533
```

The target can be a gene name or a UniProt accession. If you omit it, you'll be prompted interactively:

```powershell
python -m src.agent
Enter a gene name or UniProt ID: EGFR
```

For usage help:

```powershell
python -m src.agent --help
```

Each run saves a Markdown report and an interactive 3D HTML visualization to `data/reports/`.

## Key design decisions

- **Calibration before candidates**: known ligands are docked as a validation set, never as screening candidates — this checks whether the docking setup can actually distinguish real binders before trusting any candidate ranking.
- **Unbiased candidate sampling**: candidates are randomly sampled (fixed seed for reproducibility) from ChEMBL's approved-drug pool, rather than selected for similarity to known binders, to avoid circularity.
- **Shared state, not LLM-routed data**: large molecular data (sequences, SMILES lists, docking results) stays in a Python-side state dict; Claude only ever passes small identifiers between tool calls.

## Limitations

- AutoDock Vina scores are a coarse triage signal — useful for relative ranking within one run, not a validated measure of binding affinity.
- Candidate pool reflects ChEMBL's `max_phase=4` annotation (approved *somewhere*), not FDA approval specifically.
- ESMFold-predicted structures carry materially lower confidence than experimental ones, particularly for binding-site geometry; the report states which was used.
- Docking with no co-crystallized ligand falls back to blind docking over the whole protein, which is slower and less reliable than pocket-centered docking.

## Roadmap

- Validate against additional targets beyond EGFR
- Integrate NVIDIA BioNeMo Agent Toolkit for production-grade structure/docking tools
- Add a literature sub-agent (PubMed abstract retrieval) to enrich rationale

## Tech stack

Anthropic API (tool use) · UniProt · RCSB PDB · ESM Atlas (ESMFold) · ChEMBL (`chembl_webresource_client`) · RDKit · Open Babel · AutoDock Vina · py3Dmol
