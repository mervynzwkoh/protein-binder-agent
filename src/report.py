import os
from datetime import datetime
import py3Dmol

REPORTS_DIR = "data/reports"


def generate_report(state: dict, narrative: str) -> str:
    """Build a Markdown report from pipeline state plus the agent's own narrative.

    `state` holds the shared pipeline state (target, structure, known_ligands,
    candidates, calibration_docking_results, docking_results). `narrative` is the
    agent's free-text interpretation, included verbatim rather than regenerated here -
    the report should reflect the same reasoning the agent already did, not a second
    independent summary that could drift from it.
    """
    os.makedirs(REPORTS_DIR, exist_ok=True)
    target = state["target"]
    structure = state["structure"]
    calibration = sorted(state.get("calibration_docking_results", []), key=lambda r: r["affinity_kcal_mol"])
    candidates = sorted(state.get("docking_results", []), key=lambda r: r["affinity_kcal_mol"])

    lines = [
        f"# Target-to-Binder Triage Report: {target['gene_name']}",
        f"*Generated {datetime.now().strftime('%Y-%m-%d %H:%M')}*\n",
        f"**UniProt accession:** {target['accession']}  ",
        f"**Structure source:** {structure['source']}"
        + (f" (PDB {structure['pdb_id']})" if structure.get("pdb_id") else "") + "\n",
        "## Calibration: known ligands vs. candidates\n",
        "| Group | n | Best score (kcal/mol) | Worst score |",
        "|---|---|---|---|",
    ]
    if calibration:
        lines.append(f"| Known ChEMBL ligands | {len(calibration)} | {calibration[0]['affinity_kcal_mol']:.2f} | {calibration[-1]['affinity_kcal_mol']:.2f} |")
    if candidates:
        lines.append(f"| Repurposing candidates | {len(candidates)} | {candidates[0]['affinity_kcal_mol']:.2f} | {candidates[-1]['affinity_kcal_mol']:.2f} |")

    lines += ["", "## Agent interpretation\n", narrative, ""]

    lines += [
        "## Ranked candidate shortlist\n",
        "| Rank | ChEMBL ID | Name | Affinity (kcal/mol) |",
        "|---|---|---|---|",
    ]
    for i, c in enumerate(candidates[:10], start=1):
        lines.append(f"| {i} | {c['chembl_id']} | {c.get('name', '-')} | {c['affinity_kcal_mol']:.2f} |")

    lines += ["", "## Limitations\n",
              "- Docking scores (AutoDock Vina) are a coarse triage signal, not validated binding affinities.",
              "- Candidate pool is a random sample (fixed seed) of ChEMBL max_phase=4 molecules, Lipinski-filtered."]
    if structure["source"] != "experimental":
        lines.append("- Structure is an ESMFold prediction, not experimental - lower confidence in binding-site geometry.")

    if candidates:
        top = candidates[0]
        viz_path = generate_pose_visualization(structure["path"], top["pose_path"], target["gene_name"])
        lines += ["", f"## Top pose visualization\n", f"See `{os.path.basename(viz_path)}` (open in a browser)."]

    report_path = os.path.join(REPORTS_DIR, f"{target['gene_name']}_report.md")
    with open(report_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    return report_path


def generate_pose_visualization(receptor_path: str, pose_pdbqt_path: str, label: str) -> str:
    """Render the receptor + top docked pose as a standalone, browser-openable HTML file."""
    view = py3Dmol.view(width=800, height=600)
    with open(receptor_path) as f:
        view.addModel(f.read(), "pdb")
    view.setStyle({"model": 0}, {"cartoon": {"color": "spectrum"}})

    with open(pose_pdbqt_path) as f:
        view.addModel(f.read(), "pdbqt")
    view.setStyle({"model": 1}, {"stick": {"colorscheme": "greenCarbon"}})
    view.zoomTo()

    html_path = os.path.join(REPORTS_DIR, f"{label}_top_pose.html")
    with open(html_path, "w", encoding="utf-8") as f:
        f.write(view._make_html())
    return html_path