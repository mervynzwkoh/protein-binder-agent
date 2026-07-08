import requests
import os
import random
from chembl_webresource_client.new_client import new_client
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.MolStandardize import rdMolStandardize
import subprocess
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem
from openbabel import pybel

STRUCTURE_DIR = "data/structures"
VINA_EXE = "tools/vina_1.2.7_win.exe"
DOCKING_DIR = "data/docking"

def resolve_target(name_or_id: str) -> dict:
    """Resolve a gene/protein name or UniProt ID to accession + sequence."""
    url = "https://rest.uniprot.org/uniprotkb/search"
    params = {"query": name_or_id, "fields": "accession,gene_names,sequence,organism_id", "format": "json", "size": 50}
    resp = requests.get(url, params=params, timeout=10)
    resp.raise_for_status()
    for r in resp.json()["results"]:
        if r["organism"]["taxonId"] == 9606:
            result = r
            break
    else:
        raise ValueError(f"No human (taxon 9606) UniProt entry found for '{name_or_id}'")

    return {
        "accession": result["primaryAccession"],
        "gene_name": result["genes"][0]["geneName"]["value"],
        "sequence": result["sequence"]["value"],
    }


def get_structure(uniprot_accession: str, sequence: str) -> dict:
    """Resolve a 3D structure for a target.

    Prefers an experimental structure cross-referenced in UniProt (via RCSB PDB).
    Falls back to an ESMFold prediction from the raw sequence if none exists.
    """
    os.makedirs(STRUCTURE_DIR, exist_ok=True)

    pdb_id = _find_pdb_id(uniprot_accession)
    if pdb_id:
        path = _download_pdb(pdb_id)
        return {"path": path, "source": "experimental", "pdb_id": pdb_id}

    path = _predict_structure_esmfold(sequence, uniprot_accession)
    return {"path": path, "source": "predicted (ESMFold)", "pdb_id": None}


def _find_pdb_id(uniprot_accession: str) -> str | None:
    """Look for the first PDB cross-reference on the UniProt entry."""
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_accession}.json"
    resp = requests.get(url, timeout=10)
    resp.raise_for_status()
    data = resp.json()
    for xref in data.get("uniProtKBCrossReferences", []):
        if xref["database"] == "PDB":
            return xref["id"]
    return None


def _download_pdb(pdb_id: str) -> str:
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    resp = requests.get(url, timeout=15)
    resp.raise_for_status()
    path = os.path.join(STRUCTURE_DIR, f"{pdb_id}.pdb")
    with open(path, "w") as f:
        f.write(resp.text)
    return path


def _predict_structure_esmfold(sequence: str, label: str) -> str:
    url = "https://api.esmatlas.com/foldSequence/v1/pdb/"
    resp = requests.post(url, data=sequence, timeout=60)
    resp.raise_for_status()
    path = os.path.join(STRUCTURE_DIR, f"{label}_esmfold.pdb")
    with open(path, "w") as f:
        f.write(resp.text)
    return path


def get_known_ligands(uniprot_accession: str, max_results: int = 20) -> list[dict]:
    chembl_target_id = _find_chembl_target(uniprot_accession)
    if not chembl_target_id:
        return []
    raw = _get_chembl_activities(chembl_target_id, max_results)
    return raw


def _find_chembl_target(uniprot_accession: str) -> str | None:
    target_api = new_client.target
    target = target_api.filter(
        target_components__accession=uniprot_accession,
        target_type="SINGLE PROTEIN"
        ).only('target_chembl_id')

    if not target:
        print(f"No ChEMBL target found for UniProt Accession: {uniprot_accession}")
        return None
    
    target_chembl_id = target[0]['target_chembl_id']

    return target_chembl_id
    
        
def _get_chembl_activities(chembl_target_id: str, max_results: int) -> list[dict]:
    print("Fetching activities from server (this may take a moment)...")
    activity_api = new_client.activity
    activities = activity_api.filter(
        target_chembl_id=chembl_target_id
    ).only(
        'molecule_chembl_id', 
        'standard_type', 
        'standard_value', 
        'standard_units'
    )[:max_results]
    
    raw_activities = list(activities)
    print(f"Retrieved {len(raw_activities)} raw activity records.")
    
    valid_activities = []
    molecule_ids = set()
    
    for act in raw_activities:
        m_id = act.get('molecule_chembl_id')
        val = act.get('standard_value')
        
        if m_id and val is not None:
            valid_activities.append(act)
            molecule_ids.add(m_id)
            
    print(f"Found {len(valid_activities)} activities with valid measurements across {len(molecule_ids)} unique molecules.")
    
    if not molecule_ids:
        return []

    molecule_api = new_client.molecule
    molecule_list = list(molecule_ids)
    smiles_lookup = {}
    
    batch_size = 100
    for i in range(0, len(molecule_list), batch_size):
        batch = molecule_list[i:i + batch_size]
        mol_batch = molecule_api.filter(molecule_chembl_id__in=batch).only('molecule_chembl_id', 'molecule_structures')
        for mol in mol_batch:
            mol_id = mol.get('molecule_chembl_id')
            structures = mol.get('molecule_structures')
            if structures and 'canonical_smiles' in structures:
                smiles_lookup[mol_id] = structures['canonical_smiles']

    ligands_list = []
    for act in valid_activities:
        mol_id = act.get('molecule_chembl_id')
        
        if mol_id in smiles_lookup:
            ligand_dict = {
                'chembl_id': mol_id,
                'smiles': smiles_lookup[mol_id],
                'standard_type': act.get('standard_type'),
                'standard_value': act.get('standard_value'),
                'standard_units': act.get('standard_units')
            }
            ligands_list.append(ligand_dict)
            
    return ligands_list


def load_candidate_library(
    exclude_chembl_ids: set[str] | None = None,
    max_candidates: int = 50,
    random_seed: int = 42,
) -> list[dict]:
    """Load a fixed, target-agnostic pool of approved small-molecule drugs.

    Randomly sampled (with a fixed seed for reproducibility) from all ChEMBL
    max_phase=4 molecules, rather than taking the lowest ChEMBL IDs, so the
    pool isn't dominated by whichever molecules happen to be numbered first.
    """
    exclude_chembl_ids = exclude_chembl_ids or set()
    molecule_api = new_client.molecule
    approved = list(molecule_api.filter(
        max_phase=4, # approved for use somewhere
        molecule_structures__isnull=False,
    ).only('molecule_chembl_id', 'pref_name', 'molecule_structures')[:1000])

    random.Random(random_seed).shuffle(approved)  # avoid always drawing the same low-ID entries

    candidates = []
    for mol in approved:
        chembl_id = mol.get('molecule_chembl_id')
        if chembl_id in exclude_chembl_ids:
            continue
        structures = mol.get('molecule_structures')
        smiles = structures.get('canonical_smiles') if structures else None
        if not smiles or not _passes_lipinski(smiles):
            continue
        candidates.append({'chembl_id': chembl_id, 'name': mol.get('pref_name'), 'smiles': smiles})
        if len(candidates) >= max_candidates:
            break

    return candidates


def _passes_lipinski(smiles: str) -> bool:
    """Standard Lipinski rule of five, allowing at most one violation."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False
    violations = sum([
        Descriptors.MolWt(mol) > 500,
        Descriptors.MolLogP(mol) > 5,
        Descriptors.NumHDonors(mol) > 5,
        Descriptors.NumHAcceptors(mol) > 10,
    ])
    return violations <= 1


def dock_candidates(receptor_path: str, candidates: list[dict]) -> list[dict]:
    """Dock candidate molecules against a receptor structure.

    Returns candidates annotated with binding affinity (kcal/mol, more negative
    is better), sorted best-first. Candidates that fail to prepare or dock are
    skipped, with the reason printed rather than silently dropped.
    """
    os.makedirs(DOCKING_DIR, exist_ok=True)
    receptor_pdbqt = _prepare_receptor(receptor_path)
    center, box_size = _define_search_box(receptor_path)

    results = []
    for candidate in candidates:
        ligand_pdbqt = _prepare_ligand(candidate["smiles"], candidate["chembl_id"])
        if ligand_pdbqt is None:
            print(f"Skipping {candidate['chembl_id']}: failed to prepare ligand")
            continue

        docked = _run_vina(receptor_pdbqt, ligand_pdbqt, center, box_size)
        if docked is None:
            print(f"Skipping {candidate['chembl_id']}: docking failed")
            continue
        affinity, pose_path = docked
        results.append({**candidate, "affinity_kcal_mol": affinity, "pose_path": pose_path})

    results.sort(key=lambda r: r["affinity_kcal_mol"])
    return results


def _prepare_receptor(receptor_path: str) -> str:
    """Strip waters/heteroatoms and convert to PDBQT. Cached per receptor file."""
    receptor_name = Path(receptor_path).stem
    pdbqt_path = os.path.join(DOCKING_DIR, f"{receptor_name}_receptor.pdbqt")
    if os.path.exists(pdbqt_path):
        return pdbqt_path

    cleaned_pdb = os.path.join(DOCKING_DIR, f"{receptor_name}_clean.pdb")
    with open(receptor_path) as infile, open(cleaned_pdb, "w") as outfile:
        for line in infile:
            if line.startswith("ATOM"):
                outfile.write(line)

    mol = next(pybel.readfile("pdb", cleaned_pdb))
    mol.addh()
    mol.write("pdbqt", pdbqt_path, overwrite=True, opt={"r": None})  # rigid - no torsion tree
    return pdbqt_path


def _define_search_box(receptor_path: str) -> tuple[list[float], list[float]]:
    """Center on a co-crystallized ligand if present, else fall back to blind docking."""
    hetero_coords, protein_coords = [], []
    with open(receptor_path) as f:
        for line in f:
            if line.startswith("HETATM") and line[17:20].strip() != "HOH":
                hetero_coords.append(_parse_coords(line))
            elif line.startswith("ATOM"):
                protein_coords.append(_parse_coords(line))

    if hetero_coords:
        return _centroid(hetero_coords), [25.0, 25.0, 25.0]  # pocket-sized box

    # Blind docking fallback - no known pocket to anchor on
    span = _bounding_box_span(protein_coords)
    return _centroid(protein_coords), [dim + 10.0 for dim in span]


def _parse_coords(pdb_line: str) -> tuple[float, float, float]:
    return (float(pdb_line[30:38]), float(pdb_line[38:46]), float(pdb_line[46:54]))


def _centroid(coords: list[tuple[float, float, float]]) -> list[float]:
    n = len(coords)
    return [sum(c[i] for c in coords) / n for i in range(3)]


def _bounding_box_span(coords: list[tuple[float, float, float]]) -> list[float]:
    return [max(c[i] for c in coords) - min(c[i] for c in coords) for i in range(3)]


def _prepare_ligand(smiles: str, label: str) -> str | None:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    # Strip counter-ions/salts - a disconnected multi-fragment ligand produces
    # multiple ROOT sections in the PDBQT, which Vina rejects outright.
    if len(Chem.GetMolFrags(mol)) > 1:
        mol = rdMolStandardize.LargestFragmentChooser().choose(mol)

    mol = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol, randomSeed=42) != 0:
        return None
    AllChem.MMFFOptimizeMolecule(mol)

    sdf_path = os.path.join(DOCKING_DIR, f"{label}.sdf")
    with Chem.SDWriter(sdf_path) as writer:
        writer.write(mol)

    pdbqt_path = os.path.join(DOCKING_DIR, f"{label}_ligand.pdbqt")
    ob_mol = next(pybel.readfile("sdf", sdf_path))
    ob_mol.write("pdbqt", pdbqt_path, overwrite=True)
    return pdbqt_path


def _run_vina(receptor_pdbqt: str, ligand_pdbqt: str, center: list[float], box_size: list[float]) -> float | None:
    out_path = ligand_pdbqt.replace("_ligand.pdbqt", "_out.pdbqt")
    cmd = [
        VINA_EXE,
        "--receptor", receptor_pdbqt,
        "--ligand", ligand_pdbqt,
        "--center_x", str(center[0]), "--center_y", str(center[1]), "--center_z", str(center[2]),
        "--size_x", str(box_size[0]), "--size_y", str(box_size[1]), "--size_z", str(box_size[2]),
        "--exhaustiveness", "8",
        "--out", out_path,
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
    if result.returncode != 0:
        print(f"Vina failed: {result.stderr}")
        return None

    for line in result.stdout.splitlines():
        parts = line.split()
        if parts and parts[0] == "1":  # row "1" = best-scoring pose
            return float(parts[1]), out_path
    return None


if __name__ == "__main__":
    target = resolve_target("EGFR")
    structure = get_structure(target["accession"], target["sequence"])
    known = get_known_ligands(target["accession"])
    known_ids = {l["chembl_id"] for l in known}
    candidates = load_candidate_library(exclude_chembl_ids=known_ids, max_candidates=5)  # small for a first test

    results = dock_candidates(structure["path"], candidates)
    for r in results:
        print(f"{r['chembl_id']}: {r['affinity_kcal_mol']} kcal/mol")