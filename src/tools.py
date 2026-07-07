import requests
import os
from chembl_webresource_client.new_client import new_client

STRUCTURE_DIR = "data/structures"

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

def get_known_ligands(uniprot_accession: str, max_results: int = 50) -> list[dict]:
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


if __name__ == "__main__":
    target = resolve_target("EGFR")
    ligands = get_known_ligands(target["accession"])
    print(f"Found {len(ligands)} known ligands")
    print(ligands[0] if ligands else "none")