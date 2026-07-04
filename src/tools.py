import requests
import os

STRUCTURE_DIR = "data/structures"

def resolve_target(name_or_id: str) -> dict:
    """Resolve a gene/protein name or UniProt ID to accession + sequence."""
    url = "https://rest.uniprot.org/uniprotkb/search"
    params = {"query": name_or_id, "fields": "accession,gene_names,sequence", "format": "json", "size": 1}
    resp = requests.get(url, params=params, timeout=10)
    resp.raise_for_status()
    result = resp.json()["results"][0]
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


if __name__ == "__main__":
    # Manual smoke test — DRD2 has a known PDB structure, so this should hit the experimental path
    target = resolve_target("DRD2")
    result = get_structure(target["accession"], target["sequence"])
    print(result)