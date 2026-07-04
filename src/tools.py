import requests

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

if __name__ == "__main__":
    print(resolve_target("DRD2"))