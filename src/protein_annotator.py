import os
import requests
from tqdm import tqdm
import pandas as pd
import argparse

TAXID = 83332

def query_uniprot(uniprot_ac):
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_ac}"
    response = requests.get(url)
    if response.status_code != 200:
        return {
            "uniprot_ac": uniprot_ac,
            "gene_name": None,
            "full_name": None,
            "organism": None,
            "protein_length": None,
            "uniprot_reviewed": None,
            "protein_evidence": None,
        }
    entry = response.json()
    uniprot_ac = entry.get("primaryAccession")
    gene_symbol = None
    genes = entry.get("genes", [])
    if genes:
        g = genes[0]
        if "geneName" in g:
            gene_symbol = g["geneName"]["value"]
        elif "orderedLocusNames" in g:
            gene_symbol = g["orderedLocusNames"][0]["value"]
        elif "orfNames" in g:
            gene_symbol = g["orfNames"][0]["value"]
    full_name = (
        entry.get("proteinDescription", {})
             .get("recommendedName", {})
             .get("fullName", {})
             .get("value")
    )
    organism = entry.get("organism", {}).get("scientificName")
    protein_length = entry.get("sequence", {}).get("length")
    is_reviewed = entry.get("entryType", "").startswith("UniProtKB reviewed")
    protein_level = entry.get("proteinExistence", "").startswith("1:")
    return {
        "uniprot_ac": uniprot_ac,
        "gene_name": gene_symbol,
        "full_name": full_name,
        "organism": organism,
        "protein_length": protein_length,
        "uniprot_reviewed": is_reviewed,
        "protein_evidence": protein_level
    }


def query_pdb(uniprot_ac):
    url = f"https://www.ebi.ac.uk/pdbe/graph-api/uniprot/unipdb/{uniprot_ac}"
    r = requests.get(url).json()
    if uniprot_ac not in r:
        return {"pdb_count": 0}
    structures = r[uniprot_ac]["data"]
    num_structures = len(structures)
    return {
        "pdb_count": num_structures,
    }


def query_alphafold(uniprot_ac):
    url = f"https://alphafold.ebi.ac.uk/api/uniprot/summary/{uniprot_ac}.json"
    r = requests.get(url)
    if r.status_code != 200:
        return {
            "alphafold_conf": None,
        }
    data = r.json()
    scores = []
    for s in data.get("structures", []):
        summary = s.get("summary", {})
        score = summary.get("confidence_avg_local_score")
        if isinstance(score, (int, float)):
            scores.append(score)
    max_conf = max(scores) if scores else None
    return {
        "alphafold_conf": max_conf,
    }


def query_chembl(uniprot_ac):
    BASE = "https://www.ebi.ac.uk"
    target_url = f"{BASE}/chembl/api/data/target.json?target_components__accession={uniprot_ac}"
    r = requests.get(target_url)
    r.raise_for_status()
    data = r.json()
    if data["page_meta"]["total_count"] == 0:
        return {"chembl_count": 0}
    target_ids = [t["target_chembl_id"] for t in data["targets"]]
    all_molecules = set()
    for tid in target_ids:
        activities_url = f"{BASE}/chembl/api/data/activity.json?target_chembl_id={tid}&limit=1000"
        while activities_url:
            r = requests.get(activities_url)
            r.raise_for_status()
            acts = r.json()
            for act in acts["activities"]:
                mol_id = act.get("molecule_chembl_id")
                if mol_id:
                    all_molecules.add(mol_id)
            next_url = acts["page_meta"].get("next")
            if next_url:
                if next_url.startswith("/"):
                    activities_url = BASE + next_url
                else:
                    activities_url = next_url
            else:
                activities_url = None
    return {"chembl_count": len(all_molecules)}


def query_panther(uniprot_ac):
    url = f"https://pantherdb.org/services/oai/pantherdb/geneinfo?geneInputList={uniprot_ac}&organism={TAXID}"
    r = requests.get(url)
    data = r.json()
    if "mapped_genes" not in data["search"]:
        return {
            "panther_family_name": None,
            "panther_sf_name": None,
            "panther_annotation": None,
        }
    gene = data["search"]["mapped_genes"]["gene"]
    family_name = gene.get("family_name", None)
    sf_name = gene.get("sf_name", None)
    annotations = []
    if "annotation_type_list" not in gene:
        return {
            "panther_family_name": family_name,
            "panther_sf_name": sf_name,
            "panther_annotation": None,
        }
    for entry in gene["annotation_type_list"]["annotation_data_type"]:
        if type(entry) is not dict:
            continue
        if entry.get("content") == "ANNOT_TYPE_ID_PANTHER_PC":
            ann = entry["annotation_list"]["annotation"]["name"]
            annotations.append(ann)
    return {
        "panther_family_name": family_name,
        "panther_sf_name": sf_name,
        "panther_annotation": "; ".join(annotations),
    }


def annotate_proteins(uniprot_acs, output_file):
    columns=[
        "uniprot_ac",
        "gene_name",
        "full_name",
        "organism",
        "protein_length",
        "uniprot_reviewed",
        "protein_evidence",
        "pdb_count",
        "alphafold_conf",
        "chembl_count",
        "panther_family_name",
        "panther_sf_name",
        "panther_annotation",
    ]
    if not os.path.exists(output_file):
        df = pd.DataFrame(columns=columns)
        df.to_csv(output_file, index=False)
    for uniprot_ac in tqdm(uniprot_acs, desc="Annotating proteins"):
        df = pd.read_csv(output_file)
        if uniprot_ac in df["uniprot_ac"].values:
            continue
        uniprot_data = query_uniprot(uniprot_ac)
        pdb_data = query_pdb(uniprot_ac)
        alphafold_data = query_alphafold(uniprot_ac)
        chembl_data = query_chembl(uniprot_ac)
        panther_data = query_panther(uniprot_ac)
        data = [
            uniprot_data["uniprot_ac"],
            uniprot_data["gene_name"],
            uniprot_data["full_name"],
            uniprot_data["organism"],
            uniprot_data["protein_length"],
            uniprot_data["uniprot_reviewed"],
            uniprot_data["protein_evidence"],
            pdb_data["pdb_count"],
            alphafold_data["alphafold_conf"],
            chembl_data["chembl_count"],
            panther_data["panther_family_name"],
            panther_data["panther_sf_name"],
            panther_data["panther_annotation"],
        ]
        df_ = pd.DataFrame([data], columns=columns)
        df = pd.concat([df, df_], ignore_index=True)
        df.to_csv(output_file, index=False)


def main():
    parser = argparse.ArgumentParser(description="Annotate proteins with data from various databases.")
    parser.add_argument("--input_file", "-i", type=str, required=True, help="Input CSV file with a column 'uniprot_ac' containing UniProt accession numbers.")
    parser.add_argument("--output_file", "-o", type=str, required=True, help="Output CSV file to save the annotations.")
    parsed_args = parser.parse_args()
    input_file = parsed_args.input_file
    output_file = parsed_args.output_file
    input_df = pd.read_csv(input_file)
    uniprot_acs = input_df["uniprot_ac"].tolist()
    annotate_proteins(uniprot_acs, output_file)


if __name__ == "__main__":
    main()