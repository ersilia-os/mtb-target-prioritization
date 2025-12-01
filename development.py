import requests

TAXID = 83332

def query_panther(uniprot_ac):
    url = f"https://pantherdb.org/services/oai/pantherdb/geneinfo?geneInputList={uniprot_ac}&organism={TAXID}"
    r = requests.get(url)
    data = r.json()

    gene = data["search"]["mapped_genes"]["gene"]

    family_name = gene.get("family_name")
    sf_name = gene.get("sf_name")

    annotations = []
    for entry in gene["annotation_type_list"]["annotation_data_type"]:
        if entry.get("content") == "ANNOT_TYPE_ID_PANTHER_PC":
            ann = entry["annotation_list"]["annotation"]["name"]
            annotations.append(ann)

    return {
        "panther_family_name": family_name,
        "panther_sf_name": sf_name,
        "panther_annotation": "; ".join(annotations),
    }


uniprot_ac = "P96262"
panther_data = query_panther(uniprot_ac)
print(panther_data)