"""
Cleans and deduplicates the master UniProt dataset, resolving non-H37Rv entries
back to their canonical H37Rv UniProt accessions and standardising column names.

Inputs:
    data/processed/01_master_uniprots.csv  — master list of all MTB UniProt entries
                                             across Bosc et al. essentiality categories,
                                             ChEMBL, and known targets
                                             (produced by 01_process_bosc2021.py)

Outputs:
    data/processed/02_master_uniprots.csv  — deduplicated H37Rv-only UniProt table
                                             with in_chembl and is_known binary flags
                                             and standardised column names
"""

import os
import sys

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))

import pandas as pd

data_dir = os.path.join(root, "..", "data", "processed")
os.makedirs(data_dir, exist_ok=True)

df = pd.read_csv(os.path.join(data_dir, "01_master_uniprots.csv"))

# Identify non-H37Rv entries and propagate their in_chembl status to the H37Rv synonym
orphan_names = []
for r in df[df["strain"] != "H37Rv"].iterrows():
    orphan_names += [(r[1]["name"], r[1]["uniprot_id"])]

synonyms = {"P95276": "I6YC03", "P96830": "I6WXK4", "O69638": "I6YCQ4"}

uniprot_outs = []
chembl_uniprot = {}
for n in orphan_names:
    if str(n[0]) != "nan":
        df_ = df[df["name"] == n[0]]
        uniprot_in = df_[df_["strain"] == "H37Rv"]["uniprot_id"].tolist()
        uniprot_out = df_[df_["strain"] != "H37Rv"]["uniprot_id"].tolist()
        in_chembl = min(1, sum(df_["in_chembl"].tolist()))
        uniprot_outs += uniprot_out
        for p in uniprot_in:
            chembl_uniprot[p] = in_chembl

df = df[~df["uniprot_id"].isin(uniprot_outs)]

in_chembls = []
for c, p in df[["in_chembl", "uniprot_id"]].values:
    if p in chembl_uniprot:
        in_chembls += [chembl_uniprot[p]]
    else:
        in_chembls += [c]
df.loc[:, "in_chembl"] = in_chembls

uniprots = []
for p in df["uniprot_id"].tolist():
    if p in synonyms:
        uniprots += [synonyms[p]]
    else:
        uniprots += [p]
df.loc[:, "uniprot_id"] = uniprots

remove_orfids = ["Rv0795", "Rv2561"]
df = df[~df["ORFID"].isin(remove_orfids)]

remove_orfids = ["MT0162", "MT3771", "MT1988"]
prot2chembl = {
    "I6YC03": 1,
    "I6YCQ4": 1,
}
df = df[~df["ORFID"].isin(remove_orfids)]

chembls = []
for p, c in df[["uniprot_id", "in_chembl"]].values:
    if p in prot2chembl:
        chembls += [prot2chembl[p]]
    else:
        chembls += [c]
df.loc[:, "in_chembl"] = chembls

rename = {
    "uniprot_id": "uniprot_ac",
    "ORFID": "orf_id",
    "patric": "patric_name",
    "VI": "vi",
    "VI_lower": "vi_lower",
    "VI_higher": "vi_higher",
    "in_known": "is_known",
}

df = df.rename(columns=rename)
df = df.reset_index(drop=True)
df = df.drop(columns=["name"])

df.to_csv(os.path.join(data_dir, "02_master_uniprots.csv"), index=False)
