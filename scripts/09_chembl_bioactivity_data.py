"""
Script 09: Compile MTB bioactivity compound list from ChEMBL and PubChem

Sources:
  - ChEMBL organism-level assays: data/raw/chembl_antimicrobial_tasks/ORG_DS.csv (dose-response)
    and ORG_SP.csv (single-point)
  - PubChem target-level assays: selected AIDs from the pubchem-antimicrobial-tasks sibling
    repository (output/05_curate_bioassays/mycobacterium_tuberculosis/).
    Two AID lists are used:
      PUBCHEM_MORE — assays that also appear in pubchem_targets.csv (target-annotated)
      PUBCHEM_ONLY — additional organism-level assays not in pubchem_targets.csv

The script merges all compounds, deduplicates by InChIKey (keeping the first SMILES
encountered), tags each compound with its source(s), and writes the result to
output/09_mtb_bioactivity.csv.

Requires: pubchem-antimicrobial-tasks cloned as a sibling of this repository.
"""

import os
import sys

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))

import pandas as pd

raw_dir    = os.path.join(root, "..", "data", "raw")
output_dir = os.path.join(root, "..", "output")
os.makedirs(output_dir, exist_ok=True)

pubchem_assay_dir = os.path.join(
    root, "..", "..", "pubchem-antimicrobial-tasks",
    "output", "05_curate_bioassays", "mycobacterium_tuberculosis"
)

# PubChem AIDs to include (both lists combined)
PUBCHEM_MORE = [488890, 492952, 1626, 449762, 449764]
PUBCHEM_ONLY = [1332, 1949, 2842, 504645, 504646, 743175, 1259417]
PUBCHEM_AIDS = PUBCHEM_MORE + PUBCHEM_ONLY

# ── 1. Load ChEMBL organism-level data ───────────────────────────────────────

chembl_ds = pd.read_csv(os.path.join(raw_dir, "chembl_antimicrobial_tasks", "ORG_DS.csv"))
chembl_sp = pd.read_csv(os.path.join(raw_dir, "chembl_antimicrobial_tasks", "ORG_SP.csv"))

chembl_ds["source"] = "chembl_org_ds"
chembl_sp["source"] = "chembl_org_sp"

chembl = pd.concat([chembl_ds, chembl_sp], ignore_index=True)[["smiles", "inchikey", "source"]]
print(f"ChEMBL ORG_DS: {len(chembl_ds):,} rows")
print(f"ChEMBL ORG_SP: {len(chembl_sp):,} rows")

# ── 2. Load PubChem assay data ────────────────────────────────────────────────

if not os.path.isdir(pubchem_assay_dir):
    raise FileNotFoundError(
        f"pubchem-antimicrobial-tasks assay output not found at:\n  {pubchem_assay_dir}\n"
        "Please clone ersilia-os/pubchem-antimicrobial-tasks as a sibling of this repository."
    )

pubchem_frames = []
missing_aids = []

for aid in PUBCHEM_AIDS:
    fpath = os.path.join(pubchem_assay_dir, f"{aid}.csv")
    if not os.path.exists(fpath):
        missing_aids.append(aid)
        continue
    df = pd.read_csv(fpath)
    if "smiles" not in df.columns or "inchikey" not in df.columns:
        print(f"  Warning: AID {aid} missing expected columns, skipping")
        continue
    df = df[["smiles", "inchikey"]].copy()
    df["source"] = f"pubchem_{aid}"
    pubchem_frames.append(df)
    print(f"  AID {aid}: {len(df):,} rows")

if missing_aids:
    print(f"\nWarning: {len(missing_aids)} AID(s) not found: {missing_aids}")

pubchem = pd.concat(pubchem_frames, ignore_index=True) if pubchem_frames else pd.DataFrame(columns=["smiles", "inchikey", "source"])
print(f"\nPubChem total: {len(pubchem):,} rows across {len(pubchem_frames)} assays")

# ── 3. Merge and deduplicate by InChIKey ──────────────────────────────────────

all_cpds = pd.concat([chembl, pubchem], ignore_index=True)

# Keep the first SMILES per InChIKey
dedup = all_cpds.drop_duplicates(subset="inchikey").reset_index(drop=True)[["inchikey", "smiles"]]

print(f"\nTotal compounds before deduplication: {len(all_cpds):,}")
print(f"Unique compounds (by InChIKey):        {len(dedup):,}")

# ── 4. Save output ────────────────────────────────────────────────────────────

out_path = os.path.join(output_dir, "09_mtb_bioactivity.csv")
dedup.to_csv(out_path, index=False)
print(f"\nSaved: {out_path}")
