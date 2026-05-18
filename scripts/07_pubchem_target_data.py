"""
Script 07: PubChem target-level data for the 417 prioritised MTB targets

Requires the pubchem-antimicrobial-tasks repository to be cloned in the same
parent directory as this repository (i.e. as a sibling folder):
  ../pubchem-antimicrobial-tasks/

Pipeline:
  1. Load the manually curated assay table from data/raw/pubchem_targets.csv.
  2. Load the 417 prioritised targets from step 04.
  3. Filter the assay table to retain only rows whose uniprot_accession overlaps
     with the 417-target set.
  4. Copy the corresponding per-assay CSV files (identified by AID) from
     pubchem-antimicrobial-tasks/output/05_curate_bioassays/mycobacterium_tuberculosis/
     into data/processed/07_pubchem_target_data/.
  5. Aggregate compound data per target and save one CSV per target under
     output/07_pubchem_target_data/datasets/ with columns:
       uniprot_ac, gene_name, smiles, affinity_type, affinity_value, affinity_unit
  6. Save the filtered assay table and a per-target summary.

Outputs:
  output/07_pubchem_target_data/
    07_pubchem_assays_417.csv    — filtered assay table for the 417 targets
    07_target_summary.csv        — per-target assay and compound counts
    07_pubchem_combined.csv      — all targets combined (deduplicated rows)
    datasets/{uniprot_ac}.csv    — per-target compound lists (deduplicated rows)
  data/processed/07_pubchem_target_data/
    {aid}.csv                    — raw per-assay dataset files (copied)
"""

import os
import shutil
import sys

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))

import pandas as pd

raw_dir      = os.path.join(root, "..", "data", "raw")
proc_dir     = os.path.join(root, "..", "data", "processed")
output_dir   = os.path.join(root, "..", "output", "07_pubchem_target_data")
datasets_out = os.path.join(output_dir, "datasets")
proc_data    = os.path.join(proc_dir, "07_pubchem_target_data")
os.makedirs(output_dir, exist_ok=True)
os.makedirs(datasets_out, exist_ok=True)
os.makedirs(proc_data, exist_ok=True)

# Path to sibling pubchem-antimicrobial-tasks repository
pubchem_repo  = os.path.join(root, "..", "..", "pubchem-antimicrobial-tasks")
pubchem_assay_dir = os.path.join(
    pubchem_repo, "output", "05_curate_bioassays", "mycobacterium_tuberculosis"
)
if not os.path.isdir(pubchem_assay_dir):
    raise FileNotFoundError(
        f"pubchem-antimicrobial-tasks assay output not found at:\n  {pubchem_assay_dir}\n"
        "Please clone ersilia-os/pubchem-antimicrobial-tasks as a sibling of this repository."
    )

# ── 1. Load curated assay table ───────────────────────────────────────────────
assay_table_path = os.path.join(raw_dir, "pubchem_targets.csv")
assay_table = pd.read_csv(assay_table_path)
print(f"Curated assay table loaded: {len(assay_table):,} assays")

# ── 2. Load 417 prioritised targets ───────────────────────────────────────────
targets_417_path = os.path.join(proc_dir, "04_first_selection_417.csv")
targets_417 = pd.read_csv(targets_417_path)
target_set  = set(targets_417["uniprot_ac"].dropna().unique())
acc_to_gene = dict(zip(targets_417["uniprot_ac"], targets_417["gene_name"]))
print(f"417 prioritised targets loaded ({len(target_set)} unique UniProt accessions)")

# ── 3. Filter assay table to 417-target set ───────────────────────────────────
def _overlaps_417(accession_cell):
    """True if any semicolon-separated accession in the cell is in target_set."""
    if pd.isna(accession_cell) or not str(accession_cell).strip():
        return False
    return bool(set(str(accession_cell).split(";")).intersection(target_set))

in_417 = assay_table["uniprot_accession"].apply(_overlaps_417)
assays_417 = assay_table[in_417].copy()

print(f"\nFiltered assays (uniprot_accession overlaps 417 targets): {len(assays_417):,}")
print(f"  Unique AIDs: {assays_417['aid'].nunique()}")

out_assays = os.path.join(output_dir, "07_pubchem_assays_417.csv")
assays_417.to_csv(out_assays, index=False)
print(f"Saved: {out_assays}")

# ── 4. Copy per-assay raw files to data/processed/07_pubchem_target_data/ ─────
available_files = set(os.listdir(pubchem_assay_dir))
copied = 0
missing_aids = []

for aid in assays_417["aid"].unique():
    fname = f"{aid}.csv"
    if fname in available_files:
        src = os.path.join(pubchem_assay_dir, fname)
        dst = os.path.join(proc_data, fname)
        if not os.path.exists(dst):
            shutil.copy2(src, dst)
            copied += 1
    else:
        missing_aids.append(aid)

print(f"\nCopied {copied} new assay files to {proc_data}")
if missing_aids:
    print(f"  Warning: {len(missing_aids)} AID(s) not found in source directory: {missing_aids}")

# ── 5. Build per-target compound datasets ────────────────────────────────────
print("\nBuilding per-target compound datasets ...")

target_records = {}  # uniprot_ac -> list of dicts

for _, assay_row in assays_417.iterrows():
    aid      = assay_row["aid"]
    fname    = f"{aid}.csv"
    fpath    = os.path.join(proc_data, fname)
    if not os.path.exists(fpath):
        continue

    cpd_df = pd.read_csv(fpath)
    if "smiles" not in cpd_df.columns or "activity" not in cpd_df.columns:
        print(f"  Skipping AID {aid}: unexpected columns {list(cpd_df.columns)}")
        continue

    affinity_type = assay_row["assay_type"]
    affinity_unit = assay_row["unit"]

    accessions = [a.strip() for a in str(assay_row["uniprot_accession"]).split(";")]
    for acc in accessions:
        if acc not in target_set:
            continue
        gene = acc_to_gene.get(acc, "")
        for _, cpd in cpd_df.iterrows():
            target_records.setdefault(acc, []).append({
                "uniprot_ac":     acc,
                "gene_name":      gene,
                "smiles":         cpd["smiles"],
                "affinity_type":  affinity_type,
                "affinity_value": cpd["activity"],
                "affinity_unit":  affinity_unit,
            })

# Save one CSV per target (deduplicated by row)
saved_targets = 0
for acc, records in target_records.items():
    df_target = pd.DataFrame(records).drop_duplicates()
    out_path  = os.path.join(datasets_out, f"{acc}.csv")
    df_target.to_csv(out_path, index=False)
    saved_targets += 1

print(f"  Saved {saved_targets} per-target dataset files in {datasets_out}")

# Save combined CSV across all targets
all_records = [rec for records in target_records.values() for rec in records]
df_combined = pd.DataFrame(all_records).drop_duplicates()
out_combined = os.path.join(output_dir, "07_pubchem_combined.csv")
df_combined.to_csv(out_combined, index=False)
print(f"Saved: {out_combined} ({len(df_combined):,} records)")

# ── 6. Per-target summary ─────────────────────────────────────────────────────
summary_rows = []
for acc, records in target_records.items():
    df_t = pd.DataFrame(records)
    summary_rows.append({
        "uniprot_ac":  acc,
        "gene_name":   acc_to_gene.get(acc, ""),
        "n_assays":    assays_417[assays_417["uniprot_accession"].apply(
                           lambda x: acc in str(x).split(";"))]["aid"].nunique(),
        "n_compounds": len(df_t),
    })

summary = (
    pd.DataFrame(summary_rows)
    .sort_values("n_compounds", ascending=False)
    .reset_index(drop=True)
)
out_summary = os.path.join(output_dir, "07_target_summary.csv")
summary.to_csv(out_summary, index=False)
print(f"Saved: {out_summary}")
print(f"\nTop 10 targets by compound count:")
print(summary.head(10).to_string(index=False))
