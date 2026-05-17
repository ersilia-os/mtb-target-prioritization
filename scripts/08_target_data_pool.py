"""
Script 08: Combined compound–target dataset across all three data sources

Merges the per-target compound datasets produced by scripts 05 (ChEMBL),
06 (BindingDB), and 07 (PubChem) into a single pool, preserving the origin
of each record.

Pipeline:
  1. Load all per-target CSVs from each source's datasets/ folder.
  2. Add a `source` column (chembl / bindingdb / pubchem).
  3. Concatenate into one combined table with columns:
       uniprot_ac, gene_name, smiles, affinity_type, affinity_value,
       affinity_unit, source
  4. Compute an InChIKey for every record (NaN if SMILES cannot be parsed).
     The combined table is saved with the inchikey column included.
  5. Save the combined table and a per-target summary with raw per-source
     counts (total) and a deduplicated count by unique InChIKey
     (total_deduplicated).
  6. Save one CSV per target under output/08_target_data_pool/datasets/.

Outputs:
  output/08_target_data_pool/
    08_combined.csv              — all records from all three sources (with inchikey)
    08_target_summary.csv        — per-target counts by source, total, and
                                   total_deduplicated (unique InChIKeys)
    datasets/{uniprot_ac}.csv    — per-target compound lists
"""

import os
import sys

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))

import pandas as pd
from rdkit import Chem
from rdkit.Chem.inchi import MolToInchiKey


def smiles_to_inchikey(smiles):
    """Return the InChIKey for a SMILES string, or None if unparseable."""
    try:
        mol = Chem.MolFromSmiles(str(smiles))
        if mol is None:
            return None
        return MolToInchiKey(mol)
    except Exception:
        return None

output_dir   = os.path.join(root, "..", "output", "08_target_data_pool")
datasets_out = os.path.join(output_dir, "datasets")
os.makedirs(output_dir, exist_ok=True)
os.makedirs(datasets_out, exist_ok=True)

SOURCES = {
    "chembl":    os.path.join(root, "..", "output", "05_chembl_target_data", "datasets"),
    "bindingdb": os.path.join(root, "..", "output", "06_bindingdb"),
    "pubchem":   os.path.join(root, "..", "output", "07_pubchem_target_data", "datasets"),
}

KEEP_COLS = ["uniprot_ac", "gene_name", "smiles", "affinity_type", "affinity_value", "affinity_unit"]

# ── 1. Load all per-target files from each source ─────────────────────────────
all_dfs = []

for source, folder in SOURCES.items():
    if not os.path.isdir(folder):
        print(f"Warning: {source} datasets folder not found at {folder} — skipping.")
        continue
    files = [f for f in os.listdir(folder) if f.endswith(".csv")]
    source_dfs = []
    for fname in files:
        df = pd.read_csv(os.path.join(folder, fname))
        # Keep only the shared columns (BindingDB has an extra monomer_id column)
        cols = [c for c in KEEP_COLS if c in df.columns]
        source_dfs.append(df[cols])
    if source_dfs:
        combined_source = pd.concat(source_dfs, ignore_index=True)
        combined_source["source"] = source
        all_dfs.append(combined_source)
        print(f"  {source:10s}: {len(combined_source):>7,} records from {len(source_dfs)} targets")
    else:
        print(f"  {source:10s}: no files found")

if not all_dfs:
    raise RuntimeError("No data loaded from any source. Check that scripts 05, 06, and 07 have been run.")

# ── 2. Concatenate ─────────────────────────────────────────────────────────────
combined = pd.concat(all_dfs, ignore_index=True)
print(f"\nTotal records: {len(combined):,}")
print(f"Unique targets: {combined['uniprot_ac'].nunique()}")

# ── 3. Compute InChIKeys ───────────────────────────────────────────────────────
print("Computing InChIKeys ...")
combined["inchikey"] = combined["smiles"].apply(smiles_to_inchikey)
n_failed = combined["inchikey"].isna().sum()
if n_failed:
    print(f"  Warning: {n_failed:,} records could not be converted to InChIKey (SMILES unparseable)")

# ── 4. Save combined table ─────────────────────────────────────────────────────
out_combined = os.path.join(output_dir, "08_combined.csv")
combined.to_csv(out_combined, index=False)
print(f"Saved: {out_combined}")

# ── 5. Per-target summary ──────────────────────────────────────────────────────
summary = (
    combined
    .groupby(["uniprot_ac", "gene_name", "source"])
    .size()
    .reset_index(name="n_compounds")
    .pivot_table(
        index=["uniprot_ac", "gene_name"],
        columns="source",
        values="n_compounds",
        fill_value=0,
    )
    .reset_index()
)
summary.columns.name = None

# Ensure all source columns exist even if a source had no data
for src in ["chembl", "bindingdb", "pubchem"]:
    if src not in summary.columns:
        summary[src] = 0

summary["total"] = summary[["chembl", "bindingdb", "pubchem"]].sum(axis=1)

# Deduplicated count: unique InChIKeys per target (excluding NaN)
dedup = (
    combined.dropna(subset=["inchikey"])
    .groupby("uniprot_ac")["inchikey"]
    .nunique()
    .reset_index(name="total_deduplicated")
)
summary = summary.merge(dedup, on="uniprot_ac", how="left")
summary["total_deduplicated"] = summary["total_deduplicated"].fillna(0).astype(int)

summary = summary.sort_values("total", ascending=False).reset_index(drop=True)

out_summary = os.path.join(output_dir, "08_target_summary.csv")
summary.to_csv(out_summary, index=False)
print(f"Saved: {out_summary}")
print(f"\nTop 10 targets by total compound count:")
print(summary.head(10).to_string(index=False))

# ── 6. Save per-target CSVs ────────────────────────────────────────────────────
print(f"\nWriting per-target files ...")
saved = 0
for acc, group in combined.groupby("uniprot_ac"):
    group.to_csv(os.path.join(datasets_out, f"{acc}.csv"), index=False)
    saved += 1
print(f"  Saved {saved} per-target files in {datasets_out}")
