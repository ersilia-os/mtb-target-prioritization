"""
Script 05: ChEMBL target-level data for the 417 prioritised MTB targets

Requires the chembl-antimicrobial-tasks repository to be cloned in the same
parent directory as this repository (i.e. as a sibling folder):
  ../chembl-antimicrobial-tasks/

Pipeline:
  1. Load the assay master table and auto-fill missing UniProt accessions by
     matching target_name_curated against the MTB UniProtKB proteome.
  2. Save the enriched master table as 05_chembl_master_assay_table.csv.
  3. Load the 417 prioritised targets from step 04.
  4. From the filled master table, retain SINGLE PROTEIN, quantitative assays
     whose UniProt accession overlaps with the 417-target set and whose source
     is not BINDINGDB (BindingDB data is handled separately).
  5. Save the filtered assay table and a per-target summary.
  6. Copy the corresponding per-assay dataset files from chembl-antimicrobial-tasks
     into data/processed/05_chembl_target_data/.
  7. Aggregate compound data per target and save one CSV per target under
     output/05_chembl_target_data/datasets/ with columns:
       uniprot_ac, gene_name, smiles, affinity_type, affinity_value, affinity_unit

Outputs:
  output/05_chembl_target_data/
    05_chembl_master_assay_table.csv  — master table with UniProt gaps filled
    05_single_protein_assays_417.csv  — filtered assays for the 417 targets
    05_target_summary.csv             — per-target assay and compound counts
    datasets/{uniprot_ac}.csv         — per-target compound lists
  data/processed/05_chembl_target_data/
    {assay_id}_*.csv.gz               — raw per-assay dataset files (copied)
"""

import os
import re
import shutil
import sys

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))

import pandas as pd

raw_dir      = os.path.join(root, "..", "data", "raw")
proc_dir     = os.path.join(root, "..", "data", "processed")
output_dir   = os.path.join(root, "..", "output", "05_chembl_target_data")
datasets_out = os.path.join(output_dir, "datasets")
proc_data    = os.path.join(proc_dir, "05_chembl_target_data")
os.makedirs(output_dir, exist_ok=True)
os.makedirs(datasets_out, exist_ok=True)
os.makedirs(proc_data, exist_ok=True)

# Path to sibling chembl-antimicrobial-tasks repository
chembl_repo  = os.path.join(root, "..", "..", "chembl-antimicrobial-tasks")
datasets_qt  = os.path.join(chembl_repo, "output", "mtuberculosis", "12_datasets", "datasets_qt")
if not os.path.isdir(datasets_qt):
    raise FileNotFoundError(
        f"chembl-antimicrobial-tasks datasets not found at:\n  {datasets_qt}\n"
        "Please clone ersilia-os/chembl-antimicrobial-tasks as a sibling of this repository."
    )

# ── 1. Load master table ───────────────────────────────────────────────────────
master_path = os.path.join(raw_dir, "chembl_antimicrobial_tasks", "18_assays_master.csv")
master = pd.read_csv(master_path)
print(f"Master table loaded: {len(master):,} assays")

# ── 2. Auto-fill missing UniProt accessions ────────────────────────────────────

# ── 2a. UniProtKB proteome index ──────────────────────────────────────────────
proteome_files = sorted(
    f for f in os.listdir(raw_dir)
    if f.startswith("uniprotkb_proteome") and f.endswith(".tsv")
)
if not proteome_files:
    raise FileNotFoundError("No uniprotkb_proteome_*.tsv found in data/raw/")
proteome = pd.read_csv(os.path.join(raw_dir, proteome_files[-1]), sep="\t")

def _norm_search(s):
    """Lowercase, remove EC numbers, collapse punctuation to spaces."""
    s = str(s).lower()
    s = re.sub(r"\(ec\s[\d\.\-]+\)", "", s)
    s = re.sub(r"[^\w\s]", " ", s)
    return re.sub(r"\s+", " ", s).strip()

proteome["_prot_norm"]  = proteome["Protein names"].fillna("").apply(_norm_search)
proteome["_gene_norm"]  = proteome["Gene Names"].fillna("").apply(_norm_search)
proteome["_gene_tokens"] = proteome["_gene_norm"].apply(lambda g: set(g.split()))

RV_PAT = re.compile(r"rv\d+[a-z]?", re.I)

# Curated abbreviation → UniProt gene-name mappings.
# Keys are post-_norm_search strings.
ABBREV_TO_GENE = {
    # Reductases / oxidoreductases
    "dhfr":                              "dfra",
    "fabi":                              "inha",
    "aladh":                             "ald",
    "l aladh":                           "ald",
    "asadh":                             "asd",
    "nadh menaquinone oxidoreductase":   "ndh",
    "ndh2":                              "ndh",
    "ndh 2":                             "ndh",
    # Fatty acid / lipid metabolism
    "faal19":                            "fadd19",
    "facl28":                            "fadd28",
    "fatty acid synthase i":             "fas",
    "fatty acid synthase 1":             "fas",
    "fatty acid synthase type i":        "fas",
    "malate synthase":                   "glcb",
    # Aminoglycoside modification
    "aac 2 1c":                          "aac",
    "eis acetyltransferase":             "eis",
    # Antigen 85C
    "antigen 85c acyltransferase":       "fbpc",
    "antigen 85c mycolyltransferase":    "fbpc",
    # Protein kinases
    "protein kinase b":                  "pknb",
    "protein kinase g":                  "pkng",
    # Phosphatases
    "ptp b":                             "ptpb",
    # CYP isoforms (A1 suffix absent from UniProt gene field)
    "cyp121a1":                          "cyp121",
    "cyp124a1":                          "cyp124",
    # DNA metabolism
    "nad dependent dna ligase":          "liga",
    "dna liga":                          "liga",
    "dhq2":                              "aroq",
    "dna sliding clamp protein":         "dnan",
    # Chaperone / proteolysis
    "hsp70":                             "dnak",
    # Other enzymes
    "icl2":                              "aceaa",
    "hgprt":                             "hpt",
    "metrs":                             "mets",
    "shikimate dehydrogenase":           "aroe",
    "zinc metalloprotease 1":            "zmp1",
    "beta ketoacyl acyl carrier protein synthase mtfabh condensing enzyme": "fabh",
    "carbonic anhydrase 3":              "rv3273",
    "beta carbonic anhydrase rv3588":    "rv3588c",
    "epoxide hydrolase e":               "ephe",
    "glutamine dependent nade":          "nade",
    # Regulators
    "ethr2":                             "rv0078",
    # RNA polymerase (canonical RpoB mutants)
    "rna polymerase rpob":               "rpob",
    "rna polymerase d435v mutant":       "rpob",
    "rna polymerase h445y mutant":       "rpob",
    # Partial-length construct
    "partial length pknb":               "pknb",
}

MUTANT_PAT = re.compile(
    r"\s+([a-z]\d+[a-z](\s+mutant)?|"
    r"c171q\s+mutant|d435v\s+mutant|h445y\s+mutant|t313a|"
    r"open.gate\s+mutant|n.terminal\s+domain|"
    r"partial\s+length.*|m1\s+to\s+g279\s+residues.*)",
    re.I,
)


def _match_uniprot(norm_name, variant_names_str):
    """Return the UniProt Entry accession for a target, or '' if none found.

    Matching levels (highest to lowest confidence):
      4 — gene-token exact subset match
      3 — Rv locus tag match
      2 — protein name substring match
    Ambiguous multi-entry ties at the same level → returns ''.
    """
    variants = [norm_name] + [v.strip().lower() for v in variant_names_str.split("|")]
    variants = list(dict.fromkeys(variants))
    candidates = []

    for raw in variants:
        raw_norm = _norm_search(raw)

        # Check ABBREV_TO_GENE on full (un-stripped) name first
        pre_alias = ABBREV_TO_GENE.get(raw_norm)
        if pre_alias:
            for _, u in proteome.iterrows():
                if set(pre_alias.split()).issubset(u["_gene_tokens"]):
                    candidates.append((4, u["Entry"]))

        stripped      = MUTANT_PAT.sub("", raw).strip()
        stripped_norm = _norm_search(stripped)
        tokens        = set(stripped_norm.split()) - {""}

        alias        = ABBREV_TO_GENE.get(stripped_norm, stripped_norm)
        alias_tokens = set(alias.split()) - {""}
        for _, u in proteome.iterrows():
            if alias_tokens and alias_tokens.issubset(u["_gene_tokens"]):
                candidates.append((4, u["Entry"]))

        for tag in RV_PAT.findall(raw):
            tag_norm = _norm_search(tag)
            for _, u in proteome.iterrows():
                if tag_norm in u["_gene_tokens"]:
                    candidates.append((3, u["Entry"]))

        if len(stripped_norm) >= 6 and tokens:
            for _, u in proteome.iterrows():
                if stripped_norm in u["_prot_norm"]:
                    candidates.append((2, u["Entry"]))

    if not candidates:
        return ""

    max_score   = max(c[0] for c in candidates)
    top_entries = list(dict.fromkeys(c[1] for c in candidates if c[0] == max_score))
    return "" if len(top_entries) > 1 else top_entries[0]


def _normalise_name(name):
    """Normalise target name for grouping (strip organism prefix & parentheticals)."""
    if pd.isna(name):
        return ""
    n = str(name).lower().strip()
    for prefix in ["mycobacterium tuberculosis ", "mtb", "mt"]:
        if n.startswith(prefix):
            n = n[len(prefix):].strip()
    return re.sub(r"\s*\(.*?\)\s*$", "", n).strip()


# Identify SINGLE PROTEIN rows with missing UniProt in the master table
sp_mask    = master["target_type_curated_extra"] == "SINGLE PROTEIN"
null_mask  = master["uniprot_accession"].isna() | (master["uniprot_accession"] == "")
to_fill    = master[sp_mask & null_mask].copy()

# Group by normalised name → build variant_names string for _match_uniprot
to_fill["_norm"] = to_fill["target_name_curated"].apply(_normalise_name)
groups = (
    to_fill.groupby("_norm")["target_name_curated"]
    .apply(lambda x: " | ".join(sorted(set(x))))
    .reset_index()
    .rename(columns={"target_name_curated": "variant_names"})
)

print(f"\nAuto-filling UniProt for {len(groups)} normalised target groups "
      f"({null_mask.sum():,} rows without accession)...")

groups["auto_uniprot"] = groups.apply(
    lambda r: _match_uniprot(r["_norm"], r["variant_names"]), axis=1
)
auto_map = dict(zip(groups["_norm"], groups["auto_uniprot"]))

n_auto = sum(1 for v in auto_map.values() if v)
print(f"  Auto-filled: {n_auto} / {len(auto_map)} groups")

fill_map = auto_map

# ── 3. Build filled master table ───────────────────────────────────────────────
filled_master = master.copy()
filled_master["_norm"] = filled_master["target_name_curated"].apply(_normalise_name)

fill_mask = sp_mask & null_mask
filled_master.loc[fill_mask, "uniprot_accession"] = (
    filled_master.loc[fill_mask, "_norm"]
    .map(fill_map)
    .replace("", pd.NA)
)
filled_master = filled_master.drop(columns=["_norm"])

n_filled = (
    filled_master.loc[sp_mask & null_mask, "uniprot_accession"]
    .notna()
    .sum()
)
print(f"\nUniProt cells filled in master table: {n_filled:,} "
      f"(out of {null_mask.sum():,} SINGLE PROTEIN rows that were empty)")

out_master = os.path.join(output_dir, "05_chembl_master_assay_table.csv")
filled_master.to_csv(out_master, index=False)
print(f"Saved: {out_master}")

# ── 4. Load 417 prioritised targets ───────────────────────────────────────────
targets_417_path = os.path.join(proc_dir, "04_first_selection_417.csv")
targets_417 = pd.read_csv(targets_417_path)
target_set  = set(targets_417["uniprot_ac"].dropna().unique())
print(f"\n417 prioritised targets loaded ({len(target_set)} unique UniProt accessions)")

# ── 5. Filter to 417-target, non-BindingDB, SINGLE PROTEIN assays ─────────────
def _overlaps_417(accession_cell):
    """True if any semicolon-separated accession in the cell is in target_set."""
    if pd.isna(accession_cell) or not str(accession_cell).strip():
        return False
    return bool(set(str(accession_cell).split(";")).intersection(target_set))

sp_rows        = filled_master["target_type_curated_extra"] == "SINGLE PROTEIN"
not_bindingdb  = filled_master["source_label"] != "BINDINGDB"
in_417         = filled_master["uniprot_accession"].apply(_overlaps_417)
quantitative   = filled_master["dataset_type"] == "quantitative"

assays_417 = filled_master[sp_rows & not_bindingdb & in_417 & quantitative].copy()

print(f"\nFiltered assays (SINGLE PROTEIN, quantitative, in 417 targets, not BindingDB): "
      f"{len(assays_417):,}")
print(f"  Unique UniProt accessions covered: "
      f"{assays_417['uniprot_accession'].nunique()}")
print(f"  Unique assay IDs: {assays_417['assay_id'].nunique()}")

out_assays = os.path.join(output_dir, "05_single_protein_assays_417.csv")
assays_417.to_csv(out_assays, index=False)
print(f"Saved: {out_assays}")

# ── 6. Per-target summary ──────────────────────────────────────────────────────
summary = (
    assays_417.groupby(["uniprot_accession", "target_name_curated"])
    .agg(
        n_assays=("assay_id", "count"),
        n_compounds=("cpds", "sum"),
    )
    .reset_index()
    .sort_values("n_assays", ascending=False)
)

out_summary = os.path.join(output_dir, "05_target_summary.csv")
summary.to_csv(out_summary, index=False)
print(f"Saved: {out_summary}")
print(f"\nTop 10 targets by assay count:")
print(summary.head(10).to_string(index=False))

# ── 7. Resolve per-assay dataset file for each selected assay ─────────────────
# File naming in chembl-antimicrobial-tasks:
#   {assay_id}_{activity_type}_{unit}_qt_{cutoff}.csv.gz
# Use selected_cutoff when available, else fall back to the first evaluated_cutoff.

available_files = set(os.listdir(datasets_qt))

def _resolve_dataset_file(row):
    """Return (filename, cutoff_used) for an assay row, or (None, None)."""
    cutoffs = []
    if pd.notna(row["selected_cutoff"]):
        cutoffs = [row["selected_cutoff"]]
    elif pd.notna(row.get("evaluated_cutoffs")):
        cutoffs = [c.strip() for c in str(row["evaluated_cutoffs"]).split(";")]
    for c in cutoffs:
        fname = f"{row['assay_id']}_{row['activity_type']}_{row['unit']}_qt_{c}.csv.gz"
        if fname in available_files:
            return fname, c
    return None, None

assays_417[["_dataset_file", "_cutoff_used"]] = assays_417.apply(
    lambda r: pd.Series(_resolve_dataset_file(r)), axis=1
)

matched_files = assays_417["_dataset_file"].notna()
print(f"\nDataset files resolved: {matched_files.sum():,} / {len(assays_417):,}")

# ── 8. Copy resolved dataset files to data/processed/05_chembl_target_data/ ───
print(f"Copying dataset files to {proc_data} ...")
copied = 0
for fname in assays_417.loc[matched_files, "_dataset_file"].unique():
    src = os.path.join(datasets_qt, fname)
    dst = os.path.join(proc_data, fname)
    if not os.path.exists(dst):
        shutil.copy2(src, dst)
        copied += 1
print(f"  Copied {copied} new files ({assays_417['_dataset_file'].nunique()} total in folder)")

# ── 9. Build per-target compound lists ────────────────────────────────────────
# For multi-accession rows (e.g. "P9WPC5;P9WPC3"), assign to each individual
# accession that appears in the 417-target set.
acc_to_gene = dict(zip(targets_417["uniprot_ac"], targets_417["gene_name"]))

print(f"\nBuilding per-target compound datasets ...")

# Collect all compound data keyed by individual uniprot_ac
target_records = {}  # uniprot_ac -> list of dicts

for _, assay_row in assays_417[matched_files].iterrows():
    fname   = assay_row["_dataset_file"]
    cpd_df  = pd.read_csv(os.path.join(proc_data, fname))

    # Map to individual accessions in the 417-target set
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
                "affinity_type":  cpd["activity_type"],
                "affinity_value": cpd["value"] * 1000 if cpd["unit"] == "umol.L-1" else cpd["value"],
                "affinity_unit":  "nM" if cpd["unit"] == "umol.L-1" else cpd["unit"],
            })

# Save one CSV per target
saved_targets = 0
for acc, records in target_records.items():
    df_target = pd.DataFrame(records)
    gene      = acc_to_gene.get(acc, acc)
    out_path  = os.path.join(datasets_out, f"{acc}.csv")
    df_target.to_csv(out_path, index=False)
    saved_targets += 1

print(f"  Saved {saved_targets} per-target dataset files in {datasets_out}")
