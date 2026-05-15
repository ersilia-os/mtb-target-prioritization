# Scripts

Run in order. Each script depends on the outputs of the previous one, except `00_chembl_data_processing.py` and `01_process_bosc2021.py` which can run independently (but `00` must finish before `01` completes, as `01` reads `uniprot_ac_in_chembl.csv`).

## 00_chembl_data_processing.py
Processes the raw ChEMBL bioactivity export for *M. tuberculosis*: resolves UNCHECKED assay targets via manual curation, maps ChEMBL target IDs to UniProt accessions, builds molecule × target activity matrices, and standardises SMILES. Outputs to `data/processed/00_chembl/`.

**Manual input required:** `data/raw/chembl_unchecked_manual.csv` must exist before running. Generate the base file by running the script once (produces `chembl_unchecked.csv`), curating it by hand, then re-running.

## 01_process_bosc2021.py
Extracts all sheets from the Bosc et al. 2021 Excel file into raw CSVs, maps ORF IDs to UniProt accessions via the UniProt REST API, applies manual curation for mismatches, and builds the master UniProt table enriched with ChEMBL and known-target flags. Outputs raw extractions to `data/raw/01_bosc_individual/`, processed per-sheet CSVs to `data/processed/01_bosc2021/`, and the master table to `data/processed/01_master_uniprots.csv`.

**Note:** API calls cover all ~3900 genes — expect several minutes on first run. Sheets already present in `01_bosc2021/` are skipped on re-runs.

## 02_data_annotation.py
Cleans and deduplicates the master UniProt table: resolves non-H37Rv entries back to their canonical H37Rv accessions, applies synonym corrections, removes duplicate ORF IDs, and standardises column names. Outputs `data/processed/02_master_uniprots.csv`.

## 03_annotate_proteins.py
Annotates each protein with ChEMBL activity count, PDB structure count, AlphaFold confidence score, and PANTHER functional classification by calling `src/protein_annotator.py`. Skips execution if the output already exists. Outputs `data/processed/03_annotated_proteins.csv`.

## 04_data_merging.py
Merges the annotated protein features with the master UniProt table, derives binary flags (`is_essential`, `is_growth_defect`, `is_vulnerable`), and produces a ranked selection of 417 prioritised MTB drug targets. Also generates an upset plot and a PANTHER annotation bar chart.
**Cutoff:** `is_vulnerable` requires `high_confidence = True` and `VI < -9.17575`. This threshold was set to match the original analysis and should be reviewed before revision.
Outputs `data/processed/04_first_selection_417.csv` and plots to `output/plots/`.
