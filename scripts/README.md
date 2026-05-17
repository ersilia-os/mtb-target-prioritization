# Scripts

Run in order. Each script depends on the outputs of the previous one.

## 00_chembl_data_processing.py
Processes the raw ChEMBL bioactivity export: resolves UNCHECKED assay targets via manual curation and maps ChEMBL target IDs to UniProt accessions.
**Manual input required:** run once to generate `chembl_unchecked.csv`, curate it by hand as `chembl_unchecked_manual.csv`, then re-run.

## 01_process_bosc2021.py
Extracts the Bosc et al. 2021 target list, maps ORF IDs to UniProt accessions via the UniProt REST API, and builds the master UniProt table. API calls cover ~3900 genes; sheets already processed are skipped on re-runs.

## 02_data_annotation.py
Cleans and deduplicates the master UniProt table: resolves non-H37Rv accessions to canonical H37Rv entries and standardises column names.

## 03_annotate_proteins.py
Annotates each protein with ChEMBL activity count, PDB structure count, AlphaFold confidence score, and PANTHER functional classification.

## 04_data_merging.py
Merges annotated features, derives binary flags (`is_essential`, `is_growth_defect`, `is_vulnerable`), and produces the ranked selection of 417 prioritised MTB drug targets.
**Cutoff:** `is_vulnerable` requires `high_confidence = True` and `VI < -9.17575`. Review before revision.

## 05_chembl_target_data.py
Pulls quantitative single-protein ChEMBL assays for the 417 targets, auto-fills missing UniProt accessions against the MTB proteome, and builds per-target compound CSVs.

## 06_bindingdb_target_data.py
Queries the BindingDB REST API for all ligands tested against the 417 targets and builds per-target compound CSVs.
**Cutoff:** 100,000 nM. Raise to retrieve weaker binders.

## 07_pubchem_target_data.py
Filters PubChem bioassays from the manually curated table (`data/raw/pubchem_targets.csv`) to those covering the 417 targets, copies raw assay files from the sibling `pubchem-antimicrobial-tasks` repository, and builds per-target compound CSVs. Affinity values are binary (1 = active, 0 = inactive).

## 08_target_data_pool.py
Merges per-target compound datasets from scripts 05, 06, and 07 into a single pool with a `source` column and RDKit-computed InChIKeys. The summary includes raw per-source counts and `total_deduplicated` (unique InChIKeys per target across all sources).
