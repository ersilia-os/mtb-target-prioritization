"""
Annotates MTB proteins with ChEMBL count, PDB count, AlphaFold confidence, PANTHER
classification, and protein evidence flag by calling the protein_annotator module.

Inputs:
    data/processed/02_master_uniprots.csv   — deduplicated H37Rv UniProt table
                                              (produced by 02_data_annotation.py)

Outputs:
    data/processed/03_annotated_proteins.csv — per-protein annotation table
"""

import os
import subprocess

root = os.path.abspath(os.path.dirname(__file__))

input_file = os.path.join(root, "..", "data", "processed", "02_master_uniprots.csv")
output_file = os.path.join(root, "..", "data", "processed", "03_annotated_proteins.csv")

if os.path.exists(output_file):
    print(f"{output_file} already exists, skipping annotation")
else:
    cmd = f"python {root}/../src/protein_annotator.py --input_file {input_file} --output_file {output_file}"
    print(cmd)
    subprocess.Popen(cmd, shell=True).wait()
