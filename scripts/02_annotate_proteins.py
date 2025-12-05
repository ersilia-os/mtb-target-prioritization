import os
import subprocess

root = os.path.abspath(os.path.dirname(__file__))

input_file = os.path.join(root, "..", "data", "processed", "01_master_uniprots.csv")
output_file = os.path.join(root, "..", "data", "processed", "02_annotated_proteins.csv")

cmd = f"python {root}/../src/protein_annotator.py --input_file {input_file} --output_file {output_file}"
print(cmd)

subprocess.Popen(cmd, shell=True).wait()
