"""
Query BindingDB for all molecules tested against the 417 curated MTB targets.

Input:
    data/processed/04_first_selection_417.csv  — prioritised target list (uniprot_ac, gene_name)

Output:
    output/06_bindingdb/{uniprot_ac}.csv       — one CSV per target with activity data
    output/06_bindingdb_combined.csv           — all targets combined

Affinity cutoff: 1,000,000 nM (1000 µM). Only molecules with at least one measurement
at or below this threshold are returned by the BindingDB API.

All affinity values are in nM as returned by BindingDB.
"""

import os
import sys
import time
import json
import urllib.request
import urllib.error

import pandas as pd

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))

input_dir = os.path.join(root, "..", "data", "processed")
output_dir = os.path.join(root, "..", "output", "06_bindingdb")
os.makedirs(input_dir, exist_ok=True)
os.makedirs(output_dir, exist_ok=True)

AFFINITY_CUTOFF_NM = 10000000  # 10000 µM
SLEEP_BETWEEN_REQUESTS = 0.5  # seconds
RETRY_WAIT = 5  # seconds
BASE_URL = "https://bindingdb.org/rest/getLigandsByUniprot"


def query_bindingdb(uniprot_ac):
    """
    Query BindingDB for all ligands with affinity <= AFFINITY_CUTOFF_NM against a target.

    Returns a list of dicts with keys: monomer_id, smiles, affinity_type, affinity_value.
    Returns None on HTTP error (after one retry), and [] for empty / no-data responses.
    """
    url = f"{BASE_URL}?uniprot={uniprot_ac};{AFFINITY_CUTOFF_NM}&response=application/json"
    for attempt in range(2):
        try:
            with urllib.request.urlopen(url, timeout=30) as response:
                raw = response.read().decode("utf-8").strip()
            break
        except urllib.error.HTTPError as e:
            if attempt == 0:
                print(f"    HTTP {e.code} for {uniprot_ac}, retrying in {RETRY_WAIT}s...")
                time.sleep(RETRY_WAIT)
            else:
                print(f"    HTTP {e.code} for {uniprot_ac} on retry — skipping.")
                return None
        except Exception as e:
            if attempt == 0:
                print(f"    Error for {uniprot_ac}: {e}, retrying in {RETRY_WAIT}s...")
                time.sleep(RETRY_WAIT)
            else:
                print(f"    Error for {uniprot_ac} on retry — skipping.")
                return None

    if not raw:
        return []

    try:
        data = json.loads(raw)
    except json.JSONDecodeError:
        print(f"    JSON parse error for {uniprot_ac} — skipping.")
        return None

    # Top-level key has a typo in BindingDB: "getLindsByUniprotResponse"
    response_key = next(iter(data), None)
    if response_key is None:
        return []

    payload = data[response_key]
    affinities = payload.get("bdb.affinities", [])
    if not affinities:
        return []

    records = []
    for entry in affinities:
        monomer_id = entry.get("bdb.monomerid", "")
        smiles = entry.get("bdb.smile", "")
        affinity_type = entry.get("bdb.affinity_type", "")
        affinity_raw = entry.get("bdb.affinity", "")
        try:
            affinity_value = float(str(affinity_raw).strip())
        except (ValueError, TypeError):
            affinity_value = None
        records.append({
            "monomer_id": monomer_id,
            "smiles": smiles,
            "affinity_type": affinity_type,
            "affinity_value": affinity_value,
            "affinity_unit": "nM",
        })
    return records


def main():
    targets = pd.read_csv(os.path.join(input_dir, "04_first_selection_417.csv"))
    targets = targets[["uniprot_ac", "gene_name"]].drop_duplicates()
    print(f"Loaded {len(targets)} targets.")

    all_dfs = []
    no_data = []
    errors = []

    for i, row in targets.iterrows():
        uniprot_ac = row["uniprot_ac"]
        gene_name = row["gene_name"]
        print(f"[{i+1}/{len(targets)}] {uniprot_ac} ({gene_name})", end=" ... ")

        records = query_bindingdb(uniprot_ac)

        if records is None:
            print("ERROR")
            errors.append(uniprot_ac)
        elif len(records) == 0:
            print("no data")
            no_data.append(uniprot_ac)
        else:
            df = pd.DataFrame(records)
            df.insert(0, "gene_name", gene_name)
            df.insert(0, "uniprot_ac", uniprot_ac)
            per_target_path = os.path.join(output_dir, f"{uniprot_ac}.csv")
            df.to_csv(per_target_path, index=False)
            all_dfs.append(df)
            print(f"{len(df)} records")

        time.sleep(SLEEP_BETWEEN_REQUESTS)

    print("\n--- Summary ---")
    print(f"Targets with data:    {len(all_dfs)}")
    print(f"Targets with no data: {len(no_data)}")
    print(f"Targets with errors:  {len(errors)}")

    if all_dfs:
        combined = pd.concat(all_dfs, ignore_index=True)
        combined_path = os.path.join(root, "..", "output", "06_bindingdb_combined.csv")
        combined.to_csv(combined_path, index=False)
        print(f"Total records:        {len(combined)}")
        print(f"Combined saved to:    {combined_path}")

    if no_data:
        print(f"\nTargets with no BindingDB data ({len(no_data)}):")
        for ac in no_data:
            print(f"  {ac}")

    if errors:
        print(f"\nTargets that failed ({len(errors)}):")
        for ac in errors:
            print(f"  {ac}")


if __name__ == "__main__":
    main()
