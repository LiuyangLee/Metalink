#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
import re
import os
import glob
import argparse
import pandas as pd

parser = argparse.ArgumentParser(
    description="Merge kraken2/bracken reports into a single taxonomic profile matrix."
)
parser.add_argument("--input", "-i", required=True,
                    help="Folder containing bracken report files.")
parser.add_argument("--output", "-o", required=True,
                    help="Output file path (e.g. /path/kraken2_taxonomic_profiles.tsv).")
parser.add_argument("--suffix", "-s",
                    default=".k2_pluspf_16gb_20220908_bracken_species.report",
                    help="Suffix of bracken report files to be merged.")
args = parser.parse_args()

bracken_reports_folder_path = args.input
result_file = args.output
search_sample_suffix_name = args.suffix

re_num = re.compile(r"[0-9]")
DOMAIN_ORDER = {"U": 0, "R": 1, "B": 2, "A": 3, "E": 4, "V": 5}
LEVEL_ORDER = {"D": 0, "P": 1, "C": 2, "O": 3, "F": 4, "G": 5, "S": 6}

bracken_files_list = glob.glob(
    os.path.join(bracken_reports_folder_path, "*" + search_sample_suffix_name)
)

# Helper to write a valid empty matrix when no data are available
# (e.g. Bracken disabled for all long-read samples).
def _write_empty_matrix(out_path):
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write("Clade_name\tLevel\tDomain\n")
    print("No valid data found; wrote empty matrix to {}".format(out_path))

if not bracken_files_list:
    print("No files matched the given suffix.")
    _write_empty_matrix(result_file)
    sys.exit(0)

long_rows = []
for file_path in bracken_files_list:
    sample_name = os.path.basename(file_path).replace(search_sample_suffix_name, "")
    tax_big_pointer = ""
    with open(file_path, "r", encoding="utf-8") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            abundance, tax_marker, tax_name = parts[1], parts[3], parts[5].lstrip(" ")

            if re_num.search(tax_marker):
                continue

            if (tax_marker == "R") and (tax_name == "root"):
                tax_big_pointer = "R"
            elif (tax_marker == "U") and (tax_name == "unclassified"):
                tax_big_pointer = "U"
            elif tax_marker == "D":
                if tax_name in ("d__Bacteria", "Bacteria"):
                    tax_big_pointer = "B"
                elif tax_name in ("d__Archaea", "Archaea"):
                    tax_big_pointer = "A"
                elif tax_name in ("d__Eukaryota", "Eukaryota"):
                    tax_big_pointer = "E"
                elif tax_name in ("d__Viruses", "Viruses"):
                    tax_big_pointer = "V"


            long_rows.append([tax_name, tax_marker, tax_big_pointer, sample_name, abundance])

if not long_rows:
    _write_empty_matrix(result_file)
    sys.exit(0)

long_df = pd.DataFrame(long_rows, columns=["Clade_name", "Level", "Domain", "sample", "abundance"])
# ===== NEW: merge duplicates by summing abundance =====
dup_mask = long_df.duplicated(subset=["Clade_name", "Level", "Domain", "sample"], keep=False)
if dup_mask.any():
    print(">>> Duplicated entries detected. The following rows will be summed <<<",
          file=sys.stderr)
    to_collapse = (long_df[dup_mask]
                   .sort_values(["Clade_name", "Level", "Domain", "sample"]))
    print(to_collapse, file=sys.stderr)
    long_df["abundance"] = pd.to_numeric(long_df["abundance"], errors="coerce")
    long_df = (long_df
               .groupby(["Clade_name", "Level", "Domain", "sample"], as_index=False)
               .agg({"abundance": "sum"}))
    keys = (to_collapse.drop(columns="abundance")
                       .drop_duplicates()
                       .merge(long_df, on=["Clade_name", "Level", "Domain", "sample"]))
    print(">>> After summing:", file=sys.stderr)
    print(keys, file=sys.stderr)
else:
    print(">>> No duplicates found.", file=sys.stderr)
# ======================================================
df_final = (long_df.pivot(index=["Clade_name", "Level", "Domain"],
                           columns="sample",
                           values="abundance")
                   .fillna(0)
                   .reset_index())

df_final["Domain_rank"] = df_final["Domain"].map(DOMAIN_ORDER)
df_final["Level_rank"] = df_final["Level"].map(LEVEL_ORDER)
df_final = df_final.sort_values(
    by=["Domain_rank", "Level_rank", "Clade_name"], ascending=[True, True, True]
)
df_final = df_final.drop(columns=["Domain_rank", "Level_rank"])

os.makedirs(os.path.dirname(result_file) or ".", exist_ok=True)
df_final.to_csv(result_file, index=False, sep="\t")

print("kraken2_merge.cmd.py complete.")
