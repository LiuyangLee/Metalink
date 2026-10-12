#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Metalink.combine_genomad_vs_vb.R
"""

import argparse
import os
import pandas as pd
import csv

def parse_args():
    parser = argparse.ArgumentParser(
        description="Merge geNomad & VirSorter2 lytic/lysogenic summaries.")
    parser.add_argument("--in_genomad",
                        default="Identification/scf2lytic_or_lyso.summary.txt")
    parser.add_argument("--in_virsorter2",
                        default="Identification/scf2lytic_or_lyso.summary.txt")
    parser.add_argument("--out_3methods",
                        default="Identification/Virus_Identification.tsv")
    parser.add_argument("--out_genomad",
                        default="Identification/list.genomad.tsv")
    parser.add_argument("--out_virsorter2",
                        default="Identification/list.virsorter2.tsv")
    return parser.parse_args()

def read_tsv(path):
    if not os.path.exists(path) or os.stat(path).st_size == 0:
        return pd.DataFrame()
    return pd.read_csv(path, sep='\t', header=0, dtype=str, keep_default_na=False)


def write_tsv(df, path, header=True, quote=None):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    df.to_csv(path, sep='\t', index=False, header=header,
            quoting=quote if quote is not None else csv.QUOTE_NONE)

def process_data(in_genomad, in_virsorter2,
                 out_3methods, out_genomad, out_virsorter2):

    # 1. geNomad ----------------------------------------------------------------
    genomad = read_tsv(in_genomad)
    if genomad.empty:
        genomad = pd.DataFrame(columns=["seqname", "lytic_or_lyso_or_integrated_prophage",
                                        "integrase_presence_or_absence", "ContigName", "method"])
    else:
        genomad.rename(columns={genomad.columns[0]: "seqname"}, inplace=True)
        # separate： geNomad |provirus_
        genomad_split = genomad["seqname"].str.split(r"\|provirus_", n=1, expand=True)
        if genomad_split.shape[1] == 1:
            genomad_split[1] = pd.NA
        genomad[["ContigName", "fragment"]] = genomad_split

        genomad = genomad[["seqname", "lytic_or_lyso_or_integrated_prophage",
                        "integrase_presence_or_absence", "ContigName"]].copy()
        genomad["method"] = "genomad"

    # 2. VirSorter2 -------------------------------------------------------------
    vs2 = read_tsv(in_virsorter2)
    if vs2.empty:
        vs2 = pd.DataFrame(columns=["seqname", "lytic_or_lyso_or_integrated_prophage",
                                    "integrase_presence_or_absence", "ContigName", "method"])
    else:
        vs2.rename(columns={vs2.columns[0]: "seqname"}, inplace=True)
        # separate：VS2 片段原名带 || 后缀；经 Metalink.embed_viral_coords.py 改名的
        # 片段带 |vs2region_S_E。统一先转成 || 前缀再拆，ContigName 均为母 contig
        vs2_norm = vs2["seqname"].str.replace(r"\|vs2region_", "||vs2region_", regex=True)
        vs2_split = vs2_norm.str.split(r"\|\|", n=1, expand=True)
        if vs2_split.shape[1] == 1:
            vs2_split[1] = pd.NA
        vs2[["ContigName", "fragment"]] = vs2_split

        vs2["topology"] = vs2["fragment"].apply(lambda x: "Provirus" if x == "0_partial" else pd.NA)

        vs2 = vs2[["seqname", "lytic_or_lyso_or_integrated_prophage",
                "integrase_presence_or_absence", "ContigName"]].copy()
        vs2["method"] = "virsorter2_vibrant"

    if not genomad.empty:
        vs2 = vs2[~vs2["ContigName"].isin(genomad["ContigName"])]

    # 3. 合并 -------------------------------------------------------------------
    combined = pd.concat([genomad, vs2], ignore_index=True, sort=False)

    # 4. 写文件 -----------------------------------------------------------------
    write_tsv(combined, out_3methods, header=True)
    write_tsv(genomad[["seqname"]], out_genomad, header=False)
    write_tsv(vs2[["seqname"]], out_virsorter2, header=False)

    print("Processing completed successfully!")

def main():
    args = parse_args()
    process_data(args.in_genomad, args.in_virsorter2,
                 args.out_3methods, args.out_genomad, args.out_virsorter2)


if __name__ == "__main__":
    main()
