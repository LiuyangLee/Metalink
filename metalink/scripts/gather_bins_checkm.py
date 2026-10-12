#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
gather_bins_checkm.py <sample> <Chdb_file> <bins_dir> <genome_dir> <genome_renamed_dir> <csv_out> [min_completeness] [max_contamination] [min_quality_score]
"""

import os
import sys
import shutil
import re

if len(sys.argv) < 7:
    print("Usage: gather_bins_checkm.py <sample> <Chdb_file> <bins_dir> <genome_dir> <genome_renamed_dir> <csv_out> [min_completeness] [max_contamination] [min_quality_score]")
    sys.exit(1)

sample             = sys.argv[1]
Chdb_file          = sys.argv[2]
bins_dir           = sys.argv[3]
genome_dir         = sys.argv[4]   # bin_genome_fna_raw: 只改文件名，不改内容
genome_renamed_dir = sys.argv[5]   # bin_genome_fna_renamed: 改文件名+改header
csv_out            = sys.argv[6]
min_completeness   = float(sys.argv[7]) if len(sys.argv) > 7 else 50
max_contamination  = float(sys.argv[8]) if len(sys.argv) > 8 else 10
min_quality_score  = float(sys.argv[9]) if len(sys.argv) > 9 else 50

if not os.path.isdir(bins_dir):
    print('Bins directory not found:', bins_dir)
    sys.exit(1)

bin_passed = []
with open(Chdb_file) as f:
    next(f)
    for line in f:
        p = line.strip('\n').split('\t')
        com, con = float(p[5]), float(p[6])
        score = com - 5 * con
        if com >= min_completeness and con <= max_contamination and score >= min_quality_score:
            bin_passed.append(p[0])   # e.g. "bin.1"

print('Found {} bins passed Com>={}, Con<={}, QS>={}'.format(len(bin_passed), min_completeness, max_contamination, min_quality_score))

# --- genome_renamed_dir: rename file + rename contig headers ---
os.makedirs(genome_renamed_dir, exist_ok=True)
for bin_id in bin_passed:
    src      = os.path.join(bins_dir, f'{bin_id}.fa')          # source: bin.1.fa
    dst_name = f'{sample}_{bin_id}.fna'                          # target: SRR10325488_top2000k_bin.1.fa
    dst      = os.path.join(genome_renamed_dir, dst_name)
    prefix   = f'{sample}_{bin_id}'
    with open(src) as fin, open(dst, 'w') as fout:
        for line in fin:
            if line.startswith('>'):
                fout.write(re.sub(r'^>', f'>{prefix}---', line, count=1))
            else:
                fout.write(line)

# --- genome_dir: rename file only, keep contig headers untouched ---
os.makedirs(genome_dir, exist_ok=True)
for bin_id in bin_passed:
    src      = os.path.join(bins_dir, f'{bin_id}.fa')          # source: bin.1.fa
    dst_name = f'{sample}_{bin_id}.fa'                          # target: SRR10325488_top2000k_bin.1.fa
    dst      = os.path.join(genome_dir, dst_name)
    shutil.copy2(src, dst)                                      # pure copy, headers unchanged

# --- CSV summary ---
with open(Chdb_file) as f, open(csv_out, 'w') as fo:
    fo.write('genome,completeness,contamination,strain_heterogeneity,length,N50,QS\n')
    next(f)
    for line in f:
        p = line.strip('\n').split('\t')
        com, con = float(p[5]), float(p[6])
        score = com - 5 * con
        if com >= min_completeness and con <= max_contamination and score >= min_quality_score:
            fo.write(f'{sample}_{p[0]}.fa,{com:.2f},{con:.2f},'
                     f'{float(p[7]):.2f},{int(p[8])},{int(p[13])},{score:.2f}\n')

print('CSV summary written to', csv_out)
