#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
gather_bins_checkm.no_add_name.py <sample> <Chdb_file> <bins_dir> <genome_dir> <genome_renamed_dir> <csv_out> [min_completeness] [max_contamination] [min_quality_score]
python gather_bins_checkm.no_add_name.py L_33_42 L_33_42.Chdb.tsv bins L_33_42_genomes L_33_42_renamed_genomes L_33_42_renamed_sample_genomes L_33_42_summary.csv
"""

import os
import sys
import shutil
import re

if len(sys.argv) < 7:
    print("Usage: gather_bins_checkm.no_add_name.py <sample> <Chdb_file> <bins_dir> <genome_dir> <genome_renamed_dir> <csv_out> [min_completeness] [max_contamination] [min_quality_score]")
    sys.exit(1)

sample                      = sys.argv[1]
Chdb_file                   = sys.argv[2]
bins_dir                    = sys.argv[3]
genome_dir                  = sys.argv[4]
genome_renamed_dir          = sys.argv[5]
csv_out                     = sys.argv[6]
min_completeness            = float(sys.argv[7]) if len(sys.argv) > 7 else 50
max_contamination           = float(sys.argv[8]) if len(sys.argv) > 8 else 10
min_quality_score           = float(sys.argv[9]) if len(sys.argv) > 9 else 50

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
            fa_name = f'{p[0]}.fa'
            bin_passed.append(fa_name)

print('Found {} bins passed Com>={}, Con<={}, QS>={}'.format(len(bin_passed), min_completeness, max_contamination, min_quality_score))

# --- original renamed with bin prefix (UNCHANGED) ---
os.makedirs(genome_renamed_dir, exist_ok=True)
for fa in bin_passed:
    src = os.path.join(bins_dir, fa)
    renamed = os.path.splitext(fa)[0] + '.fna'
    dst = os.path.join(genome_renamed_dir, renamed)
    prefix = os.path.splitext(fa)[0]
    with open(src) as fin, open(dst, 'w') as fout:
        for line in fin:
            if line.startswith('>'):
                fout.write(re.sub(r'^>', f'>{prefix}---', line, count=1))
            else:
                fout.write(line)

os.makedirs(genome_dir, exist_ok=True)
for fa in bin_passed:
    shutil.copy2(os.path.join(bins_dir, fa),
                 os.path.join(genome_dir, fa))

with open(Chdb_file) as f, open(csv_out, 'w') as fo:
    fo.write('genome,completeness,contamination,strain_heterogeneity,length,N50,QS\n')
    next(f)
    for line in f:
        p = line.strip('\n').split('\t')
        com, con = float(p[5]), float(p[6])
        score = com - 5 * con
        if com >= min_completeness and con <= max_contamination and score >= min_quality_score:
            fo.write(f'{p[0]}.fa,{com:.2f},{con:.2f},'
                     f'{float(p[7]):.2f},{int(p[8])},{int(p[13])},{score:.2f}\n')

print('CSV summary written to', csv_out)
