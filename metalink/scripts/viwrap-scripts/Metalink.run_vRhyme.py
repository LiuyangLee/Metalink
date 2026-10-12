#!/usr/bin/env python3

try:
    import warnings
    import sys
    import os
    import re
    warnings.filterwarnings("ignore")
except Exception as e:
    sys.stderr.write(str(e) + "\n\n")
    exit(1)

def run_vrhyme(viral_scaffold, vRhyme_outdir, mapping_bam, threads):
    viral_scaffold_ffn = viral_scaffold.rsplit(".", 1)[0] + ".ffn"
    viral_scaffold_faa = viral_scaffold.rsplit(".", 1)[0] + ".faa"

    cmd = f'vRhyme -i {viral_scaffold} -g {viral_scaffold_ffn} -p {viral_scaffold_faa} -b {mapping_bam} -t {int(threads)} -o {vRhyme_outdir} --red 5 1> /dev/null'
    os.system(cmd)

viral_scaffold, vRhyme_outdir, mapping_bam, threads = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
run_vrhyme(viral_scaffold, vRhyme_outdir, mapping_bam, threads)
