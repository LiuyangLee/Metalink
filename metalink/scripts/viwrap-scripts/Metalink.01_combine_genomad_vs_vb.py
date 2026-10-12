#!/usr/bin/env python3
# metalink modified version

import sys
import os
import argparse
#import scripts
#from scripts import module
from datetime import datetime
from pathlib import Path
from glob import glob
import shutil

# modified from ViWrap

try:
    import warnings
    import sys
    import os
    import re
    import pandas as pd
    from statistics import mean
    from collections import defaultdict
    warnings.filterwarnings("ignore")
    from pathlib import Path
    from glob import glob
except Exception as e:
    sys.stderr.write(str(e) + "\n\n")
    exit(1)

def prepare_output_dir(dir_path):
    if os.path.exists(dir_path):
        shutil.rmtree(dir_path)
    os.makedirs(dir_path, exist_ok=True)


def fetch_arguments(parser,root_dir,db_path_default):
    parser.set_defaults(func=main)
    parser.set_defaults(program="run")
    parser.add_argument('--input_metagenome','-i', dest='input_metagenome', required=True, default='none', help=r'(required) input metagenome assembly. It can be a metagenome or entire virome assembly. The extension of the input nucleotide sequence should be ".fna"')
    parser.add_argument('--out_dir','-o', dest='out_dir', required=False, default='./Metalink_MetaG_virus_sample_outdir', help=r'(required) output directory to deposit all results (default = ./Metalink_MetaG_virus_sample_outdir) output folder to deposit all results.')
    parser.add_argument('--threads','-t', dest='threads', required=False, default=10, help=r'number of threads (default = 10)')
    parser.add_argument('--conda_env_python', dest='conda_env_python', required=True, default='none', help=r'(required) path to the python conda environment used by viwrap helper scripts')
    parser.add_argument('--conda_env_seqkit', dest='conda_env_seqkit', required=True, default='none', help=r'(required) path to the sra conda environment used by viwrap helper scripts')
    parser.add_argument('--root_dir', dest='root_dir', required=False, default=root_dir,help=argparse.SUPPRESS)

def set_defaults(args):
    print(args['out_dir'])
    print(args['input_metagenome'])
    ## Store outdirs
    args['genomad_outdir'] = os.path.join(args['out_dir'],f"01_identification/{Path(args['input_metagenome']).stem}/genomad_result")
    args['virsorter_outdir'] = os.path.join(args['out_dir'],f"01_identification/{Path(args['input_metagenome']).stem}/vs_vb_result")
    args['genomad_vs_vb_outdir'] = os.path.join(args['out_dir'],f"01_identification/{Path(args['input_metagenome']).stem}/genomad_vs_vb_result")
    print(args['genomad_vs_vb_outdir'])

def main(args):
    # Welcome and logger
    print("### Welcome to Metalink.master_run.01_virus_identify ###\n")
	## Set up the logger
    ## Set the default args:
    set_defaults(args)
    prepare_output_dir(args['genomad_vs_vb_outdir'])

    # Step 3 combine_genomad_vs_vb
    # Call the python script to combine genomad and virsorter2 results
    cmd_args = [
        '--in_genomad', os.path.join(args['genomad_outdir'], 'scf2lytic_or_lyso.summary.txt'),
        '--in_virsorter2', os.path.join(args['virsorter_outdir'], 'scf2lytic_or_lyso.summary.txt'),
        '--out_3methods', os.path.join(args['genomad_vs_vb_outdir'], 'scf2lytic_or_lyso.summary.txt'),
        '--out_genomad', os.path.join(args['genomad_vs_vb_outdir'], 'list.genomad.tsv'),
        '--out_virsorter2', os.path.join(args['genomad_vs_vb_outdir'], 'list.virsorter2_vibrant.tsv')
    ]
    args_str = ' '.join(cmd_args)
    #combine_genomad_vs_vb_cmd = f"conda run -n {os.path.join(args['conda_env_dir'], 'R4.2')} Rscript {os.path.join(args['root_dir'], 'viwrap-scripts/combine_genomad_vs2_vb.cmd.R')} {args_str}"
    combine_genomad_vs_vb_cmd = f"conda run -n {args['conda_env_python']} python {os.path.join(args['root_dir'],'viwrap-scripts/Metalink.combine_genomad_vs_vb.py')} {args_str}"
    os.system(combine_genomad_vs_vb_cmd)
    print(combine_genomad_vs_vb_cmd)

    # seqkit command to extract the added virus from virsorter2
    seqkit_input = os.path.join(args['genomad_vs_vb_outdir'], 'list.virsorter2_vibrant.tsv')
    vs2_fasta = os.path.join(args['virsorter_outdir'], 'final_vs2_virus.fna')
    vs2_output = os.path.join(args['virsorter_outdir'], 'add.virsorter2.fna')
    seqkit_cmd = f"conda run -n {args['conda_env_seqkit']} seqkit grep --threads {args['threads']} -f {seqkit_input} {vs2_fasta} > {vs2_output}"
    os.system(seqkit_cmd)
    print(seqkit_cmd)

    # cat sequences of genomad and virsorter2
    genomad_fasta = os.path.join(args['genomad_outdir'], 'final_genomad_virus.fna')
    final_output = os.path.join(args['genomad_vs_vb_outdir'], 'final_virus.fna')
    cat_cmd = f"cat {genomad_fasta} {vs2_output} > {final_output}"
    os.system(cat_cmd)
    print(cat_cmd)

    genomad_annotation = os.path.join(args['genomad_outdir'], 'final_genomad_virus.annotation.txt')
    vs2_annotation     = os.path.join(args['virsorter_outdir'], 'final_vs2_virus.annotation.txt')
    seqkit_input       = os.path.join(args['genomad_vs_vb_outdir'], 'list.virsorter2_vibrant.tsv')
    final_annotation   = os.path.join(args['genomad_vs_vb_outdir'], 'final_virus.annotation.txt')

    with open(genomad_annotation) as src, open(final_annotation, 'w') as dst:
        dst.write(src.read())
    with open(seqkit_input) as f:
        wanted = set(line.rstrip('\n') for line in f if line.strip())
    with open(vs2_annotation) as fin, open(final_annotation, 'a') as fout:
        for line in fin:
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 2 and parts[1] in wanted:
                fout.write(line)

    # 1. ORF ID no header
    final_annotation = os.path.join(args['genomad_vs_vb_outdir'], 'final_virus.annotation.txt')
    orf_id_file      = os.path.join(args['genomad_vs_vb_outdir'], 'list.final_virus.orf_ids.tsv')
    extract_ORFID_cmd = f"awk 'NR>1 {{print $1}}' {final_annotation} > {orf_id_file}"
    os.system(extract_ORFID_cmd)
    print(extract_ORFID_cmd)

    # 2. faa
    genomad_faa = os.path.join(args['genomad_outdir'], 'final_genomad_virus.faa')
    vs2_faa     = os.path.join(args['virsorter_outdir'], 'final_vs2_virus.faa')
    final_faa   = os.path.join(args['genomad_vs_vb_outdir'], 'final_virus.faa')
    extract_faa_cmd = f"conda run -n {args['conda_env_seqkit']} cat {genomad_faa} {vs2_faa} | seqkit grep -f {orf_id_file} -o {final_faa}"
    os.system(extract_faa_cmd)
    print(extract_faa_cmd)

    # 2. ffn
    genomad_ffn = os.path.join(args['genomad_outdir'], 'final_genomad_virus.ffn')
    vs2_ffn     = os.path.join(args['virsorter_outdir'], 'final_vs2_virus.ffn')
    final_ffn   = os.path.join(args['genomad_vs_vb_outdir'], 'final_virus.ffn')
    extract_ffn_cmd = f"conda run -n {args['conda_env_seqkit']} cat {genomad_ffn} {vs2_ffn} | seqkit grep -f {orf_id_file} -o {final_ffn}"
    os.system(extract_ffn_cmd)
    print(extract_ffn_cmd)


def cli():
    root_dir = os.path.dirname(os.path.abspath(__file__)) # The dir that contains the Metalink-virus script
    current_dir = os.getcwd() # The current dir
    db_path_default = os.path.join(current_dir, "Metalink_virus_db") # Set the default Metalink_virus_db dir under the current directory

    parser = argparse.ArgumentParser(usage = argparse.SUPPRESS, formatter_class=argparse.RawDescriptionHelpFormatter,
        description = """Analyzing wrapper for microbiome

Usage: python Metalink.01_combine_genomad_vs_vb.py [options]

Task: Combine genomad and vs_vb virus to get the union results
        """,
	)

    fetch_arguments(parser, root_dir, db_path_default)

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = vars(parser.parse_args())
    args["func"](args)


if __name__ == "__main__":
    output = cli()
