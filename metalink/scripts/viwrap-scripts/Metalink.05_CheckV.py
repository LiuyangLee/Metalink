#!/usr/bin/env python3

import sys
import os
import argparse
import shutil
#import scripts
#from scripts import module
from datetime import datetime
from pathlib import Path
from glob import glob

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

def store_seq(input_seq_file): # The input sequence file should be a file with full path
    head = "" # Store the header line
    seq_dict = {} # Store the sequence dict

    with open(input_seq_file, "r") as seq_lines:
        for line in seq_lines:
            line = line.rstrip("\n") # Remove "\n" in the end
            if ">" in line:
                head = line.split(None, 1)[0] # Cut at the first " " or "\t", use the first part
                seq_dict[head] = ""
            else:
                seq_dict[head] += line

    return seq_dict

def Nlinker(infolder, outdir, extension, n):
# Copied from vRhyme auxiliary scripts by Kristopher Kieft, UW-Madison
# Using n of Ns to link scaffolds
    files = os.listdir(infolder)
    len_ext = len(extension)
    files = [i for i in files if i[-len_ext:] == extension]

    N_string = ''.join("N" * n)
    for f in files:
        file = infolder + "/" + f
        base = f.rsplit(".",1)[0]
        with open(file, 'r') as fasta:
            seq = ''
            for line in fasta:
                if not line.startswith(">"):
                    seq += line.strip("\n")
                    seq += N_string

        with open(outdir + "/" + base + '.linked.' + extension, "w") as outfile:
            outfile.write(">" + base + "\n" + seq[:-n] + "\n") # -n to remove last n added characters

def parse_checkv_result(input_dir, outfile):
    header = ""
    checkv_summary = []
    for entry in os.scandir(input_dir):
        if not entry.is_dir():
            continue
        if entry.name == 'Nlinked_viral_gn':
            continue
        summary_file = os.path.join(entry.path, 'quality_summary.tsv')
        if not os.path.isfile(summary_file):
            continue
        with open(summary_file) as fh:
            for line in fh:
                line = line.rstrip('\n')
                if line.startswith('contig_id'):
                    header = line
                else:
                    checkv_summary.append(line)
    with open(outfile, 'w') as out:
        out.write(header + '\n')
        for line in checkv_summary:
            out.write(line + '\n')

def write_down_seq(seq_dict, path_to_file):
    # Two inputs are required:
    # (1) The dict of the sequence
    # (2) The path that you want to write your sequence down

    seq_file = open(path_to_file,"w")
    for head in seq_dict:
        seq_file.write(head + "\n")
        seq_file.write(seq_dict[head] + "\n")
    seq_file.close()

def combine_all_vRhyme_fasta(vRhyme_best_bin_dir, vRhyme_unbinned_viral_gn_dir, all_vRhyme_fasta):
    walk = os.walk(vRhyme_best_bin_dir)
    walk2 = ''
    if vRhyme_unbinned_viral_gn_dir:
        walk2 = os.walk(vRhyme_unbinned_viral_gn_dir)
    all_vRhyme_fasta_dict = {}

    for path, dir_list, file_list in walk:
        for file_name in file_list:
            if "fna" in file_name:
                file_name_with_path = os.path.join(path, file_name)
                each_vRhyme_faa_dict = store_seq(file_name_with_path)
                all_vRhyme_fasta_dict.update(each_vRhyme_faa_dict) # Add dict from one to another one

    if walk2:
        for path, dir_list, file_list in walk2:
            for file_name in file_list:
                if "fna" in file_name:
                    file_name_with_path = os.path.join(path, file_name)
                    each_vRhyme_faa_dict = store_seq(file_name_with_path)
                    all_vRhyme_fasta_dict.update(each_vRhyme_faa_dict) # Add dict from one to another one

    write_down_seq(all_vRhyme_fasta_dict, all_vRhyme_fasta)

def fetch_arguments(parser,root_dir,db_path_default):
    parser.set_defaults(func=main)
    parser.set_defaults(program="run")
    parser.add_argument('--indir_vrhyme','-i', dest='indir_vrhyme', required=False, default='./indir_vrhyme', help=r'(required) output directory to deposit all results (default = ./indir_vrhyme) output folder to deposit all results.')
    parser.add_argument('--outdir_checkv','-o', dest='outdir_checkv', required=False, default='./outdir_checkv', help=r'(required) output directory to deposit all results (default = ./outdir_checkv) output folder to deposit all results.')
    parser.add_argument('--checkv_db', dest='checkv_db', required=True, default='none', help=r'(required) path to the CheckV database')
    parser.add_argument('--threads','-t', dest='threads', required=False, default=10, help=r'number of threads (default = 10)')
    parser.add_argument('--conda_env_checkv', dest='conda_env_checkv', required=True, default='none', help=r'(required) path to the checkm conda environment used by viwrap helper scripts')
    parser.add_argument('--root_dir', dest='root_dir', required=False, default=root_dir,help=argparse.SUPPRESS)

def set_defaults(args):
    ## Store databases
    args['CheckV_db'] = args['checkv_db']

    print(args['outdir_checkv'])
    print(args['CheckV_db'])

    ## Store outdirs
    args['nlinked_viral_gn_dir'] = os.path.join(args['outdir_checkv'],f"Nlinked_viral_gn/")

def main(args):
    # Welcome and logger
    print("### Welcome to Metalink.checkv ###\n")
	## Set up the logger
    set_defaults(args)
    prepare_output_dir(args['outdir_checkv'])
    # Step 6 Run CheckV
    time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
    print(f"{time_current} | Run CheckV to evaluate virus genome quality. In processing...")

    ## Step 6.1 Link multiple scaffolds within a bin
    vRhyme_best_bin_dir_modified = os.path.join(args['indir_vrhyme'], 'vRhyme_best_bins_fasta_modified')
    vRhyme_unbinned_viral_gn_dir = os.path.join(args['indir_vrhyme'], 'vRhyme_unbinned_viral_gn_fasta')
    prepare_output_dir(args['nlinked_viral_gn_dir'])
    Nlinker(vRhyme_best_bin_dir_modified, args['nlinked_viral_gn_dir'], 'fna', 1000)
    Nlinker(vRhyme_unbinned_viral_gn_dir, args['nlinked_viral_gn_dir'], 'fna', 1000)

    ## Step 6.2 Run CheckV in parallel and parse the result
    os.system(f"conda run -n {args['conda_env_checkv']} python {os.path.join(args['root_dir'],'viwrap-scripts/run_CheckV.py')} {args['nlinked_viral_gn_dir']} {args['outdir_checkv']} {args['threads']} {args['CheckV_db']} >/dev/null 2>&1")
    CheckV_quality_summary = os.path.join(args['outdir_checkv'], 'CheckV_quality_summary.txt')
    parse_checkv_result(args['outdir_checkv'], CheckV_quality_summary)

    all_vRhyme_fasta_Nlinked = os.path.join(args['outdir_checkv'], 'all_vRhyme_fasta.Nlinked_viral_gn.fna')
    combine_all_vRhyme_fasta(args['nlinked_viral_gn_dir'], '', all_vRhyme_fasta_Nlinked)

    time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
    print(f"{time_current} | Run CheckV to evaluate virus genome quality. Finished")



def cli():
    root_dir = os.path.dirname(os.path.abspath(__file__)) # The dir that contains the Metalink-virus script
    current_dir = os.getcwd() # The current dir
    db_path_default = os.path.join(current_dir, "Metalink_virus_db") # Set the default Metalink_virus_db dir under the current directory

    parser = argparse.ArgumentParser(usage = argparse.SUPPRESS, formatter_class=argparse.RawDescriptionHelpFormatter,
        description = """Analyzing wrapper for microbiome

Usage: python Metalink.05_CheckV.py [options]

Task: Check virus genomes quality
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
