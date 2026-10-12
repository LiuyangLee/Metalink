#!/usr/bin/env python3
# metalink modified version
import sys
import os
import argparse
import logging
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


def write_down_seq(seq_dict, path_to_file):
    # Two inputs are required:
    # (1) The dict of the sequence
    # (2) The path that you want to write your sequence down

    seq_file = open(path_to_file,"w")
    for head in seq_dict:
        seq_file.write(head + "\n")
        seq_file.write(seq_dict[head] + "\n")
    seq_file.close()

def store_seq_with_full_head(input_seq_file): # The input sequence file should be a file with full path
    head = "" # Store the header line
    seq_dict = {} # Store the sequence dict

    with open(input_seq_file, "r") as seq_lines:
        for line in seq_lines:
            line = line.replace("\n","") # Remove "\n"
            if line.startswith('>'):
                head = line
                seq_dict[head] = ""
            else:
                seq_dict[head] += line

    seq_lines.close()

    return seq_dict


def make_all_custom_MAGs_combined_fasta_file_and_scaffold2MAG_map(custom_MAGs_dir, all_custom_MAGs_combined_fasta_file):
    scaffold2MAG_map = {} # scaffold => [MAG, 'not_viral']
    fasta_files = glob(os.path.join(custom_MAGs_dir, '*.fna'))
    if not fasta_files:
        sys.exit(f'Please make sure there are input MAGs in {custom_MAGs_dir} and all of them end with ".fna"')

    all_custom_MAGs_combined_fasta_seq = {} # Store the seqs of MAGs
    for fasta_file in fasta_files:
        fasta_seq = store_seq_with_full_head(fasta_file)
        MAG = Path(fasta_file).stem
        for header, seq in fasta_seq.items():
            scaffold = header[1:]
            scaffold2MAG_map[scaffold] = [MAG, 'not_viral']
            all_custom_MAGs_combined_fasta_seq[header] = seq
    write_down_seq(all_custom_MAGs_combined_fasta_seq, all_custom_MAGs_combined_fasta_file)
    return scaffold2MAG_map

def update_scaffold2MAG_map_according_to_geNomad_result_and_make_filtered_MAGs(outdir_iphop, scaffold2MAG_map):
    # Step 1: Efficiently create scaffold2length dict
    all_custom_MAGs_combined_seq = store_seq_with_full_head(os.path.join(outdir_iphop, 'all_custom_MAGs_combined.fna'))
    scaffold2length = {header[1:]: len(seq) for header, seq in all_custom_MAGs_combined_seq.items()}

    # Step 2: Efficiently parse the geNomad result
    viral_scaffold_set = set()
    scaffold_containing_provirus2provirus_length = defaultdict(int)
    virus_summary_file = os.path.join(outdir_iphop, 'genomad_output', 'all_custom_MAGs_combined_summary/all_custom_MAGs_combined_virus_summary.tsv')
    with open(virus_summary_file, 'r') as file:
        next(file)  # Skip the header line
        for line in file:
            seq, length = line.rstrip('\n').split('\t')[:2]
            length = int(length)
            if '|provirus_' in seq:
                scaffold = seq.split('|provirus_', 1)[0]
                scaffold_containing_provirus2provirus_length[scaffold] += length
            else:
                viral_scaffold_set.add(seq)

    # Step 3: Efficiently update scaffold2MAG_map
    scaffold2MAG_map_new = {}
    for scaffold, (MAG, _) in scaffold2MAG_map.items():
        viral_tag = 'viral' if scaffold in viral_scaffold_set or scaffold in scaffold_containing_provirus2provirus_length and float(scaffold_containing_provirus2provirus_length[scaffold]) / scaffold2length[scaffold] >= 0.85 else 'not_viral'
        scaffold2MAG_map_new[scaffold] = [MAG, viral_tag]

    # Step 4: Efficiently make filtered MAGs
    MAG2non_viral_scaffolds = defaultdict(list)
    for scaffold, (MAG, viral_tag) in scaffold2MAG_map_new.items():
        if viral_tag == 'not_viral':
            MAG2non_viral_scaffolds[MAG].append(scaffold)

    # Creating directory if it doesn't exist
    filtered_dir = os.path.join(outdir_iphop, 'custom_MAGs_filtered_dir')
    prepare_output_dir(filtered_dir)

    for MAG, scaffolds in MAG2non_viral_scaffolds.items():
        MAG_seq = {f'>{s}': all_custom_MAGs_combined_seq[f'>{s}'] for s in scaffolds}
        write_down_seq(MAG_seq, os.path.join(filtered_dir, f"{MAG}.fna"))

    return scaffold2MAG_map_new

def write_down_scaffold2MAG_map_file(scaffold2MAG_map_file, scaffold2MAG_map):
    with open(scaffold2MAG_map_file, 'w') as f:
        f.write("#Scaffold\tMAG\tViral tag\n")
        for scaffold in scaffold2MAG_map:
            MAG, viral_tag = scaffold2MAG_map[scaffold]
            f.write(f"{scaffold}\t{MAG}\t{viral_tag}\n")

def fetch_arguments(parser,root_dir,db_path_default):
    parser.set_defaults(func=main)
    parser.set_defaults(program="run")
    parser.add_argument('--input_Nlinked_viral_gn','-i', dest='input_Nlinked_viral_gn', required=False, default='./input_Nlinked_viral_gn', help=r'(required) input fasta of Nlinked vRhyme virus (default = ./input_Nlinked_viral_gn).')
    parser.add_argument('--custom_MAGs_dir', dest='custom_MAGs_dir', required=False, default='none', help=r'custom MAGs dir that contains only *.fna files for MAGs reconstructed from the same metagenome, this will be used in iPHoP for host prediction; note that it should be the absolute address path')
    parser.add_argument('--iPHoP_db_custom', dest='iPHoP_db_custom', required=False, default='none', help=r'custom iPHoP db that will be generated by your input custom MAGs; note that it should be the absolute address path')
    parser.add_argument('--iPHoP_db_custom_pre', dest='iPHoP_db_custom_pre', required=False, default='none', help=r'custom iPHoP db that has been made from the previous run, this will be used in iPHoP for host prediction by custom db; note that it should be the absolute address path')
    parser.add_argument('--outdir_iphop','-o', dest='outdir_iphop', required=False, default='./Metalink_MetaG_virus_sample_outdir', help=r'(required) output directory to deposit all results (default = ./Metalink_MetaG_virus_sample_outdir) output folder to deposit all results.')
    parser.add_argument('--iphop_db', dest='iphop_db', required=True, default='none', help=r'(required) path to the iPHoP database')
    parser.add_argument('--genomad_db', dest='genomad_db', required=True, default='none', help=r'(required) path to the geNomad database')
    parser.add_argument('--threads','-t', dest='threads', required=False, default=10, help=r'number of threads (default = 10)')
    parser.add_argument('--conda_env_iphop', dest='conda_env_iphop', required=True, default='none', help=r'(required) path to the iPHoP conda environment used by viwrap helper scripts')
    parser.add_argument('--conda_env_genomad', dest='conda_env_genomad', required=True, default='none', help=r'(required) path to the genomad conda environment used by viwrap helper scripts')
    parser.add_argument('--conda_env_gtdbtk', dest='conda_env_gtdbtk', required=True, default='none', help=r'(required) path to the gtdbtk conda environment used by viwrap helper scripts')
    parser.add_argument('--gtdbtk_db', dest='gtdbtk_db', required=True, default='none', help=r'(required) path to the GTDB-Tk database')
    parser.add_argument('--root_dir', dest='root_dir', required=False, default=root_dir,help=argparse.SUPPRESS)

def set_defaults(args):
    ## Store databases
    args['iPHoP_db'] = args['iphop_db']
    ## Store outdirs
    print(args['outdir_iphop'])
    args['outdir_iphop_custom'] = os.path.join(args['outdir_iphop'],f"outdir_iphop_custom_MAGs")

def main(args):
    # Welcome and logger
    print("### Welcome to Metalink.07_iPHoP ###\n")
    set_defaults(args)

    if os.path.exists(args['outdir_iphop']):
        shutil.rmtree(args['outdir_iphop'])

    logging.basicConfig(
        level=logging.INFO,
        format="%(message)s",
        handlers=[
            logging.StreamHandler(sys.stdout)
        ]
    )
    logger = logging.getLogger(__name__)

    # Give requirements for the usage of options "custom_MAGs_dir" and "iPHoP_db_custom"
    if args['custom_MAGs_dir'] == 'none':
        if not os.environ.get('METALINK_SILENT'):
            print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M')} INFO] custom_MAGs_dir is set to 'none'. iPHoP will run without custom MAGs database.", file=sys.stderr)
    elif args['custom_MAGs_dir'] != 'none' and not os.path.exists(args['custom_MAGs_dir']):
        #print(args['custom_MAGs_dir'])
        sys.exit(f"Could not find custom MAGs directory {args['custom_MAGs_dir']}. Maybe the directory is not correct")
    elif args['custom_MAGs_dir'] != 'none' and not os.path.isabs(args['custom_MAGs_dir']):
        #print(args['custom_MAGs_dir'])
        sys.exit(f"The input custom MAGs directory {args['custom_MAGs_dir']} is not an absolute address, please change it to an absolute address")
    elif args['custom_MAGs_dir'] != 'none' and os.path.exists(args['custom_MAGs_dir']):
        #print(args['custom_MAGs_dir'])
        for file in glob(os.path.join(args['custom_MAGs_dir'],'*.fna')):
            if '.fna' not in file:
                sys.exit(f"Make sure all MAGs in custom MAGs directory {args['custom_MAGs_dir']} end with '.fna', and no additional files within the directory")
        if args['iPHoP_db_custom'] == 'none':
            sys.exit(f"Make sure you have assigned a proper address to iPHoP_db_custom by using the option '--iPHoP_db_custom'")
        elif args['iPHoP_db_custom'] != 'none' and os.path.exists(args['iPHoP_db_custom']):
            #print(args['iPHoP_db_custom'])
            sys.exit(f"Please make sure that {args['iPHoP_db_custom']} is not present before Metalink-virus run. If present, please remove the folder. If you want to use the iPHoP_db_custom generated by your previous run, please use '--iPHoP_db_custom_pre' instead")
        elif args['iPHoP_db_custom'] != 'none' and not os.path.exists(args['iPHoP_db_custom']) and not os.path.isabs(args['iPHoP_db_custom']):
            #print(args['iPHoP_db_custom'])
            sys.exit(f"Please make sure that {args['iPHoP_db_custom']} is an absolute address")

    # Give requirements for the usage of the option of "iPHoP_db_custom_pre"
    if args['iPHoP_db_custom_pre'] != 'none':
        print(args['iPHoP_db_custom_pre'])
        if not os.path.exists(args['iPHoP_db_custom_pre']):
            sys.exit(f"Could not find 'iPHoP_db_custom_pre' directory {args['iPHoP_db_custom_pre']}. Maybe the directory is not correct")
        elif os.path.exists(args['iPHoP_db_custom_pre']) and not os.path.isabs(args['iPHoP_db_custom_pre']):
            sys.exit(f"The input 'iPHoP_db_custom_pre' directory {args['iPHoP_db_custom_pre']} is not an absolute address, please change it to an absolute address")
        elif args['custom_MAGs_dir'] != 'none' or args['iPHoP_db_custom'] != 'none':
            sys.exit(f"If you use --iPHoP_db_custom_pre, you should not use --custom_MAGs_dir or --iPHoP_db_custom. Simply assign the custom iPHoP db that has been made from the previous run")



    # Step 9 Host prediction
    time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
    logger.info(f"{time_current} | Conduct Host prediction by iPHoP. In processing...")

    ## Step 9.1 Host prediction by iPHoP

    iphop_cmd = f"conda run -n {args['conda_env_iphop']} python {os.path.join(args['root_dir'],'viwrap-scripts/run_iPHoP.py')} {args['input_Nlinked_viral_gn']} {args['outdir_iphop']} {args['iPHoP_db']} {args['threads']}"
    print(iphop_cmd)
    os.system(iphop_cmd)

    time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
    logger.info(f"{time_current} | Conduct Host prediction by iPHoP. Finished")

    ## Step 9.2 Host prediction by iPHoP by adding custom MAGs to host db
    if args['custom_MAGs_dir'] != 'none' and args['iPHoP_db_custom'] != 'none' and args['iPHoP_db_custom_pre'] == 'none':
        time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
        logger.info(f"{time_current} | Conduct Host prediction by iPHoP using custom MAGs. In processing...")

        all_custom_MAGs_combined_fasta_file = os.path.join(args['outdir_iphop'], 'all_custom_MAGs_combined.fna')
        scaffold2MAG_map = make_all_custom_MAGs_combined_fasta_file_and_scaffold2MAG_map(args['custom_MAGs_dir'], all_custom_MAGs_combined_fasta_file)
        geNomad_cmd = f"conda run -n {args['conda_env_genomad']} python {os.path.join(args['root_dir'],'viwrap-scripts/Metalink.run_geNomad.py')} {all_custom_MAGs_combined_fasta_file} {args['outdir_iphop']} {args['threads']} {args['genomad_db']}"

        os.system(geNomad_cmd)
        print(geNomad_cmd)

        scaffold2MAG_map = update_scaffold2MAG_map_according_to_geNomad_result_and_make_filtered_MAGs(args['outdir_iphop'], scaffold2MAG_map)
        custom_MAGs_filtered_dir  = os.path.join(args['outdir_iphop'], 'custom_MAGs_filtered_dir')
        args['custom_MAGs_GTDB'] = os.path.join(args['outdir_iphop'],'custom_MAGs_GTDB-tk_results')

        python_cmd1 = f"conda run -n {args['conda_env_gtdbtk']} python {os.path.join(args['root_dir'],'viwrap-scripts/Metalink.add_custom_MAGs_to_host_db__make_gtdbtk_results.py')} {args['custom_MAGs_GTDB']} {custom_MAGs_filtered_dir} {args['threads']} {args['gtdbtk_db']}"
        os.system(python_cmd1)
        print(python_cmd1)

        python_cmd2 = f"conda run -n {args['conda_env_iphop']} python {os.path.join(args['root_dir'],'viwrap-scripts/Metalink.add_custom_MAGs_to_host_db__add_to_db.py')} {args['custom_MAGs_GTDB']} {custom_MAGs_filtered_dir} {args['iPHoP_db']} {args['iPHoP_db_custom']}"
        os.system(python_cmd2)
        print(python_cmd2)

        iPHoP_cmd = f"conda run -n {args['conda_env_iphop']} python {os.path.join(args['root_dir'],'viwrap-scripts/run_iPHoP.py')} {args['input_Nlinked_viral_gn']} {args['outdir_iphop_custom']} {args['iPHoP_db_custom']} {args['threads']}"
        os.system(iPHoP_cmd)
        print(iPHoP_cmd)

        scaffold2MAG_map_file = os.path.join(args['outdir_iphop'], 'custom_MAGs_scaffold_filtering_map.txt')
        write_down_scaffold2MAG_map_file(scaffold2MAG_map_file, scaffold2MAG_map)
        #os.system(f"rm -rf {custom_MAGs_filtered_dir} {os.path.join(args['outdir_iphop'], 'genomad_output')}; rm {all_custom_MAGs_combined_fasta_file}")

        time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
        logger.info(f"{time_current} | Conduct Host prediction by iPHoP using custom MAGs. Finished")
    elif args['custom_MAGs_dir'] == 'none' and args['iPHoP_db_custom'] == 'none' and args['iPHoP_db_custom_pre'] != 'none': # iPHoP db custom was provided by the previous run
        time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
        logger.info(f"{time_current} | Conduct Host prediction by iPHoP using custom MAGs and using iPHoP db custom provided by the previous run. In processing...")

        os.system(f"conda run -n {args['conda_env_iphop']} python {os.path.join(args['root_dir'],'viwrap-scripts/run_iPHoP.py')} {args['input_Nlinked_viral_gn']} {args['outdir_iphop_custom']} {args['iPHoP_db_custom_pre']} {args['threads']}")

        time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
        logger.info(f"{time_current} | Conduct Host prediction by iPHoP using custom MAGs and using iPHoP db custom provided by the previous run. Finished")


def cli():
    root_dir = os.path.dirname(os.path.abspath(__file__)) # The dir that contains the Metalink-virus script
    current_dir = os.getcwd() # The current dir
    db_path_default = os.path.join(current_dir, "Metalink_virus_db") # Set the default Metalink_virus_db dir under the current directory

    parser = argparse.ArgumentParser(usage = argparse.SUPPRESS, formatter_class=argparse.RawDescriptionHelpFormatter,
        description = """Analyzing wrapper for microbiome

Usage: python Metalink.07_iPHoP.py [options]

Task: Predict virus and host interactions from metagenomes
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
