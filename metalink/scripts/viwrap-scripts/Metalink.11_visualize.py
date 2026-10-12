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

def generate_result_visualization_inputs(Metalink_virus_visualization_outdir, Metalink_virus_summary_outdir, VIBRANT_db):
    Result_visualization_inputs_folder = os.path.join(Metalink_virus_visualization_outdir, 'Result_visualization_inputs')
    prepare_output_dir(Result_visualization_inputs_folder)
    # Step 1 Make input for bar-plot 1 - virus statistics
    virus_statistics_file = os.path.join(Result_visualization_inputs_folder, 'Virus_statistics.txt')
    f = open(virus_statistics_file, 'w')
    f.write('viral scaffold no.\tvirus no.\tspecies cluster no.\tvirus taxonomy info no.\tvirus with host prediction no.\n')
    virus_summary_info_df = pd.read_csv(os.path.join(Metalink_virus_summary_outdir, 'vMAG_summary_info.tsv'), sep = '\t', index_col = 0) # Column 0 used as the row labels of the dataframe
    viral_scaffold_no = virus_summary_info_df['scaffold_num'].sum() # Give the viral scaffold number
    virus_no = virus_summary_info_df.shape[0] # Give the number of rows
    species_cluster_no = len(open(os.path.join(Metalink_virus_summary_outdir, 'vMAG_species_summary_info.tsv')).readlines(  )) - 1

    no_of_virus_taxonomy_info = virus_summary_info_df[virus_summary_info_df['superkingdom'].str.startswith('d__Viruses:', na=False)].shape[0]
    no_of_virus_with_host_prediction = len(open(os.path.join(Metalink_virus_summary_outdir, 'Host_prediction_to_genus_m90.csv')).readlines(  )) - 1
    ## Convert all numbers into strings
    viral_scaffold_no = str(viral_scaffold_no)
    virus_no = str(virus_no)
    species_cluster_no = str(species_cluster_no)

    no_of_virus_taxonomy_info = str(no_of_virus_taxonomy_info)
    no_of_virus_with_host_prediction = str(no_of_virus_with_host_prediction)

    f.write(f"{viral_scaffold_no}\t{virus_no}\t{species_cluster_no}\t{no_of_virus_taxonomy_info}\t{no_of_virus_with_host_prediction}\n")
    f.close()

    # Step 2.1 Make input for pie-chart 1 - virus phylum relative abundance
    phylum2rel_abun = {}


    # for project level only summary the infor of best species
    #vMAG_species_summary_info_df = pd.read_csv(os.path.join(Metalink_virus_summary_outdir, 'vMAG_species_summary_info.tsv'),sep='\t', index_col=0)
    #virus2rel_abun = vMAG_species_summary_info_df['relative_abundance'].to_dict()
    # for sample level
    #virus_summary_info_df = pd.read_csv(os.path.join(Metalink_virus_summary_outdir, 'vMAG_summary_info.tsv'), sep = '\t', index_col = 0)
    virus2rel_abun = virus_summary_info_df['relative_abundance'].to_dict()

    virus2phylum = {}
    for virus, phylum in virus_summary_info_df['phylum'].items():
        if pd.isna(phylum) or phylum == 'no support':
            phylum_clean = 'unassigned'
        else:
            phylum_clean = phylum.split(':')[0].strip()
            if not phylum_clean:
                phylum_clean = 'unassigned'
        virus2phylum[virus] = phylum_clean

    phylum2virus = defaultdict(list)
    for virus in virus2phylum:
        phylum2virus[virus2phylum[virus]].append(virus)
    phylum2rel_abun = {
        f: sum(virus2rel_abun[v] for v in vs)
        for f, vs in phylum2virus.items()
    }
    with open(os.path.join(Result_visualization_inputs_folder, 'Virus_phylum_relative_abundance.txt'), 'w') as f:
        f.write('phylum\trelative abundance\n')
        for fam, abun in phylum2rel_abun.items():
            f.write(f"{fam}\t{abun}\n")

    # Step 2.2 Make input for pie-chart 1 - virus family relative abundance
    family2rel_abun = {}
    virus2family = {}
    for virus, family in virus_summary_info_df['family'].items():
        if pd.isna(family) or family == 'no support':
            family_clean = 'unassigned'
        else:
            family_clean = family.split(':')[0].strip()
            if not family_clean:
                family_clean = 'unassigned'
        virus2family[virus] = family_clean

    family2virus = defaultdict(list)
    for virus in virus2family:
        family2virus[virus2family[virus]].append(virus)
    family2rel_abun = {
        f: sum(virus2rel_abun[v] for v in vs)
        for f, vs in family2virus.items()
    }
    with open(os.path.join(Result_visualization_inputs_folder, 'Virus_family_relative_abundance.txt'), 'w') as f:
        f.write('family\trelative abundance\n')
        for fam, abun in family2rel_abun.items():
            f.write(f"{fam}\t{abun}\n")

    # Step 3 Make input for bar-plot 2 - KO ID relative abundance
    ko2rel_abun = {} # ko => rel_abun

    ko2rel_abun = defaultdict(float)
    sum_rel_abun_of_ko = 0


    #for _, row in vMAG_species_summary_info_df.iterrows(): # for project level only summary the infor of best species
    for _, row in virus_summary_info_df.iterrows():
        virus = row.name                # 行索引就是 Virus 名
        amg_str = row['AMG_KOs']
        if pd.isna(amg_str) or amg_str == '':
            continue
        for item in amg_str.split(';'):
            if not item.strip():
                continue
            # 解析 K00018(3) 格式
            ko, copy_str = item.split('(')[0], item.split('(')[1].replace(')', '')
            copy = int(copy_str)
            contrib = copy * virus2rel_abun[virus]
            ko2rel_abun[ko] += contrib
            sum_rel_abun_of_ko += contrib

    if sum_rel_abun_of_ko != 0:
        for ko in ko2rel_abun:
            ko2rel_abun[ko] = ko2rel_abun[ko] / sum_rel_abun_of_ko

    f = open(os.path.join(Result_visualization_inputs_folder, 'KO_ID_relative_abundance.txt'), 'w')
    f.write('KO\trelative abundance\n')
    for ko in ko2rel_abun:
        f.write(f"{ko}\t{ko2rel_abun[ko]}\n")
    f.close()

    # Step 4 Make input for pie-chart 2 - KO metabolism relative abundance
    metabolism2kos = {}  # metabolism => [kos]
    ko2metabolism = {}   # ko => metabolism
    KEGG_pathway_file = os.path.join(VIBRANT_db, 'files/VIBRANT_KEGG_pathways_summary.tsv')

    # Read the file and process lines in one go
    with open(KEGG_pathway_file, 'r') as file:
        for line in file:
            line = line.rstrip('\n')
            if 'map' in line:
                _, metabolism, _, ko_array = line.split('\t', 3)  # Direct unpacking for efficiency
                metabolism = metabolism.strip('"')  # More efficient stripping
                kos = ko_array.split('~')
                metabolism2kos[metabolism] = kos
                for ko in kos:
                    ko2metabolism[ko] = metabolism

    # Collect all the existing metabolisms from the KO found in this run
    metabolism_of_this_run2kos = defaultdict(list)
    for ko in ko2rel_abun:
        metabolism = ko2metabolism.get(ko, 'Unclassified metabolism')
        metabolism_of_this_run2kos[metabolism].append(ko)

    # Calculate relative abundance
    metabolism_of_this_run2rel_abun = {}
    sum_rel_abun_of_metabolism_of_this_run = 0
    for metabolism, kos in metabolism_of_this_run2kos.items():
        rel_abun = sum(ko2rel_abun.get(ko, 0) for ko in kos)  # Efficient summing with generator expression
        metabolism_of_this_run2rel_abun[metabolism] = rel_abun
        sum_rel_abun_of_metabolism_of_this_run += rel_abun

    # Normalize relative abundance
    if sum_rel_abun_of_metabolism_of_this_run != 0:
        for metabolism in metabolism_of_this_run2rel_abun:
            metabolism_of_this_run2rel_abun[metabolism] /= sum_rel_abun_of_metabolism_of_this_run

    # Write to file
    output_file_path = os.path.join(Result_visualization_inputs_folder, 'KO_metabolism_relative_abundance.txt')
    with open(output_file_path, 'w') as f:
        f.write('KO metabolism\trelative abundance\n')
        for metabolism, rel_abun in metabolism_of_this_run2rel_abun.items():
            f.write(f"{metabolism}\t{rel_abun}\n")

def fetch_arguments(parser,root_dir,db_path_default):
    parser.set_defaults(func=main)
    parser.set_defaults(program="run")
    parser.add_argument('--indir_summary','-i', dest='indir_summary', required=True, default='./indir_summary', help=r'(required) input directory of summary infor.')
    parser.add_argument('--outdir_visualize','-o', dest='outdir_visualize', required=False, default='./outdir_visualize', help=r'(required) output directory to deposit all results (default = ./outdir_visualize) output folder to deposit all fisualization results.')
    parser.add_argument('--threads','-t', dest='threads', required=False, default=10, help=r'number of threads (default = 10)')
    parser.add_argument('--vibrant_db', dest='vibrant_db', required=True, default='none', help=r'(required) path to the VIBRANT database')
    parser.add_argument('--root_dir', dest='root_dir', required=False, default=root_dir,help=argparse.SUPPRESS)

def set_defaults(args):
    ## Store databases
    args['VIBRANT_db'] = args['vibrant_db']
    print(args['outdir_visualize'])

def main(args):
    # Welcome and logger
    print("### Welcome to Metalink.master_run.visualize ###\n")
	## Set up the logger
    set_defaults(args)
    prepare_output_dir(args['outdir_visualize'])
    time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
    print(f"{time_current} | Visualizing...")
    prepare_output_dir(args['outdir_visualize'])
    summary_file = os.path.join(args['indir_summary'], 'vMAG_species_summary_info.tsv')

    summary_EMPTY = False
    if not os.path.exists(summary_file) or os.path.getsize(summary_file) == 0:
        summary_EMPTY = True
        print(f"{time_current} | Input vMAG_species_summary_info.tsv has no information, skipping all visualization steps.")

    if not summary_EMPTY:
        # Step 12 Visualize the result
        generate_result_visualization_inputs(args['outdir_visualize'], args['indir_summary'], args['VIBRANT_db'])
        visualization_input_dir = os.path.join(args['outdir_visualize'],'Result_visualization_inputs')
        visualization_output_dir = os.path.join(args['outdir_visualize'],'Result_visualization_outputs')

        os.system(f"python {os.path.join(args['root_dir'],'viwrap-scripts/Metalink.run_Visualization.py')} -i {visualization_input_dir} -r {args['outdir_visualize']} -o {visualization_output_dir}")

        time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
        print(f"{time_current} | Visualize the results. Finished")

def cli():
    root_dir = os.path.dirname(os.path.abspath(__file__)) # The dir that contains the Metalink-virus script
    current_dir = os.getcwd() # The current dir
    db_path_default = os.path.join(current_dir, "Metalink_virus_db") # Set the default Metalink_virus_db dir under the current directory

    parser = argparse.ArgumentParser(usage = argparse.SUPPRESS, formatter_class=argparse.RawDescriptionHelpFormatter,
        description = """Analyzing wrapper for microbiome

Usage: python Metalink.11_visualize.py [options]

Task: Visualize Metalink outputs for virus genomes from metagenomes
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
