#!/usr/bin/env python3

import sys
import os
import argparse
from datetime import datetime
from pathlib import Path
from glob import glob
import shutil
import warnings

warnings.filterwarnings("ignore")

def combine_iphop_results(args, combined_host_pred_to_genome_result, combined_host_pred_to_genus_result):
    host_pred_to_genome_header = ''
    host_pred_to_genome_m90_result = {}
    host_pred_to_genus_header = ''
    host_pred_to_genus_m90_result = {}

    if args['custom_MAGs_dir'] == 'none':
        host_pred_to_genome_m90 = os.path.join(args['indir_iphop'], "Host_prediction_to_genome_m90.csv")
        host_pred_to_genus_m90 = os.path.join(args['indir_iphop'], "Host_prediction_to_genus_m90.csv")

        with open(host_pred_to_genome_m90, 'r') as lines:
            for line in lines:
                line = line.strip('\n')
                if line.startswith('Virus,'):
                    host_pred_to_genome_header = line
                else:
                    tmp = line.split(',')
                    Virus, Confidence_score = tmp[0], tmp[4]
                    if Virus not in host_pred_to_genome_m90_result:
                        host_pred_to_genome_m90_result[Virus] = tmp
                    elif float(Confidence_score) > float(host_pred_to_genome_m90_result[Virus][4]):
                        host_pred_to_genome_m90_result[Virus] = tmp

        with open(host_pred_to_genus_m90, 'r') as lines:
            for line in lines:
                line = line.strip('\n')
                if line.startswith('Virus,'):
                    host_pred_to_genus_header = line
                else:
                    tmp = line.split(',')
                    Virus, Confidence_score = tmp[0], tmp[3]
                    if Virus not in host_pred_to_genus_m90_result:
                        host_pred_to_genus_m90_result[Virus] = tmp
                    elif float(Confidence_score) > float(host_pred_to_genus_m90_result[Virus][3]):
                        host_pred_to_genus_m90_result[Virus] = tmp

    if args['custom_MAGs_dir'] != 'none' or args['iPHoP_db_custom_pre'] != 'none':
        host_pred_to_genome_m90_custom = os.path.join(args['iphop_custom_outdir'], "Host_prediction_to_genome_m90.csv")
        host_pred_to_genus_m90_custom = os.path.join(args['iphop_custom_outdir'], "Host_prediction_to_genus_m90.csv")
        with open(host_pred_to_genome_m90_custom, 'r') as lines:
            for line in lines:
                line = line.strip('\n')
                if line.startswith('Virus,'):
                    host_pred_to_genome_header = line
                else:
                    tmp = line.split(',')
                    Virus, Confidence_score = tmp[0], tmp[4]
                    if Virus not in host_pred_to_genome_m90_result:
                        host_pred_to_genome_m90_result[Virus] = tmp
                    elif float(Confidence_score) > float(host_pred_to_genome_m90_result[Virus][4]):
                        host_pred_to_genome_m90_result[Virus] = tmp

        with open(host_pred_to_genus_m90_custom, 'r') as lines:
            for line in lines:
                line = line.strip('\n')
                if line.startswith('Virus,'):
                    host_pred_to_genus_header = line
                else:
                    tmp = line.split(',')
                    Virus, Confidence_score = tmp[0], tmp[3]
                    if Virus not in host_pred_to_genus_m90_result:
                        host_pred_to_genus_m90_result[Virus] = tmp
                    elif float(Confidence_score) > float(host_pred_to_genus_m90_result[Virus][3]):
                        host_pred_to_genus_m90_result[Virus] = tmp

    with open(combined_host_pred_to_genome_result, 'w') as f:
        f.write(host_pred_to_genome_header + '\n')
        for virus in sorted(host_pred_to_genome_m90_result):
            f.write(','.join(host_pred_to_genome_m90_result[virus]) + '\n')

    with open(combined_host_pred_to_genus_result, 'w') as f:
        f.write(host_pred_to_genus_header + '\n')
        for virus in sorted(host_pred_to_genus_m90_result):
            f.write(','.join(host_pred_to_genus_m90_result[virus]) + '\n')


def fetch_arguments(parser, root_dir, db_path_default):
    parser.set_defaults(func=main)
    parser.add_argument('--outdir_summary', '-o', dest='outdir_summary', required=False,
                        default='./Metalink_MetaG_virus_sample_outdir',
                        help='output directory to deposit all results.')
    parser.add_argument('--indir_iphop','-i', dest='indir_iphop', required=False, default='./indir_iphop', help=r'(required) output directory to deposit all results (default = ./indir_iphop) input folder of iphop results.')
    parser.add_argument('--root_dir', dest='root_dir', required=False, default=root_dir,
                        help=argparse.SUPPRESS)
    parser.add_argument('--custom_MAGs_dir', dest='custom_MAGs_dir', required=False, default='none',
                        help='custom MAGs dir that contains only *.fna files for MAGs reconstructed from the same metagenome, this will be used in iPHoP for host prediction; note that it should be the absolute address path')
    parser.add_argument('--iPHoP_db_custom', dest='iPHoP_db_custom', required=False, default='none',
                        help='custom iPHoP db that will be generated by your input custom MAGs; note that it should be the absolute address path')
    parser.add_argument('--iPHoP_db_custom_pre', dest='iPHoP_db_custom_pre', required=False, default='none',
                        help='custom iPHoP db that has been made from the previous run, this will be used in iPHoP for host prediction by custom db; note that it should be the absolute address path')


def set_defaults(args):
    args['iphop_custom_outdir'] = os.path.join(args['indir_iphop'], f"outdir_iphop_custom_MAGs")


def main(args):
    print("### Welcome to Metalink.10_summary_project ###\n")

    set_defaults(args)
    os.makedirs(args['outdir_summary'], exist_ok=True)

    # Step 11 merge host prediction
    combined_host_pred_to_genome_result = os.path.join(args['outdir_summary'],
                                                       'Host_prediction_to_genome_m90.csv')
    combined_host_pred_to_genus_result = os.path.join(args['outdir_summary'],
                                                      'Host_prediction_to_genus_m90.csv')
    combine_iphop_results(args, combined_host_pred_to_genome_result, combined_host_pred_to_genus_result)

    time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
    print(f"{time_current} | All summary tasks finished.")


def cli():
    root_dir = os.path.dirname(os.path.abspath(__file__))
    current_dir = os.getcwd()
    db_path_default = os.path.join(current_dir, "Metalink_virus_db")

    parser = argparse.ArgumentParser(
        usage=argparse.SUPPRESS,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""Analyzing wrapper for microbiome
Usage: python Metalink.10_summary_project.py [options]
Task: Summary Metalink outputs for virus genomes from metagenomes
        """
    )

    fetch_arguments(parser, root_dir, db_path_default)

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = vars(parser.parse_args())
    args["func"](args)


if __name__ == "__main__":
    cli()
