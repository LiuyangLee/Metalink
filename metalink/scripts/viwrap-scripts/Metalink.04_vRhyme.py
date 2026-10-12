#!/usr/bin/env python3
# metalink modified version
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

def write_down_seq(seq_dict, path_to_file):
    # Two inputs are required:
    # (1) The dict of the sequence
    # (2) The path that you want to write your sequence down

    seq_file = open(path_to_file,"w")
    for head in seq_dict:
        seq_file.write(head + "\n")
        seq_file.write(seq_dict[head] + "\n")
    seq_file.close()

def get_vRhyme_best_bin_lytic_and_lysogenic_info(vRhyme_best_bin_dir, vrhyme_outdir, scf2lytic_or_lyso_summary):
    # Step 1 Get vRhyme_bin2scf dict
    vRhyme_bin2scf = {} # vRhyme_bin => [scfs]; scf here is the short scf name
    vRhyme_bin_addrs = glob(os.path.join(vRhyme_best_bin_dir, '*.fna'))
    for vRhyme_bin_addr in vRhyme_bin_addrs:
        vRhyme_bin = Path(vRhyme_bin_addr).stem
        vRhyme_bin_seq = store_seq(vRhyme_bin_addr)
        scfs = []
        for header in vRhyme_bin_seq:
            scf = header.replace('>', '', 1).split('__', 1)[1]
            scfs.append(scf)
        vRhyme_bin2scf[vRhyme_bin] = scfs

    # Step 2 Store scf2lytic_or_lyso dict
    scf2lytic_or_lyso = {} # scf => [lytic_or_lyso_or_integrated_prophage, integrase_presence_or_absence]
    if scf2lytic_or_lyso_summary:
        with open(scf2lytic_or_lyso_summary, 'r') as lines:
            for line in lines:
                line = line.rstrip('\n')
                tmp = line.split('\t')
                if tmp[0] != 'scaffold':
                    scf, lytic_or_lyso_or_integrated_prophage, integrase_presence_or_absence = tmp[0], tmp[1], tmp[2]
                    if ' (' in lytic_or_lyso_or_integrated_prophage:
                        lytic_or_lyso_or_integrated_prophage = lytic_or_lyso_or_integrated_prophage.split(' (', 1)[0]
                    scf2lytic_or_lyso[scf] = [lytic_or_lyso_or_integrated_prophage, integrase_presence_or_absence]

    # Step 3 Get vRhyme_bin2lytic_and_lyso_info dict
    vRhyme_bin2lytic_and_lyso_info = {} # vRhyme_bin => lytic_and_lyso_info ([pattern, assignment], for example, ['1 lysogenic_scaffold + N lytic_scaffold', 'lysogenic_virus'])

    # Rules
    # A bin with one or more lytic members and one lysogenic member should not cause concern.
    # A bin with one or more lytic members and one integrated prophage should be examined.
    # A bin with two or more lysogenic members, each encoding an integrase, is likely contamination.
    # A bin with two or more integrated prophages, regardless of integrases, from multiple parent sequences is likely contamination.

    # Formulas:
    # Case 1: 1 lysogenic_scaffold + N lytic_scaffold => lysogenic_virus
    # Case 2: 1 integrated_prophage + N lytic_scaffold => split into scaffolds
    # Case 3: N (N >= 2) of lysogenic_scaffold or integrated_prophage + N lytic_scaffold => split into scaffolds
    # Case 4: N lytic_scaffold => lytic_virus
    for vRhyme_bin in vRhyme_bin2scf:
        scfs = vRhyme_bin2scf[vRhyme_bin]
        no_of_lysogenic_scaffold = 0
        no_of_integrated_prophage = 0
        no_of_lytic_scaffold = 0
        for scf in scfs:
            if scf2lytic_or_lyso:
                if scf2lytic_or_lyso[scf][0] == 'lysogenic_scaffold':
                    no_of_lysogenic_scaffold = no_of_lysogenic_scaffold + 1
                elif scf2lytic_or_lyso[scf][0] == 'integrated_prophage':
                    no_of_integrated_prophage = no_of_integrated_prophage + 1
                elif scf2lytic_or_lyso[scf][0] == 'lytic_scaffold':
                    no_of_lytic_scaffold = no_of_lytic_scaffold + 1

        if no_of_lysogenic_scaffold == 1 and no_of_integrated_prophage == 0 and no_of_lytic_scaffold >= 1:
            vRhyme_bin2lytic_and_lyso_info[vRhyme_bin] = ['1 lysogenic_scaffold + N lytic_scaffold', 'lysogenic_virus']
        elif no_of_integrated_prophage == 1 and no_of_lysogenic_scaffold == 0 and no_of_lytic_scaffold >= 1:
            vRhyme_bin2lytic_and_lyso_info[vRhyme_bin] = ['1 integrated_prophage + N lytic_scaffold', 'split into scaffolds']
        elif (no_of_lysogenic_scaffold + no_of_integrated_prophage) >= 2 and no_of_lytic_scaffold >= 0:
            vRhyme_bin2lytic_and_lyso_info[vRhyme_bin] = ['N (N > 1) of lysogenic_scaffold or integrated_prophage + N lytic_scaffold', 'split into scaffolds']
        elif no_of_lysogenic_scaffold == 0 and no_of_integrated_prophage == 0 and  no_of_lytic_scaffold >= 2:
            vRhyme_bin2lytic_and_lyso_info[vRhyme_bin] = ['N lytic_scaffold', 'lytic_virus']

    # Step 4 Write down vRhyme_best_bin_lytic_and_lysogenic_info
    f = open(os.path.join(vrhyme_outdir, 'vRhyme_best_bin_lytic_and_lysogenic_info.txt'),'w')
    f.write("vRhyme_bin\tpattern\tassignment\n")
    if vRhyme_bin2lytic_and_lyso_info:
        for vRhyme_bin in vRhyme_bin2lytic_and_lyso_info:
            line = vRhyme_bin + '\t' + '\t'.join(vRhyme_bin2lytic_and_lyso_info[vRhyme_bin]) + '\n'
            f.write(line)
    f.close()

def get_vRhyme_best_bin_scaffold_complete_info(CheckV_quality_summary, vRhyme_best_bin_scaffold_complete_info):
    vRhyme_best_bin_scaffold_complete_info_dict = {} # scf => [vRhyme_bin, Complete or Not-complete]
    with open(CheckV_quality_summary,'r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            tmp = line.split('\t')
            if tmp[0] != 'contig_id':
                scf, complete_info = tmp[0], tmp[7]
                scf = scf.split('__', 1)[1]  # scf here is the short scf name
                vRhyme_bin = scf.split('__', 1)[0].replace('vRhyme_', 'vRhyme_bin_', 1)
                if complete_info == 'Complete':
                    vRhyme_best_bin_scaffold_complete_info_dict[scf] = [vRhyme_bin, 'Complete']
                elif complete_info != 'Complete':
                    vRhyme_best_bin_scaffold_complete_info_dict[scf] = [vRhyme_bin, 'Not-complete']
    lines.close()

    f = open(vRhyme_best_bin_scaffold_complete_info, 'w')
    f.write("scaffold\tcomplete_info\n")
    for scf in vRhyme_best_bin_scaffold_complete_info_dict:
        line = scf + '\t' + '\t'.join(vRhyme_best_bin_scaffold_complete_info_dict[scf]) + '\n'
        f.write(line)
    f.close()

def parse_checkv_result(input_dir, outfile):
    walk = os.walk(input_dir)
    header = "" # Store the header line of quality summary file
    checkv_summary = [] # Store all checkv_summary_result
    for path, dir_list, file_list in walk:
        for dir_name in dir_list:
            dir_name_with_path = os.path.join(path, dir_name)

            # Get only the first level dir
            if 'tmp' not in dir_name_with_path:
                each_quality_summary_file = f'{dir_name_with_path}/quality_summary.tsv'
                with open(each_quality_summary_file, "r") as each_file:
                    for line in each_file:
                        line = line.rstrip("\n")
                        if line.startswith("contig_id"):
                            header = line
                        else:
                            checkv_summary.append(line)
                each_file.close()

    f = open(outfile, "w")
    f.write(header + "\n")
    for line in checkv_summary:
        f.write(line + "\n")
    f.close()

def get_pro2viral_gn_map(vRhyme_best_bin_dir_modified, vRhyme_unbinned_viral_gn_dir, pro2viral_gn_map):
    pro2viral_gn_dict = {}
    walk = os.walk(vRhyme_best_bin_dir_modified)
    for path,dir_list,file_list in walk:
        for file_name in file_list:
            if "faa" in file_name:
                viral_gn = Path(file_name).stem
                file_name_with_path = os.path.join(path, file_name)
                with open(file_name_with_path, "r") as faa_lines:
                    for line in faa_lines:
                        if ">" in line:
                            line = line.rstrip("\n")
                            pro = line.replace(">", "", 1).split("\t", 1)[0]
                            pro2viral_gn_dict[pro] = viral_gn
                faa_lines.close()

    walk2 = os.walk(vRhyme_unbinned_viral_gn_dir)
    for path, dir_list, file_list in walk2:
        for file_name in file_list:
            if "faa" in file_name:
                viral_gn = Path(file_name).stem
                file_name_with_path = os.path.join(path, file_name)
                with open(file_name_with_path,"r") as faa_lines2:
                    for line in faa_lines2:
                        if ">" in line:
                            line = line.rstrip("\n")
                            pro = line.replace(">", "", 1).split("\t", 1)[0]
                            pro2viral_gn_dict[pro] = viral_gn
                faa_lines2.close()

    file = open(pro2viral_gn_map,"w")
    file.write('protein_id,contig_id\n')
    for pro in pro2viral_gn_dict:
        file.write(f'{pro},{pro2viral_gn_dict[pro]}\n')
    file.close()

def get_contig2viral_gn_map(vRhyme_best_bin_dir_modified, vRhyme_unbinned_viral_gn_dir, contig2viral_gn_map):
    contig2viral_gn_dict = {}
    walk = os.walk(vRhyme_best_bin_dir_modified)
    for path,dir_list,file_list in walk:
        for file_name in file_list:
            if "fna" in file_name:
                viral_gn = Path(file_name).stem
                file_name_with_path = os.path.join(path, file_name)
                with open(file_name_with_path, "r") as fna_lines:
                    for line in fna_lines:
                        if ">" in line:
                            line = line.rstrip("\n")
                            contig = line.replace(">", "", 1).split("\t", 1)[0]
                            contig2viral_gn_dict[contig] = viral_gn
                fna_lines.close()

    walk2 = os.walk(vRhyme_unbinned_viral_gn_dir)
    for path, dir_list, file_list in walk2:
        for file_name in file_list:
            if "fna" in file_name:
                viral_gn = Path(file_name).stem
                file_name_with_path = os.path.join(path, file_name)
                with open(file_name_with_path,"r") as fna_lines2:
                    for line in fna_lines2:
                        if ">" in line:
                            line = line.rstrip("\n")
                            contig = line.replace(">", "", 1).split("\t", 1)[0]
                            contig2viral_gn_dict[contig] = viral_gn
                fna_lines2.close()

    file = open(contig2viral_gn_map,"w")
    file.write('contig_newid,contig_id\n')
    for contig in contig2viral_gn_dict:
        file.write(f'{contig},{contig2viral_gn_dict[contig]}\n')
    file.close()

def make_vRhyme_best_bins_fasta_modified(vRhyme_best_bin_dir,
                                         vRhyme_best_bin_dir_modified,
                                         vRhyme_best_bin_lytic_and_lysogenic_info,
                                         vRhyme_best_bin_scaffold_complete_info,
                                         args):
    prefix = Path(args['input_metagenome']).stem
    vRhyme_bin_to_split = []
    with open(vRhyme_best_bin_lytic_and_lysogenic_info, 'r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            tmp = line.split('\t')
            if tmp[0] != 'vRhyme_bin' and tmp[2] == 'split into scaffolds':
                vRhyme_bin_to_split.append(tmp[0])
    with open(vRhyme_best_bin_scaffold_complete_info, 'r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            tmp = line.split('\t')
            if tmp[0] != 'scaffold' and tmp[2] == 'Complete':
                if tmp[1] not in vRhyme_bin_to_split:
                    vRhyme_bin_to_split.append(tmp[1])

    fasta_addrs = glob(os.path.join(vRhyme_best_bin_dir, '*.fna'))
    for fasta_addr in fasta_addrs:
        fasta_stem = Path(fasta_addr).stem
        new_stem   = f"{prefix}_{fasta_stem}"
        fasta_addr_new = os.path.join(vRhyme_best_bin_dir_modified, f"{new_stem}.fna")
        faa_addr       = os.path.join(vRhyme_best_bin_dir, f"{fasta_stem}.faa")
        faa_addr_new   = os.path.join(vRhyme_best_bin_dir_modified, f"{new_stem}.faa")
        ffn_addr       = os.path.join(vRhyme_best_bin_dir, f"{fasta_stem}.ffn")
        ffn_addr_new   = os.path.join(vRhyme_best_bin_dir_modified, f"{new_stem}.ffn")

        if fasta_stem not in vRhyme_bin_to_split:
            fasta_seq = store_seq(fasta_addr)
            faa_seq   = store_seq(faa_addr)
            ffn_seq   = store_seq(ffn_addr)
            def remap(d):
                out = {}
                for old_head, seq in d.items():
                    body = old_head.lstrip('>')
                    if '__' not in body:
                        out[old_head] = seq
                        continue
                    bin_part, scf_part = body.split('__', 1)
                    bin_id = bin_part.replace('vRhyme', "vRhyme_bin", 1)
                    # drop the duplicated sample prefix:
                    # {prefix}_{bin}__{prefix}---{cid} -> {prefix}_{bin}---{cid}
                    cid = scf_part.removeprefix(prefix + '---')
                    new_head = f'>{prefix}_{bin_id}---{cid}'
                    out[new_head] = seq
                return out
            write_down_seq(remap(fasta_seq), fasta_addr_new)
            write_down_seq(remap(faa_seq),   faa_addr_new)
            write_down_seq(remap(ffn_seq),   ffn_addr_new)


def make_unbinned_viral_gn(viral_scaffold,
                           vRhyme_best_bin_dir,
                           vRhyme_unbinned_viral_gn_dir,
                           args=None):          # 保持兼容，可接 args
    # --- 0. 路径与前缀 -------------------------------------------------------------
    viral_scaffold_faa = viral_scaffold.rsplit(".", 1)[0] + ".faa"
    viral_scaffold_ffn = viral_scaffold.rsplit(".", 1)[0] + ".ffn"

    prefix = Path(args['input_metagenome']).stem if args else ""   # 取样品前缀

    # --- 1. 读入序列字典 ----------------------------------------------------------
    viral_scaffold_fasta_dict = store_seq(viral_scaffold)
    viral_scaffold_faa_dict   = store_seq(viral_scaffold_faa)
    viral_scaffold_ffn_dict   = store_seq(viral_scaffold_ffn)

    # --- 2. 建立“已 binned”的 scaffold 集合 ---------------------------------------
    viral_scaffold_fasta_binned_dict = {}   # scaffold_id -> bin_name
    for path_, _, file_list in os.walk(vRhyme_best_bin_dir):
        for file_name in file_list:
            if file_name.endswith(".fna"):
                file_name_with_path = os.path.join(path_, file_name)
                for h in store_seq(file_name_with_path):
                    # modified bins carry {prefix}_{bin}---{cid}; rebuild the original
                    # contig id {prefix}---{cid} for membership checks
                    scf = f"{prefix}---" + h.lstrip('>').split('---', 1)[1]   # 取 scaffold 部分
                    viral_scaffold_fasta_binned_dict[scf] = 1

    # --- 3. 生成 unbinned scaffold 编号与名称 -------------------------------------
    viral_scaffold_fasta_unbinned_dict = {}   # scaffold_id -> unbinned_gn_name
    cnt = 1
    for header in viral_scaffold_fasta_dict:
        scaffold_id = header.lstrip('>')
        if scaffold_id not in viral_scaffold_fasta_binned_dict:
            name = f"{prefix}_vRhyme_unbinned_{cnt}" if prefix else f"vRhyme_unbinned_{cnt}"
            viral_scaffold_fasta_unbinned_dict[scaffold_id] = name
            cnt += 1

    # --- 4. 写 .fna --------------------------------------------------------------
    for scf, name in viral_scaffold_fasta_unbinned_dict.items():
        out_fna = os.path.join(vRhyme_unbinned_viral_gn_dir, f"{name}.fna")
        with open(out_fna, 'w') as f:
            f.write(f'>{name}---{scf.removeprefix(prefix + "---")}\n')
            f.write(f'{viral_scaffold_fasta_dict[">"+scf]}\n')

    # --- 5. 写 .faa --------------------------------------------------------------
    # 先建立 protein -> unbinned_gn_name 映射
    pro2gn = {}
    for pro_head in viral_scaffold_faa_dict:
        pro_id = pro_head.lstrip('>')
        scf_id = pro_id.rsplit('_', 1)[0]
        if scf_id in viral_scaffold_fasta_unbinned_dict:
            pro2gn[pro_id] = viral_scaffold_fasta_unbinned_dict[scf_id]

    # 按 scaffold 写文件
    for scf, name in viral_scaffold_fasta_unbinned_dict.items():
        out_faa = os.path.join(vRhyme_unbinned_viral_gn_dir, f"{name}.faa")
        with open(out_faa, 'w') as f:
            for pro_id, gn in pro2gn.items():
                if pro_id.rsplit('_', 1)[0] == scf:
                    f.write(f'>{gn}---{pro_id.removeprefix(prefix + "---")}\n')
                    f.write(f'{viral_scaffold_faa_dict[">"+pro_id]}\n')

    # --- 6. 写 .ffn  （与 .faa 完全对称） ----------------------------------------
    pro2gn_ffn = {}
    for pro_head in viral_scaffold_ffn_dict:
        pro_id = pro_head.lstrip('>')
        scf_id = pro_id.rsplit('_', 1)[0]
        if scf_id in viral_scaffold_fasta_unbinned_dict:
            pro2gn_ffn[pro_id] = viral_scaffold_fasta_unbinned_dict[scf_id]

    for scf, name in viral_scaffold_fasta_unbinned_dict.items():
        out_ffn = os.path.join(vRhyme_unbinned_viral_gn_dir, f"{name}.ffn")
        with open(out_ffn, 'w') as f:
            for pro_id, gn in pro2gn_ffn.items():
                if pro_id.rsplit('_', 1)[0] == scf:
                    f.write(f'>{gn}---{pro_id.removeprefix(prefix + "---")}\n')
                    f.write(f'{viral_scaffold_ffn_dict[">"+pro_id]}\n')

def combine_all_vRhyme_faa(vRhyme_best_bin_dir_modified,
                           vRhyme_unbinned_viral_gn_dir,
                           args):
    prefix = Path(args['input_metagenome']).stem
    out_faa = os.path.join(args['outdir_vrhyme'], f"all_vRhyme_faa.faa")

    merged = {}
    for d in (vRhyme_best_bin_dir_modified, vRhyme_unbinned_viral_gn_dir):
        for faa in glob(os.path.join(d, '*.faa')):
            merged.update(store_seq(faa))

    write_down_seq(merged, out_faa)
    return out_faa


def fetch_arguments(parser,root_dir,db_path_default):
    parser.set_defaults(func=main)
    parser.set_defaults(program="run")
    parser.add_argument('--viral_scaffold','-v', dest='viral_scaffold', required=True, default='none', help=r'(required) input viral scaffold. The extension of the input nucleotide sequence should be ".fna"')
    parser.add_argument('--scf2lytic_or_lyso_summary','-s', dest='scf2lytic_or_lyso_summary', required=True, default='none', help=r'(required) input viral scf2lytic_or_lyso_summary infor.')
    parser.add_argument('--input_metagenome','-i', dest='input_metagenome', required=True, default='none', help=r'(required) input metagenome assembly. It can be a metagenome or entire virome assembly. The extension of the input nucleotide sequence should be ".fna"')
    parser.add_argument('--outdir_vrhyme','-o', dest='outdir_vrhyme', required=False, default='./outdir_vrhyme', help=r'(required) output directory to deposit all results (default = ./outdir_vrhyme) output folder to deposit all results.')
    parser.add_argument('--threads','-t', dest='threads', required=False, default=10, help=r'number of threads (default = 10)')
    parser.add_argument('--conda_env_vrhyme', dest='conda_env_vrhyme', required=True, default='none', help=r'(required) path to the vRhyme conda environment used by viwrap helper scripts')
    parser.add_argument('--conda_env_checkv', dest='conda_env_checkv', required=True, default='none', help=r'(required) path to the checkm conda environment used by viwrap helper scripts')
    parser.add_argument('--root_dir', dest='root_dir', required=False, default=root_dir,help=argparse.SUPPRESS)
    parser.add_argument('--input_bam', '-bam', dest = 'input_bam', required=False, default='none', help=r'input bam file. When this option is used, --input_reads should not be used. This option will intake a pre-made coverage file that was generated by short reads metagenomic mapping')
    parser.add_argument('--checkv_db', dest='checkv_db', required=True, default='none', help=r'(required) path to the CheckV database')


def set_defaults(args):

    print(args['outdir_vrhyme'])
    print(args['input_metagenome'])
    print(args['viral_scaffold'])

    ## Store outdirs
    args['CheckV_db'] = args['checkv_db']


def main(args):
    # Welcome and logger
    print("### Welcome to Metalink.master_run.vRhyme ###\n")
	## Set up the logger
    ## Set the default args:
    set_defaults(args)

    if os.path.exists(args['outdir_vrhyme']):
        shutil.rmtree(args['outdir_vrhyme'])

    # Step 4 Run vRhyme
    time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
    print(f"{time_current} | Run vRhyme to bin viral scaffolds. In processing...")

    BIN_EMPTY = False
    VS_EMPTY = False
    if not os.path.exists(args['viral_scaffold']) or os.path.getsize(args['viral_scaffold']) == 0:
        VS_EMPTY = True


    if not VS_EMPTY:
        ## Step 4.1 Run vRhyme to get the original vRhyme_best_bins
        os.system(f"conda run -n {args['conda_env_vrhyme']} python {os.path.join(args['root_dir'],'viwrap-scripts/Metalink.run_vRhyme.py')} {args['viral_scaffold']} {args['outdir_vrhyme']} {args['input_bam']} {args['threads']}")
        vRhyme_best_bin_dir = os.path.join(args['outdir_vrhyme'], 'vRhyme_best_bins_fasta')
        if not os.path.isdir(vRhyme_best_bin_dir) or not glob(os.path.join(vRhyme_best_bin_dir, '*.fasta')):
            print(f"{time_current} | No valid bins, creating empty placeholder files for bin directory.")
            BIN_EMPTY = True

        ## Step 4.4 Get modified vRhyme_best_bins acccording to both lytic and lysogenic and scaffold complete information
        vRhyme_best_bin_dir_modified = os.path.join(args['outdir_vrhyme'], 'vRhyme_best_bins_fasta_modified')
        prepare_output_dir(vRhyme_best_bin_dir_modified)
        vRhyme_best_bin_lytic_and_lysogenic_info = os.path.join(args['outdir_vrhyme'], 'vRhyme_best_bin_lytic_and_lysogenic_info.txt')
        vRhyme_best_bin_scaffold_complete_info = os.path.join(args['outdir_vrhyme'], 'vRhyme_best_bin_scaffold_complete_info.txt')
        # Rename .fasta to .fna
        if not BIN_EMPTY:
            with os.scandir(vRhyme_best_bin_dir) as it:
                for entry in it:
                    if entry.is_file() and entry.name.endswith('.fasta'):
                        new_path = entry.path[:-6] + '.fna'
                        os.rename(entry.path, new_path)

            ## Step 4.2 Get the lytic and lysogenic information for vRhyme_best_bins
            get_vRhyme_best_bin_lytic_and_lysogenic_info(vRhyme_best_bin_dir, args['outdir_vrhyme'], args['scf2lytic_or_lyso_summary'])
            ## Step 4.3 Get the scaffold complete information for vRhyme_best_bins
            vRhyme_best_bin_CheckV_result = os.path.join(args['outdir_vrhyme'], 'vRhyme_best_bins_fasta_CheckV_result')
            CheckV_cmd = f"conda run -n {args['conda_env_checkv']} python {os.path.join(args['root_dir'],'viwrap-scripts/run_CheckV.py')} {vRhyme_best_bin_dir} {vRhyme_best_bin_CheckV_result} {args['threads']} {args['CheckV_db']} >/dev/null 2>&1"
            os.system(CheckV_cmd)
            print(CheckV_cmd)
            CheckV_quality_summary = os.path.join(vRhyme_best_bin_CheckV_result, 'CheckV_quality_summary.txt')
            parse_checkv_result(vRhyme_best_bin_CheckV_result, CheckV_quality_summary)
            get_vRhyme_best_bin_scaffold_complete_info(CheckV_quality_summary, vRhyme_best_bin_scaffold_complete_info)
            os.system(f"rm -rf {vRhyme_best_bin_CheckV_result}")
            time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
            print(f"{time_current} | Run vRhyme to bin viral scaffolds. Finished")
            make_vRhyme_best_bins_fasta_modified(vRhyme_best_bin_dir,vRhyme_best_bin_dir_modified,vRhyme_best_bin_lytic_and_lysogenic_info,vRhyme_best_bin_scaffold_complete_info,args)
        else:
            with open(os.path.join(args['outdir_vrhyme'], 'vRhyme_best_bin_lytic_and_lysogenic_info.txt'), 'w') as fh:
                fh.write('vRhyme_bin\tpattern\tassignment\n')
            with open(os.path.join(args['outdir_vrhyme'], 'vRhyme_best_bin_scaffold_complete_info.txt'), 'w') as fh:
                fh.write('scaffold\tcomplete_info\n')

        ## Step 5.1 Make unbinned viral gn folder
        vRhyme_unbinned_viral_gn_dir = os.path.join(args['outdir_vrhyme'], 'vRhyme_unbinned_viral_gn_fasta')
        prepare_output_dir(vRhyme_unbinned_viral_gn_dir)
        ## Step 5.2 Make all vRhyme viral gn combined faa file
        make_unbinned_viral_gn(args['viral_scaffold'],vRhyme_best_bin_dir_modified,vRhyme_unbinned_viral_gn_dir,args)

        ## Step 5.3 Prepare pro2viral_gn map file
        pro2viral_gn_map = os.path.join(args['outdir_vrhyme'], 'pro2viral_gn_map.csv')
        get_pro2viral_gn_map(vRhyme_best_bin_dir_modified, vRhyme_unbinned_viral_gn_dir, pro2viral_gn_map)
        contig2viral_gn_map = os.path.join(args['outdir_vrhyme'], 'contig2viral_gn_map.csv')
        get_contig2viral_gn_map(vRhyme_best_bin_dir_modified, vRhyme_unbinned_viral_gn_dir, contig2viral_gn_map)

        all_vRhyme_faa = os.path.join(args['outdir_vrhyme'], 'all_vRhyme_faa.faa')
        all_vRhyme_faa = combine_all_vRhyme_faa(vRhyme_best_bin_dir_modified,vRhyme_unbinned_viral_gn_dir,args)

        both_dir         = os.path.join(args['outdir_vrhyme'], 'vRhyme_both_binned_unbinned')
        both_raw_dir     = os.path.join(args['outdir_vrhyme'], 'vRhyme_both_binned_unbinned_rawname')

        prepare_output_dir(both_dir)
        prepare_output_dir(both_raw_dir)

        for src_dir in [vRhyme_best_bin_dir_modified, vRhyme_unbinned_viral_gn_dir]:
            for fna in glob(os.path.join(src_dir, '*.fna')):
                shutil.copy2(fna, both_dir)
        pattern = re.compile(r'^>(.+?_vRhyme_(?:bin|unbinned)_\d+)---(.*)$')

        for fna in glob(os.path.join(both_dir, '*.fna')):
            out_file = os.path.join(both_raw_dir, os.path.basename(fna))
            with open(fna) as fin, open(out_file, 'w') as fout:
                for line in fin:
                    if line.startswith('>'):
                        m = pattern.match(line.rstrip('\n'))
                        if m:
                            # {sample}_vRhyme_bin_N---{cid} -> {sample}---{cid}
                            sample = re.sub(r'_vRhyme_(?:bin|unbinned)_\d+$', '', m.group(1))
                            new = f'>{sample}---{m.group(2)}\n'
                        else:
                            new = line
                        fout.write(new)
                    else:
                        fout.write(line)

        time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
        print(f"{time_current} | Run vRhyme to bin viral contigs from input contigs. Finished")

    else:
        print(f"{time_current} | Input viral scaffold file has no sequences, skipping all binning steps.")
        for placeholder in ('vRhyme_both_binned_unbinned','vRhyme_best_bins_fasta_modified',
                            'vRhyme_unbinned_viral_gn_fasta',
                            'vRhyme_coverage_files','vRhyme_both_binned_unbinned_rawname'
                            ):
            prepare_output_dir(os.path.join(args['outdir_vrhyme'], placeholder))

        with open(os.path.join(args['outdir_vrhyme'], 'vRhyme_best_bin_lytic_and_lysogenic_info.txt'), 'w') as fh:
            fh.write('vRhyme_bin\tpattern\tassignment\n')
        with open(os.path.join(args['outdir_vrhyme'], 'pro2viral_gn_map.csv'), 'w') as fh:
            fh.write('protein_id,contig_id\n')
        with open(os.path.join(args['outdir_vrhyme'], 'contig2viral_gn_map.csv'), 'w') as fh:
            fh.write('contig_newid,contig_id\n')
        with open(os.path.join(args['outdir_vrhyme'], 'vRhyme_best_bin_scaffold_complete_info.txt'), 'w') as fh:
            fh.write('scaffold\tcomplete_info\n')

        open(os.path.join(args['outdir_vrhyme'], 'all_vRhyme_faa.faa'), 'w').close()

        return

def cli():
    root_dir = os.path.dirname(os.path.abspath(__file__)) # The dir that contains the Metalink-virus script
    current_dir = os.getcwd() # The current dir
    db_path_default = os.path.join(current_dir, "Metalink_virus_db") # Set the default Metalink_virus_db dir under the current directory

    parser = argparse.ArgumentParser(usage = argparse.SUPPRESS, formatter_class=argparse.RawDescriptionHelpFormatter,
        description = """Analyzing wrapper for microbiome

Usage: python Metalink.04_vRhyme.py [options]

Task: Binning virus genomes from metagenomes
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
