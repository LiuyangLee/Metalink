#!/usr/bin/env python3

import sys
import os
import argparse
#import scripts
#from scripts import module
from datetime import datetime
from pathlib import Path
from glob import glob
import shutil
from Bio import SeqIO
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

def has_long_contig(input_fasta: str, length_limit) -> bool:
    try:
        length_limit = int(length_limit)
    except (ValueError, TypeError):
        raise ValueError(f"length_limit must be int, current: {length_limit!r}")
    with open(input_fasta, 'rt', encoding='utf-8', errors='ignore') as fh:
        for rec in SeqIO.parse(fh, 'fasta'):
            if len(rec) >= length_limit:
                return True
    return False

def make_empty_virus_outputs(out_dir, prefix='final_vs2_virus'):
    """
    Create 0-byte fasta + header-only annotation/summary for any virus caller.
    prefix examples:
        'final_vs2_virus'        -> vs2
        'final_genomad_virus'    -> genomad
        'final_vibrant_virus'    -> vibrant
        'final_virus'            -> genomad_vs_vb
    """
    ann_header = 'protein\tscaffold\tKO\tAMG\tKO name\tKO evalue\tKO score\tKO v-score\tPfam\tPfam name\tPfam evalue\tPfam score\tPfam v-score\tVOG\tVOG name\tVOG evalue\tVOG score\tVOG v-score\n'
    with open(os.path.join(out_dir, f'{prefix}.annotation.txt'), 'w') as fh:
        fh.write(ann_header)
    # header-only lytic/lysogenic summary
    #sum_header = 'scaffold\tlytic_or_lyso_or_integrated_prophage\tintegrase_presence_or_absence\n'
    sum_header = ''
    with open(os.path.join(out_dir, f'scf2lytic_or_lyso.summary.txt'), 'w') as fh:
        fh.write(sum_header)
    # 0-byte fasta
    for ext in ('fna', 'ffn', 'faa'):
        open(os.path.join(out_dir, f'{prefix}.{ext}'), 'w').close()

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

def get_keep2_mc_seq(virsorter_outdir, checkv_outdir, keep2_list_file, manual_check_list_file, keep2_fasta, manual_check_fasta):
    # Step 1 Store keep2_list, manual_check_list
    keep2_list = {} # seq => [length, score, hallmark, viral_gene, host_gene]
    with open(keep2_list_file, 'r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            if not line.startswith('#seq\t'):
                tmp = line.split('\t')
                seq, length, score, hallmark, viral_gene, host_gene = tmp[0], float(tmp[1]), float(tmp[2]), float(tmp[3]), float(tmp[4]), float(tmp[5])
                keep2_list[seq] = [length, score, hallmark, viral_gene, host_gene]
    lines.close()

    manual_check_list = {} # seq => [length, score, hallmark, viral_gene, host_gene]
    with open(manual_check_list_file, 'r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            if not line.startswith('#seq\t'):
                tmp = line.split('\t')
                seq, length, score, hallmark, viral_gene, host_gene = tmp[0], float(tmp[1]), float(tmp[2]), float(tmp[3]), float(tmp[4]), float(tmp[5])
                manual_check_list[seq] = [length, score, hallmark, viral_gene, host_gene]
    lines.close()

    # Step 2 Make keep2_fasta, manual_check_fasta
    all_seq = store_seq(os.path.join(checkv_outdir, 'combined.fna'))

    all_seq_keep2 = {}
    for header in all_seq:
        header_wo_array = header.replace('>', '', 1)
        if header_wo_array in keep2_list:
            all_seq_keep2[header] = all_seq[header]

    all_seq_manual_check = {}
    for header in all_seq:
        header_wo_array = header.replace('>', '', 1)
        if header_wo_array in manual_check_list:
            all_seq_manual_check[header] = all_seq[header]

    write_down_seq(all_seq_keep2, keep2_fasta)
    write_down_seq(all_seq_manual_check, manual_check_fasta)

def get_keep2_vb_passed_list(virsorter_outdir, keep2_vb_result, keep2_list_vb_passed_file):
    # Step 1 Store keep2_list
    keep2_list_file = os.path.join(virsorter_outdir, 'keep2_list.txt')
    keep2_list = {} # seq => [length, score, hallmark, viral_gene, host_gene]
    with open(keep2_list_file, 'r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            if not line.startswith('#seq\t'):
                tmp = line.split('\t')
                seq, length, score, hallmark, viral_gene, host_gene = tmp[0], float(tmp[1]), float(tmp[2]), float(tmp[3]), float(tmp[4]), float(tmp[5])
                keep2_list[seq] = [length, score, hallmark, viral_gene, host_gene]
    lines.close()
    keep2_list_vb_passed = {} # seq => [length, score, hallmark, viral_gene, host_gene]
    keep2_vb_result_seq = store_seq(keep2_vb_result)
    for header in keep2_vb_result_seq:
        header_wo_array = header.replace('>', '', 1)
        if '_fragment' in header_wo_array:
            header_wo_array = header_wo_array.rsplit('_fragment', 1)[0]
        keep2_list_vb_passed[header_wo_array] = keep2_list[header_wo_array]
    f = open(keep2_list_vb_passed_file, 'w')
    f.write('#seq\tlength\tscore\thallmark\tviral_gene\thost_gene' + '\n')
    for seq in keep2_list_vb_passed:
        f.write(seq + '\t' + '\t'.join(str(item) for item in keep2_list_vb_passed[seq]) + '\n')
    f.close()


def get_manual_check_vb_passed_list(virsorter_outdir, manual_check_vb_result, manual_check_list_vb_passed_file):
    # Step 1 Store manual_check_list
    manual_check_list_file = os.path.join(virsorter_outdir, 'manual_check_list.txt')

    manual_check_list = {} # seq => [length, score, hallmark, viral_gene, host_gene]
    with open(manual_check_list_file, 'r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            if not line.startswith('#seq\t'):
                tmp = line.split('\t')
                seq, length, score, hallmark, viral_gene, host_gene = tmp[0], float(tmp[1]), float(tmp[2]), float(tmp[3]), float(tmp[4]), float(tmp[5])
                manual_check_list[seq] = [length, score, hallmark, viral_gene, host_gene]
    lines.close()

    manual_check_list_vb_passed = {} # seq => [length, score, hallmark, viral_gene, host_gene]
    manual_check_vb_result_seq = store_seq(manual_check_vb_result)
    for header in manual_check_vb_result_seq:
        header_wo_array = header.replace('>', '', 1)
        if '_fragment' in header_wo_array:
            header_wo_array = header_wo_array.rsplit('_fragment', 1)[0]
        manual_check_list_vb_passed[header_wo_array] = manual_check_list[header_wo_array]

    f = open(manual_check_list_vb_passed_file, 'w')
    f.write('#seq\tlength\tscore\thallmark\tviral_gene\thost_gene' + '\n')
    for seq in manual_check_list_vb_passed:
        f.write(seq + '\t' + '\t'.join(str(item) for item in manual_check_list_vb_passed[seq]) + '\n')
    f.close()

def get_final_vs2_virus(virsorter_outdir, checkv_outdir, keep1_list_file, keep2_list_vb_passed_file, manual_check_list_vb_passed_file, final_vs2_virus_fasta_file):
    # Step 1 Store keep1_list,  keep2_list_vb_passed, manual_check_list_vb_passed
    keep1_list = {} # seq => [length, score, hallmark, viral_gene, host_gene]
    with open(keep1_list_file, 'r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            if not line.startswith('#seq\t'):
                tmp = line.split('\t')
                seq, length, score, hallmark, viral_gene, host_gene = tmp[0], float(tmp[1]), float(tmp[2]), float(tmp[3]), float(tmp[4]), float(tmp[5])
                keep1_list[seq] = [length, score, hallmark, viral_gene, host_gene]
    lines.close()
    keep2_list_vb_passed = {} # seq => [length, score, hallmark, viral_gene, host_gene]
    if os.path.exists(keep2_list_vb_passed_file):
        with open(keep2_list_vb_passed_file, 'r') as lines:
            for line in lines:
                line = line.rstrip('\n')
                if not line.startswith('#seq\t'):
                    tmp = line.split('\t')
                    seq, length, score, hallmark, viral_gene, host_gene = tmp[0], float(tmp[1]), float(tmp[2]), float(tmp[3]), float(tmp[4]), float(tmp[5])
                    keep2_list_vb_passed[seq] = [length, score, hallmark, viral_gene, host_gene]
        lines.close()
    manual_check_list_vb_passed = {} # seq => [length, score, hallmark, viral_gene, host_gene]
    if os.path.exists(manual_check_list_vb_passed_file):
        with open(manual_check_list_vb_passed_file, 'r') as lines:
            for line in lines:
                line = line.rstrip('\n')
                if not line.startswith('#seq\t'):
                    tmp = line.split('\t')
                    seq, length, score, hallmark, viral_gene, host_gene = tmp[0], float(tmp[1]), float(tmp[2]), float(tmp[3]), float(tmp[4]), float(tmp[5])
                    manual_check_list_vb_passed[seq] = [length, score, hallmark, viral_gene, host_gene]
        lines.close()

    # Step 2 Make final_vs2_virus.fna
    all_seq = store_seq(os.path.join(checkv_outdir, 'combined.fna'))

    all_seq_final = {}
    for header in all_seq:
        header_wo_array = header.replace('>', '', 1)
        if header_wo_array in keep1_list or header_wo_array in keep2_list_vb_passed or header_wo_array in manual_check_list_vb_passed:
            all_seq_final[header] = all_seq[header]
    write_down_seq(all_seq_final, final_vs2_virus_fasta_file)

def make_short_headers_fasta(input_fasta, output_fasta): # Remove spaces or '\t' in fasta headers and make a new fasta file
    input_fasta_seq = store_seq(input_fasta)

    output_fasta_seq = {}
    for header in input_fasta_seq:
        seq = input_fasta_seq[header]
        output_fasta_seq[header] = seq

    write_down_seq(output_fasta_seq, output_fasta)

def filter_fasta_by_length(input_fasta, output_fasta, length_limit): # Filter sequences in a FASTA file based on sequence length
    try:
        length_limit = int(length_limit)
    except (ValueError, TypeError):
        raise ValueError(f"length_limit must be int, current: {length_limit!r}")

    input_fasta_seq = store_seq(input_fasta)
    output_fasta_seq = {}
    for header in input_fasta_seq:
        seq = input_fasta_seq[header]
        if len(seq) >= int(length_limit):
            output_fasta_seq[header] = seq

    write_down_seq(output_fasta_seq, output_fasta)


def screen_virsorter2_result(virsorter_outdir, checkv_outdir, keep1_list_file, keep2_list_file, discard_list_file, manual_check_list_file):
    # ---------- Early exit if upstream results contain only headers ----------
    final_viral_score = os.path.join(virsorter_outdir, 'final-viral-score.tsv')
    quality_summary   = os.path.join(checkv_outdir, 'quality_summary.tsv')
    header = '#seq\tlength\tscore\thallmark\tviral_gene\thost_gene\n'

    # If either source file has ≤1 line (header only), export empty lists and return
    if (sum(1 for _ in open(final_viral_score)) <= 1 or
        sum(1 for _ in open(quality_summary)) <= 1):
        for f in (keep1_list_file, keep2_list_file,
                  discard_list_file, manual_check_list_file):
            with open(f, 'w') as fh:
                fh.write(header)
        return

    seq2info = defaultdict(list) # seq => [length, score, hallmark, viral_gene, host_gene]
    # Step 1 Parse final_viral_score file from VirSorter2 pass folder
    with open(final_viral_score, 'r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            if not line.startswith('seqname\t'):
                tmp = line.split('\t')
                seq, length, score, hallmark = tmp[0], float(tmp[5]), float(tmp[3]), float(tmp[6])
                seq2info[seq].append(length)
                seq2info[seq].append(score)
                seq2info[seq].append(hallmark)
    lines.close()

    # Step 2 Parse quality_summary file from VirSorter2 CheckV_result folder
    with open(quality_summary, 'r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            if not line.startswith('contig_id\t'):
                tmp = line.split('\t')
                seq, viral_gene, host_gene = tmp[0], float(tmp[5]), float(tmp[6])
                seq2info[seq].append(viral_gene)
                seq2info[seq].append(host_gene)
    lines.close()

    # Step 3 Make keep_list_file, discard_list_file, manual_check_list_file
    keep1_list = {} # seq => [length, score, hallmark, viral_gene, host_gene]
    keep2_list = {} # seq => [length, score, hallmark, viral_gene, host_gene]
    discard_list = {} # seq => [length, score, hallmark, viral_gene, host_gene]
    manual_check_list = {} # seq => [length, score, hallmark, viral_gene, host_gene]

    for seq in seq2info:
        length, score, hallmark, viral_gene, host_gene = seq2info[seq][0], seq2info[seq][1], seq2info[seq][2], seq2info[seq][3], seq2info[seq][4]
        if viral_gene > 0:
            keep1_list[seq] = [length, score, hallmark, viral_gene, host_gene]
        elif viral_gene == 0 and (score >= 0.95 or hallmark > 2 or host_gene == 0):
            keep2_list[seq] = [length, score, hallmark, viral_gene, host_gene]
        elif viral_gene == 0 and host_gene == 1 and length >= 10000:
            manual_check_list[seq] = [length, score, hallmark, viral_gene, host_gene]

    for seq in manual_check_list:
        if seq in keep1_list:
            del keep1_list[seq]
        if seq in keep2_list:
            del keep2_list[seq]

    for seq in seq2info:
        if seq not in keep1_list and seq not in keep2_list and seq not in manual_check_list:
            discard_list[seq] = seq2info[seq]

    # Replace the seq name after CheckV trimming
    cat_cmd = f"cat {checkv_outdir}/proviruses.fna {checkv_outdir}/viruses.fna > {checkv_outdir}/combined.fna"
    os.system(cat_cmd)
    all_seq = store_seq(os.path.join(checkv_outdir, 'combined.fna')) # all the seq after CheckV trimming
    seq_map = {} # seq1 (seq name before CheckV trimming) => seq2 (seq name after CheckV trimming; can contain "_{i}" at the end)
    for header in all_seq:
        seq2 = header.replace('>', '', 1)
        pattern = r'_\d+$' # Delete the _{i} pattern at the end of the scaffold
        seq1 = seq2
        if re.search(pattern, seq1):
            seq1 = re.sub(pattern, '', seq1)
        seq_map[seq1] = seq2

    f = open(keep1_list_file, 'w')
    f.write('#seq\tlength\tscore\thallmark\tviral_gene\thost_gene' + '\n')
    for seq in keep1_list:
        seq2 = seq_map[seq]
        f.write(seq2 + '\t' + '\t'.join(str(item) for item in keep1_list[seq]) + '\n')
    f.close()

    f = open(keep2_list_file, 'w')
    f.write('#seq\tlength\tscore\thallmark\tviral_gene\thost_gene' + '\n')
    for seq in keep2_list:
        seq2 = seq_map[seq]
        f.write(seq2 + '\t' + '\t'.join(str(item) for item in keep2_list[seq]) + '\n')
    f.close()

    f = open(discard_list_file, 'w')
    f.write('#seq\tlength\tscore\thallmark\tviral_gene\thost_gene' + '\n')
    for seq in discard_list:
        seq2 = seq_map[seq]
        f.write(seq2 + '\t' + '\t'.join(str(item) for item in discard_list[seq]) + '\n')
    f.close()

    f = open(manual_check_list_file, 'w')
    f.write('#seq\tlength\tscore\thallmark\tviral_gene\thost_gene' + '\n')
    for seq in manual_check_list:
        seq2 = seq_map[seq]
        f.write(seq2 + '\t' + '\t'.join(str(item) for item in manual_check_list[seq]) + '\n')
    f.close()

def parse_virsorter_lytic_and_lysogenic_info(virsorter_outdir):
    # Step 1 Get scf 2 lytic or lysogenic dict
    virsorter_annotation = os.path.join(virsorter_outdir, 'final_vs2_virus.annotation.txt')
    vs2_virus_file = os.path.join(virsorter_outdir, 'final_vs2_virus.fna')

    scf2lytic_or_lyso = {} # scf => [lytic_or_lyso_or_integrated_prophage, integrase_presence_or_absence]
    # lytic_or_lyso_or_integrated_prophage can contain: lytic_scaffold, integrated_prophage (parent scaffold), and lysogenic_scaffold
    # integrase_presence_or_absence can contain: integrase_present and integrase_absent
    vs2_virus_fasta_seq = store_seq(vs2_virus_file)

    scf2lytic_or_integrated_prophage = {} # scf => 'lytic_scaffold' or 'integrated_prophage (parent scaffold)'
    for header in vs2_virus_fasta_seq:
        scf = header.replace('>', '', 1)
        if '||full' in scf or '||lt2gene' in scf:
            scf2lytic_or_integrated_prophage[scf] = 'lytic_scaffold'
        elif re.search(r'\|\|\d+_partial', scf):
            parent_scaffold = scf.rsplit('||', 1)[0]
            scf2lytic_or_integrated_prophage[scf] = f"integrated_prophage ({parent_scaffold})"

    scf2integrase_presence = {} # scf => integrase_presence_or_absence
    integrase = ['VOG00041', 'VOG15133', 'VOG20969', 'VOG02658', 'VOG04024', 'VOG01778', 'VOG02371']
    with open(virsorter_annotation, 'r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            tmp = line.split('\t')
            if tmp[0] != 'protein':
                scf, VOG = tmp[1], tmp[13]
                scf2integrase_presence[scf] = 'integrase_absent'
                if VOG in integrase:
                    scf2integrase_presence[scf] = 'integrase_present'

    for scf in scf2integrase_presence:
        lytic_or_lyso_or_integrated_prophage = ''

        if scf in scf2lytic_or_integrated_prophage and 'integrated_prophage' in scf2lytic_or_integrated_prophage[scf]:
            lytic_or_lyso_or_integrated_prophage = scf2lytic_or_integrated_prophage[scf]
        elif scf in scf2lytic_or_integrated_prophage and scf2lytic_or_integrated_prophage[scf] == 'lytic_scaffold':
            if scf2integrase_presence[scf] == 'integrase_absent':
                lytic_or_lyso_or_integrated_prophage = 'lytic_scaffold'
            elif scf2integrase_presence[scf] == 'integrase_present':
                lytic_or_lyso_or_integrated_prophage = 'lysogenic_scaffold'

        integrase_presence_or_absence = scf2integrase_presence[scf]
        scf2lytic_or_lyso[scf] = [lytic_or_lyso_or_integrated_prophage, integrase_presence_or_absence]

    # Step 2 Write down scf2lytic_or_lyso result
    f = open(os.path.join(virsorter_outdir, 'scf2lytic_or_lyso.summary.txt'),'w')
    f.write('scaffold\tlytic_or_lyso_or_integrated_prophage\tintegrase_presence_or_absence\n')
    for scf in scf2lytic_or_lyso:
        line = scf + '\t' + '\t'.join(scf2lytic_or_lyso[scf]) + '\n'
        f.write(line)
    f.close()

def parse_genomad_lytic_and_lysogenic_info(genomad_outdir, metagenomic_scaffold_stem_name):
    # Step 1 Get scf 2 lytic or lysogenic dict
    virus_summary_file_addr = f"{genomad_outdir}/{metagenomic_scaffold_stem_name}_summary/{metagenomic_scaffold_stem_name}_virus_summary.tsv"
    integrase_result_file_addr = f"{genomad_outdir}/integrase_result_dir/{metagenomic_scaffold_stem_name}_virus_proteins_integrase_mmseqs2.tsv"

    scf2lytic_or_lyso = {} # scf => [lytic_or_lyso_or_integrated_prophage, integrase_presence_or_absence]
    # lytic_or_lyso_or_integrated_prophage can contain: lytic_scaffold, integrated_prophage (parent scaffold), and lysogenic_scaffold
    # integrase_presence_or_absence can contain: integrase_present and integrase_absent
    scf2topology = {} # scf => Provirus or other
    with open(virus_summary_file_addr, 'r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            if not line.startswith('seq_name'):
                scf, topology = line.split('\t')[0], line.split('\t')[2]
                if topology == 'Provirus':
                    scf2topology[scf] = 'Provirus'
                else:
                    scf2topology[scf] = 'other'

    scf2integrase_presence = {} # scf => integrase_present; only store scfs that have integrases
    with open(integrase_result_file_addr, 'r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            integrase_hit = line.split('\t')[0].split(' ', 1)[0]
            scf = integrase_hit.rsplit('_', 1)[0]
            scf2integrase_presence[scf] = 'integrase_present'

    for scf in scf2topology:
        lytic_or_lyso_or_integrated_prophage = ''
        if scf2topology[scf] == 'other' and scf not in scf2integrase_presence:
            lytic_or_lyso_or_integrated_prophage = 'lytic_scaffold'
        if scf2topology[scf] == 'other' and scf in scf2integrase_presence:
            lytic_or_lyso_or_integrated_prophage = 'lysogenic_scaffold'
        if scf2topology[scf] == 'Provirus':
            parent_scaffold = scf.split('|', 1)[0]
            lytic_or_lyso_or_integrated_prophage = f"integrated_prophage ({parent_scaffold})"
        integrase_presence_or_absence = ''
        if scf in scf2integrase_presence:
            integrase_presence_or_absence = 'integrase_present'
        else:
            integrase_presence_or_absence = 'integrase_absent'
        scf2lytic_or_lyso[scf] = [lytic_or_lyso_or_integrated_prophage, integrase_presence_or_absence]

    # Step 2 Write down scf2lytic_or_lyso result
    f = open(os.path.join(genomad_outdir, 'scf2lytic_or_lyso.summary.txt'),'w')
    f.write('scaffold\tlytic_or_lyso_or_integrated_prophage\tintegrase_presence_or_absence\n')
    for scf in scf2lytic_or_lyso:
        line = scf + '\t' + '\t'.join(scf2lytic_or_lyso[scf]) + '\n'
        f.write(line)
    f.close()



def fetch_arguments(parser,root_dir,db_path_default):
    parser.set_defaults(func=main)
    parser.set_defaults(program="run")
    parser.add_argument('--input_metagenome','-i', dest='input_metagenome', required=True, default='none', help=r'(required) input metagenome assembly. It can be a metagenome or entire virome assembly. The extension of the input nucleotide sequence should be ".fna"')
    parser.add_argument('--out_dir','-o', dest='out_dir', required=False, default='./Metalink_MetaG_virus_sample_outdir', help=r'(required) output directory to deposit all results (default = ./Metalink_MetaG_virus_sample_outdir) output folder to deposit all results.')
    parser.add_argument('--identify_method', dest='identify_method', required=False, default='genomad',help=r'(required) the virus identifying method to choose: genomad - Use geNomad; vs_vb - Use VirSorter2 and VIBRANT to get the overlapped viruses; genomad_vs_vb - Combine genomad and vs_vb virus to get the union results (default)')
    parser.add_argument('--threads','-t', dest='threads', required=False, default=10, help=r'number of threads (default = 10)')
    parser.add_argument('--conda_env_python', dest='conda_env_python', required=True, default='none', help=r'(required) path to the python conda environment used by viwrap helper scripts')
    parser.add_argument('--conda_env_seqkit', dest='conda_env_seqkit', required=True, default='none', help=r'(required) path to the sra conda environment used by viwrap helper scripts')
    parser.add_argument('--conda_env_checkv', dest='conda_env_checkv', required=True, default='none', help=r'(required) path to the checkm conda environment used by viwrap helper scripts')
    parser.add_argument('--conda_env_vibrant', dest='conda_env_vibrant', required=True, default='none', help=r'(required) path to the vibrant conda environment used by viwrap helper scripts')
    parser.add_argument('--conda_env_virsorter2', dest='conda_env_virsorter2', required=True, default='none', help=r'(required) path to the virsorter2 conda environment used by viwrap helper scripts')
    parser.add_argument('--conda_env_genomad', dest='conda_env_genomad', required=True, default='none', help=r'(required) path to the genomad conda environment used by viwrap helper scripts')
    parser.add_argument('--checkv_db', dest='checkv_db', required=True, default='none', help=r'(required) path to the CheckV database')
    parser.add_argument('--virsorter2_db', dest='virsorter2_db', required=True, default='none', help=r'(required) path to the VirSorter2 database')
    parser.add_argument('--genomad_db', dest='genomad_db', required=True, default='none', help=r'(required) path to the geNomad database')
    parser.add_argument('--vibrant_db', dest='vibrant_db', required=True, default='none', help=r'(required) path to the VIBRANT database')
    parser.add_argument('--root_dir', dest='root_dir', required=False, default=root_dir,help=argparse.SUPPRESS)
    parser.add_argument('--input_length_limit', dest='input_length_limit', required=False, default=2000, help=r'length in basepairs to limit input sequences (default=2000, can increase but not decrease); 2000 at least suggested for VIBRANT (vb)-based pipeline, 5000 at least suggested for VirSorter2 (vs)-based pipeline')
    parser.add_argument('--virome','-v', dest='virome', action='store_true', required=False, default=False, help=r"edit VIBRANT's sensitivity if the input dataset is a virome. It is suggested to use it if you know that the input assembly is virome or metagenome")

def set_defaults(args):
    ## Store databases
    args['VIBRANT_db'] = args['vibrant_db']
    args['VirSorter2_db'] = args['virsorter2_db']
    args['geNomad_db'] = args['genomad_db']
    args['CheckV_db'] = args['checkv_db']

    print(args['out_dir'])
    print(args['input_metagenome'])

    ## Store outdirs
    args['genomad_outdir'] = os.path.join(args['out_dir'],f"01_identification/{Path(args['input_metagenome']).stem}/genomad_result")
    args['virsorter_outdir'] = os.path.join(args['out_dir'],f"01_identification/{Path(args['input_metagenome']).stem}/vs_vb_result")
    args['genomad_vs_vb_outdir'] = os.path.join(args['out_dir'],f"01_identification/{Path(args['input_metagenome']).stem}/genomad_vs_vb_result")
    args['01_identification_outdir'] = os.path.join(args['out_dir'],f"01_identification/{Path(args['input_metagenome']).stem}")
    args['inner_vs_outdir'] = os.path.join(args['virsorter_outdir'],"vs_vb_result")
    print(args['genomad_vs_vb_outdir'])

def main(args):
    # Welcome and logger
    print("### Welcome to Metalink.master_run.01_virus_identify ###\n")
    ## Set the default args:
    set_defaults(args)
    if os.path.exists(args['out_dir']):
        #sys.exit(f"Error: The output directory '{args['out_dir']}' already exists. Please specify a different directory or remove the existing one.")
        pass
    else:
        # Create the output directory since it does not exist
        os.makedirs(args['out_dir'], exist_ok=True)

    os.makedirs(args['01_identification_outdir'], exist_ok=True)

    time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
    print(f"{time_current} | Check contig length. In processing...")

    if has_long_contig(args['input_metagenome'], args['input_length_limit']):
        print(f"{time_current} | Input metagenome of {args['input_metagenome']} contains contigs >= {args['input_length_limit']} bp.")
        # Step 2 Run VIBRANT or VirSorter2 or geNomad
        if args['identify_method'] == 'vs_vb':

            prepare_output_dir(args['virsorter_outdir'])
            time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
            print(f"{time_current} | Run VirSorter2 to identify viruses from input metagenome. Also plus CheckV to QC and trim, and KEGG, Pfam, and VOG HMMs to annotate viruses. In processing...")


            header_replaced_input_metagenome = os.path.join(args['virsorter_outdir'], f"{Path(args['input_metagenome']).stem}.fna")
            make_short_headers_fasta(args['input_metagenome'], header_replaced_input_metagenome)
            run_VirSorter2 = f"conda run -n {args['conda_env_virsorter2']} python {os.path.join(args['root_dir'],'viwrap-scripts/Metalink.run_VirSorter2.py')} {header_replaced_input_metagenome} {args['virsorter_outdir']} {args['threads']} {args['input_length_limit']}"
            print(run_VirSorter2)
            os.system(run_VirSorter2)
            os.system(f"rm {header_replaced_input_metagenome}")
            time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
            print(f"{time_current} | Run VirSorter2 to identify viruses from input metagenome. Finished")

            # *******  NEW: decide whether any virus was predicted  ********
            final_viral_combined = os.path.join(args['inner_vs_outdir'],'final-viral-combined.fa')
            VS_EMPTY = False
            if not os.path.exists(final_viral_combined) or os.path.getsize(final_viral_combined) == 0:
                VS_EMPTY = True
                print(f"{time_current} | VirSorter2 produced no viral sequences, skipping all VirSorter2-dependent steps.")

            # ----------  If empty -> create only header / 0-byte files  ----------
            if VS_EMPTY:
                make_empty_virus_outputs(args['virsorter_outdir'], prefix='final_vs2_virus')
            else:
            # ------------------  normal pipeline when viruses exist  ------------------
                # CheckV QC
                run_VirSorter2_CheckV = f"conda run -n {args['conda_env_checkv']} python {os.path.join(args['root_dir'],'viwrap-scripts/Metalink.run_VirSorter2_CheckV.py')} {args['virsorter_outdir']} {args['threads']} {args['CheckV_db']} >/dev/null 2>&1"
                print(run_VirSorter2_CheckV)
                os.system(run_VirSorter2_CheckV)
                time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
                print(f"{time_current} | Run CheckV to QC and trim viruses identified from VirSorter2. Finished")

                keep1_list_file = os.path.join(args['inner_vs_outdir'], 'keep1_list.txt')
                keep2_list_file = os.path.join(args['inner_vs_outdir'], 'keep2_list.txt')
                discard_list_file = os.path.join(args['inner_vs_outdir'], 'discard_list.txt')
                manual_check_list_file = os.path.join(args['inner_vs_outdir'], 'manual_check_list.txt')
                args['inner_checkv_outdir'] = os.path.join(args['virsorter_outdir'], 'CheckV_result')
                screen_virsorter2_result(args['inner_vs_outdir'], args['inner_checkv_outdir'], keep1_list_file, keep2_list_file, discard_list_file, manual_check_list_file)

                keep2_fasta = os.path.join(args['inner_vs_outdir'], 'keep2.fna')
                manual_check_fasta = os.path.join(args['inner_vs_outdir'], 'manual_check.fna')
                get_keep2_mc_seq(args['inner_vs_outdir'], args['inner_checkv_outdir'], keep2_list_file, manual_check_list_file, keep2_fasta, manual_check_fasta)

                if os.path.exists(keep2_fasta) and os.path.getsize(keep2_fasta) != 0:
                    keep2_fasta_VIBRANT_cmd = f"conda run -n {args['conda_env_vibrant']} python {os.path.join(args['root_dir'],'viwrap-scripts/run_VIBRANT.py')} {keep2_fasta} {args['inner_vs_outdir']} {args['threads']} {args['virome']} {args['input_length_limit']} {args['vibrant_db']}"
                    print(keep2_fasta_VIBRANT_cmd)
                    os.system(keep2_fasta_VIBRANT_cmd)
                    keep2_vb_result = os.path.join(args['inner_vs_outdir'], 'VIBRANT_keep2/VIBRANT_phages_keep2/keep2.phages_combined.fna')
                    keep2_list_vb_passed_file = os.path.join(args['inner_vs_outdir'], 'keep2_list_vb_passed.txt')
                    if not os.path.exists(keep2_vb_result) or os.path.getsize(keep2_vb_result) == 0:
                        open(keep2_list_vb_passed_file, 'a').close()
                    else:
                        get_keep2_vb_passed_list(args['inner_vs_outdir'], keep2_vb_result, keep2_list_vb_passed_file)
                    shutil.rmtree(os.path.join(args['inner_vs_outdir'], 'VIBRANT_keep2'), ignore_errors=True)
                if os.path.exists(manual_check_fasta) and os.path.getsize(manual_check_fasta) != 0:
                    manual_check_fasta_VIBRANT_cmd = f"conda run -n {args['conda_env_vibrant']} python {os.path.join(args['root_dir'],'viwrap-scripts/run_VIBRANT.py')} {manual_check_fasta} {args['inner_vs_outdir']} {args['threads']} {args['virome']} {args['input_length_limit']} {args['vibrant_db']}"
                    print(manual_check_fasta_VIBRANT_cmd)
                    os.system(manual_check_fasta_VIBRANT_cmd)
                    manual_check_vb_result = os.path.join(args['inner_vs_outdir'], 'VIBRANT_manual_check/VIBRANT_phages_manual_check/manual_check.phages_combined.fna')
                    manual_check_list_vb_passed_file = os.path.join(args['inner_vs_outdir'], 'manual_check_list_vb_passed.txt')
                    if not os.path.exists(manual_check_vb_result) or os.path.getsize(manual_check_vb_result) == 0:
                        open(manual_check_list_vb_passed_file, 'a').close()
                    else:
                        get_manual_check_vb_passed_list(args['inner_vs_outdir'], manual_check_vb_result, manual_check_list_vb_passed_file)
                    shutil.rmtree(os.path.join(args['inner_vs_outdir'], 'VIBRANT_manual_check'), ignore_errors=True)

                keep2_list_vb_passed_file = os.path.join(args['inner_vs_outdir'], 'keep2_list_vb_passed.txt')
                manual_check_list_vb_passed_file = os.path.join(args['inner_vs_outdir'], 'manual_check_list_vb_passed.txt')
                final_vs2_virus_fasta_file = os.path.join(args['virsorter_outdir'], 'final_vs2_virus.fna')
                get_final_vs2_virus(args['inner_vs_outdir'], args['inner_checkv_outdir'], keep1_list_file, keep2_list_vb_passed_file, manual_check_list_vb_passed_file, final_vs2_virus_fasta_file)


                FINAL_VS_EMPTY = False
                if not os.path.exists(final_vs2_virus_fasta_file) or os.path.getsize(final_vs2_virus_fasta_file) == 0:
                    FINAL_VS_EMPTY = True
                    print(f"{time_current} | VirSorter2 produced no viral sequences that passed VIBRANT manual check, skipping all VirSorter2-dependent steps.")
                # ----------  If empty -> create only header / 0-byte files  ----------
                if FINAL_VS_EMPTY:
                    make_empty_virus_outputs(args['virsorter_outdir'], prefix='final_vs2_virus')
                else:
                    args['tmp_annotate_dir'] = os.path.join(args['virsorter_outdir'],'tmp_annotate_dir')
                    run_annotate_by_VIBRANT_db = f"conda run -n {args['conda_env_vibrant']} python {os.path.join(args['root_dir'],'viwrap-scripts/Metalink.run_annotate_by_VIBRANT_db.py')} {args['VIBRANT_db']} vs_vb {args['virsorter_outdir']} {args['genomad_outdir']} {args['tmp_annotate_dir']} {args['threads']}"
                    print(run_annotate_by_VIBRANT_db)
                    os.system(run_annotate_by_VIBRANT_db)
                    parse_virsorter_lytic_and_lysogenic_info(args['virsorter_outdir'])

                    if os.path.exists(args['tmp_annotate_dir']):
                        shutil.rmtree(args['tmp_annotate_dir'])

                    time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
                    print(f"{time_current} | Use KEGG, Pfam, and VOG HMMs to annotate viruses. Finished")

                    # Embed parent-contig coords of VS2/CheckV fragments into seq IDs (|vs2region_S_E)
                    # run with the current interpreter (the script is stdlib-only): avoids
                    # conda-run lock stalls and buffered-output loss, and aborts loudly
                    embed_viral_coords_cmd = (
                        f"{sys.executable} "
                        f"{os.path.join(args['root_dir'],'viwrap-scripts/Metalink.embed_viral_coords.py')} "
                        f"--input_metagenome {args['input_metagenome']} "
                        f"--virsorter_outdir {args['virsorter_outdir']} --fail_loud")
                    print(embed_viral_coords_cmd)
                    if os.system(embed_viral_coords_cmd) != 0:
                        sys.exit('[embed_viral_coords] failed; aborting 01_identify to avoid inconsistent downstream results')

        elif args['identify_method'] == 'genomad':

            # Step 1 genomad part
            prepare_output_dir(args['genomad_outdir'])

            time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
            print(f"{time_current} | Run geNomad to identify and annotate viruses from input metagenome. In processing...")

            length_filtered_input_metagenome = os.path.join(args['genomad_outdir'], f"{Path(args['input_metagenome']).stem}.fna")
            filter_fasta_by_length(args['input_metagenome'], length_filtered_input_metagenome, args['input_length_limit'])

            run_geNomad = f"conda run -n {args['conda_env_genomad']} python {os.path.join(args['root_dir'],'viwrap-scripts/Metalink.run_geNomad.py')} {length_filtered_input_metagenome} {args['genomad_outdir']} {args['threads']} {args['genomad_db']}"
            print(run_geNomad)
            os.system(run_geNomad)

            default_genomad_outdir = os.path.join(args['genomad_outdir'], 'genomad_output')
            for name in os.listdir(default_genomad_outdir):
                os.rename(os.path.join(default_genomad_outdir, name),
                        os.path.join(args['genomad_outdir'], name))
            os.rmdir(default_genomad_outdir)
            os.system(f"rm {length_filtered_input_metagenome}")
            parse_genomad_lytic_and_lysogenic_info(args['genomad_outdir'], Path(args['input_metagenome']).stem)

            # Copy resulted fna and faa files
            genomad_virus_fna = f"{args['genomad_outdir']}/{Path(args['input_metagenome']).stem}_summary/{Path(args['input_metagenome']).stem}_virus.fna"
            genomad_virus_faa = f"{args['genomad_outdir']}/{Path(args['input_metagenome']).stem}_summary/{Path(args['input_metagenome']).stem}_virus_proteins.faa"
            os.system(f"cp {genomad_virus_fna} {args['genomad_outdir']}/final_genomad_virus.fna; cp {genomad_virus_faa} {args['genomad_outdir']}/final_genomad_virus.faa")

            FINAL_Genomad_EMPTY = False
            if not os.path.exists(genomad_virus_fna) or os.path.getsize(genomad_virus_fna) == 0:
                FINAL_Genomad_EMPTY = True
                print(f"{time_current} | geNomad produced no viral sequences, skipping all geNomad-dependent steps.")
            # ----------  If empty -> create only header / 0-byte files  ----------
            if FINAL_Genomad_EMPTY:
                make_empty_virus_outputs(args['genomad_outdir'], prefix='final_genomad_virus')
            else:
                args['tmp_annotate_dir'] = os.path.join(args['genomad_outdir'],'tmp_annotate_dir')
                run_annotate_by_VIBRANT_db = f"conda run -n {args['conda_env_vibrant']} python {os.path.join(args['root_dir'],'viwrap-scripts/Metalink.run_annotate_by_VIBRANT_db.py')} {args['VIBRANT_db']} genomad {args['inner_vs_outdir']} {args['genomad_outdir']} {args['tmp_annotate_dir']} {args['threads']}"
                print(run_annotate_by_VIBRANT_db)
                os.system(run_annotate_by_VIBRANT_db)

                if os.path.exists(args['tmp_annotate_dir']):
                    shutil.rmtree(args['tmp_annotate_dir'])

            time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
            print(f"{time_current} | Run geNomad to identify and annotate viruses from input metagenome. Finished")

        elif args['identify_method'] == 'genomad_vs_vb':
            prepare_output_dir(args['genomad_vs_vb_outdir'])

            # Step 1 genomad part
            prepare_output_dir(args['genomad_outdir'])

            time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
            print(f"{time_current} | Run geNomad to identify and annotate viruses from input metagenome. In processing...")


            length_filtered_input_metagenome = os.path.join(args['genomad_outdir'], f"{Path(args['input_metagenome']).stem}.fna")
            filter_fasta_by_length(args['input_metagenome'], length_filtered_input_metagenome, args['input_length_limit'])
            run_geNomad = f"conda run -n {args['conda_env_genomad']} python {os.path.join(args['root_dir'],'viwrap-scripts/Metalink.run_geNomad.py')} {length_filtered_input_metagenome} {args['genomad_outdir']} {args['threads']} {args['genomad_db']}"
            print(run_geNomad)
            os.system(run_geNomad)

            default_genomad_outdir = os.path.join(args['genomad_outdir'], 'genomad_output')
            for name in os.listdir(default_genomad_outdir):
                os.rename(os.path.join(default_genomad_outdir, name),
                        os.path.join(args['genomad_outdir'], name))
            os.rmdir(default_genomad_outdir)
            os.system(f"rm {length_filtered_input_metagenome}")
            parse_genomad_lytic_and_lysogenic_info(args['genomad_outdir'], Path(args['input_metagenome']).stem)

            # Copy resulted fna and faa files
            genomad_virus_fna = f"{args['genomad_outdir']}/{Path(args['input_metagenome']).stem}_summary/{Path(args['input_metagenome']).stem}_virus.fna"
            genomad_virus_faa = f"{args['genomad_outdir']}/{Path(args['input_metagenome']).stem}_summary/{Path(args['input_metagenome']).stem}_virus_proteins.faa"
            os.system(f"cp {genomad_virus_fna} {args['genomad_outdir']}/final_genomad_virus.fna; cp {genomad_virus_faa} {args['genomad_outdir']}/final_genomad_virus.faa")

            FINAL_Genomad_EMPTY = False
            if not os.path.exists(genomad_virus_fna) or os.path.getsize(genomad_virus_fna) == 0:
                FINAL_Genomad_EMPTY = True
                print(f"{time_current} | geNomad produced no viral sequences, skipping all geNomad-dependent steps.")
            # ----------  If empty -> create only header / 0-byte files  ----------
            if FINAL_Genomad_EMPTY:
                make_empty_virus_outputs(args['genomad_outdir'], prefix='final_genomad_virus')
            else:
                args['tmp_annotate_dir'] = os.path.join(args['genomad_outdir'],'tmp_annotate_dir')
                run_annotate_by_VIBRANT_db = f"conda run -n {args['conda_env_vibrant']} python {os.path.join(args['root_dir'],'viwrap-scripts/Metalink.run_annotate_by_VIBRANT_db.py')} {args['VIBRANT_db']} genomad {args['inner_vs_outdir']} {args['genomad_outdir']} {args['tmp_annotate_dir']} {args['threads']}"
                print(run_annotate_by_VIBRANT_db)
                os.system(run_annotate_by_VIBRANT_db)

                if os.path.exists(args['tmp_annotate_dir']):
                    shutil.rmtree(args['tmp_annotate_dir'])

            time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
            print(f"{time_current} | Run geNomad to identify and annotate viruses from input metagenome. Finished")

            # Step 2 VirSorter2 part verified by vibrant
            time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
            print(f"{time_current} | Run VirSorter2 to identify viruses from input metagenome. Also plus CheckV to QC and trim, and KEGG, Pfam, and VOG HMMs to annotate viruses. In processing...")

            prepare_output_dir(args['virsorter_outdir'])
            header_replaced_input_metagenome = os.path.join(args['virsorter_outdir'], f"{Path(args['input_metagenome']).stem}.fna")
            make_short_headers_fasta(args['input_metagenome'], header_replaced_input_metagenome)
            run_VirSorter2 = (f"conda run -n {args['conda_env_virsorter2']} python {os.path.join(args['root_dir'],'viwrap-scripts/Metalink.run_VirSorter2.py')} {header_replaced_input_metagenome} {args['virsorter_outdir']} {args['threads']} {args['input_length_limit']}")
            print(run_VirSorter2)
            os.system(run_VirSorter2)
            os.system(f"rm {header_replaced_input_metagenome}")
            time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
            print(f"{time_current} | Run VirSorter2 to identify viruses from input metagenome. Finished")

            # *******  NEW: decide whether any virus was predicted  ********
            final_viral_combined = os.path.join(args['inner_vs_outdir'],'final-viral-combined.fa')
            VS_EMPTY = False
            if not os.path.exists(final_viral_combined) or os.path.getsize(final_viral_combined) == 0:
                VS_EMPTY = True
                print(f"{time_current} | VirSorter2 produced no viral sequences, skipping all VirSorter2-dependent steps.")

            # ----------  If empty -> create only header / 0-byte files  ----------
            if VS_EMPTY:
                make_empty_virus_outputs(args['virsorter_outdir'], prefix='final_vs2_virus')
            else:
            # ------------------  normal pipeline when viruses exist  ------------------
                # CheckV QC
                run_VirSorter2_CheckV = (f"conda run -n {args['conda_env_checkv']} python {os.path.join(args['root_dir'],'viwrap-scripts/Metalink.run_VirSorter2_CheckV.py')} {args['virsorter_outdir']} {args['threads']} {args['CheckV_db']} >/dev/null 2>&1")
                print(run_VirSorter2_CheckV)
                os.system(run_VirSorter2_CheckV)
                time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
                print(f"{time_current} | Run CheckV to QC and trim viruses identified from VirSorter2. Finished")

                keep1_list_file = os.path.join(args['inner_vs_outdir'], 'keep1_list.txt')
                keep2_list_file = os.path.join(args['inner_vs_outdir'], 'keep2_list.txt')
                discard_list_file = os.path.join(args['inner_vs_outdir'], 'discard_list.txt')
                manual_check_list_file = os.path.join(args['inner_vs_outdir'], 'manual_check_list.txt')
                args['inner_checkv_outdir'] = os.path.join(args['virsorter_outdir'], 'CheckV_result')
                screen_virsorter2_result(args['inner_vs_outdir'], args['inner_checkv_outdir'], keep1_list_file, keep2_list_file, discard_list_file, manual_check_list_file)

                keep2_fasta = os.path.join(args['inner_vs_outdir'], 'keep2.fna')
                manual_check_fasta = os.path.join(args['inner_vs_outdir'], 'manual_check.fna')
                get_keep2_mc_seq(args['inner_vs_outdir'], args['inner_checkv_outdir'], keep2_list_file, manual_check_list_file, keep2_fasta, manual_check_fasta)

                if os.path.exists(keep2_fasta) and os.path.getsize(keep2_fasta) != 0:
                    keep2_fasta_VIBRANT_cmd = f"conda run -n {args['conda_env_vibrant']} python {os.path.join(args['root_dir'],'viwrap-scripts/run_VIBRANT.py')} {keep2_fasta} {args['inner_vs_outdir']} {args['threads']} {args['virome']} {args['input_length_limit']} {args['vibrant_db']}"
                    print(keep2_fasta_VIBRANT_cmd)
                    os.system(keep2_fasta_VIBRANT_cmd)
                    keep2_vb_result = os.path.join(args['inner_vs_outdir'], 'VIBRANT_keep2/VIBRANT_phages_keep2/keep2.phages_combined.fna')
                    keep2_list_vb_passed_file = os.path.join(args['inner_vs_outdir'], 'keep2_list_vb_passed.txt')
                    if not os.path.exists(keep2_vb_result) or os.path.getsize(keep2_vb_result) == 0:
                        open(keep2_list_vb_passed_file, 'a').close()
                    else:
                        get_keep2_vb_passed_list(args['inner_vs_outdir'], keep2_vb_result, keep2_list_vb_passed_file)
                    shutil.rmtree(os.path.join(args['inner_vs_outdir'], 'VIBRANT_keep2'), ignore_errors=True)
                if os.path.exists(manual_check_fasta) and os.path.getsize(manual_check_fasta) != 0:
                    manual_check_fasta_VIBRANT_cmd = f"conda run -n {args['conda_env_vibrant']} python {os.path.join(args['root_dir'],'viwrap-scripts/run_VIBRANT.py')} {manual_check_fasta} {args['inner_vs_outdir']} {args['threads']} {args['virome']} {args['input_length_limit']} {args['vibrant_db']}"
                    print(manual_check_fasta_VIBRANT_cmd)
                    os.system(manual_check_fasta_VIBRANT_cmd)
                    manual_check_vb_result = os.path.join(args['inner_vs_outdir'], 'VIBRANT_manual_check/VIBRANT_phages_manual_check/manual_check.phages_combined.fna')
                    manual_check_list_vb_passed_file = os.path.join(args['inner_vs_outdir'], 'manual_check_list_vb_passed.txt')
                    if not os.path.exists(manual_check_vb_result) or os.path.getsize(manual_check_vb_result) == 0:
                        open(manual_check_list_vb_passed_file, 'a').close()
                    else:
                        get_manual_check_vb_passed_list(args['inner_vs_outdir'], manual_check_vb_result, manual_check_list_vb_passed_file)
                    shutil.rmtree(os.path.join(args['inner_vs_outdir'], 'VIBRANT_manual_check'), ignore_errors=True)

                keep2_list_vb_passed_file = os.path.join(args['inner_vs_outdir'], 'keep2_list_vb_passed.txt')
                manual_check_list_vb_passed_file = os.path.join(args['inner_vs_outdir'], 'manual_check_list_vb_passed.txt')
                final_vs2_virus_fasta_file = os.path.join(args['virsorter_outdir'], 'final_vs2_virus.fna')
                get_final_vs2_virus(args['inner_vs_outdir'], args['inner_checkv_outdir'], keep1_list_file, keep2_list_vb_passed_file, manual_check_list_vb_passed_file, final_vs2_virus_fasta_file)

                FINAL_VS_EMPTY = False
                if not os.path.exists(final_vs2_virus_fasta_file) or os.path.getsize(final_vs2_virus_fasta_file) == 0:
                    FINAL_VS_EMPTY = True
                    print(f"{time_current} | VirSorter2 produced no viral sequences that passed VIBRANT manual check, skipping all VirSorter2-dependent steps.")

                # ----------  If empty -> create only header / 0-byte files  ----------
                if FINAL_VS_EMPTY:
                    make_empty_virus_outputs(args['virsorter_outdir'], prefix='final_vs2_virus')
                else:
                    args['tmp_annotate_dir'] = os.path.join(args['virsorter_outdir'],'tmp_annotate_dir')
                    run_annotate_by_VIBRANT_db = f"conda run -n {args['conda_env_vibrant']} python {os.path.join(args['root_dir'],'viwrap-scripts/Metalink.run_annotate_by_VIBRANT_db.py')} {args['VIBRANT_db']} vs_vb {args['virsorter_outdir']} {args['genomad_outdir']} {args['tmp_annotate_dir']} {args['threads']}"
                    print(run_annotate_by_VIBRANT_db)
                    os.system(run_annotate_by_VIBRANT_db)

                    parse_virsorter_lytic_and_lysogenic_info(args['virsorter_outdir'])

                    if os.path.exists(args['tmp_annotate_dir']):
                        shutil.rmtree(args['tmp_annotate_dir'])

                    time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
                    print(f"{time_current} | Use KEGG, Pfam, and VOG HMMs to annotate viruses. Finished")

                    # Embed parent-contig coords of VS2/CheckV fragments into seq IDs (|vs2region_S_E)
                    # run with the current interpreter (the script is stdlib-only): avoids
                    # conda-run lock stalls and buffered-output loss, and aborts loudly
                    embed_viral_coords_cmd = (
                        f"{sys.executable} "
                        f"{os.path.join(args['root_dir'],'viwrap-scripts/Metalink.embed_viral_coords.py')} "
                        f"--input_metagenome {args['input_metagenome']} "
                        f"--virsorter_outdir {args['virsorter_outdir']} --fail_loud")
                    print(embed_viral_coords_cmd)
                    if os.system(embed_viral_coords_cmd) != 0:
                        sys.exit('[embed_viral_coords] failed; aborting 01_identify to avoid inconsistent downstream results')


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
            combine_genomad_vs_vb_cmd = f"conda run -n {args['conda_env_python']} python {os.path.join(args['root_dir'],'viwrap-scripts/Metalink.combine_genomad_vs_vb.py')} {args_str}"
            os.system(combine_genomad_vs_vb_cmd)
            print(combine_genomad_vs_vb_cmd)

            # seqkit command to extract the added virus from virsorter2
            seqkit_input = os.path.join(args['genomad_vs_vb_outdir'], 'list.virsorter2_vibrant.tsv')
            vs2_fasta = os.path.join(args['virsorter_outdir'], 'final_vs2_virus.fna')
            vs2_output = os.path.join(args['virsorter_outdir'], 'add.virsorter2.fna')
            seqkit_cmd = f"conda run -n {args['conda_env_seqkit']} seqkit grep -f {seqkit_input} {vs2_fasta} > {vs2_output}"
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

        else:
            sys.exit(f"Please make sure your input for --identify_method option is one of these: \"genomad\", \"vs_vb\", and \"genomad_vs_vb\"; you can also omit this in the command line, the default is \"genomad_vs_vb\"")
    else:
        print(
            f"Warning: {time_current} | No contig >= {args['input_length_limit']} bp in {args['input_metagenome']}; "
            "virus identification will yield empty results."
        )
        if args['identify_method'] == 'vs_vb':
            prepare_output_dir(args['virsorter_outdir'])
            make_empty_virus_outputs(args['virsorter_outdir'], prefix='final_vs2_virus')
        elif args['identify_method'] == 'genomad':
            prepare_output_dir(args['genomad_outdir'])
            make_empty_virus_outputs(args['genomad_outdir'], prefix='final_genomad_virus')
        elif args['identify_method'] == 'genomad_vs_vb':
            prepare_output_dir(args['genomad_vs_vb_outdir'])
            make_empty_virus_outputs(args['genomad_vs_vb_outdir'], prefix='final_virus')


def cli():
    root_dir = os.path.dirname(os.path.abspath(__file__)) # The dir that contains the Metalink-virus script
    current_dir = os.getcwd() # The current dir
    db_path_default = os.path.join(current_dir, "Metalink_virus_db") # Set the default Metalink_virus_db dir under the current directory

    parser = argparse.ArgumentParser(usage = argparse.SUPPRESS, formatter_class=argparse.RawDescriptionHelpFormatter,
        description = """Analyzing wrapper for microbiome

Usage: python Metalink.01_identify.py [options]

Task: Identify virus genomes from metagenomes
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
