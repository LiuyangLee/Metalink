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

def clean_contig(s: str, A: str) -> str:
    # bin-member ids are {A}---{cid} (A = {sample}_vRhyme_bin_N / _vRhyme_unbinned_N);
    # restore the original contig id {sample}---{cid}
    sample = re.sub(r'_vRhyme_(bin|unbinned)_\d+$', '', A)
    s = re.sub(rf'^{re.escape(A)}---', f'{sample}---', s)
    s = re.sub(rf',{re.escape(A)}---', f',{sample}---', s)
    return s

def merge_faa4_into_annotation(annotation_tsv: str, faa_file: str):
    pat = re.compile(r'^>([^\s#]+)\s.*#\s*(\d+)\s*#\s*(\d+)\s*#\s*([-+]?1)')
    faa4 = {}                       # seqid -> (start, end, strand)
    with open(faa_file) as fh:
        for line in fh:
            if line.startswith('>'):
                m = pat.search(line)
                if m:
                    seqid, s, e, st = m.groups()
                    faa4[seqid] = (s, e, st)
    tmp = annotation_tsv + '.tmp'
    with open(annotation_tsv) as fin, open(tmp, 'w') as fout:
        header = fin.readline().rstrip('\n')
        fout.write(f'{header}\tORF start\tORF end\tORF strand\n')
        for line in fin:
            row = line.rstrip('\n').split('\t')
            genome = row[0]
            for idx in (1, 2):
                if idx < len(row):
                    row[idx] = clean_contig(row[idx], genome)
            seqid = row[1]
            s, e, st = faa4.get(seqid, ('', '', ''))
            fout.write('\t'.join(row) + f'\t{s}\t{e}\t{st}\n')
    os.replace(tmp, annotation_tsv)

def combine_iphop_results(args, combined_host_pred_to_genome_result, combined_host_pred_to_genus_result):
    host_pred_to_genome_header = ''
    host_pred_to_genome_m90_result = {} # Virus => [Virus, Host_genome, Host_taxonomy, Main_method, Confidence_score, Additional_methods]
    host_pred_to_genus_header = ''
    host_pred_to_genus_m90_result = {} # Virus => [Virus, AAI_to_closest_RaFAH_reference, Host_genus, Confidence_score, List_of_methods]
    if args['custom_MAGs_dir'] == 'none':
        host_pred_to_genome_m90 = os.path.join(args['iphop_outdir'], "Host_prediction_to_genome_m90.csv")
        host_pred_to_genus_m90 = os.path.join(args['iphop_outdir'], "Host_prediction_to_genus_m90.csv")

        host_pred_EMPTY = False
        if not os.path.exists(host_pred_to_genome_m90) or os.path.getsize(host_pred_to_genome_m90) == 0:
            host_pred_EMPTY = True

        if not host_pred_EMPTY:
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

        host_pred_EMPTY = False
        if not os.path.exists(host_pred_to_genome_m90_custom) or os.path.getsize(host_pred_to_genome_m90_custom) == 0:
            host_pred_EMPTY = True

        if not host_pred_EMPTY:
            with open(host_pred_to_genome_m90_custom, 'r') as lines:
                for line in lines:
                    line = line.strip('\n')
                    if not line.startswith('Virus,'):
                        tmp = line.split(',')
                        Virus, Confidence_score = tmp[0], tmp[4]
                        if Virus not in host_pred_to_genome_m90_result:
                            host_pred_to_genome_m90_result[Virus] = tmp
                        elif float(Confidence_score) > float(host_pred_to_genome_m90_result[Virus][4]):
                            host_pred_to_genome_m90_result[Virus] = tmp

            with open(host_pred_to_genus_m90_custom, 'r') as lines:
                for line in lines:
                    line = line.strip('\n')
                    if not line.startswith('Virus,'):
                        tmp = line.split(',')
                        Virus, Confidence_score = tmp[0], tmp[3]
                        if Virus not in host_pred_to_genus_m90_result:
                            host_pred_to_genus_m90_result[Virus] = tmp
                        elif float(Confidence_score) > float(host_pred_to_genus_m90_result[Virus][3]):
                            host_pred_to_genus_m90_result[Virus] = tmp

    if not host_pred_EMPTY:
        f = open(combined_host_pred_to_genome_result, 'w')
        f.write(host_pred_to_genome_header + '\n')
        for virus in sorted(host_pred_to_genome_m90_result):
            f.write(','.join(host_pred_to_genome_m90_result[virus]) + '\n')
        f.close()

        f = open(combined_host_pred_to_genus_result , 'w')
        f.write(host_pred_to_genus_header + '\n')
        for virus in sorted(host_pred_to_genus_m90_result):
            f.write(','.join(host_pred_to_genus_m90_result[virus]) + '\n')
        f.close()
    else:
        f = open(combined_host_pred_to_genome_result, 'w')
        f.write(host_pred_to_genome_header + '\n')
        f.close()

        f = open(combined_host_pred_to_genus_result , 'w')
        f.write(host_pred_to_genus_header + '\n')
        f.close()

def get_virus_genome_annotation_result(args):
    if 'vs' in args['identify_method']:
        # Step 1 Store annotation result
        annotation_result_file = ''
        if args['identify_method'] == 'vs_vb':
            annotation_result_file = os.path.join(args['virsorter_outdir'], 'final_vs2_virus.annotation.txt')
        elif args['identify_method'] == 'genomad_vs_vb':
            annotation_result_file = os.path.join(args['genomad_vs_vb_outdir'], 'final_virus.annotation.txt')

        annotation_result = {} # protein => [items in each line]
        annotation_result_header = ''
        with open (annotation_result_file, 'r') as lines:
            for line in lines:
                line = line.rstrip('\n')
                if line.startswith('protein\t'):
                    annotation_result_header = line
                else:
                    tmp = line.split('\t')
                    tmp = [item.strip('"') if item.startswith('"') and item.endswith('"') else item for item in tmp]
                    protein = tmp[0]
                    annotation_result[protein] = tmp
        lines.close()

        # Step 2 Store gn2long_proteins and long_protein2gn dict
        gn2long_proteins = {} # gn => [long_proteins]
        long_protein2gn = {} # long_protein => gn

        all_faa_addrs = glob(os.path.join(args['Metalink_virus_summary_outdir'],'Virus_genomes_files/*.faa'))
        for faa_addr in all_faa_addrs:
            gn = Path(faa_addr).stem
            faa_seqs = store_seq(faa_addr)
            long_proteins = [x.replace('>', '', 1) for x in faa_seqs]
            gn2long_proteins[gn] = long_proteins

            for long_protein in long_proteins:
                long_protein2gn[long_protein] = gn

        # Step 3 Get new VIBRANT annotation result
        annotation_result_new = {} # long_protein => [items in each line]
        for long_protein in long_protein2gn:
            #print(long_protein)
            # {gn}---{cid_ord} -> {sample}---{cid_ord} (annotation keys keep the
            # original sample-prefixed contig ids)
            gn_sample = re.sub(r'_vRhyme_(bin|unbinned)_\d+$', '', long_protein2gn[long_protein])
            protein = gn_sample + '---' + long_protein.split('---', 1)[1]
            #print(protein)
            items = annotation_result[protein]
            #print(items)
            items[0] = long_protein
            items[1] = long_protein.rsplit('_', 1)[0]
            #print(items[1])
            items.insert(0, long_protein2gn[long_protein])
            if len(items[4]) > 0 and items[4] == 'AMG': # Replace AMG to AMG (unfiltered)
                items[4] = 'AMG (unfiltered)'
            annotation_result_new[long_protein] = items

        # Step 4 Write down the result
        result = os.path.join(args['Metalink_virus_summary_outdir'],'Virus_annotation_results.tsv')
        f = open(result, 'w')
        f.write('viral genome\t' + annotation_result_header + '\n')
        for long_protein in annotation_result_new:
            line = '\t'.join(annotation_result_new[long_protein])
            f.write(line + '\n')
        f.close()

    elif args['identify_method'] == 'genomad':
        # Step 1 Store annotation result
        ## Store the VIBRANT db annotating result
        annotation_result_file1 = os.path.join(args['genomad_outdir'], 'final_genomad_virus.annotation.txt')
        annotation_result = {} # protein => [items in each line]
        annotation_result_header = ''
        with open (annotation_result_file1, 'r') as lines:
            for line in lines:
                line = line.rstrip('\n')
                if line.startswith('protein\t'):
                    annotation_result_header = line
                else:
                    tmp = line.split('\t')
                    tmp = [item.strip('"') if item.startswith('"') and item.endswith('"') else item for item in tmp]
                    protein = tmp[0]
                    annotation_result[protein] = tmp
        lines.close()

        ## Store the geNomad-self annotating result
        annotation_result_file2 = os.path.join(args['genomad_outdir'], f"{Path(args['input_metagenome']).stem}_annotate/{Path(args['input_metagenome']).stem}_genes.tsv")
        annotation_result_header = annotation_result_header + "\tgeNomad marker\tgeNomad evalue\tgeNomad score\tgeNomad marker annotation\tgeNomad marker annotation description"
        with open (annotation_result_file2, 'r') as lines:
            for line in lines:
                line = line.rstrip('\n')
                if not line.startswith('gene'):
                    tmp = line.split('\t')
                    protein, geNomad_marker, geNomad_evalue, geNomad_score, geNomad_marker_annotation, geNomad_marker_annotation_description = tmp[0], tmp[8], tmp[9], tmp[10], tmp[18], tmp[19]
                    if protein in annotation_result:
                        annotation_result[protein].extend([geNomad_marker, geNomad_evalue, geNomad_score, geNomad_marker_annotation, geNomad_marker_annotation_description])
                    else:
                        annotation_result[protein] = ['' for _ in range(18)] + [geNomad_marker, geNomad_evalue, geNomad_score, geNomad_marker_annotation, geNomad_marker_annotation_description]

        # Step 2 Store gn2long_proteins and long_protein2gn dict
        gn2long_proteins = {} # gn => [long_proteins]
        long_protein2gn = {} # long_protein => gn

        all_faa_addrs = glob(os.path.join(args['Metalink_virus_summary_outdir'],'Virus_genomes_files/*.faa'))
        for faa_addr in all_faa_addrs:
            gn = Path(faa_addr).stem
            faa_seqs = store_seq(faa_addr)
            long_proteins = [x.replace('>', '', 1) for x in faa_seqs]
            gn2long_proteins[gn] = long_proteins

            for long_protein in long_proteins:
                long_protein2gn[long_protein] = gn

        # Step 3 Get new VIBRANT annotation result
        annotation_result_new = {} # long_protein => [items in each line]
        for long_protein in long_protein2gn:
            # {gn}---{cid_ord} -> {sample}---{cid_ord} (annotation keys keep the
            # original sample-prefixed contig ids)
            gn_sample = re.sub(r'_vRhyme_(bin|unbinned)_\d+$', '', long_protein2gn[long_protein])
            protein = gn_sample + '---' + long_protein.split('---', 1)[1]
            items = annotation_result[protein]
            items[0] = long_protein
            items[1] = long_protein.rsplit('_', 1)[0]
            items.insert(0, long_protein2gn[long_protein])
            if len(items[4]) > 0 and items[4] == 'AMG': # Replace AMG to AMG (unfiltered)
                items[4] = 'AMG (unfiltered)'
            annotation_result_new[long_protein] = items

        # Step 4 Write down the result
        result = os.path.join(args['Metalink_virus_summary_outdir'],'Virus_annotation_results.tsv')
        f = open(result, 'w')
        f.write('viral genome\t' + annotation_result_header + '\n')
        for long_protein in annotation_result_new:
            line = '\t'.join(annotation_result_new[long_protein])
            f.write(line + '\n')
        f.close()

def check_left(AMG, pros_ordered, AMG_list):
    AMG_index = pros_ordered.index(AMG) # The index of the given AMG
    logic = True
    for i in range(0, AMG_index):
        if pros_ordered[i] not in AMG_list:
            logic = False
    return logic

def check_right(AMG, pros_ordered, AMG_list):
    AMG_index = pros_ordered.index(AMG) # The index of the given AMG
    logic = True
    for i in range((AMG_index + 1), len(pros_ordered)):
        if pros_ordered[i] not in AMG_list:
            logic = False
    return logic

def get_amg_pro_info(AMG_dir, virus_annotation_result_file, viral_gn_dir, VIBRANT_db):
    amg_pro2info = {} # amg_pro => [long_scf, ko, ko_name, metabolisms, pathways]

    # Step 1 Store the metabolism and pathway kos
    metabolism2kos = {} # metabolism => [kos]
    pathway2kos = {} # pathway => [kos]
    KEGG_pathway_file = os.path.join(VIBRANT_db, 'files/VIBRANT_KEGG_pathways_summary.tsv')
    with open(KEGG_pathway_file,  'r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            tmp = line.split('\t')
            if 'map' in tmp[0]:
                metabolism, pathway, ko_arrary = tmp[1], tmp[2], tmp[3]
                if metabolism[0] == '"' and metabolism[-1] == '"':
                    metabolism = metabolism.strip('"').rstrip('"')
                if pathway[0] == '"' and pathway[-1] == '"':
                    pathway = pathway.strip('"').rstrip('"')
                kos = ko_arrary.split('~')
                metabolism2kos[metabolism] = kos
                pathway2kos[pathway] = kos
    lines.close()

    all_amg_kos = [] # Store all the amg kos
    VIBRANT_AMGs_file = os.path.join(VIBRANT_db, 'files/VIBRANT_AMGs.tsv')
    with open(VIBRANT_AMGs_file,  'r') as lines: # 'U' to deal with either dos or unix file format
        for line in lines:
            line = line.rstrip('\n')
            if line != 'KO' and line.startswith('K'):
                all_amg_kos.append(line)
    lines.close()

    # Step 2 Store ko2metablisms and ko2pathways dicts
    ko2metablisms = {} # ko => metabolisms
    ko2pathways = {} # ko => pathways
    for ko in all_amg_kos:
        metabolism_hits = set()
        pathway_hits = set()
        for metabolism in metabolism2kos:
            if ko in metabolism2kos[metabolism]:
                metabolism_hits.add(metabolism)
        for pathway in pathway2kos:
            if ko in pathway2kos[pathway]:
                pathway_hits.add(pathway)
        ko2metablisms[ko] = ' | '.join(list(metabolism_hits))
        ko2pathways[ko] = ' | '.join(list(pathway_hits))

    # Step 3 Parse virus_annotation_result to get useful information for AMGs
    pro2v_scores = {} # pro => [KEGG_v_score, Pfam_v_score, VOG_v_score]
    with open(virus_annotation_result_file, 'r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            if not line.startswith('viral genome'):
                tmp = line.split('\t')
                if tmp[4] == 'AMG (unfiltered)':
                    amg_pro, long_scf, ko, ko_name = tmp[1], tmp[2], tmp[3], tmp[5]
                    if ko_name[0] == '"' and ko_name[-1] == '"':
                        ko_name = ko_name.strip('"').rstrip('"')
                    metabolisms = ko2metablisms[ko]
                    pathways = ko2pathways[ko]
                    amg_pro2info[amg_pro] = [long_scf, ko, ko_name, metabolisms, pathways]
                pro, KEGG_v_score, Pfam_v_score, VOG_v_score = tmp[1], tmp[8], tmp[13], tmp[18]
                pro2v_scores[pro] = [KEGG_v_score, Pfam_v_score, VOG_v_score]
    lines.close()

    # Step 4 Filter tail AMGs (AMGs that are located in either ends of a given scaffold)
    ## Step 4.1 Store the scf2amg_list dict
    scf2amg_list = defaultdict(list)
    for amg in amg_pro2info.keys():
        scf = amg.rsplit('_', 1)[0]
        scf2amg_list[scf].append(amg)

    ## Step 4.2 Store scaffold to protein list dict
    scf2pros_ordered = {} # scf => [pros]; The proteins are re-ordered
    scf2pros = defaultdict(list) # scf => [pros]

    all_virus_protein_seq = {} # The seq dict for all virus proteins
    all_viral_faa_addrs = glob(f"{viral_gn_dir}/*.faa")
    for each_viral_faa_addr in all_viral_faa_addrs:
        each_viral_faa_seq = store_seq(each_viral_faa_addr)
        all_virus_protein_seq.update(each_viral_faa_seq)

    headers = [header.replace('>', '', 1) for header in all_virus_protein_seq]
    for header in headers:
        scf = header.rsplit('_', 1)[0]
        pro = header
        scf2pros[scf].append(pro)

    for scf in scf2pros:
        pros = scf2pros[scf] # pros is now a list
        # Make dict from [pros] as: pro => order_num
        pro2order_num = {pro:int(pro.rsplit('_', 1)[1]) for pro in pros}

        pros_ordered = list(sorted(pro2order_num, key = pro2order_num.get)) # pros_ordered is now a re-ordered list
        scf2pros_ordered[scf] = pros_ordered

    ## Step 4.3 Get the tail AMG label results
    amg2label = {} # AMG => [position, tail_AMG (or not_tail_AMG)]; position can be '1 in 11'
    for scf in scf2amg_list:
        pros_ordered = scf2pros_ordered[scf]
        amg_list = scf2amg_list[scf]

        ### First, start from left to right
        for amg in amg_list:
            amg_index = pros_ordered.index(amg) # The index of the given AMG
            if amg_index == 0:
                amg2label[amg] = [f'{(amg_index + 1)} in {len(pros_ordered)}', 'tail_AMG']
            if amg_index > 0:
                if check_left(amg, pros_ordered, amg_list): # If the left genes are all AMGs
                    amg2label[amg] = [f'{(amg_index + 1)} in {len(pros_ordered)}', 'tail_AMG']

        ### Second, start from right to left
        for amg in amg_list:
            amg_index = pros_ordered.index(amg) # The index of the given AMG
            if amg_index == (len(pros_ordered) - 1) :
                amg2label[amg] = [f'{(amg_index + 1)} in {len(pros_ordered)}', 'tail_AMG']
            if amg_index < (len(pros_ordered) - 1):
                if check_right(amg, pros_ordered, amg_list): # If the right genes are all AMGs
                    amg2label[amg] = [f'{(amg_index + 1)} in {len(pros_ordered)}', 'tail_AMG']

        ### Label the rest as 'not_tail_AMG'
        for amg in amg_list:
            amg_index = pros_ordered.index(amg) # The index of the given AMG
            if amg not in amg2label:
                amg2label[amg] = [f'{(amg_index + 1)} in {len(pros_ordered)}', 'not_tail_AMG']

    ## Step 4.4 Get tail_amg
    tail_amg = set()
    for amg in amg2label:
        if amg2label[amg][1] == 'tail_AMG':
            tail_amg.add(amg)

    # Step 5 Filter AMGs that have any v-scores (KEGG and Pfam v-scores) >= 1
    amgs_w_high_v_scores = [] # Store the list of AMGs that have any v-scores (KEGG and Pfam v-scores) >= 1
    for amg in amg_pro2info.keys():
        KEGG_v_score, Pfam_v_score = pro2v_scores[amg][0], pro2v_scores[amg][1]
        logic = True
        if KEGG_v_score and float(KEGG_v_score) >= 1:
            logic = False
        if Pfam_v_score and float(Pfam_v_score) >= 1:
            logic = False
        if not logic:
            amgs_w_high_v_scores.append(amg)

    # Step 6 Filter AMGs with flanking genes of v-scores < 0.25
    ## Step 6.1 Get AMG to flanking genes dict
    flanking_num = 2 # The number of flanking gene on both sides to be taken into consideration
    amg2flanking_genes = {} # amg => [flanking genes]

    for scf in scf2amg_list:
        pros_ordered = scf2pros_ordered[scf]
        amg_list = scf2amg_list[scf]

        for amg in amg_list:
            if amg2label[amg][1] == 'not_tail_AMG': # Only focusing not_tail_AMG
                flanking_genes_left = []
                flanking_genes_right = []
                amg_index = pros_ordered.index(amg) # The index of the given AMG
                # Find left side flanking genes
                for i in range((amg_index - 1), -1, -1):
                    if pros_ordered[i] not in amg_list:
                        flanking_genes_left.append(pros_ordered[i])
                        if len(flanking_genes_left) == flanking_num:
                            break
                # Find right side flanking genes:
                for i in range((amg_index + 1), len(pros_ordered), 1):
                    if pros_ordered[i] not in amg_list:
                        flanking_genes_right.append(pros_ordered[i])
                        if len(flanking_genes_right) == flanking_num:
                            break
                flanking_genes = flanking_genes_left + flanking_genes_right
                amg2flanking_genes[amg] = flanking_genes

    ## Step 6.2 Filter AMGs with flanking genes of v-scores < 0.25
    amgs_w_flanking_genes_of_low_v_score = []
    for amg in amg2flanking_genes:
        flanking_genes = amg2flanking_genes[amg]
        low_v_score_gene_count = 0
        for gene in flanking_genes:
            if gene in pro2v_scores:
                KEGG_v_score = pro2v_scores[gene][0]
                if (KEGG_v_score and float(KEGG_v_score) < 0.25):
                    low_v_score_gene_count = low_v_score_gene_count + 1

        if low_v_score_gene_count == len(flanking_genes): # If all the flanking genes are low v-score genes
            amgs_w_flanking_genes_of_low_v_score.append(amg)

    # Step 7 Filter AMGs based on COG category
    ## The number of KOs that are corresponding to 'T' and 'B' COG categories: 598
    amgs_belong_to_not_correct_COG = []
    KOs_of_T_or_B =  ['K07782', 'K03688', 'K07655', 'K20976', 'K10682', 'K11521', 'K00575', 'K21828', 'K06269', 'K19704', 'K20977', 'K14981', 'K18352', 'K04409', 'K03407', 'K03741', 'K02768', 'K07734', 'K03114', 'K20263', 'K02661', 'K07648', 'K08968', 'K05527', 'K02688', 'K07649', 'K06413', 'K07177', 'K20972', 'K07715', 'K11240', 'K07665', 'K10961', 'K11183', 'K19734', 'K01139', 'K20973', 'K11230', 'K07180', 'K11913', 'K21023', 'K21558', 'K00870', 'K01090', 'K02204', 'K11638', 'K02821', 'K17611', 'K14983', 'K04464', 'K02490', 'K07774', 'K18701', 'K19077', 'K10036', 'K21090', 'K06041', 'K04441', 'K11524', 'K08482', 'K07775', 'K01420', 'K02798', 'K07680', 'K02476', 'K11522', 'K11691', 'K20339', 'K14979', 'K07695', 'K07700', 'K12771', 'K20540', 'K16326', 'K04460', 'K14989', 'K16511', 'K08082', 'K21024', 'K07165', 'K15836', 'K09158', 'K18941', 'K11383', 'K19694', 'K07067', 'K02475', 'K07651', 'K14414', 'K11917', 'K00888', 'K18350', 'K07670', 'K05792', 'K20340', 'K15562', 'K18219', 'K20201', 'K07768', 'K20249', 'K04752', 'K11354', 'K14758', 'K02487', 'K07169', 'K18143', 'K17612', 'K11232', 'K03410', 'K01525', 'K12294', 'K18679', 'K08866', 'K18073', 'K00916', 'K07179', 'K07131', 'K07770', 'K10943', 'K08083', 'K07659', 'K07689', 'K21025', 'K03414', 'K03415', 'K07773', 'K07772', 'K11231', 'K18099', 'K02478', 'K11330', 'K03651', 'K02216', 'K06889', 'K14023', 'K03408', 'K06379', 'K18940', 'K09975', 'K11640', 'K11211', 'K19801', 'K11226', 'K02486', 'K11615', 'K05877', 'K13486', 'K13599', 'K13040', 'K18444', 'K02488', 'K19705', 'K05795', 'K14064', 'K11916', 'K04382', 'K07645', 'K07712', 'K10914', 'K07769', 'K18672', 'K07710', 'K01524', 'K20954', 'K19812', 'K02484', 'K11357', 'K11912', 'K02178', 'K06149', 'K07641', 'K07683', 'K05875', 'K20112', 'K00924', 'K07203', 'K03803', 'K07183', 'K20252', 'K20334', 'K07639', 'K02647', 'K12766', 'K19575', 'K15012', 'K00088', 'K11355', 'K20253', 'K07685', 'K05873', 'K07097', 'K07643', 'K14949', 'K07176', 'K11356', 'K10715', 'K20248', 'K04749', 'K11692', 'K03776', 'K06382', 'K02590', 'K07200', 'K14980', 'K21020', 'K15861', 'K21562', 'K11894', 'K18987', 'K07702', 'K21696', 'K21088', 'K14575', 'K14803', 'K20956', 'K02483', 'K07637', 'K18144', 'K20966', 'K00990', 'K07706', 'K20958', 'K21560', 'K11233', 'K19661', 'K07662', 'K18349', 'K11623', 'K07638', 'K02481', 'K10697', 'K11520', 'K11523', 'K08475', 'K02479', 'K06595', 'K14978', 'K05791', 'K02850', 'K07684', 'K16957', 'K10039', 'K11633', 'K08296', 'K07654', 'K08884', 'K15474', 'K13490', 'K07716', 'K15852', 'K06023', 'K02658', 'K21564', 'K14065', 'K20975', 'K07714', 'K11624', 'K07709', 'K02482', 'K06639', 'K02784', 'K07166', 'K12776', 'K08293', 'K13069', 'K07686', 'K21021', 'K03607', 'K07644', 'K04757', 'K11641', 'K12266', 'K10001', 'K18044', 'K13527', 'K11932', 'K21556', 'K07813', 'K11526', 'K18680', 'K07678', 'K03598', 'K12765', 'K00951', 'K17073', 'K11229', 'K11629', 'K07675', 'K07167', 'K17762', 'K07718', 'K08485', 'K21085', 'K15427', 'K18669', 'K13589', 'K07168', 'K10942', 'K19623', 'K19641', 'K11915', 'K06597', 'K07640', 'K07126', 'K09969', 'K18045', 'K07694', 'K20957', 'K06276', 'K08484', 'K11620', 'K07212', 'K02202', 'K19733', 'K11444', 'K07636', 'K03413', 'K01768', 'K21884', 'K14571', 'K06714', 'K02218', 'K13924', 'K19692', 'K09996', 'K08476', 'K21158', 'K18967', 'K07697', 'K02677', 'K02667', 'K13060', 'K07778', 'K05770', 'K13061', 'K18072', 'K13525', 'K08479', 'K13487', 'K01994', 'K18304', 'K01356', 'K14051', 'K19617', 'K04371', 'K01104', 'K06027', 'K07673', 'K11198', 'K08294', 'K09773', 'K16314', 'K19732', 'K07660', 'K20330', 'K13593', 'K07650', 'K17060', 'K06378', 'K11634', 'K07699', 'K10005', 'K03969', 'K13339', 'K02206', 'K05876', 'K13491', 'K02030', 'K13488', 'K07705', 'K07690', 'K07692', 'K07668', 'K20959', 'K19621', 'K18968', 'K07693', 'K07661', 'K21086', 'K18046', 'K06867', 'K10022', 'K18344', 'K07216', 'K19833', 'K09137', 'K11712', 'K03563', 'K18043', 'K21157', 'K18670', 'K02489', 'K05518', 'K19690', 'K10125', 'K19622', 'K07711', 'K07814', 'K02773', 'K19666', 'K07720', 'K07704', 'K10018', 'K11914', 'K20918', 'K07663', 'K11443', 'K11189', 'K06217', 'K07681', 'K04751', 'K11384', 'K11908', 'K21685', 'K07664', 'K19693', 'K07181', 'K07182', 'K19830', 'K12146', 'K20250', 'K05874', 'K06641', 'K07667', 'K19081', 'K11630', 'K03597', 'K07154', 'K07674', 'K06207', 'K07708', 'K04563', 'K15423', 'K21561', 'K11617', 'K11328', 'K19731', 'K07717', 'K11525', 'K13591', 'K07178', 'K02589', 'K20488', 'K14982', 'K07777', 'K11228', 'K11332', 'K13598', 'K15011', 'K13590', 'K11201', 'K14055', 'K11750', 'K14061', 'K17752', 'K18351', 'K21397', 'K07658', 'K20264', 'K07172', 'K10014', 'K18842', 'K10912', 'K07816', 'K07653', 'K13489', 'K04348', 'K21217', 'K07671', 'K02659', 'K20945', 'K03411', 'K07313', 'K13243', 'K07688', 'K07666', 'K21019', 'K02806', 'K07669', 'K21084', 'K01697', 'K11637', 'K02584', 'K17610', 'K07687', 'K17508', 'K04095', 'K20487', 'K07771', 'K07698', 'K03974', 'K19609', 'K03412', 'K12767', 'K13587', 'K11618', 'K10941', 'K19078', 'K02491', 'K09997', 'K17763', 'K08269', 'K07656', 'K18765', 'K13245', 'K19800', 'K20960', 'K20978', 'K13041', 'K07707', 'K02657', 'K14987', 'K10013', 'K07739', 'K03503', 'K07314', 'K17061', 'K19292', 'K13642', 'K06200', 'K03406', 'K05805', 'K12295', 'K14986', 'K21901', 'K04768', 'K07652', 'K07701', 'K07672', 'K02480', 'K04767', 'K02660', 'K11227', 'K18345', 'K07703', 'K18098', 'K02424', 'K06268', 'K07315', 'K18096', 'K19082', 'K04333', 'K06958', 'K20974', 'K06596', 'K07646', 'K06206', 'K19708', 'K21555', 'K07781', 'K19852', 'K13338', 'K20955', 'K07713', 'K19806', 'K10126', 'K00906', 'K09933', 'K07657', 'K20961', 'K03666', 'K20962', 'K18691', 'K02214', 'K18986', 'K07776', 'K02668', 'K20074', 'K11184', 'K21603', 'K19610', 'K20963', 'K13246', 'K19691', 'K03973', 'K02831', 'K19616', 'K02515', 'K07173', 'K16956', 'K10681', 'K20919', 'K11614', 'K07642', 'K00914', 'K07696', 'K02467', 'K07153', 'K21022', 'K11329', 'K21563', 'K16961', 'K13584', 'K02477', 'K13815', 'K20965', 'K15475', 'K02485']
    for amg in amg_pro2info.keys():
        if amg_pro2info[amg][1] in KOs_of_T_or_B:
            amgs_belong_to_not_correct_COG.append(amg)

    # Step 8 Get flitered AMGs
    amg_pro_filtered2info = {} # amg_pro => [long_scf, ko, ko_name, metabolisms, pathways]
    for amg in amg_pro2info:
        if amg not in tail_amg and amg not in amgs_w_high_v_scores and amg not in amgs_w_flanking_genes_of_low_v_score and amg not in amgs_belong_to_not_correct_COG:
            amg_pro_filtered2info[amg] = amg_pro2info[amg]

    return amg_pro_filtered2info

def get_amg_pro_info_for_wo_reads(AMG_dir, virus_annotation_result_file, Metalink_virus_summary_outdir, VIBRANT_db):
    amg_pro2info = {} # amg_pro => [scf, ko, ko_name, metabolisms, pathways]

    # Step 1 Store the metabolism and pathway kos
    metabolism2kos = {} # metabolism => [kos]
    pathway2kos = {} # pathway => [kos]
    KEGG_pathway_file = os.path.join(VIBRANT_db, 'files/VIBRANT_KEGG_pathways_summary.tsv')
    with open(KEGG_pathway_file,  'r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            tmp = line.split('\t')
            if 'map' in tmp[0]:
                metabolism, pathway, ko_arrary = tmp[1], tmp[2], tmp[3]
                if metabolism[0] == '"' and metabolism[-1] == '"':
                    metabolism = metabolism.strip('"').rstrip('"')
                if pathway[0] == '"' and pathway[-1] == '"':
                    pathway = pathway.strip('"').rstrip('"')
                kos = ko_arrary.split('~')
                metabolism2kos[metabolism] = kos
                pathway2kos[pathway] = kos
    lines.close()

    all_amg_kos = [] # Store all the amg kos
    VIBRANT_AMGs_file = os.path.join(VIBRANT_db, 'files/VIBRANT_AMGs.tsv')
    with open(VIBRANT_AMGs_file,  'r') as lines: # 'U' to deal with either dos or unix file format
        for line in lines:
            line = line.rstrip('\n')
            if line != 'KO' and line.startswith('K'):
                all_amg_kos.append(line)
    lines.close()

    # Step 2 Store ko2metablisms and ko2pathways dicts
    ko2metablisms = {} # ko => metabolisms
    ko2pathways = {} # ko => pathways
    for ko in all_amg_kos:
        metabolism_hits = set()
        pathway_hits = set()
        for metabolism in metabolism2kos:
            if ko in metabolism2kos[metabolism]:
                metabolism_hits.add(metabolism)
        for pathway in pathway2kos:
            if ko in pathway2kos[pathway]:
                pathway_hits.add(pathway)
        ko2metablisms[ko] = ' | '.join(list(metabolism_hits))
        ko2pathways[ko] = ' | '.join(list(pathway_hits))

    # Step 3 Parse virus_annotation_result to get useful information for AMGs
    pro2v_scores = {} # pro => [KEGG_v_score, Pfam_v_score, VOG_v_score]
    with open(virus_annotation_result_file, 'r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            if not line.startswith('protein'):
                tmp = line.split('\t')
                if tmp[3] == 'AMG (unfiltered)':
                    amg_pro, scf, ko, ko_name = tmp[0], tmp[1], tmp[2], tmp[4]
                    if ko_name[0] == '"' and ko_name[-1] == '"':
                        ko_name = ko_name.strip('"').rstrip('"')
                    metabolisms = ko2metablisms[ko]
                    pathways = ko2pathways[ko]
                    amg_pro2info[amg_pro] = [scf, ko, ko_name, metabolisms, pathways]
                pro, KEGG_v_score, Pfam_v_score, VOG_v_score = tmp[0], tmp[7], tmp[12], tmp[17]
                pro2v_scores[pro] = [KEGG_v_score, Pfam_v_score, VOG_v_score]
    lines.close()

    # Step 4 Filter tail AMGs (AMGs that are located in either ends of a given scaffold)
    ## Step 4.1 Store the scf2amg_list dict
    scf2amg_list = defaultdict(list)
    for amg in amg_pro2info.keys():
        scf = amg.rsplit('_', 1)[0]
        scf2amg_list[scf].append(amg)

    ## Step 4.2 Store scaffold to protein list dict
    scf2pros_ordered = {} # scf => [pros]; The proteins are re-ordered
    scf2pros = defaultdict(list) # scf => [pros]
    all_virus_protein_seq = store_seq(os.path.join(Metalink_virus_summary_outdir,'final_virus.faa')) # The seq dict for all virus proteins

    headers = [header.replace('>', '', 1) for header in all_virus_protein_seq]
    for header in headers:
        scf = header.rsplit('_', 1)[0]
        pro = header
        scf2pros[scf].append(pro)

    for scf in scf2pros:
        pros = scf2pros[scf] # pros is now a list
        # Make dict from [pros] as: pro => order_num
        pro2order_num = {pro:int(pro.rsplit('_', 1)[1]) for pro in pros}

        pros_ordered = list(sorted(pro2order_num, key = pro2order_num.get)) # pros_ordered is now a re-ordered list
        scf2pros_ordered[scf] = pros_ordered

    ## Step 4.3 Get the tail AMG label results
    amg2label = {} # AMG => [position, tail_AMG (or not_tail_AMG)]; position can be '1 in 11'
    for scf in scf2amg_list:
        pros_ordered = scf2pros_ordered[scf]
        amg_list = scf2amg_list[scf]

        ### First, start from left to right
        for amg in amg_list:
            amg_index = pros_ordered.index(amg) # The index of the given AMG
            if amg_index == 0:
                amg2label[amg] = [f'{(amg_index + 1)} in {len(pros_ordered)}', 'tail_AMG']
            if amg_index > 0:
                if check_left(amg, pros_ordered, amg_list): # If the left genes are all AMGs
                    amg2label[amg] = [f'{(amg_index + 1)} in {len(pros_ordered)}', 'tail_AMG']

        ### Second, start from right to left
        for amg in amg_list:
            amg_index = pros_ordered.index(amg) # The index of the given AMG
            if amg_index == (len(pros_ordered) - 1) :
                amg2label[amg] = [f'{(amg_index + 1)} in {len(pros_ordered)}', 'tail_AMG']
            if amg_index < (len(pros_ordered) - 1):
                if check_right(amg, pros_ordered, amg_list): # If the right genes are all AMGs
                    amg2label[amg] = [f'{(amg_index + 1)} in {len(pros_ordered)}', 'tail_AMG']

        ### Label the rest as 'not_tail_AMG'
        for amg in amg_list:
            amg_index = pros_ordered.index(amg) # The index of the given AMG
            if amg not in amg2label:
                amg2label[amg] = [f'{(amg_index + 1)} in {len(pros_ordered)}', 'not_tail_AMG']

    ## Step 4.4 Get tail_amg
    tail_amg = set()
    for amg in amg2label:
        if amg2label[amg][1] == 'tail_AMG':
            tail_amg.add(amg)

    # Step 5 Filter AMGs that have any v-scores (KEGG and Pfam v-scores) >= 1
    amgs_w_high_v_scores = [] # Store the list of AMGs that have any v-scores (KEGG and Pfam v-scores) >= 1
    for amg in amg_pro2info.keys():
        KEGG_v_score, Pfam_v_score = pro2v_scores[amg][0], pro2v_scores[amg][1]
        logic = True
        if KEGG_v_score and float(KEGG_v_score) >= 1:
            logic = False
        if Pfam_v_score and float(Pfam_v_score) >= 1:
            logic = False
        if not logic:
            amgs_w_high_v_scores.append(amg)

    # Step 6 Filter AMGs with flanking genes of v-scores < 0.25
    ## Step 6.1 Get AMG to flanking genes dict
    flanking_num = 2 # The number of flanking gene on both sides to be taken into consideration
    amg2flanking_genes = {} # amg => [flanking genes]

    for scf in scf2amg_list:
        pros_ordered = scf2pros_ordered[scf]
        amg_list = scf2amg_list[scf]

        for amg in amg_list:
            if amg2label[amg][1] == 'not_tail_AMG': # Only focusing not_tail_AMG
                flanking_genes_left = []
                flanking_genes_right = []
                amg_index = pros_ordered.index(amg) # The index of the given AMG
                # Find left side flanking genes
                for i in range((amg_index - 1), -1, -1):
                    if pros_ordered[i] not in amg_list:
                        flanking_genes_left.append(pros_ordered[i])
                        if len(flanking_genes_left) == flanking_num:
                            break
                # Find right side flanking genes:
                for i in range((amg_index + 1), len(pros_ordered), 1):
                    if pros_ordered[i] not in amg_list:
                        flanking_genes_right.append(pros_ordered[i])
                        if len(flanking_genes_right) == flanking_num:
                            break
                flanking_genes = flanking_genes_left + flanking_genes_right
                amg2flanking_genes[amg] = flanking_genes

    ## Step 6.2 Filter AMGs with flanking genes of v-scores < 0.25
    amgs_w_flanking_genes_of_low_v_score = []
    for amg in amg2flanking_genes:
        flanking_genes = amg2flanking_genes[amg]
        low_v_score_gene_count = 0
        for gene in flanking_genes:
            if gene in pro2v_scores:
                KEGG_v_score = pro2v_scores[gene][0]
                if (KEGG_v_score and float(KEGG_v_score) < 0.25):
                    low_v_score_gene_count = low_v_score_gene_count + 1

        if low_v_score_gene_count == len(flanking_genes): # If all the flanking genes are low v-score genes
            amgs_w_flanking_genes_of_low_v_score.append(amg)

    # Step 7 Filter AMGs based on COG category
    ## The number of KOs that are corresponding to 'T' and 'B' COG categories: 598
    amgs_belong_to_not_correct_COG = []
    KOs_of_T_or_B =  ['K07782', 'K03688', 'K07655', 'K20976', 'K10682', 'K11521', 'K00575', 'K21828', 'K06269', 'K19704', 'K20977', 'K14981', 'K18352', 'K04409', 'K03407', 'K03741', 'K02768', 'K07734', 'K03114', 'K20263', 'K02661', 'K07648', 'K08968', 'K05527', 'K02688', 'K07649', 'K06413', 'K07177', 'K20972', 'K07715', 'K11240', 'K07665', 'K10961', 'K11183', 'K19734', 'K01139', 'K20973', 'K11230', 'K07180', 'K11913', 'K21023', 'K21558', 'K00870', 'K01090', 'K02204', 'K11638', 'K02821', 'K17611', 'K14983', 'K04464', 'K02490', 'K07774', 'K18701', 'K19077', 'K10036', 'K21090', 'K06041', 'K04441', 'K11524', 'K08482', 'K07775', 'K01420', 'K02798', 'K07680', 'K02476', 'K11522', 'K11691', 'K20339', 'K14979', 'K07695', 'K07700', 'K12771', 'K20540', 'K16326', 'K04460', 'K14989', 'K16511', 'K08082', 'K21024', 'K07165', 'K15836', 'K09158', 'K18941', 'K11383', 'K19694', 'K07067', 'K02475', 'K07651', 'K14414', 'K11917', 'K00888', 'K18350', 'K07670', 'K05792', 'K20340', 'K15562', 'K18219', 'K20201', 'K07768', 'K20249', 'K04752', 'K11354', 'K14758', 'K02487', 'K07169', 'K18143', 'K17612', 'K11232', 'K03410', 'K01525', 'K12294', 'K18679', 'K08866', 'K18073', 'K00916', 'K07179', 'K07131', 'K07770', 'K10943', 'K08083', 'K07659', 'K07689', 'K21025', 'K03414', 'K03415', 'K07773', 'K07772', 'K11231', 'K18099', 'K02478', 'K11330', 'K03651', 'K02216', 'K06889', 'K14023', 'K03408', 'K06379', 'K18940', 'K09975', 'K11640', 'K11211', 'K19801', 'K11226', 'K02486', 'K11615', 'K05877', 'K13486', 'K13599', 'K13040', 'K18444', 'K02488', 'K19705', 'K05795', 'K14064', 'K11916', 'K04382', 'K07645', 'K07712', 'K10914', 'K07769', 'K18672', 'K07710', 'K01524', 'K20954', 'K19812', 'K02484', 'K11357', 'K11912', 'K02178', 'K06149', 'K07641', 'K07683', 'K05875', 'K20112', 'K00924', 'K07203', 'K03803', 'K07183', 'K20252', 'K20334', 'K07639', 'K02647', 'K12766', 'K19575', 'K15012', 'K00088', 'K11355', 'K20253', 'K07685', 'K05873', 'K07097', 'K07643', 'K14949', 'K07176', 'K11356', 'K10715', 'K20248', 'K04749', 'K11692', 'K03776', 'K06382', 'K02590', 'K07200', 'K14980', 'K21020', 'K15861', 'K21562', 'K11894', 'K18987', 'K07702', 'K21696', 'K21088', 'K14575', 'K14803', 'K20956', 'K02483', 'K07637', 'K18144', 'K20966', 'K00990', 'K07706', 'K20958', 'K21560', 'K11233', 'K19661', 'K07662', 'K18349', 'K11623', 'K07638', 'K02481', 'K10697', 'K11520', 'K11523', 'K08475', 'K02479', 'K06595', 'K14978', 'K05791', 'K02850', 'K07684', 'K16957', 'K10039', 'K11633', 'K08296', 'K07654', 'K08884', 'K15474', 'K13490', 'K07716', 'K15852', 'K06023', 'K02658', 'K21564', 'K14065', 'K20975', 'K07714', 'K11624', 'K07709', 'K02482', 'K06639', 'K02784', 'K07166', 'K12776', 'K08293', 'K13069', 'K07686', 'K21021', 'K03607', 'K07644', 'K04757', 'K11641', 'K12266', 'K10001', 'K18044', 'K13527', 'K11932', 'K21556', 'K07813', 'K11526', 'K18680', 'K07678', 'K03598', 'K12765', 'K00951', 'K17073', 'K11229', 'K11629', 'K07675', 'K07167', 'K17762', 'K07718', 'K08485', 'K21085', 'K15427', 'K18669', 'K13589', 'K07168', 'K10942', 'K19623', 'K19641', 'K11915', 'K06597', 'K07640', 'K07126', 'K09969', 'K18045', 'K07694', 'K20957', 'K06276', 'K08484', 'K11620', 'K07212', 'K02202', 'K19733', 'K11444', 'K07636', 'K03413', 'K01768', 'K21884', 'K14571', 'K06714', 'K02218', 'K13924', 'K19692', 'K09996', 'K08476', 'K21158', 'K18967', 'K07697', 'K02677', 'K02667', 'K13060', 'K07778', 'K05770', 'K13061', 'K18072', 'K13525', 'K08479', 'K13487', 'K01994', 'K18304', 'K01356', 'K14051', 'K19617', 'K04371', 'K01104', 'K06027', 'K07673', 'K11198', 'K08294', 'K09773', 'K16314', 'K19732', 'K07660', 'K20330', 'K13593', 'K07650', 'K17060', 'K06378', 'K11634', 'K07699', 'K10005', 'K03969', 'K13339', 'K02206', 'K05876', 'K13491', 'K02030', 'K13488', 'K07705', 'K07690', 'K07692', 'K07668', 'K20959', 'K19621', 'K18968', 'K07693', 'K07661', 'K21086', 'K18046', 'K06867', 'K10022', 'K18344', 'K07216', 'K19833', 'K09137', 'K11712', 'K03563', 'K18043', 'K21157', 'K18670', 'K02489', 'K05518', 'K19690', 'K10125', 'K19622', 'K07711', 'K07814', 'K02773', 'K19666', 'K07720', 'K07704', 'K10018', 'K11914', 'K20918', 'K07663', 'K11443', 'K11189', 'K06217', 'K07681', 'K04751', 'K11384', 'K11908', 'K21685', 'K07664', 'K19693', 'K07181', 'K07182', 'K19830', 'K12146', 'K20250', 'K05874', 'K06641', 'K07667', 'K19081', 'K11630', 'K03597', 'K07154', 'K07674', 'K06207', 'K07708', 'K04563', 'K15423', 'K21561', 'K11617', 'K11328', 'K19731', 'K07717', 'K11525', 'K13591', 'K07178', 'K02589', 'K20488', 'K14982', 'K07777', 'K11228', 'K11332', 'K13598', 'K15011', 'K13590', 'K11201', 'K14055', 'K11750', 'K14061', 'K17752', 'K18351', 'K21397', 'K07658', 'K20264', 'K07172', 'K10014', 'K18842', 'K10912', 'K07816', 'K07653', 'K13489', 'K04348', 'K21217', 'K07671', 'K02659', 'K20945', 'K03411', 'K07313', 'K13243', 'K07688', 'K07666', 'K21019', 'K02806', 'K07669', 'K21084', 'K01697', 'K11637', 'K02584', 'K17610', 'K07687', 'K17508', 'K04095', 'K20487', 'K07771', 'K07698', 'K03974', 'K19609', 'K03412', 'K12767', 'K13587', 'K11618', 'K10941', 'K19078', 'K02491', 'K09997', 'K17763', 'K08269', 'K07656', 'K18765', 'K13245', 'K19800', 'K20960', 'K20978', 'K13041', 'K07707', 'K02657', 'K14987', 'K10013', 'K07739', 'K03503', 'K07314', 'K17061', 'K19292', 'K13642', 'K06200', 'K03406', 'K05805', 'K12295', 'K14986', 'K21901', 'K04768', 'K07652', 'K07701', 'K07672', 'K02480', 'K04767', 'K02660', 'K11227', 'K18345', 'K07703', 'K18098', 'K02424', 'K06268', 'K07315', 'K18096', 'K19082', 'K04333', 'K06958', 'K20974', 'K06596', 'K07646', 'K06206', 'K19708', 'K21555', 'K07781', 'K19852', 'K13338', 'K20955', 'K07713', 'K19806', 'K10126', 'K00906', 'K09933', 'K07657', 'K20961', 'K03666', 'K20962', 'K18691', 'K02214', 'K18986', 'K07776', 'K02668', 'K20074', 'K11184', 'K21603', 'K19610', 'K20963', 'K13246', 'K19691', 'K03973', 'K02831', 'K19616', 'K02515', 'K07173', 'K16956', 'K10681', 'K20919', 'K11614', 'K07642', 'K00914', 'K07696', 'K02467', 'K07153', 'K21022', 'K11329', 'K21563', 'K16961', 'K13584', 'K02477', 'K13815', 'K20965', 'K15475', 'K02485']
    for amg in amg_pro2info.keys():
        if amg_pro2info[amg][1] in KOs_of_T_or_B:
            amgs_belong_to_not_correct_COG.append(amg)

    # Step 8 Get flitered AMGs
    amg_pro_filtered2info = {} # amg_pro => [scf, ko, ko_name, metabolisms, pathways]
    for amg in amg_pro2info:
        if amg not in tail_amg and amg not in amgs_w_high_v_scores and amg not in amgs_w_flanking_genes_of_low_v_score and amg not in amgs_belong_to_not_correct_COG:
            amg_pro_filtered2info[amg] = amg_pro2info[amg]

    return amg_pro_filtered2info

def write_down_amg_pro2info(AMG_dir, amg_pro2info):
    f = open(os.path.join(AMG_dir,'AMG_pro2info.txt'), 'w')  # The name of the output file
    header = 'Genome\tAMG\tScaffold\tAMG_KO\tAMG_KO_name\tMetabolism\tPathway\n'
    f.write(header)
    amg_pro2info = dict(sorted(amg_pro2info.items())) # Store the dict by the keys
    for amg_pro in amg_pro2info:
        genome = amg_pro.split('---', 1)[0]
        info = '\t'.join(amg_pro2info[amg_pro])
        line = f"{genome}\t{amg_pro}\t{info}\n"
        f.write(line)
    f.close()

def write_down_amg_pro2info_for_wo_reads(AMG_dir, amg_pro2info):
    f = open(os.path.join(AMG_dir,'AMG_pro2info.txt'), 'w')  # The name of the output file
    header = 'AMG\tScaffold\tAMG_KO\tAMG_KO_name\tMetabolism\tPathway\n'
    f.write(header)
    amg_pro2info = dict(sorted(amg_pro2info.items())) # Store the dict by the keys
    for amg_pro in amg_pro2info:
        info = '\t'.join(amg_pro2info[amg_pro])
        line = f"{amg_pro}\t{info}\n"
        f.write(line)
    f.close()

def pick_amg_pro(AMG_dir, amg_pro2info, viral_gn_dir):
    amg_pro_seq_file = os.path.join(AMG_dir,'AMG_pros.faa')  # sequence file for AMG proteins

    all_amg_pro_seq = {} # Store all the AMG proteins
    all_faa_addrs = glob(f'{viral_gn_dir}/*.faa')
    for faa_addr in all_faa_addrs:
        faa_seq = store_seq(faa_addr)
        for header_w_array in faa_seq:
            header = header_w_array.replace('>', '', 1)
            if header in amg_pro2info:
                all_amg_pro_seq[header_w_array] = faa_seq[header_w_array]
    write_down_seq(all_amg_pro_seq, amg_pro_seq_file)

def pick_amg_pro_for_wo_reads(AMG_dir, amg_pro2info, final_virus_faa_file):
    amg_pro_seq_file = os.path.join(AMG_dir,'AMG_pros.faa')  # sequence file for AMG proteins

    all_amg_pro_seq = {} # Store all the AMG proteins
    faa_seq = store_seq(final_virus_faa_file)
    for header_w_array in faa_seq:
        header = header_w_array.replace('>', '', 1)
        if header in amg_pro2info:
            all_amg_pro_seq[header_w_array] = faa_seq[header_w_array]
    write_down_seq(all_amg_pro_seq, amg_pro_seq_file)


def get_amg_statistics(amg_pro2info):
    gn2amg_statistics = {} # gn => amg_statistics; for example, K00018(3);K01953(4)
    amg_list = amg_pro2info.keys()
    gn2amgs = defaultdict(list)
    for amg in amg_list:
        gn = amg.split('---', 1)[0]
        gn2amgs[gn].append(amg)

    for gn in gn2amgs:
        amg_statistics_list = []
        ko2hit_num = {} # ko => hit_num
        amgs = gn2amgs[gn] # The list containing all the AMGs in this genome
        for amg in amgs:
            ko = amg_pro2info[amg][1]
            ko2hit_num[ko] = ko2hit_num.get(ko, 0) + 1
        for ko in ko2hit_num:
            hit_num = ko2hit_num[ko]
            amg_statistics_list.append(f'{ko}({hit_num})')
        gn2amg_statistics[gn] = ';'.join(amg_statistics_list)

    return gn2amg_statistics

def get_amg_statistics_for_wo_reads(amg_pro2info):
    gn2amg_statistics = {} # gn => amg_statistics; for example, K00018(3);K01953(4)
    amg_list = amg_pro2info.keys()
    gn2amgs = defaultdict(list)
    for amg in amg_list:
        gn = amg_pro2info[amg][0] # gn is just the scf name here for single-scaffold genome
        gn2amgs[gn].append(amg)

    for gn in gn2amgs:
        amg_statistics_list = []
        ko2hit_num = {} # ko => hit_num
        amgs = gn2amgs[gn] # The list containing all the AMGs in this genome
        for amg in amgs:
            ko = amg_pro2info[amg][1]
            ko2hit_num[ko] = ko2hit_num.get(ko, 0) + 1
        for ko in ko2hit_num:
            hit_num = ko2hit_num[ko]
            amg_statistics_list.append(f'{ko}({hit_num})')
        gn2amg_statistics[gn] = ';'.join(amg_statistics_list)

    return gn2amg_statistics

def write_down_gn2amg_statistics(AMG_dir, gn2amg_statistics):
    f = open(os.path.join(AMG_dir,'Gn2amg_statistics.txt'), 'w')  # The name of the output file
    header = 'Gn\tAMG_KOs\n'
    f.write(header)
    for gn in gn2amg_statistics:
        amg_statistics = gn2amg_statistics.get(gn, '')
        line = f"{gn}\t{amg_statistics}\n"
        f.write(line)
    f.close()

def get_checkv_useful_info(CheckV_quality_summary):
    checkv_table = pd.read_csv(CheckV_quality_summary, sep = '\t', index_col = 0)
    checkv_table = checkv_table[['checkv_quality', 'miuvig_quality', 'completeness', 'completeness_method']]
    checkv_dict = checkv_table.to_dict() # col => row => value
    return checkv_dict

def get_viral_gn_size_and_scf_no_and_pro_count(viral_gn_dir):
    gn2size_and_scf_no_and_pro_count = {} # gn => [size, scf_no, pro_count]
    all_gn_addrs = glob(f'{viral_gn_dir}/*.fna')
    for gn_addr in all_gn_addrs:
        gn = Path(gn_addr).stem
        gn_seq = store_seq(gn_addr)

        size = 0
        for header in gn_seq:
            size += len(gn_seq[header])

        scf_no = len(gn_seq)

        gn_faa_addr = gn_addr.replace('.fna', '.faa', 1)
        gn_faa_seq = store_seq(gn_faa_addr)
        pro_count = len(gn_faa_seq)

        gn2size_and_scf_no_and_pro_count[gn] = [size, scf_no, pro_count]
    return gn2size_and_scf_no_and_pro_count

def get_viral_gn_size_and_scf_no_and_pro_count_for_wo_reads(final_virus_fasta_file):
    gn2size_and_scf_no_and_pro_count = {} # gn => [size, scf_no, pro_count]
    final_virus_faa_file = final_virus_fasta_file.replace('.fna', '.faa', 1)
    final_virus_fasta_file_seq = store_seq(final_virus_fasta_file)
    final_virus_faa_file_seq = store_seq(final_virus_faa_file)
    for header in final_virus_fasta_file_seq:
        header_wo_array = header.replace('>', '', 1)
        gn = header_wo_array
        size = len(final_virus_fasta_file_seq[header])
        scf_no = 1
        pro_count = 0
        for pro_header in final_virus_faa_file_seq:
            gn_from_pro_header = pro_header.replace('>', '', 1).rsplit('_', 1)[0]
            if gn_from_pro_header == gn:
                pro_count = pro_count + 1
        gn2size_and_scf_no_and_pro_count[gn] = [size, scf_no, pro_count]
    return gn2size_and_scf_no_and_pro_count

def get_gn_lyso_lytic_result(scf2lytic_or_lyso_summary,
                             vRhyme_best_bin_lytic_and_lysogenic_info,
                             viral_gn_dir):
    gn2lyso_lytic_result = {}  # gn => lyso_lytic_property
    prefix = Path(vRhyme_best_bin_lytic_and_lysogenic_info).parent.name   # 样品名

    # Step 1 同原逻辑
    scf2lytic_or_lyso = {}
    with open(scf2lytic_or_lyso_summary, 'r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            tmp = line.split('\t')
            if tmp[0] != 'scaffold':
                scf, lyso = tmp[0], tmp[1]
                if ' (' in lyso:
                    lyso = lyso.split(' (', 1)[0]
                scf2lytic_or_lyso[scf] = [lyso, tmp[2]]

    # Step 2 同原逻辑
    vRhyme_bin2lytic_and_lysogenic_info = {}
    with open(vRhyme_best_bin_lytic_and_lysogenic_info, 'r') as lines:
        for line in lines:
            line = line.rstrip('\n')
            tmp = line.split('\t')
            if tmp[0] != 'vRhyme_bin' and tmp[2] != 'split into scaffolds':
                vRhyme_bin2lytic_and_lysogenic_info[tmp[0]] = tmp[2]

    # Step 3  bin
    for fna in glob(os.path.join(viral_gn_dir, f'{prefix}_vRhyme_bin_*.fna')):
        gn = Path(fna).stem
        bin_key = gn[gn.rfind('vRhyme_'):]          # 稳妥去前缀
        gn2lyso_lytic_result[gn] = vRhyme_bin2lytic_and_lysogenic_info[bin_key]

    # Step 4  unbinned
    for fna in glob(os.path.join(viral_gn_dir, f'{prefix}_vRhyme_unbinned_*.fna')):
        gn = Path(fna).stem
        bin_key = gn[gn.rfind('vRhyme_unbinned_'):]  # 稳妥去前缀
        scf = ''
        for h in store_seq(fna):
            # {name}---{cid} -> {sample}---{cid} (scf2lytic_or_lyso keys keep the
            # original sample-prefixed contig ids)
            gn_sample = re.sub(r'_vRhyme_(bin|unbinned)_\d+$', '', gn)
            scf = gn_sample + '---' + h.lstrip('>').split('---', 1)[1]
        gn2lyso_lytic_result[gn] = scf2lytic_or_lyso[scf][0]

    return gn2lyso_lytic_result

def get_virus_summary_info(checkv_dict, gn2lyso_lytic_result,
                           gn2size_and_scf_no_and_pro_count,
                           gn2amg_statistics, virus_summary_info, virus_species_summary_info,
                           taxa_tsv=None, clu_tsv=None, cov_tsv=None,
                           contig2gn_csv=None):
    gns = list(gn2size_and_scf_no_and_pro_count.keys())

    virus_summary_info_dict = defaultdict(dict)
    virus_summary_info_dict.update(checkv_dict)

    gn2tax = get_gtdb_style_taxa(taxa_tsv) if taxa_tsv and os.path.exists(taxa_tsv) else {}
    gn2cluster = get_cluster_map(clu_tsv) if clu_tsv and os.path.exists(clu_tsv) else {}
    gn2cov = get_coverage_map(cov_tsv) if cov_tsv and os.path.exists(cov_tsv) else {}

    gn2contig_newid = {}
    if contig2gn_csv and os.path.exists(contig2gn_csv):
        contig_df = pd.read_csv(contig2gn_csv, dtype=str)
        for row in contig_df.itertuples(index=False):
            contig_newid = row[0]
            gn  = row[1]
            if gn not in gn2contig_newid:
                gn2contig_newid[gn] = []
            gn2contig_newid[gn].append(contig_newid)
        gn2contig_newid = {gn: ','.join(lst) for gn, lst in gn2contig_newid.items()}

    for gn in gns:
        lytic_state = gn2lyso_lytic_result.get(gn, '')
        virus_summary_info_dict['lytic_state'][gn] = lytic_state
        virus_summary_info_dict['genome_size'][gn] = gn2size_and_scf_no_and_pro_count[gn][0]
        virus_summary_info_dict['scaffold_num'][gn] = gn2size_and_scf_no_and_pro_count[gn][1]
        virus_summary_info_dict['protein_count'][gn] = gn2size_and_scf_no_and_pro_count[gn][2]
        virus_summary_info_dict['AMG_KOs'][gn] = gn2amg_statistics.get(gn, '')

    virus_summary_info_df = pd.DataFrame(virus_summary_info_dict)

    tax_df = pd.DataFrame.from_dict(gn2tax, orient='index')
    virus_summary_info_df = virus_summary_info_df.join(tax_df)

    virus_summary_info_df['cluster'] = virus_summary_info_df.index.map(gn2cluster).fillna('')
    virus_summary_info_df['coverm_abundance'] = virus_summary_info_df.index.map(gn2cov)

    virus_summary_info_df['coverm_abundance'] = pd.to_numeric(virus_summary_info_df['coverm_abundance'], errors='coerce')
    total_abundance = virus_summary_info_df['coverm_abundance'].sum(skipna=True) or 1.0
    virus_summary_info_df['relative_abundance'] = (virus_summary_info_df['coverm_abundance'] / total_abundance)

    virus_summary_info_df.insert(
            1, 'scaffold',
            virus_summary_info_df.index.map(gn2contig_newid).fillna('')
        )

    col_order = ['cluster', 'scaffold', 'lytic_state', 'checkv_quality', 'miuvig_quality',
                 'genome_size', 'scaffold_num', 'protein_count', 'AMG_KOs',
                 'completeness','completeness_method','coverm_abundance','relative_abundance','lineage',
                 'superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']

    virus_summary_info_df['scaffold'] = virus_summary_info_df.apply(
        lambda r: clean_contig(r['scaffold'], str(r.name)),
        axis=1
    )
    virus_summary_info_df = (virus_summary_info_df[col_order].sort_index())
    virus_summary_info_df.to_csv(virus_summary_info, sep='\t', index_label='Virus')

    def agg_cluster(g):
        out = g.head(1).copy()
        out['coverm_abundance'] = g['coverm_abundance'].sum()
        return out
    species_df = (virus_summary_info_df
                .groupby('cluster', group_keys=False)
                .apply(agg_cluster)
                .assign(relative_abundance=lambda d: d['coverm_abundance'] /
                                                    d['coverm_abundance'].sum())
                .sort_values('cluster'))
    species_df.to_csv(virus_species_summary_info, sep='\t', index_label='Virus')

def agg_cluster(g):
    if len(g) == 1:
        return g.assign(relative_abundance=0.0)  # 先占位，后面统一归一化
    # 多行：加权平均 coverm_abundance
    weighted = (g['genome_size'] * g['coverm_abundance']).sum() / g['genome_size'].sum()
    return g.assign(coverm_abundance=weighted, relative_abundance=0.0)

def get_gtdb_style_taxa(taxa_tsv):
    tax = pd.read_csv(taxa_tsv, sep='\t', comment='#', names=[
        'bin', 'classification', 'reason', 'lineage',
        'lineage_scores', 'superkingdom', 'phylum', 'class',
        'order', 'family', 'genus', 'species'
    ])
    tax['bin'] = tax['bin'].apply(lambda x: Path(x).stem)
    return tax.set_index('bin').to_dict('index')

    def parse_one(lin: str):
        out = {r: '' for r in ['superkingdom', 'phylum', 'class', 'order',
                               'family', 'genus', 'species']}
        out['lineage'] = '' if pd.isna(lin) else lin
        if pd.isna(lin) or lin == '' or lin.startswith('no taxid'):
            return out
        ranks = ['d__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']
        for item in lin.split(';'):
            for idx, prefix in enumerate(ranks):
                if item.startswith(prefix):
                    key = ['superkingdom', 'phylum', 'class', 'order',
                           'family', 'genus', 'species'][idx]
                    out[key] = item.replace(prefix, '').strip()
        return out
    parsed = tax['lineage'].apply(parse_one).apply(pd.Series)
    parsed.index = tax['bin']
    return parsed.to_dict('index')



def get_cluster_map(clu_tsv):
    clu = pd.read_csv(clu_tsv, sep='\t', header=None)
    clu.columns = ['object', 'cluster']
    return dict(zip(clu['object'], clu['cluster']))

def get_coverage_map(cov_tsv):
    cov = pd.read_csv(cov_tsv, sep='\t', header=None)
    cov.columns = ['Genome', 'coverm_abundance']
    return dict(zip(cov['Genome'], cov['coverm_abundance']))

def fetch_arguments(parser,root_dir,db_path_default):
    parser.set_defaults(func=main)
    parser.set_defaults(program="run")
    parser.add_argument('--input_metagenome','-i', dest='input_metagenome', required=True, default='none', help=r'(required) input metagenome assembly. It can be a metagenome or entire virome assembly. The extension of the input nucleotide sequence should be ".fna"')
    parser.add_argument('--out_dir','-o', dest='out_dir', required=False, default='./Metalink_MetaG_virus_sample_outdir', help=r'(required) output directory to deposit all results (default = ./Metalink_MetaG_virus_sample_outdir) output folder to deposit all results.')
    parser.add_argument('--identify_method', dest='identify_method', required=False, default='genomad',help=r'(required) the virus identifying method to choose: genomad - Use geNomad; vs_vb - VirSorter2 and CheckV, then use VIBRANT to check the low quality viruses from VirSorter2; genomad_vs_vb - Combine genomad and vs_vb virus to get the union results (default)')
    parser.add_argument('--threads','-t', dest='threads', required=False, default=10, help=r'number of threads (default = 10)')
    parser.add_argument('--conda_env_seqkit', dest='conda_env_seqkit', required=True, default='none', help=r'(required) path to the sra conda environment used by viwrap helper scripts')
    parser.add_argument('--vibrant_db', dest='vibrant_db', required=True, default='none', help=r'(required) path to the VIBRANT database')
    parser.add_argument('--root_dir', dest='root_dir', required=False, default=root_dir,help=argparse.SUPPRESS)
    parser.add_argument('--custom_MAGs_dir', dest='custom_MAGs_dir', required=False, default='none', help=r'custom MAGs dir that contains only *.fna files for MAGs reconstructed from the same metagenome, this will be used in iPHoP for host prediction; note that it should be the absolute address path')
    parser.add_argument('--iPHoP_db_custom', dest='iPHoP_db_custom', required=False, default='none', help=r'custom iPHoP db that will be generated by your input custom MAGs; note that it should be the absolute address path')
    parser.add_argument('--iPHoP_db_custom_pre', dest='iPHoP_db_custom_pre', required=False, default='none', help=r'custom iPHoP db that has been made from the previous run, this will be used in iPHoP for host prediction by custom db; note that it should be the absolute address path')


def set_defaults(args):
    ## Store databases
    args['VIBRANT_db'] = args['vibrant_db']

    print(args['out_dir'])
    print(args['input_metagenome'])

    ## Store outdirs
    args['genomad_outdir'] = os.path.join(args['out_dir'],f"01_identification/{Path(args['input_metagenome']).stem}/genomad_result")
    args['virsorter_outdir'] = os.path.join(args['out_dir'],f"01_identification/{Path(args['input_metagenome']).stem}/vs_vb_result")
    args['genomad_vs_vb_outdir'] = os.path.join(args['out_dir'],f"01_identification/{Path(args['input_metagenome']).stem}/genomad_vs_vb_result")
    args['vrhyme_outdir'] = os.path.join(args['out_dir'],f"04_vRhyme/{Path(args['input_metagenome']).stem}")
    args['iphop_outdir'] = os.path.join(args['out_dir'],f"07_iPHoP/{Path(args['input_metagenome']).stem}")
    args['iphop_custom_outdir'] = os.path.join(args['out_dir'],f"07_iPHoP/iPHoP_outdir_custom_MAGs/{Path(args['input_metagenome']).stem}")
    args['checkv_outdir'] = os.path.join(args['out_dir'],f"05_CheckV/{Path(args['input_metagenome']).stem}")
    args['vclust_outdir'] = os.path.join(args['out_dir'],f"06_vclust/{Path(args['input_metagenome']).stem}")
    args['coverm_outdir'] = os.path.join(args['out_dir'],f"09_coverm/{Path(args['input_metagenome']).stem}")
    args['bin_taxonomy_outdir'] = os.path.join(args['out_dir'],f"08_taxonomy/{Path(args['input_metagenome']).stem}")
    args['Metalink_virus_summary_outdir'] = os.path.join(args['out_dir'],f"10_summary/{Path(args['input_metagenome']).stem}")
    args['final_virus_faa'] = os.path.join(args['out_dir'],f"01_identification/{Path(args['input_metagenome']).stem}/args['identify_method']/final_virus.faa")

    if args['identify_method'] == 'vs_vb':
        args['final_virus_faa'] = os.path.join(args['virsorter_outdir'], 'final_virus.faa')
    elif args['identify_method'] == 'genomad':
        args['final_virus_faa'] = os.path.join(args['genomad_outdir'], 'final_virus.faa')
    elif args['identify_method'] == 'genomad_vs_vb':
        args['final_virus_faa'] = os.path.join(args['genomad_vs_vb_outdir'], 'final_virus.faa')


def main(args):
    # Welcome and logger
    print("### Welcome to Metalink.master_run.vRhyme ###\n")
	## Set up the logger
    ## Set the default args:
    set_defaults(args)
    prepare_output_dir(args['Metalink_virus_summary_outdir'])
    time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
    print(f"{time_current} | Summary all steps...")

    AMG_dir = os.path.join(args['Metalink_virus_summary_outdir'],'AMG_results')
    prepare_output_dir(AMG_dir)
    viral_gn_dir = os.path.join(args['Metalink_virus_summary_outdir'],'Virus_genomes_files')
    prepare_output_dir(viral_gn_dir)
    virus_summary_info = os.path.join(args['Metalink_virus_summary_outdir'], 'vMAG_summary_info.tsv')
    virus_species_summary_info = os.path.join(args['Metalink_virus_summary_outdir'], 'vMAG_species_summary_info.tsv')
    combined_host_pred_to_genome_result = os.path.join(args['Metalink_virus_summary_outdir'],'Host_prediction_to_genome_m90.csv')
    combined_host_pred_to_genus_result = os.path.join(args['Metalink_virus_summary_outdir'],'Host_prediction_to_genus_m90.csv')
    virus_annotation_result_file = os.path.join(args['Metalink_virus_summary_outdir'],'Virus_annotation_results.tsv')
    taxa_tsv = os.path.join(args['bin_taxonomy_outdir'],f"BAT.bin2classification.tsv")
    AMG_pro2info_file = os.path.join(AMG_dir,'AMG_pro2info.txt')
    AMG_pros_faa_file = os.path.join(AMG_dir,'AMG_pros.faa')
    AMG_Gn2amg_statistics_file = os.path.join(AMG_dir,'Gn2amg_statistics.txt')
    contig2gn_csv=os.path.join(args['vrhyme_outdir'],f"contig2viral_gn_map.csv")
    contig_length   = os.path.join(args['Metalink_virus_summary_outdir'], 'final_virus.contig_length.txt')

    scf2lytic_or_lyso_summary = ''
    virus_contigs_fna = ''
    if args['identify_method'] == 'vs_vb':
        scf2lytic_or_lyso_summary = os.path.join(args['virsorter_outdir'], 'scf2lytic_or_lyso.summary.txt')
        virus_contigs_fna = os.path.join(args['virsorter_outdir'], 'final_virus.fna')
    elif args['identify_method'] == 'genomad':
        scf2lytic_or_lyso_summary = os.path.join(args['genomad_outdir'], 'scf2lytic_or_lyso.summary.txt')
        virus_contigs_fna = os.path.join(args['genomad_outdir'], 'final_virus.fna')
    elif args['identify_method'] == 'genomad_vs_vb':
        scf2lytic_or_lyso_summary = os.path.join(args['genomad_vs_vb_outdir'], 'scf2lytic_or_lyso.summary.txt')
        virus_contigs_fna = os.path.join(args['genomad_vs_vb_outdir'], 'final_virus.fna')

    def _only_header(path):
        if not os.path.exists(path) or os.path.getsize(path) == 0:
            return True
        with open(path, 'rb') as f:
            f.readline()
            return f.readline() == b''
    scf_EMPTY = (
        not os.path.exists(scf2lytic_or_lyso_summary) or os.path.getsize(scf2lytic_or_lyso_summary) == 0 or
        not os.path.exists(taxa_tsv) or os.path.getsize(taxa_tsv) == 0 or
        _only_header(contig2gn_csv)
    )
    if scf_EMPTY:
        print(f"{time_current} | Some input tables, like scf2lytic_or_lyso.summary.txt, BAT.bin2classification.tsv, or contig2viral_gn_map.csv, have no information, skipping all summary steps.")


    if not scf_EMPTY:
        # Step 11 Get all virus sequence information
        ## Step 11.1 Move all virus genome fasta, ffn, and faa files
        vRhyme_best_bin_dir_modified = os.path.join(args['vrhyme_outdir'], 'vRhyme_best_bins_fasta_modified')
        vRhyme_unbinned_viral_gn_dir = os.path.join(args['vrhyme_outdir'], 'vRhyme_unbinned_viral_gn_fasta')

        os.system(f'find {vRhyme_best_bin_dir_modified} -type f -exec cp {{}} {viral_gn_dir} \\;')
        os.system(f'find {vRhyme_unbinned_viral_gn_dir} -type f -exec cp {{}} {viral_gn_dir} \\;')

        ## Step 11.2 Combine host prediction result
        combine_iphop_results(args, combined_host_pred_to_genome_result, combined_host_pred_to_genus_result)

        ## Step 11.3 Get virus genome annotation result
        get_virus_genome_annotation_result(args)

        ## Step 11.4 Get AMG results
        amg_pro2info = get_amg_pro_info(AMG_dir, virus_annotation_result_file, viral_gn_dir, args['VIBRANT_db']) # Get the amg_pro2info dict (amg here is filtered)
        write_down_amg_pro2info(AMG_dir, amg_pro2info) # Write down the amg_pro2info dict
        pick_amg_pro(AMG_dir, amg_pro2info, viral_gn_dir) # Pick the AMG proteins and write down the AMG proteins
        gn2amg_statistics = get_amg_statistics(amg_pro2info) # Get the gn2amg_statistics dict
        write_down_gn2amg_statistics(AMG_dir, gn2amg_statistics) # Write down the gn2amg_statistics dict

        ## Step 11.5 Get VIBRANT lytic and lysogenic information and genome information
        vRhyme_best_bin_lytic_and_lysogenic_info = os.path.join(args['vrhyme_outdir'], 'vRhyme_best_bin_lytic_and_lysogenic_info.txt')
        checkv_dict = get_checkv_useful_info(os.path.join(args['checkv_outdir'], 'CheckV_quality_summary.txt'))
        gn2lyso_lytic_result = {}
        if args['identify_method'] == 'genomad' or args['identify_method'] == 'vs_vb' or args['identify_method'] == 'genomad_vs_vb':
            gn2lyso_lytic_result = get_gn_lyso_lytic_result(scf2lytic_or_lyso_summary, vRhyme_best_bin_lytic_and_lysogenic_info, viral_gn_dir)
        gn2size_and_scf_no_and_pro_count = get_viral_gn_size_and_scf_no_and_pro_count(viral_gn_dir)

        cluster_tsv = os.path.join(args['vclust_outdir'], 'clusters.tsv')
        coverage_tsv = os.path.join(args['coverm_outdir'], 'bin_coverage.tsv')
        get_virus_summary_info(checkv_dict, gn2lyso_lytic_result,
                            gn2size_and_scf_no_and_pro_count,
                            gn2amg_statistics, virus_summary_info,virus_species_summary_info,
                            taxa_tsv, cluster_tsv, coverage_tsv, contig2gn_csv)

        if os.path.exists(args['final_virus_faa']):
            merge_faa4_into_annotation(virus_annotation_result_file, args['final_virus_faa'])

        ## Step 11.6 contig length of virus
        calculate_length_cmd = f"conda run -n {args['conda_env_seqkit']} seqkit fx2tab -l -n -i -H {virus_contigs_fna} > {contig_length}"
        os.system(calculate_length_cmd)
        print(calculate_length_cmd)


        time_current = f"[{str(datetime.now().replace(microsecond=0))}]"
        print(f"{time_current} | Get virus sequence information. Finished")
    else:
        with open(virus_summary_info, 'w') as fh:
            fh.write('')
        with open(virus_species_summary_info, 'w') as fh:
            fh.write('')
        with open(combined_host_pred_to_genome_result, 'w') as fh:
            fh.write('')
        with open(combined_host_pred_to_genus_result, 'w') as fh:
            fh.write('')
        with open(virus_annotation_result_file, 'w') as fh:
            fh.write('')
        with open(AMG_pro2info_file, 'w') as fh:
            fh.write('')
        with open(AMG_pros_faa_file, 'w') as fh:
            fh.write('')
        with open(AMG_Gn2amg_statistics_file, 'w') as fh:
            fh.write('')
        with open(contig_length, 'w') as fh:
            fh.write('')

def cli():
    root_dir = os.path.dirname(os.path.abspath(__file__)) # The dir that contains the Metalink-virus script
    current_dir = os.getcwd() # The current dir
    db_path_default = os.path.join(current_dir, "Metalink_virus_db") # Set the default Metalink_virus_db dir under the current directory

    parser = argparse.ArgumentParser(usage = argparse.SUPPRESS, formatter_class=argparse.RawDescriptionHelpFormatter,
        description = """Metalink: Analyzing wrapper for microbiome

Usage: python Metalink.10_summary.py [options]

Task: Summary Metalink outputs for virus genomes from metagenomes
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
