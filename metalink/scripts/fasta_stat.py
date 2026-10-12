# -*- coding: utf-8 -*-
"""
 * @Date: 2023-05-11 19:33:53
 * @LastEditors: LiuyangLi
 * @LastEditTime: 2023-05-12 09:38:11
 * @FilePath: /liliuyang/pylib/genome/fasta_stat.py
 * @Description:
"""
# Import modules, initialize command arguments and variables
import argparse
import re
import pandas as pd
import os

# Statistics for a single fasta file
def fasta_stat(fasta_dir, fasta_file, extension):
    # Regular expression for fasta file name, remove "." and extension
    re_genome_name = re.compile(r"\.{}".format(extension))
    genome_name = re_genome_name.sub("", fasta_file)
    dict = {}; dict_refer = {}
    # Initialize L50 = 0; L75 = 0; L90 = 0
    seqs_sum = 0
    base_sum = 0
    GC_sum = 0
    N_sum = 0
    Ln = 0; base_sum_n = 0
    L50 = 0; L75 = 0; L90 = 0
    N50 = 0; N75 = 0; N90 = 0
    n1000000 = 0; n1000000_len = 0
    n100000 = 0; n100000_len = 0
    n10000 = 0; n10000_len = 0
    n1000 = 0; n1000_len = 0
    # Read fasta file, count scaffolds (each + total) number, length, GC content and other basic information
    with open(os.path.join(fasta_dir, fasta_file), 'r') as read_fas:
        for line in read_fas:
            if line[0] == '>':
                key = line.strip('[ >\n]')
                dict[key] = 0
                seqs_sum += 1
            else:
                value = line.strip()
                seqs_len = len(value)
                base_sum += seqs_len
                dict[key] += seqs_len
                GC_sum += len(re.findall('[gcGC]', value))
                N_sum += len(re.findall('[nN]', value))
    read_fas.close()

    # Check if empty file, if yes return all-zero result directly
    if seqs_sum == 0:
        stat_file = [genome_name, 0.0000,
                    0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0,
                    0, 0, 0, 0]
        return(stat_file)

    # Calculate N50, N75, N90 and scaffold counts/lengths for each length interval
    dict_sort = sorted(dict.items(), key = lambda x:x[1], reverse = True)
    base_sum_n50 = 5 * base_sum/10
    base_sum_n75 = 7.5 * base_sum/10
    base_sum_n90 = 9 * base_sum/10
    for value in dict_sort:
        base_sum_n += value[1]
        Ln += 1
        if N50:
            if N75:
                if base_sum_n >= base_sum_n90:
                    N90 = value[1]
                    L90 = Ln
                    break
            else:
                if base_sum_n >= base_sum_n75:
                    N75 = value[1]
                    L75 = Ln
        else:
            if base_sum_n >= base_sum_n50:
                N50 = value[1]
                L50 = Ln
    for value in dict_sort:
        if value[1] >= 1000000:
            n1000000 += 1
            n1000000_len += value[1]
        if value[1] >= 100000:
            n100000 += 1
            n100000_len += value[1]
        if value[1] >= 10000:
            n10000 += 1
            n10000_len += value[1]
        if value[1] >= 1000:
            n1000 += 1
            n1000_len += value[1]
        if value[1] < 1000:
            break
    gc = round(100 * GC_sum / base_sum, 4)
    ave_len = int(base_sum/seqs_sum)
    l_scaf = dict_sort[0][1]
    s_scaf = dict_sort[-1][1]
    stat_file = [genome_name, gc,
                seqs_sum, base_sum, N_sum, ave_len, l_scaf, s_scaf,
                L50, N50, L75, N75, L90, N90,
                n1000000, n100000, n10000, n1000,
                n1000000_len, n100000_len, n10000_len, n1000_len]
    return(stat_file)

# Statistics for all fasta files in directory
def fasta_sum(fasta_dir, fasta_extension, tsv_file):
    INDEX = ['ID',
            'GC (%)',
            'Total scaffolds',
            'Total base (bp)',
            'Total N (bp)',
            'Average length (bp)',
            'Longest scaffold (bp)',
            'Shortest scaffold (bp)',
            'L50',
            'N50',
            'L75',
            'N75',
            'L90',
            'N90',
            'Total scaffolds (>= 1000000 bp)',
            'Total scaffolds (>= 100000 bp)',
            'Total scaffolds (>= 10000 bp)',
            'Total scaffolds (>= 1000 bp)',
            'Total length (>= 1000000 bp)',
            'Total length (>= 100000 bp)',
            'Total length (>= 10000 bp)',
            'Total length (>= 1000 bp)']
    # Search for all files with specified extension
    fasta_list = []
    for each_fasta in os.listdir(fasta_dir):
        if each_fasta.endswith(f".{fasta_extension}"):
            fasta_list.append(each_fasta)

    # Create an empty DataFrame
    df = pd.DataFrame(columns = INDEX)
    # Iterate through all genome files for statistics
    for fasta in fasta_list:
        genome_stat = fasta_stat(fasta_dir, fasta, fasta_extension)
        genome_stat_s = pd.Series(genome_stat, index = INDEX, name = fasta)
        df = pd.concat([df, genome_stat_s.to_frame().T], ignore_index=False)
    # Set first column to 4 decimal places
    df.iloc[:, 1] = df.iloc[:, 1].apply(lambda x: '{:.4f}'.format(x))
    # Set all columns except first to integers
    df.iloc[:, 2:] = df.iloc[:, 2:].map(lambda x: '{:.0f}'.format(x))
    df.to_csv(tsv_file, sep = '\t', index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = 'This script is used to count basic information of sequences in all fasta files in a folder', add_help = False, usage = '\npython fasta_stat.py -i [fasta.dir] -x [fasta.extension] -o [output.file]\npython fasta_stat.py --input [fasta.dir] --extension [fasta.extension] --output [output.file]')
    parser.add_argument('-i', '--input', metavar = '[fasta.dir]', help = 'Input fasta directory', required = True)
    parser.add_argument('-x', '--extension', metavar = '[fasta.extension]', help = 'Input fasta file extension (e.g., fa, fna, faa, ffn, etc.)', required = True)
    parser.add_argument('-o', '--output', metavar = '[output.file]', help = 'Output statistics result', required = True)
    parser.add_argument('-h', '--help', action = 'help', help = 'Help information')
    args=parser.parse_args()
    input_dir = args.input
    file_extension = args.extension
    output_file = args.output
    fasta_sum(input_dir, file_extension, output_file)
    print("Genome statistics finished!")
