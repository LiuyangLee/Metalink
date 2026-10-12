# -*- coding: utf-8 -*-
"""
 * @Date: 2023-05-11 19:33:53
 * @LastEditors: LiuyangLi
 * @LastEditTime: 2023-05-12 09:38:11
 * @FilePath: /liliuyang/pylib/genome/fasta_stat_faa.py
 * @Description:
"""
# Import modules, initialize command arguments and variables
import argparse
import re
import pandas as pd
import os

# Statistics for a single fasta file
def fasta_stat(fasta_dir, fasta_file, extension):
    # Regex for fasta file name, remove "." and extension
    re_genome_name = re.compile(r"\.{}".format(extension))
    genome_name = re_genome_name.sub("", fasta_file)
    dict = {}; dict_refer = {}
    seqs_sum = 0
    base_sum = 0

    # Read fasta file, count protein_count (predicted genes) and coding_bases
    with open(os.path.join(fasta_dir, fasta_file), 'r') as read_fas:
        # Check if file is empty
        file_content = read_fas.read().strip()
        if not file_content:
            # Return 0 values for empty file
            return [genome_name, 0, 0]

        # Reset file pointer to beginning and process line by line
        read_fas.seek(0)
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
    read_fas.close()
    # Convert amino-acid count to approximate coding nucleotide bases (aa * 3)
    stat_file = [genome_name, seqs_sum, base_sum * 3]
    return(stat_file)

# Statistics for all fasta files in a directory
def fasta_sum(fasta_dir, fasta_extension, tsv_file):
    INDEX = ['ID',
            'Protein count',
            'Coding bases']
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
    # Set first column to 0 decimal places
    df.iloc[:, 1] = df.iloc[:, 1].apply(lambda x: '{:.0f}'.format(x))
    # Set all columns except first to integers
    df.iloc[:, 2:] = df.iloc[:, 2:].map(lambda x: '{:.0f}'.format(x))
    df.to_csv(tsv_file, sep = '\t', index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = 'This script is used to count basic information of sequences in all fasta files in a folder', add_help = False, usage = '\npython fasta_stat_faa.py -i [fasta.dir] -x [fasta.extension] -o [output.file]\npython fasta_stat.py --input [fasta.dir] --extension [fasta.extension] --output [output.file]')
    parser.add_argument('-i', '--input', metavar = '[fasta.dir]', help = 'Input fasta file directory', required = True)
    parser.add_argument('-x', '--extension', metavar = '[fasta.extension]', help = 'Input fasta file extension (e.g., fa, fna, faa, ffn, etc.)', required = True)
    parser.add_argument('-o', '--output', metavar = '[output.file]', help = 'Output statistics file', required = True)
    parser.add_argument('-h', '--help', action = 'help', help = 'Help information')
    args=parser.parse_args()
    input_dir = args.input
    file_extension = args.extension
    output_file = args.output
    fasta_sum(input_dir, file_extension, output_file)
    print("Genome statistics finished!")
