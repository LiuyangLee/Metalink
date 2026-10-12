'''
/*
 * @Date: 2022-12-04 17:15:47
 * @LastEditors: LiuyangLi
 * @LastEditTime: 2024-10-25 16:57:02
 * @FilePath: b-pass_trim.smk
 * @Description: base mode (pass fastp-based reads trimming)
 */
'''
import os
import sys

Dir_script_metalink = config["Directory"]["script"]["metalink"]
kegg_modules_ko_info_path = config["Directory"]["script"]["kegg_modules_ko_info_path"]
kraken2_taxonomy_db = config["Directory"]["database"]["kraken2_taxonomy_db"]

folder_checks = [
                path + 'results/',
                path + 'results/b-sample/',
                path + 'results/b-sample/00_raw_stats/',
                path + 'results/b-sample/01_trim/',
                path + 'results/b-project/'
                ]

for item in folder_checks:
    if not os.path.isdir(item): os.mkdir(item)

rule b_sample_00_compute_raw_stats:
    input:
        r1 = lambda wildcards: SAMPLE_PATHS[wildcards.sample]["read_1"],
        r2 = lambda wildcards: SAMPLE_PATHS[wildcards.sample]["read_2"]
    output:
        r1 = path + 'results/b-sample/00_raw_stats/{sample}_raw_1.fastq.stats',
        r2 = path + 'results/b-sample/00_raw_stats/{sample}_raw_2.fastq.stats',
    conda: config["conda_envs"]["metalink-main"]
    threads: config["Params_b_sample"]["00_raw_stats"]["threads"]
    log:
        path + 'logs/b-sample/compute_raw_stats/{sample}.compute_raw_stats.log'
    benchmark:
        path + 'benchmarks/b-sample/compute_raw_stats/{sample}.compute_raw_stats.benchmark'
    shell:
        """
        seqkit stat -T -b -j {threads} {input.r1} \
        | awk 'BEGIN{{OFS="\\t"}} NR==1{{print; next}} {{$1="{wildcards.sample}"; print}}' \
        > {output.r1} 2> {log}
        seqkit stat -T -b -j {threads} {input.r2} \
        | awk 'BEGIN{{OFS="\\t"}} NR==1{{print; next}} {{$1="{wildcards.sample}"; print}}' \
        > {output.r2} 2>> {log}
        """

rule b_sample_01_link_reads:
    input:
        r1 = lambda wildcards: SAMPLE_PATHS[wildcards.sample]["read_1"],
        r2 = lambda wildcards: SAMPLE_PATHS[wildcards.sample]["read_2"],
    output:
        r1 = path + 'results/b-sample/01_trim/{sample}_trimmed_1.fastq.gz',
        r2 = path + 'results/b-sample/01_trim/{sample}_trimmed_2.fastq.gz',
    threads: 1
    log:
        path + 'logs/b-sample/01_link_reads/{sample}.link_reads.log'
    benchmark:
        path + 'benchmarks/b-sample/01_link_reads/{sample}.link_reads.benchmark'
    shell:
        """
        ln -sf $(realpath {input.r1}) {output.r1} 2> {log}
        ln -sf $(realpath {input.r2}) {output.r2} 2>> {log}
        """

rule b_sample_01_compute_trimmed_stats:
    input:
        r1 = path + 'results/b-sample/00_raw_stats/{sample}_raw_1.fastq.stats',
        r2 = path + 'results/b-sample/00_raw_stats/{sample}_raw_2.fastq.stats',
    output:
        r1 = path + 'results/b-sample/01_trim_stats/{sample}_trimmed_1.fastq.stats',
        r2 = path + 'results/b-sample/01_trim_stats/{sample}_trimmed_2.fastq.stats',
    conda: config["conda_envs"]["metalink-main"]
    threads: config["Params_b_sample"]["01_trim_stats"]["threads"]
    log:
        path + 'logs/b-sample/01_trim_stats/{sample}.compute_trimmed_stats.log'
    benchmark:
        path + 'benchmarks/b-sample/01_trim_stats/{sample}.compute_trimmed_stats.benchmark'
    shell:
        """
        cp {input.r1} {output.r1} 2> {log}
        cp {input.r2} {output.r2} 2>> {log}
        """

rule b_project_01_cat_raw_stats:
    input:
        r1 = expand(path + 'results/b-sample/00_raw_stats/{sample}_raw_1.fastq.stats', sample=SAMPLE_INDEX),
        r2 = expand(path + 'results/b-sample/00_raw_stats/{sample}_raw_2.fastq.stats', sample=SAMPLE_INDEX),
    output:
        r1 = temp(path + 'results/b-project/01_cat_raw_stats/{project}/01_raw_reads_1_stats.tsv'),
        r2 = temp(path + 'results/b-project/01_cat_raw_stats/{project}/01_raw_reads_2_stats.tsv'),
        raw_reads_stats = path + 'results/b-project/01_cat_raw_stats/{project}/01_raw_reads_stats.tsv',
    conda: config["conda_envs"]["metalink-main"]
    threads: config["Params_b_project"]["01_cat_raw_stats"]["threads"]
    log:
        r1 = path + 'logs/b-project/01_cat_raw_stats/{project}/01_cat_raw_stats.log',
    benchmark:
        path + 'benchmarks/b-project/01_cat_raw_stats/{project}/01_cat_raw_stats.benchmark'
    shell:
        """
        for file1 in {input.r1}; do
            # Check if the output_r1 file exists and has a size greater than 0
            if [ -s {output.r1} ]; then
                # Execute sed and awk commands to process the file1 content
                echo "pass"
            else
                # Create a header line in the output_r1 file
                echo -e "ID\tformat_raw_r1\ttype_raw_r1\tnum_seqs_raw_r1\tsum_len_raw_r1\tmin_len_raw_r1\tavg_len_raw_r1\tmax_len_raw_r1" > {output.r1}
            fi
            sed 1d $file1 | awk -F '\t' 'BEGIN {{OFS="\t"}} {{gsub(\"_trimmed_1.fastq.gz\", \"\", $1); print}}' >> {output.r1}
            #sed 1d $file1 | awk -F '\t' 'BEGIN {{OFS="\t"}} {{gsub(\"_1.fq.gz\", \"\", $1); print}}' >> {output.r1}
        done 2> {log}

        for file1 in {input.r2}; do
            # Check if the output_r2 file exists and has a size greater than 0
            if [ -s {output.r2} ]; then
                # Execute sed and awk commands to process the file1 content
                echo "pass"
            else
                # Create a header line in the output_r2 file
                echo -e "ID\tformat_raw_r2\ttype_raw_r2\tnum_seqs_raw_r2\tsum_len_raw_r2\tmin_len_raw_r2\tavg_len_raw_r2\tmax_len_raw_r2" > {output.r2}
            fi
            sed 1d $file1 | awk -F '\t' 'BEGIN {{OFS="\t"}} {{gsub(\"_trimmed_2.fastq.gz\", \"\", $1); print}}' >> {output.r2}
            #sed 1d $file1 | awk -F '\t' 'BEGIN {{OFS="\t"}} {{gsub(\"_2.fq.gz\", \"\", $1); print}}' >> {output.r2}
        done 2>> {log}

        python -c 'import pandas as pd; df1 = pd.read_csv("{output.r1}", sep=\"\t\", header=0); df2 = pd.read_csv("{output.r2}", sep=\"\t\", header=0); merged_df = pd.merge(df1, df2, on="ID", how="outer"); merged_df.to_csv("{output.raw_reads_stats}", index=False, sep=\"\t\")' 2>> {log}
        """

rule b_project_02_cat_trimmed_stats:
    input:
        r1 = expand(path + 'results/b-sample/01_trim_stats/{sample}_trimmed_1.fastq.stats', sample=SAMPLE_INDEX),
        r2 = expand(path + 'results/b-sample/01_trim_stats/{sample}_trimmed_2.fastq.stats', sample=SAMPLE_INDEX),
    output:
        r1 = temp(path + 'results/b-project/02_cat_trimmed_stats/{project}/02_clean_reads_1_stats.tsv'),
        r2 = temp(path + 'results/b-project/02_cat_trimmed_stats/{project}/02_clean_reads_2_stats.tsv'),
        clean_reads_stats = path + 'results/b-project/02_cat_trimmed_stats/{project}/02_clean_reads_stats.tsv',
    conda: config["conda_envs"]["metalink-main"]
    threads: config["Params_b_project"]["02_cat_trimmed_stats"]["threads"]
    log:
        r1 = path + 'logs/b-project/02_cat_trimmed_stats/{project}/02_cat_trimmed_stats.log',
    benchmark:
        path + 'benchmarks/b-project/02_cat_trimmed_stats/{project}/02_cat_trimmed_stats.benchmark'
    shell:
        """
        for file1 in {input.r1}; do
            # Check if the output_r1 file exists and has a size greater than 0
            if [ -s {output.r1} ]; then
                # Execute sed and awk commands to process the file1 content
                echo "pass"
            else
                # Create a header line in the output_r1 file
                echo -e "ID\tformat_trimmed_r1\ttype_trimmed_r1\tnum_seqs_trimmed_r1\tsum_len_trimmed_r1\tmin_len_trimmed_r1\tavg_len_trimmed_r1\tmax_len_trimmed_r1" > {output.r1}
            fi
            #sed 1d $file1 | awk -F '\t' 'BEGIN {{OFS="\t"}} {{gsub(\"/data3/MASH/data/01_trim/\", \"\", $1); gsub(\"_trimmed_1.fastq.gz\", \"\", $1); print}}' >> {output.r1}
            sed 1d $file1 | awk -F '\t' 'BEGIN {{OFS="\t"}} {{gsub(\"_trimmed_1.fastq.gz\", \"\", $1); print}}' >> {output.r1}
        done 2> {log}

        for file1 in {input.r2}; do
            # Check if the output_r2 file exists and has a size greater than 0
            if [ -s {output.r2} ]; then
                # Execute sed and awk commands to process the file1 content
                echo "pass"
            else
                # Create a header line in the output_r2 file
                echo -e "ID\tformat_trimmed_r2\ttype_trimmed_r2\tnum_seqs_trimmed_r2\tsum_len_trimmed_r2\tmin_len_trimmed_r2\tavg_len_trimmed_r2\tmax_len_trimmed_r2" > {output.r2}
            fi
            #sed 1d $file1 | awk -F '\t' 'BEGIN {{OFS="\t"}} {{gsub(\"/data3/MASH/data/01_trim/\", \"\", $1); gsub(\"_trimmed_2.fastq.gz\", \"\", $1); print}}' >> {output.r2}
            sed 1d $file1 | awk -F '\t' 'BEGIN {{OFS="\t"}} {{gsub(\"_trimmed_2.fastq.gz\", \"\", $1); print}}' >> {output.r2}
        done 2>> {log}

        python -c 'import pandas as pd; df1 = pd.read_csv("{output.r1}", sep=\"\t\", header=0); df2 = pd.read_csv("{output.r2}", sep=\"\t\", header=0); merged_df = pd.merge(df1, df2, on="ID", how="outer"); merged_df.to_csv("{output.clean_reads_stats}", index=False, sep=\"\t\")' 2>> {log}
        """
