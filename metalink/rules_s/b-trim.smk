'''
/*
 * @Date: 2022-12-04 17:15:47
 * @LastEditors: LiuyangLi
 * @LastEditTime: 2024-10-25 16:57:02
 * @FilePath: /metalink/b-sample.smk
 * @Description: merged reads to assembly and short reads annotation
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

base_sample_read_type = config["read_type"]
#base_sample_read_type=paired,single,pacbio-raw,pacbio-corr,pacbio-hifi,nano-raw,nano-corr,nano-hq

rule b_sample_00_compute_raw_stats:
    input:
        r_s = lambda wildcards: SAMPLE_PATHS[wildcards.sample]["read_s"],
    output:
        r_s = path + 'results/b-sample/00_raw_stats/{sample}_raw_s.fastq.stats',
    conda:
        config["conda_envs"]["metalink-main"]
    threads: config["Params_b_sample"]["00_raw_stats"]["threads"]
    log:
        path + 'logs/b-sample/compute_raw_stats/{sample}.compute_raw_stats.log'
    benchmark:
        path + 'benchmarks/b-sample/compute_raw_stats/{sample}.compute_raw_stats.benchmark'
    shell:
        """
        seqkit stat -T -b -j {threads} {input.r_s} \
        | awk 'BEGIN{{OFS="\\t"}} NR==1{{print; next}} {{$1="{wildcards.sample}"; print}}' \
        > {output.r_s} 2> {log}
        """

if base_sample_read_type in ("single", "single-rna"):
    rule b_sample_01_trim_fastp:
        input:
            r_s = lambda wildcards: SAMPLE_PATHS[wildcards.sample]["read_s"],
        output:
            r_s = path + 'results/b-sample/01_trim/{sample}_trimmed_s.fastq.gz',
            json = path + 'results/b-sample/01_trim_json/{sample}_fastp.json',
            html = path + 'results/b-sample/01_trim_html/{sample}_fastp.html',
        params:
            length_required = config["Params_b_sample"]["01_trim"]["length_required"],
            n_base_limit = config["Params_b_sample"]["01_trim"]["n_base_limit"],
            compression = config["Params_b_sample"]["01_trim"]["compression"],
            fastp_extra_args = config["Params_b_sample"]["01_trim"]["extra_args"],
        conda:
            config["conda_envs"]["metalink-main"]
        priority: 1
        threads: config["Params_b_sample"]["01_trim"]["threads"]
        log:
            path + 'logs/b-sample/01_trim/{sample}.trim_fastp.log'
        benchmark:
            path + 'benchmarks/b-sample/01_trim/{sample}.trim_fastp.benchmark'
        shell:
            """
            fastp \
                -i {input.r_s} \
                -o {output.r_s} \
                -j {output.json} \
                -h {output.html} \
                --thread={threads} \
                --length_required={params.length_required} \
                --n_base_limit={params.n_base_limit} \
                --compression={params.compression} {params.fastp_extra_args} > {log} 2>&1
            """
else:
    rule b_sample_01_trim_fastplong:
        input:
            r_s = lambda wildcards: SAMPLE_PATHS[wildcards.sample]["read_s"],
        output:
            r_s = path + 'results/b-sample/01_trim/{sample}_trimmed_s.fastq.gz',
            json = path + 'results/b-sample/01_trim_json/{sample}_fastp.json',
            html = path + 'results/b-sample/01_trim_html/{sample}_fastp.html',
        params:
            length_required = config["Params_b_sample"]["01_trim"].get("length_required", 15),
            n_base_limit = config["Params_b_sample"]["01_trim"].get("n_base_limit", 5),
            compression = config["Params_b_sample"]["01_trim"].get("compression", 4),
            fastplong_extra_args = config["Params_b_sample"]["01_trim"].get("extra_args", ""),
        conda:
            config["conda_envs"]["metalink-main"]
        priority: 1
        threads: config["Params_b_sample"]["01_trim"]["threads"]
        log:
            path + 'logs/b-sample/01_trim/{sample}.trim_fastplong.log'
        benchmark:
            path + 'benchmarks/b-sample/01_trim/{sample}.trim_fastplong.benchmark'
        shell:
            """
            fastplong \
                -i {input.r_s} \
                -o {output.r_s} \
                -j {output.json} \
                -h {output.html} \
                --thread={threads} \
                --length_required={params.length_required} \
                --n_base_limit={params.n_base_limit} \
                --compression={params.compression} {params.fastplong_extra_args} > {log} 2>&1
            """

rule b_sample_01_compute_trimmed_stats:
    input:
        r_s = path + 'results/b-sample/01_trim/{sample}_trimmed_s.fastq.gz',
    output:
        r_s = path + 'results/b-sample/01_trim_stats/{sample}_trimmed_s.fastq.stats',
    conda:
        config["conda_envs"]["metalink-main"]
    threads: config["Params_b_sample"]["01_trim_stats"]["threads"]
    log:
        path + 'logs/b-sample/01_trim_stats/{sample}.compute_trimmed_stats.log'
    benchmark:
        path + 'benchmarks/b-sample/01_trim_stats/{sample}.compute_trimmed_stats.benchmark'
    shell:
        """
        seqkit stat -T -b -j {threads} {input.r_s} > {output.r_s} 2> {log}
        """

rule b_project_01_cat_raw_stats:
    input:
        r_s = expand(path + 'results/b-sample/00_raw_stats/{sample}_raw_s.fastq.stats', sample=SAMPLE_INDEX),
    output:
        r_s = path + 'results/b-project/01_cat_raw_stats/{project}/01_raw_reads_stats.tsv',
    conda:
        config["conda_envs"]["metalink-main"]
    threads: config["Params_b_project"]["01_cat_raw_stats"]["threads"]
    log:
        r_s = path + 'logs/b-project/01_cat_raw_stats/{project}/01_cat_raw_stats.log',
    benchmark:
        path + 'benchmarks/b-project/01_cat_raw_stats/{project}/01_cat_raw_stats.benchmark'
    shell:
        """
        for file1 in {input.r_s}; do
            # Check if the output_r_s file exists and has a size greater than 0
            if [ -s {output.r_s} ]; then
                # Execute sed and awk commands to process the file1 content
                echo "pass"
            else
                # Create a header line in the output_r_s file
                echo -e "ID\tformat_raw_s\ttype_raw_s\tnum_seqs_raw_s\tsum_len_raw_s\tmin_len_raw_s\tavg_len_raw_s\tmax_len_raw_s" > {output.r_s}
            fi
            sed 1d $file1 | awk -F '\t' 'BEGIN {{OFS="\t"}} {{gsub(\"_raw_s.fastq.gz\", \"\", $1); print}}' >> {output.r_s}
        done 2> {log}
        """

rule b_project_02_cat_trimmed_stats:
    input:
        r_s = expand(path + 'results/b-sample/01_trim_stats/{sample}_trimmed_s.fastq.stats', sample=SAMPLE_INDEX),
    output:
        r_s = path + 'results/b-project/02_cat_trimmed_stats/{project}/02_clean_reads_stats.tsv',
    conda:
        config["conda_envs"]["metalink-main"]
    threads: config["Params_b_project"]["02_cat_trimmed_stats"]["threads"]
    log:
        r_s = path + 'logs/b-project/02_cat_trimmed_stats/{project}/02_cat_trimmed_stats.log',
    benchmark:
        path + 'benchmarks/b-project/02_cat_trimmed_stats/{project}/02_cat_trimmed_stats.benchmark'
    shell:
        """
        for file1 in {input.r_s}; do
            # Check if the output_r_s file exists and has a size greater than 0
            if [ -s {output.r_s} ]; then
                # Execute sed and awk commands to process the file1 content
                echo "pass"
            else
                # Create a header line in the output_r_s file
                echo -e "ID\tformat_trimmed_s\ttype_trimmed_s\tnum_seqs_trimmed_s\tsum_len_trimmed_s\tmin_len_trimmed_s\tavg_len_trimmed_s\tmax_len_trimmed_s" > {output.r_s}
            fi
            sed 1d $file1 | awk -F '\t' 'BEGIN {{OFS="\t"}} {{gsub(\"_trimmed_s.fastq.gz\", \"\", $1); print}}' >> {output.r_s}
        done 2> {log}
        """
