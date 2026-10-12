'''
/*
 * @Date: 2022-12-04 17:15:47
 * @LastEditors: LiuyangLi
 * @LastEditTime: 2024-10-25 16:57:02
 * @FilePath: b-pass_assembly.smk
 * @Description: base mode (pass megahit-based assembly)
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
                path + 'results/b-sample/02_megahit/',
                path + 'results/b-sample/03_assembly_filter/',
                path + 'results/b-sample/03_index/',
                path + 'results/b-sample/03_mapping/',
                path + 'results/b-sample/04_bam_summary/',
                path + 'results/b-sample/07_contigs/',
                path + 'results/b-project/'
                ]

for item in folder_checks:
    if not os.path.isdir(item): os.mkdir(item)


base_sample_mapping_method = config["Params_b_sample"]["03_mapping"]["method"]
#base_sample_mapping_method=bwamem,bwamem2,strobealign,minibwa
def _base_get_mk(wildcards, rule_type):
    if base_sample_mapping_method == "strobealign":
        if rule_type == "index":
            return rules.b_sample_03_index_strobealign.output.mk
        else:
            return rules.b_sample_03_mapping_strobealign.output.mk
    elif base_sample_mapping_method == "bwamem":
        if rule_type == "index":
            return rules.b_sample_03_index_bwamem.output.mk
        else:
            return rules.b_sample_03_mapping_bwamem.output.mk
    elif base_sample_mapping_method == "bwamem2":
        if rule_type == "index":
            return rules.b_sample_03_index_bwamem2.output.mk
        else:
            return rules.b_sample_03_mapping_bwamem2.output.mk
    elif base_sample_mapping_method == "minibwa":
        if rule_type == "index":
            return rules.b_sample_03_index_minibwa.output.mk
        else:
            return rules.b_sample_03_mapping_minibwa.output.mk
    else:
        raise ValueError(f"Unknown mapping method: {base_sample_mapping_method}")

base_get_03_mapping_rule_mk = lambda wc: _base_get_mk(wc, "mapping")

rule b_sample_03_assembly_filter:
    input:
        path = lambda wildcards: SAMPLE_PATHS[wildcards.sample]["contig_file"],
    output:
        mk = path + 'markers/b-sample/03_assembly_filter/{sample}.mk',
        fa = path + 'results/b-sample/03_assembly_filter/{sample}/{sample}.assembly_filter.fa',
    threads: config["Params_b_sample"]["03_mapping"]["threads"]
    conda: config["conda_envs"]["metalink-main"]
    priority: 80
    params:
        path = path + 'results/b-sample/02_megahit/{sample}/',
        min_contig_length = config["Params_b_sample"]["07_contigs_depth"]["min_contig_length"],
    log:
        path + 'logs/b-sample/03_assembly_filter/{sample}.assembly_filter.log'
    benchmark:
        path + 'benchmarks/b-sample/03_assembly_filter/{sample}.assembly_filter.benchmark'
    shell:
        """
        seqkit seq --remove-gaps --threads {threads} -m {params.min_contig_length} {input.path} \
            | awk -F' ' '{{print $1}}' > {output.fa} 2> {log}
        touch {output.mk}
        """

rule b_sample_03_index_bwamem:
    input:
        assembly_filter = path + 'results/b-sample/03_assembly_filter/{sample}/{sample}.assembly_filter.fa',
    output:
        mk = path + 'markers/b-sample/03_index/bwamem/{sample}.mk',
    params:
        contigs_map = path + 'results/b-sample/03_index/bwamem/{sample}.assembly_filter.fna',
        prefix      = path + 'results/b-sample/03_index/bwamem/{sample}.assembly_filter',
    conda: config["conda_envs"]["metalink-main"]
    priority: 60
    threads: config["Params_b_sample"]["03_mapping"]["threads"]
    log:
        path + 'logs/b-sample/03_index/bwamem/{sample}.03_index.log'
    benchmark:
        path + 'benchmarks/b-sample/03_index/bwamem/{sample}.03_index.benchmark'
    shell:
        """
        mkdir -p $(dirname {params.contigs_map}) 2> {log}
        ln -sf {input.assembly_filter} {params.contigs_map} 2>> {log}
        bwa index -p {params.prefix} {params.contigs_map} 2>> {log}
        touch {output.mk}
        """

rule b_sample_03_mapping_bwamem:
    input:
        r1 = path + 'results/b-sample/01_trim/{sample}_trimmed_1.fastq.gz',
        r2 = path + 'results/b-sample/01_trim/{sample}_trimmed_2.fastq.gz',
        mk = path + 'markers/b-sample/03_index/bwamem/{sample}.mk',
    output:
        bam = path + 'results/b-sample/03_mapping/bwamem/{sample}.bam',
        mk = path + 'markers/b-sample/03_mapping/bwamem/{sample}.mk',
    params:
        prefix = path + 'results/b-sample/03_index/bwamem/{sample}.assembly_filter',
        sort_mem = config["Params_b_sample"]["03_mapping"].get("samtools_sort_mem", "2G"),
        mapping_extra_args = config["Params_b_sample"]["03_mapping"]["extra_args"],
    conda: config["conda_envs"]["metalink-main"]
    priority: 60
    threads: config["Params_b_sample"]["03_mapping"]["threads"]
    log:
        path + 'logs/b-sample/03_mapping/bwamem/{sample}.03_mapping.log'
    benchmark:
        path + 'benchmarks/b-sample/03_mapping/bwamem/{sample}.03_mapping.benchmark'
    shell:
        """
        (bwa mem -t {threads} {params.mapping_extra_args} {params.prefix} {input.r1} {input.r2} \
            2> {log} \
            | samtools view -bS -h -@ {threads} \
            2>> {log} \
            | samtools sort -@ {threads} -m {params.sort_mem} \
            -o {output.bam}) >> {log} 2>&1
        touch {output.mk}
        """


rule b_sample_03_index_bwamem2:
    input:
        assembly_filter = path + 'results/b-sample/03_assembly_filter/{sample}/{sample}.assembly_filter.fa',
    output:
        mk = path + 'markers/b-sample/03_index/bwamem2/{sample}.mk',
    params:
        contigs_map = path + 'results/b-sample/03_index/bwamem2/{sample}.assembly_filter.fna',
        prefix      = path + 'results/b-sample/03_index/bwamem2/{sample}.assembly_filter',
    conda: config["conda_envs"]["metalink-main"]
    priority: 60
    threads: config["Params_b_sample"]["03_mapping"]["threads"]
    log:
        path + 'logs/b-sample/03_index/bwamem2/{sample}.03_index.log'
    benchmark:
        path + 'benchmarks/b-sample/03_index/bwamem2/{sample}.03_index.benchmark'
    shell:
        """
        mkdir -p $(dirname {params.contigs_map}) 2> {log}
        ln -sf {input.assembly_filter} {params.contigs_map} 2>> {log}
        bwa-mem2 index -p {params.prefix} {params.contigs_map} >> {log} 2>&1
        touch {output.mk}
        """

rule b_sample_03_mapping_bwamem2:
    input:
        r1 = path + 'results/b-sample/01_trim/{sample}_trimmed_1.fastq.gz',
        r2 = path + 'results/b-sample/01_trim/{sample}_trimmed_2.fastq.gz',
        mk = path + 'markers/b-sample/03_index/bwamem2/{sample}.mk',
    output:
        bam = path + 'results/b-sample/03_mapping/bwamem2/{sample}.bam',
        mk = path + 'markers/b-sample/03_mapping/bwamem2/{sample}.mk',
    params:
        prefix = path + 'results/b-sample/03_index/bwamem2/{sample}.assembly_filter',
        sort_mem = config["Params_b_sample"]["03_mapping"].get("samtools_sort_mem", "2G"),
        mapping_extra_args = config["Params_b_sample"]["03_mapping"]["extra_args"],
    conda: config["conda_envs"]["metalink-main"]
    priority: 60
    threads: config["Params_b_sample"]["03_mapping"]["threads"]
    log:
        path + 'logs/b-sample/03_mapping/bwamem2/{sample}.03_mapping.log'
    benchmark:
        path + 'benchmarks/b-sample/03_mapping/bwamem2/{sample}.03_mapping.benchmark'
    shell:
        """
        (bwa-mem2 mem -t {threads} {params.mapping_extra_args} {params.prefix} {input.r1} {input.r2} \
            2> {log} \
            | samtools view -bS -h -@ {threads} \
            2>> {log} \
            | samtools sort -@ {threads} -m {params.sort_mem} \
            -o {output.bam}) >> {log} 2>&1
        touch {output.mk}
        """

rule b_sample_03_index_minibwa:
    input:
        assembly_filter = path + 'results/b-sample/03_assembly_filter/{sample}/{sample}.assembly_filter.fa',
    output:
        mk = path + 'markers/b-sample/03_index/minibwa/{sample}.mk',
    params:
        contigs_map = path + 'results/b-sample/03_index/minibwa/{sample}.assembly_filter.fna',
        prefix      = path + 'results/b-sample/03_index/minibwa/{sample}.assembly_filter',
    conda: config["conda_envs"]["metalink-main"]
    priority: 60
    threads: config["Params_b_sample"]["03_mapping"]["threads"]
    log:
        path + 'logs/b-sample/03_index/minibwa/{sample}.03_index.log'
    benchmark:
        path + 'benchmarks/b-sample/03_index/minibwa/{sample}.03_index.benchmark'
    shell:
        """
        mkdir -p $(dirname {params.contigs_map}) 2> {log}
        ln -sf {input.assembly_filter} {params.contigs_map} 2>> {log}
        minibwa index {params.contigs_map} {params.prefix} -t {threads} 2>> {log}
        touch {output.mk}
        """

rule b_sample_03_mapping_minibwa:
    input:
        r1 = path + 'results/b-sample/01_trim/{sample}_trimmed_1.fastq.gz',
        r2 = path + 'results/b-sample/01_trim/{sample}_trimmed_2.fastq.gz',
        mk = path + 'markers/b-sample/03_index/minibwa/{sample}.mk',
    output:
        bam = path + 'results/b-sample/03_mapping/minibwa/{sample}.bam',
        mk = path + 'markers/b-sample/03_mapping/minibwa/{sample}.mk',
    params:
        prefix = path + 'results/b-sample/03_index/minibwa/{sample}.assembly_filter',
        sort_mem = config["Params_b_sample"]["03_mapping"].get("samtools_sort_mem", "2G"),
        mapping_extra_args = config["Params_b_sample"]["03_mapping"]["extra_args"],
    conda: config["conda_envs"]["metalink-main"]
    priority: 60
    threads: config["Params_b_sample"]["03_mapping"]["threads"]
    log:
        path + 'logs/b-sample/03_mapping/minibwa/{sample}.03_mapping.log'
    benchmark:
        path + 'benchmarks/b-sample/03_mapping/minibwa/{sample}.03_mapping.benchmark'
    shell:
        """
        (minibwa map -t {threads} {params.mapping_extra_args} \
            {params.prefix} \
            {input.r1} {input.r2} \
            2> {log} \
            | samtools view -bS -h -@ {threads} \
            2>> {log} \
            | samtools sort -@ {threads} -m {params.sort_mem} \
            -o {output.bam}) >> {log} 2>&1
        touch {output.mk}
        """

rule b_sample_03_index_strobealign:
    priority: 60
    input:
        assembly_filter = path + 'results/b-sample/03_assembly_filter/{sample}/{sample}.assembly_filter.fa',
        r1 = expand(path + 'results/b-sample/01_trim/{sample}_trimmed_1.fastq.gz', sample=SAMPLE_INDEX),
    output:
        mk = path + 'markers/b-sample/03_index/strobealign/{sample}.mk',
    params:
        contigs_map = path + 'results/b-sample/03_index/strobealign/{sample}.assembly_filter.fna',
        len_lst     = path + 'results/b-sample/03_index/strobealign/{sample}.unique_scaled_length.txt',
    conda: config["conda_envs"]["metalink-main"]
    threads: config["Params_b_sample"]["03_mapping"]["threads"]
    log:
        path + 'logs/b-sample/03_index/strobealign/{sample}.03_index.log'
    benchmark:
        path + 'benchmarks/b-sample/03_index/strobealign/{sample}.03_index.benchmark'
    shell:
        """
        mkdir -p $(dirname {params.contigs_map}) 2> {log}
        ln -sf {input.assembly_filter} {params.contigs_map} 2>> {log}

        function process_number() {{
            local num=${{1:-0}}
            if   (( num <= 70 )); then echo 50
            elif (( num <= 90 )); then echo 75
            elif (( num <= 110 )); then echo 100
            elif (( num <= 135 )); then echo 125
            elif (( num <= 175 )); then echo 150
            elif (( num <= 375 )); then echo 250
            else echo 400
            fi
        }}

        for fq in {input.r1}
        do
            Length=`seqkit head -n 1000 $fq \
                |seqkit fx2tab -l | cut -f 2 | awk '{{ if (length > max) max = length }} END {{ print max }}'`
            ScaledLength=$(process_number $Length)
            echo $ScaledLength
        done > {params.len_lst} 2>> {log}

        sort -nu {params.len_lst} -o {params.len_lst} 2>> {log}

        for L in $(cat {params.len_lst}); do
            strobealign --create-index -t {threads} -r $L {params.contigs_map}
        done >> {log} 2>&1

        touch {output.mk}
        """

rule b_sample_03_mapping_strobealign:
    input:
        r1 = path + 'results/b-sample/01_trim/{sample}_trimmed_1.fastq.gz',
        r2 = path + 'results/b-sample/01_trim/{sample}_trimmed_2.fastq.gz',
        mk = path + 'markers/b-sample/03_index/strobealign/{sample}.mk',
    output:
        bam = path + 'results/b-sample/03_mapping/strobealign/{sample}.bam',
        mk = path + 'markers/b-sample/03_mapping/strobealign/{sample}.mk',
    params:
        contigs_map = path + 'results/b-sample/03_index/strobealign/{sample}.assembly_filter.fna',
        sort_mem = config["Params_b_sample"]["03_mapping"].get("samtools_sort_mem", "2G"),
        mapping_extra_args = config["Params_b_sample"]["03_mapping"]["extra_args"],
    conda: config["conda_envs"]["metalink-main"]
    priority: 60
    threads: config["Params_b_sample"]["03_mapping"]["threads"]
    log:
        path + 'logs/b-sample/03_mapping/strobealign/{sample}.03_mapping.log'
    benchmark:
        path + 'benchmarks/b-sample/03_mapping/strobealign/{sample}.03_mapping.benchmark'
    shell:
        """
        function process_number() {{
            local num=${{1:-0}}
            if   (( num <= 70 )); then echo 50
            elif (( num <= 90 )); then echo 75
            elif (( num <= 110 )); then echo 100
            elif (( num <= 135 )); then echo 125
            elif (( num <= 175 )); then echo 150
            elif (( num <= 375 )); then echo 250
            else echo 400
            fi
        }}
        Length=`seqkit head -n 1000 {input.r1} | seqkit fx2tab -l | cut -f 2 | awk '{{ if (length > max) max = length }} END {{ print max }}'`
        ScaledLength=$(process_number $Length)
        (strobealign --use-index -t {threads} -r $ScaledLength {params.mapping_extra_args} {params.contigs_map} {input.r1} {input.r2} \
            2> {log} \
            | samtools sort -@ {threads} -m {params.sort_mem} -o {output.bam}) >> {log} 2>&1
        touch {output.mk}
        """

rule b_sample_04_bam_summary:
    priority: 80
    input:
        bam_raw = path + 'results/b-sample/03_mapping/' + base_sample_mapping_method + '/{sample}.bam',
    output:
        stat_raw = path + 'results/b-sample/04_bam_summary/{sample}.contigs_bam_stats.tsv',
    threads: config["Params_b_sample"]["03_mapping"]["threads"]
    conda: config["conda_envs"]["metalink-main"]
    log:
        path + 'logs/b-sample/04_bam_summary/{sample}.bam_summary.log'
    benchmark:
        path + 'benchmarks/b-sample/04_bam_summary/{sample}.bam_summary.benchmark'
    shell:
        """
        samtools flagstat {input.bam_raw} > {output.stat_raw} 2> {log}
        """

rule b_sample_07_contigs_prodigal:
    priority: 80
    input:
        mk1 = base_get_03_mapping_rule_mk,
        bam = path + 'results/b-sample/03_mapping/' + base_sample_mapping_method + '/{sample}.bam',
    output:
        contigs_fna = path + 'results/b-sample/07_contigs/contigs_fna/{sample}.fna',
        protein_ffn = path + 'results/b-sample/07_contigs/contigs_protein_ffn/{sample}.ffn',
        protein_faa = path + 'results/b-sample/07_contigs/contigs_protein_faa/{sample}.faa',
        length = path + 'results/b-sample/07_contigs/contigs_ffn_length/{sample}',
        assembly_filter_depth = path + 'results/b-sample/07_contigs/contigs_filter_depth/{sample}',
        ORF_loci = path + 'results/b-sample/07_contigs/contigs_ORF_loci/{sample}',
    conda: config["conda_envs"]["metalink-main"]
    threads: config["Params_b_sample"]["07_contigs_depth"]["threads"]
    params:
        sample='{sample}',
        ctg = path + 'results/b-sample/03_index/' + base_sample_mapping_method + '/{sample}.assembly_filter.fna',
        assembly_filter_depth = path + 'results/b-sample/07_contigs/contigs_filter_depth/{sample}.tmp',
        min_read_percent_identity = config["Params_b_sample"]["07_contigs_depth"]["min_read_percent_identity"],
        contig_end_exclusion = config["Params_b_sample"]["07_contigs_depth"]["contig_end_exclusion"],
        min_contig_length = config["Params_b_sample"]["07_contigs_depth"]["min_contig_length"],
        extra_args = config["Params_b_sample"]["07_contigs_depth"]["extra_args"],
    log:
        path + 'logs/b-sample/07_contigs_prodigal/{sample}.log'
    shell:
        """
        if [ -s {params.ctg} ]; then
            rm -f {output.contigs_fna} {output.protein_ffn} {output.protein_faa} {output.length} {output.assembly_filter_depth} {output.ORF_loci} {params.assembly_filter_depth}
            sed -e "s|^>|>{params.sample}---|" {params.ctg} > {output.contigs_fna}
            prodigal {params.extra_args} \
                -p meta \
                -q \
                -i {output.contigs_fna} \
                -d {output.protein_ffn} \
                -a {output.protein_faa} -o /dev/null
            seqkit fx2tab -l -n -i -H -j {threads} {output.protein_ffn} > {output.length}
            jgi_summarize_bam_contig_depths {input.bam} --outputDepth {params.assembly_filter_depth} --maxEdgeBases {params.contig_end_exclusion} --percentIdentity {params.min_read_percent_identity} --minContigLength {params.min_contig_length}
            awk -F'\t' '{{print $1"\t"$2"\t"$3}}' {params.assembly_filter_depth} | sed -e "s|^|{params.sample}---|" > {output.assembly_filter_depth}
            awk 'BEGIN {{OFS="\t"; print "ORF_name", "start", "end", "strand"}}
            />/ {{
                split($0, parts, " # ")
                orf_name = parts[1]
                sub(/^>/, "", orf_name)
                print orf_name, parts[2], parts[3], parts[4]
            }}' {output.protein_faa} > {output.ORF_loci}
            rm -f {params.assembly_filter_depth}
        else
            echo "The assembly_filter is empty, pass."
            touch {output.contigs_fna} {output.protein_ffn} {output.protein_faa} {output.length} {output.assembly_filter_depth} {output.ORF_loci}
        fi >{log} 2>&1
        """

rule b_project_03_contigs_stats:
    input:
        contigs_fna = expand(path + 'results/b-sample/07_contigs/contigs_fna/{sample}.fna', sample=SAMPLE_INDEX),
    output:
        mk = path + 'markers/b-project/03_contigs_stats/{project}/03_contigs_stats.mk',
    threads: config["Params_b_project"]["03_contigs_stats"]["threads"]
    conda: config["conda_envs"]["metalink-main"]
    params:
        out_path = path + 'results/b-project/03_contigs_stats/{project}/',
        genometype_1 = 'fna',
        genometype_2 = 'ffn',
        in_path_1 = path + 'results/b-sample/07_contigs/contigs_fna/',
        in_path_2 = path + 'results/b-sample/07_contigs/contigs_protein_ffn/',
        output_1 = path +'results/b-project/03_contigs_stats/{project}/03_contigs_stats.tsv',
        output_2 = path + 'results/b-project/03_contigs_stats/{project}/03_contigs_ffn_stats.tsv',
    priority: 80
    log:
        path + 'logs/b-project/03_contigs_stats/{project}/03_contigs_stats.log'
    benchmark:
        path + 'benchmarks/b-project/03_contigs_stats/{project}/03_contigs_stats.benchmark'
    shell:
        """
        mkdir -p {params.out_path}
        python {Dir_script_metalink}/fasta_stat.py -i {params.in_path_1} -x {params.genometype_1} -o {params.output_1} > {log} 2>&1
        python {Dir_script_metalink}/fasta_stat_faa.py -i {params.in_path_2} -x {params.genometype_2} -o {params.output_2} >> {log} 2>&1
        touch {output.mk}
        """

rule b_project_04_merge_reads_contigs_stats:
    input:
        raw   = path + "results/b-project/01_cat_raw_stats/{project}/01_raw_reads_stats.tsv",
        clean = path + "results/b-project/02_cat_trimmed_stats/{project}/02_clean_reads_stats.tsv",
        ctg   = path + "markers/b-project/03_contigs_stats/{project}/03_contigs_stats.mk",
    output:
        reads              = temp(path + "results/b-project/04_merge_reads_contigs_stats/{project}/04_reads_stats.tsv"),
        reads_contigs      = temp(path + "results/b-project/04_merge_reads_contigs_stats/{project}/04_reads_contigs_stats.tsv"),
        reads_contigs_ffn  = path + "results/b-project/04_merge_reads_contigs_stats/{project}/04_reads_contigs_ffn_stats.tsv",
    params:
        ctg_tsv = path +'results/b-project/03_contigs_stats/{project}/03_contigs_stats.tsv',
        ctg_ffn = path + 'results/b-project/03_contigs_stats/{project}/03_contigs_ffn_stats.tsv',
    threads: config["Params_b_project"]["04_merge_reads_contigs_stats"]["threads"]
    log:
        path + "logs/b-project/04_merge_reads_contigs_stats/{project}/04_merge_reads_contigs_stats.log"
    benchmark:
        path + "benchmarks/b-project/04_merge_reads_contigs_stats/{project}/04_merge_reads_contigs_stats.benchmark"
    conda: config["conda_envs"]["metalink-main"]
    shell:
        """
        python3 -c '
import pandas as pd
def merge_and_save(left_path, right_path, out_path, on="ID"):
    left = pd.read_csv(left_path, sep="\t")
    right = pd.read_csv(right_path, sep="\t")
    merged = left.merge(right, on=on, how="outer", suffixes=("", "_drop"))
    merged = merged.loc[:, ~merged.columns.str.endswith("_drop")]
    merged.to_csv(out_path, sep="\t", index=False)
merge_and_save("{input.raw}", "{input.clean}", "{output.reads}")
merge_and_save("{output.reads}", "{params.ctg_tsv}", "{output.reads_contigs}")
merge_and_save("{output.reads_contigs}", "{params.ctg_ffn}", "{output.reads_contigs_ffn}")
' > {log} 2>&1
        """
