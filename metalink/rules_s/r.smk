'''
/*
 * @Date: 2023-11-16 16:06:02
 * @LastEditors: LiuyangLi
 * @LastEditTime: 2024-07-04 16:22:32
 * @FilePath: r.smk
 * @Description:
 */
'''

import os
import sys

reads_function_classify_method = config["Params_r_sample"]["04_function"]["reads_function_classify_method"]
taxid_tax_path = config["Directory"]["database"]["taxid_tax_path"]
taxid_fun_path = config["Directory"]["database"]["taxid_fun_path"]
KO_length_ref_path = config["Directory"]["script"]["KO_length_ref_path"]
Dir_script_metalink = config["Directory"]["script"]["metalink"]
kegg_modules_ko_info_path = config["Directory"]["script"]["kegg_modules_ko_info_path"]

# Third-generation (long-read) sequencing types for which Bracken should be skipped
_long_read_types = ['pacbio-raw', 'pacbio-corr', 'pacbio-hifi', 'nano-raw', 'nano-corr', 'nano-hq']


folder_checks = [
                path + 'results/',
                path + 'results/r-sample/',
                path + 'results/r-sample/02_FragGeneScanRs/',
                path + 'results/r-sample/03_taxonomy/',
                path + 'results/r-sample/03_taxonomy/kraken2/',
                path + 'results/r-sample/03_taxonomy/bracken/',
                path + 'results/r-sample/04_function/',
                path + 'results/r-sample/04_function/kraken2/',
                path + 'results/r-sample/04_function/kraken2/annotations/',
                path + 'results/r-sample/04_function/emapper/',
                path + 'results/r-sample/04_function/emapper/annotations/',
                path + 'results/r-project/',
                ]
for item in folder_checks:
    if not os.path.isdir(item): os.mkdir(item)

#reads_function_classify_method=kraken2,emapper
def get_reads_sample_04_function_classify_rule_mk(wildcards):
    if reads_function_classify_method == "kraken2":
        return rules.r_sample_04_function_kraken2.output.mk
    elif reads_function_classify_method == "emapper":
        return rules.r_sample_04_function_emapper.output.mk
    else:
        raise ValueError(f"Unknown reads_function_classify_method: {reads_function_classify_method}")


_sample_cfg = config.get("Params_r_sample", {})
TARGET_READS = _sample_cfg.get("00_sample_reads", {}).get("target_reads")
if TARGET_READS is None:
    TARGET_READS = _sample_cfg.get("target_reads", 0)


rule r_sample_00_sample_reads:
    input:
        r_s = path + 'results/b-sample/01_trim/{sample}_trimmed_s.fastq.gz',
    output:
        r_s = path + 'results/r-sample/00_sampled_trim/{sample}_sampled_s.fastq.gz',
        stats = path + 'results/r-sample/00_sampled_trim/{sample}_sampled.stats'
    params:
        target_reads = TARGET_READS,
        seed = config["Params_r_sample"]["00_sample_reads"]["sample_seed"],
    conda:
        config["conda_envs"]["metalink-main"]
    threads: config["Params_r_sample"]["01_sampled_reads_trim_merge"]["threads"]
    log:
        path + 'logs/r-sample/00_sample_reads/{sample}.sample_reads.log'
    benchmark:
        path + 'benchmarks/r-sample/00_sample_reads/{sample}.sample_reads.benchmark'
    shell:
        """
        raw_count=$(seqkit stat -T {input.r_s} | tail -1 | cut -f 4 | tr -d ',')
        echo "Original reads: $raw_count" > {output.stats}

        if [ {params.target_reads} -gt 0 ] && [ $raw_count -gt {params.target_reads} ]; then
            echo "Downsampling to {params.target_reads} reads (seed={params.seed})" >> {output.stats}
            seqkit sample -s {params.seed} -n {params.target_reads} {input.r_s} -o {output.r_s} 2> {log}
            final_count=$(seqkit stat -T {output.r_s} | tail -1 | cut -f 4 | tr -d ',')
            echo "Sampled reads: $final_count" >> {output.stats}
        else
            echo "No sampling needed (target={params.target_reads}, raw=$raw_count)" >> {output.stats}
            ln -sf $(realpath {input.r_s}) {output.r_s}
        fi
        """

rule r_sample_01_sampled_fastq_to_fasta:
    input:
        r_s = path + 'results/r-sample/00_sampled_trim/{sample}_sampled_s.fastq.gz',
    output:
        fna = temp(path + 'results/r-sample/02_FragGeneScanRs/{sample}_read_s.fna'),
    conda:
        config["conda_envs"]["metalink-main"]
    priority: 1
    threads: config["Params_r_sample"]["01_sampled_reads_trim_merge"]["threads"]
    log:
        path + 'logs/r-sample/01_sampled_fastq_to_fasta/{sample}.log'
    benchmark:
        path + 'benchmarks/r-sample/01_sampled_fastq_to_fasta/{sample}.benchmark'
    shell:
        """
        seqkit fq2fa -j {threads} {input.r_s} > {output.fna} 2> {log}
        """

rule r_sample_02_FragGeneScanRs:
    input:
        fna = path + 'results/r-sample/02_FragGeneScanRs/{sample}_read_s.fna',
    output:
        fna = path + 'results/r-sample/02_FragGeneScanRs/{sample}_read_s.fna.gz',
        ffn = path + 'results/r-sample/02_FragGeneScanRs/{sample}_read_s.ffn.gz',
        faa = path + 'results/r-sample/02_FragGeneScanRs/{sample}_read_s.faa.gz',
        ffn_length = path + 'results/r-sample/02_FragGeneScanRs/{sample}_read_s.ffn_length.gz',
    threads: config["Params_r_sample"]["02_FragGeneScanRs"]["threads"]
    priority: 2
    params:
        model = config["Params_r_sample"]["02_FragGeneScanRs"]["model"],
        ffn = path + 'results/r-sample/02_FragGeneScanRs/{sample}_read_s.ffn',
        faa = path + 'results/r-sample/02_FragGeneScanRs/{sample}_read_s.faa',
        ffn_length = path + 'results/r-sample/02_FragGeneScanRs/{sample}_read_s.ffn_length',
    conda:
        config["conda_envs"]["metalink-main"]
    log:
        path + 'logs/r-sample/02_FragGeneScanRs/{sample}.log'
    benchmark:
        path + 'benchmarks/r-sample/02_FragGeneScanRs/{sample}.benchmark'
    shell:
        """
        FragGeneScanRs -p {threads} -w 0 -t {params.model} \
            -s {input.fna} \
            -n {params.ffn} \
            -a {params.faa} 2> {log}
        seqkit fx2tab -l -n -i -H -j {threads} {params.ffn} > {params.ffn_length} 2>> {log}
        pigz -p {threads} -f {input.fna} 2>> {log}
        pigz -p {threads} -f {params.faa} 2>> {log}
        pigz -p {threads} -f {params.ffn} 2>> {log}
        pigz -p {threads} -f {params.ffn_length} 2>> {log}
        """

rule r_sample_02_compute_FragGeneScanRs_stats:
    input:
        fna = path + 'results/r-sample/02_FragGeneScanRs/{sample}_read_s.fna.gz',
        ffn = path + 'results/r-sample/02_FragGeneScanRs/{sample}_read_s.ffn.gz',
        faa = path + 'results/r-sample/02_FragGeneScanRs/{sample}_read_s.faa.gz',
    output:
        fna = path + 'results/r-sample/02_FragGeneScanRs_stats/{sample}_read_s.fna.stats',
        ffn = path + 'results/r-sample/02_FragGeneScanRs_stats/{sample}_read_s.ffn.stats',
        faa = path + 'results/r-sample/02_FragGeneScanRs_stats/{sample}_read_s.faa.stats',
    conda:
        config["conda_envs"]["metalink-main"]
    threads: config["Params_r_sample"]["02_FragGeneScanRs_stats"]["threads"]
    log:
        path + 'logs/r-sample/02_FragGeneScanRs_stats/{sample}.compute_trimmed_stats.log'
    benchmark:
        path + 'benchmarks/r-sample/02_FragGeneScanRs_stats/{sample}.compute_trimmed_stats.benchmark'
    shell:
        """
        seqkit stat -T -b -j {threads} {input.fna} > {output.fna} 2> {log}
        seqkit stat -T -b -j {threads} {input.ffn} > {output.ffn} 2>> {log}
        seqkit stat -T -b -j {threads} {input.faa} > {output.faa} 2>> {log}
        """

rule r_sample_03_taxonomy_kraken2:
    input:
        r_s = path + 'results/r-sample/00_sampled_trim/{sample}_sampled_s.fastq.gz',
    output:
        kraken_report = path + 'results/r-sample/03_taxonomy/kraken2/{sample}.kraken2.report.tsv',
        detail = path + 'results/r-sample/03_taxonomy/kraken2/{sample}.kraken2.read2classification.tsv.gz',
    threads: config["Params_r_sample"]["03_taxonomy"]["threads"]
    params:
        databases_path = config["Directory"]["database"]["kraken2_taxonomy_db"],
        confidence = config["Params_r_sample"]["03_taxonomy_kraken2"]["confidence"],
        detail = path + 'results/r-sample/03_taxonomy/kraken2/{sample}.kraken2.read2classification.tsv',
        extra_args = config["Params_r_sample"]["03_taxonomy_kraken2"]["extra_args"],
    conda:
        config["conda_envs"]["metalink-main"]
    resources: kraken2_taxonomy=1
    log:
        path + 'logs/r-sample/03_taxonomy/kraken2/{sample}.log'
    benchmark:
        path + 'benchmarks/r-sample/03_taxonomy/kraken2/{sample}.benchmark'
    shell:
        """
        kraken2 {params.extra_args} --threads {threads} \
            --db {params.databases_path} \
            --use-names \
            --confidence {params.confidence} \
            --report {output.kraken_report} \
            {input.r_s} \
            1> {params.detail} 2> {log}
        pigz -f -p {threads} {params.detail} > {output.detail} 2>> {log}
        """

rule r_sample_03_taxonomy_bracken:
    input:
        kraken_report = path + 'results/r-sample/03_taxonomy/kraken2/{sample}.kraken2.report.tsv',
        r_s = path + 'results/r-sample/00_sampled_trim/{sample}_sampled_s.fastq.gz',
    output:
        bracken_report = path + 'results/r-sample/03_taxonomy/bracken/{sample}.bracken.report.tsv',
        bracken_kraken = path + 'results/r-sample/03_taxonomy/bracken/{sample}.bracken_kraken.report.tsv',
    threads: config["Params_r_sample"]["03_taxonomy"]["threads"]
    params:
        databases_path = config["Directory"]["database"]["kraken2_taxonomy_db"],
        read_type = READ_TYPE,
    conda:
        config["conda_envs"]["metalink-main"]
    log:
        path + 'logs/r-sample/03_taxonomy/bracken/{sample}.log'
    benchmark:
        path + 'benchmarks/r-sample/03_taxonomy/bracken/{sample}.benchmark'
    shell:
        """
        if [ "{params.read_type}" = "pacbio-raw" ] || [ "{params.read_type}" = "pacbio-corr" ] || [ "{params.read_type}" = "pacbio-hifi" ] || [ "{params.read_type}" = "nano-raw" ] || [ "{params.read_type}" = "nano-corr" ] || [ "{params.read_type}" = "nano-hq" ]; then
            echo "[Bracken skip] Long-read type {params.read_type}; producing empty Bracken outputs." > {log}
            touch {output.bracken_report} {output.bracken_kraken}
        else
            Length=`seqkit head -n 1000 {input.r_s} \
              |seqkit fx2tab -l | cut -f 2 | awk '{{ if (length > max) max = length }} END {{ print max }}'` 2>> {log}

            function process_number() {{
                local num=${{1:-0}}
                if   (( num <= 62 )); then echo 50
                elif (( num <= 87 )); then echo 75
                elif (( num <= 112 )); then echo 100
                elif (( num <= 137 )); then echo 125
                elif (( num <= 175 )); then echo 150
                elif (( num <= 225 )); then echo 200
                elif (( num <= 275 )); then echo 250
                else echo 300
                fi
            }}
            ScaledLength=$(process_number $Length) 2>> {log}

            bracken \
                -r ${{ScaledLength}} \
                -d {params.databases_path} \
                -t {threads} \
                -i {input.kraken_report} \
                -w {output.bracken_kraken} \
                -o {output.bracken_report} >> {log} 2>&1
        fi
        """

rule r_sample_04_function_kraken2:
    priority: 20
    input:
        ffn = path + 'results/r-sample/02_FragGeneScanRs/{sample}_read_s.ffn.gz',
    output:
        kraken_report = path + 'results/r-sample/04_function/kraken2/annotations/{sample}.kraken2.report.tsv',
        detail = path + 'results/r-sample/04_function/kraken2/annotations/{sample}.kraken2.annotations.gz',
        mk = path + 'markers/r-sample/04_function/kraken2/annotations/{sample}.mk',
    threads: config["Params_r_sample"]["04_function"]["threads"]
    params:
        databases_path = config["Directory"]["database"]["kraken2_function_db"],
        confidence = config["Params_r_sample"]["04_function_kraken2"]["confidence"],
        detail = path + 'results/r-sample/04_function/kraken2/annotations/{sample}.kraken2.read2classification.tsv',
        tmp = path + 'results/r-sample/04_function/kraken2/annotations/{sample}.kraken2.read2classification.tsv.gz',
        extra_args = config["Params_r_sample"]["04_function_kraken2"]["extra_args"],
    conda:
        config["conda_envs"]["metalink-main"]
    log:
        path + 'logs/r-sample/04_function/kraken2/annotations/{sample}.log'
    benchmark:
        path + 'benchmarks/r-sample/04_function/kraken2/annotations/{sample}.benchmark'
    shell:
        """
        kraken2 {params.extra_args} --threads {threads} \
            --db {params.databases_path} \
            --use-names \
            --confidence {params.confidence} \
            --report {output.kraken_report} \
            {input.ffn} 1> {params.detail} 2> {log}
        pigz -f -p {threads} {params.detail} > {params.tmp} 2>> {log}
        mv {params.tmp} {output.detail}
        touch {output.mk}
        """

rule r_sample_04_function_emapper:
    priority: 20
    input:
        faa = path + 'results/r-sample/02_FragGeneScanRs/{sample}_read_s.faa.gz',
    output:
        mk = path + 'markers/r-sample/04_function/emapper/{sample}.mk',
    threads: config["Params_r_sample"]["04_function"]["threads"]
    conda:
        config["conda_envs"]["metalink-annotation"]
    params:
        faa = path + 'results/r-sample/02_FragGeneScanRs/{sample}_read_s.faa',
        seed_ortholog_evalue = config["Params_r_sample"]["04_function_emapper"]["seed_ortholog_evalue"],
        method = config["Params_r_sample"]["04_function_emapper"]["method"],
        query_cover = config["Params_r_sample"]["04_function_emapper"]["query_cover"],
        data_dir = config["Directory"]["database"]["eggnog_db"],
        out = path + 'results/r-sample/04_function/emapper/annotations/{sample}',
        extra_args = config["Params_r_sample"]["04_function_emapper"]["extra_args"],
    log:
        path + 'logs/r-sample/04_function/emapper/{sample}.log'
    benchmark:
        path + 'benchmarks/r-sample/04_function/emapper/{sample}.benchmark'
    shell:
        """
        pigz -d --keep -p {threads} -f {input.faa} > {params.faa} 2> {log}
        emapper.py \
            -i {params.faa} \
            --cpu {threads} \
            --seed_ortholog_evalue {params.seed_ortholog_evalue} \
            -m {params.method} \
            --query_cover {params.query_cover} \
            --data_dir {params.data_dir} \
            -o {params.out} --override {params.extra_args} >> {log} 2>&1
        pigz -p {threads} -r -f {params.out}.emapper.annotations {params.out}.emapper.hits {params.out}.emapper.seed_orthologs 2>> {log}
        rm -f {params.faa}
        touch {output.mk}
        """

rule r_sample_04_function_done:
    priority: 20
    input:
        mk = get_reads_sample_04_function_classify_rule_mk,
    output:
        mk = path + 'markers/r-sample/04_function/classify_done/{sample}.mk'
    shell:
        "cp {input.mk} {output.mk}"

rule r_sample_05_merge_tax_fun_depth:
    priority: 30
    '''
    Use --read_type single as only one kraken2 result would be input to combine_taxonomy_function_depth.r.cmd.R
    '''
    input:
        read_tax_detail = path + 'results/r-sample/03_taxonomy/kraken2/{sample}.kraken2.read2classification.tsv.gz',
        read_ffn_length = path + 'results/r-sample/02_FragGeneScanRs/{sample}_read_s.ffn_length.gz',
        reads_fun_mk = path + 'markers/r-sample/04_function/classify_done/{sample}.mk',
    output:
        summary_Tax = path + 'results/r-sample/05_merge_tax_fun_depth/summary/{sample}.summary_Tax.tsv',
        summary_ORF_KO = path + 'results/r-sample/05_merge_tax_fun_depth/summary/{sample}.summary_ORF_KO.tsv',
        summary_ORF_Tax_KO = path + 'results/r-sample/05_merge_tax_fun_depth/summary/{sample}.summary_ORF_Tax_KO.tsv',
        summary_ORF_Source = path + 'results/r-sample/05_merge_tax_fun_depth/summary/{sample}.summary_ORF_Source.tsv',
        summary_ORF_Tax_Source = path + 'results/r-sample/05_merge_tax_fun_depth/summary/{sample}.summary_ORF_Tax_Source.tsv',
        KO_depth = path + 'results/r-sample/05_merge_tax_fun_depth/KO_depth/{sample}.KO_depth.tsv',
        Taxa_abundance = path + 'results/r-sample/05_merge_tax_fun_depth/Taxa_abundance/{sample}.Taxa_abundance.tsv',
        KO_Taxa_SumCoverage = path + 'results/r-sample/05_merge_tax_fun_depth/KO_Taxa_SumCoverage/{sample}.KO_Taxa_SumCoverage.Taxa.tsv',
        KO_Taxa_RelativeAbundance = path + 'results/r-sample/05_merge_tax_fun_depth/KO_Taxa_RelativeAbundance/{sample}.KO_Taxa_RelativeAbundance.Taxa.tsv',
    params:
        reads_function_classify_method = reads_function_classify_method,
        taxid_tax = taxid_tax_path,
        taxid_fun = taxid_fun_path,
        in_KO_length_ref = KO_length_ref_path,
        read_fun_detail = path + 'results/r-sample/04_function/' + reads_function_classify_method + '/annotations/{sample}.' + reads_function_classify_method + '.annotations.gz',
        merge_detail = path + 'results/r-sample/05_merge_tax_fun_depth/detail_for_each_seq/{sample}.detail_for_each_seq.tsv',
        read_type = "single",
        minORFLen = config["Params_r_sample"]["05_merge_tax_fun_depth"]["minORFLen"],
        KO_mode = config["Params_r_sample"]["05_merge_tax_fun_depth"]["KO_mode"],
        min_KO_num_in_one_taxon = config["Params_r_sample"]["05_merge_tax_fun_depth"]["min_KO_num_in_one_taxon"],
        rank = ",".join(RANKS),
        out_CAMI_format = path + 'results/r-sample/05_merge_tax_fun_depth/Taxa_abundance/{sample}.Taxa_abundance.CAMI_format.tsv',
        CAMI_format_output = config["Params_r_sample"]["05_merge_tax_fun_depth"]["CAMI_format_output"],
        output_detail = config["Params_r_sample"]["05_merge_tax_fun_depth"]["output_detail"],
    threads: config["Params_r_sample"]["05_merge_tax_fun_depth"]["threads"]
    conda:
        config["conda_envs"]["metalink-main"]
    log:
        path + 'logs/r-sample/05_merge_tax_fun_depth/{sample}.log'
    benchmark:
        path + 'benchmarks/r-sample/05_merge_tax_fun_depth/{sample}.benchmark'
    shell:
        """
        echo "Start reads_sample_05_merge_tax_fun_depth for {wildcards.sample}" > {log}
        mkdir -p $(dirname {output.summary_Tax}) 2>{log}
        reads_miss=""
        for f in {input.read_tax_detail} {params.read_fun_detail} {input.read_ffn_length}; do
            [ -s "$f" ] || reads_miss="$reads_miss $f"
        done
        if [ -n "$reads_miss" ]; then
            echo "ERROR: reads inputs missing/empty for {wildcards.sample}; cannot merge taxonomy/function/depth." >> {log}
            echo "  missing/empty reads inputs:$reads_miss" >> {log}
            echo "ERROR [r_sample_05_merge_tax_fun_depth] sample {wildcards.sample}: reads inputs missing/empty; check upstream results. See log: {log}" >&2
            exit 1
        else
            echo "All reads inputs non-empty, running merge script." >> {log}
            Rscript {Dir_script_metalink}/combine_taxonomy_function_depth.r.cmd.R \
                --reads_function_classify_method {params.reads_function_classify_method} \
                --in_read_singleton_tax_detail {input.read_tax_detail} \
                --in_taxid_tax {params.taxid_tax} \
                --in_taxid_fun {params.taxid_fun} \
                --in_KO_length_ref {params.in_KO_length_ref} \
                --in_read_fun_detail {params.read_fun_detail} \
                --in_read_ffn_length {input.read_ffn_length} \
                --sample_name {wildcards.sample} \
                --read_type {params.read_type} \
                --KO_mode {params.KO_mode} \
                --min_KO_num_in_one_taxon {params.min_KO_num_in_one_taxon} \
                --rank {params.rank} \
                --minORFLen {params.minORFLen} \
                --summary_Tax {output.summary_Tax} \
                --summary_ORF_KO {output.summary_ORF_KO} \
                --summary_ORF_Tax_KO {output.summary_ORF_Tax_KO} \
                --summary_ORF_Source {output.summary_ORF_Source} \
                --summary_ORF_Tax_Source {output.summary_ORF_Tax_Source} \
                --out_Taxa_abundance {output.Taxa_abundance} \
                --out_KO_eggnog_depth {output.KO_depth} \
                --out_KO_Taxa_SumCoverage_wide {output.KO_Taxa_SumCoverage} \
                --out_KO_Taxa_RelativeAbundance_wide {output.KO_Taxa_RelativeAbundance} \
                --CAMI_format_output {params.CAMI_format_output} \
                --out_CAMI_format {params.out_CAMI_format} \
                --output_detail {params.output_detail} \
                --detail {params.merge_detail} >>{log} 2>&1
            if [ "{params.output_detail}" = "True" ]; then
                pigz -p {threads} -r -f {params.merge_detail} > {params.merge_detail}.gz 2>> {log}
            fi
        fi
        """

rule r_sample_05_visualize:
    priority: 40
    input:
        summary_Tax = path + 'results/r-sample/05_merge_tax_fun_depth/summary/{sample}.summary_Tax.tsv',
        summary_ORF_KO = path + 'results/r-sample/05_merge_tax_fun_depth/summary/{sample}.summary_ORF_KO.tsv',
        summary_ORF_Tax_KO = path + 'results/r-sample/05_merge_tax_fun_depth/summary/{sample}.summary_ORF_Tax_KO.tsv',
    output:
        mk = path + 'markers/r-sample/05_visualize/{sample}.mk',
    params:
        prefix = path + 'results/r-sample/05_visualization/{sample}',
    threads: 1
    conda:
        config["conda_envs"]["metalink-main"]
    log:
        path + 'logs/r-sample/05_visualize/{sample}.log'
    benchmark:
        path + 'benchmarks/r-sample/05_visualize/{sample}.benchmark'
    shell:
        """
        mkdir -p $(dirname {params.prefix}) 2>>{log}
        missing=""
        for f in {input.summary_Tax} {input.summary_ORF_KO} {input.summary_ORF_Tax_KO}; do
            [ -s "$f" ] || missing="$missing $f"
        done
        if [ -n "$missing" ]; then
            echo "Missing/empty summary inputs:$missing; skip visualization for {wildcards.sample}" >>{log}
        else
            echo "All summary files non-empty, running reads visualization." >>{log}
            Rscript {Dir_script_metalink}/plot_summary.cmd.R \
                --TaxFormat read \
                -t {input.summary_Tax} \
                -k {input.summary_ORF_KO} \
                -a {input.summary_ORF_Tax_KO} \
                -o {params.prefix} \
                -f pdf \
                -s 15 \
                -w 10 -h 6 >>{log} 2>&1
        fi
        touch {output.mk}
        """

rule r_sample_06_module_completeness:
    priority: 40
    input:
        KO_Taxa_RelativeAbundance = path + 'results/r-sample/05_merge_tax_fun_depth/KO_Taxa_RelativeAbundance/{sample}.KO_Taxa_RelativeAbundance.Taxa.tsv',
    output:
        mk = path + 'markers/r-sample/06_module_completeness/{sample}.{rank}.mk',
    params:
        KO_Taxa_RelativeAbundance_ranked = path + 'results/r-sample/05_merge_tax_fun_depth/KO_Taxa_RelativeAbundance/{sample}.KO_Taxa_RelativeAbundance.Taxa.{rank}.tsv',
        KO_Taxa_SumCoverage_ranked = path + 'results/r-sample/05_merge_tax_fun_depth/KO_Taxa_SumCoverage/{sample}.KO_Taxa_SumCoverage.Taxa.{rank}.tsv',
        KO_Taxa_SumCoverage = path + 'results/r-sample/05_merge_tax_fun_depth/KO_Taxa_SumCoverage/{sample}.KO_Taxa_SumCoverage.Taxa.tsv',
        Module_Taxa_RA = path + 'results/r-sample/06_module_completeness/{sample}/{rank}/',
        pathway_split = "TRUE",
        min_KO_num_in_one_taxon = config["Params_r_sample"]["06_mclink_taxon"]["min_KO_num_in_one_taxon"],
        min_mean_KO_RA_in_one_taxon = config["Params_r_sample"]["06_mclink_taxon"]["min_mean_KO_RA_in_one_taxon"],
        completeness_threshold = config["Params_r_sample"]["06_mclink_taxon"]["completeness_threshold"],
        plus_scale_method = config["Params_r_sample"]["06_mclink_taxon"]["plus_scale_method"],
        comma_scale_method = config["Params_r_sample"]["06_mclink_taxon"]["comma_scale_method"],
        TopN = config["Params_r_sample"]["06_mclink_taxon"]["TopN"],
    threads: config["Params_r_sample"]["05_merge_tax_fun_depth"]["threads"]
    conda:
        config["conda_envs"]["metalink-main"]
    log:
        path + 'logs/r-sample/06_module_completeness/{sample}.{rank}.log'
    benchmark:
        path + 'benchmarks/r-sample/06_module_completeness/{sample}.{rank}.benchmark'
    shell:
        """
        mkdir -p {params.Module_Taxa_RA} 2> {log}
        if [ -s {input.KO_Taxa_RelativeAbundance} ]; then
            echo "KO_Taxa_RelativeAbundance non-empty, running mclink." >> {log}
            if [ -f "{params.KO_Taxa_RelativeAbundance_ranked}" ]; then
                inputfile="{params.KO_Taxa_RelativeAbundance_ranked}"
            else
                inputfile="{input.KO_Taxa_RelativeAbundance}"
            fi
            covfile=""
            if [ -f "{params.KO_Taxa_SumCoverage_ranked}" ]; then
                covfile="{params.KO_Taxa_SumCoverage_ranked}"
            elif [ -s "{params.KO_Taxa_SumCoverage}" ]; then
                covfile="{params.KO_Taxa_SumCoverage}"
            fi
            f_arg=""
            [ -n "$covfile" ] && f_arg="-f $covfile"
            bash {Dir_script_metalink}/mclink-scripts/mclink.taxa.sh \
                -c {threads} \
                -i $inputfile \
                $f_arg \
                -k {kegg_modules_ko_info_path} \
                -s {Dir_script_metalink}/mclink-scripts \
                -p {params.pathway_split} \
                -n {params.min_KO_num_in_one_taxon} \
                -a {params.min_mean_KO_RA_in_one_taxon} \
                -t {params.TopN} \
                -e {params.completeness_threshold} \
                -l {params.plus_scale_method} \
                -m {params.comma_scale_method} \
                -o {params.Module_Taxa_RA} >>{log} 2>&1
        else
            echo "KO_Taxa_RelativeAbundance empty/missing for {wildcards.sample}; skip reads_sample_03_module_completeness." >>{log}
        fi
        touch {output.mk}
        """

rule r_project_bracken_taxonomic_merge:
    input:
        bracken_kraken = expand(path + 'results/r-sample/03_taxonomy/bracken/{sample}.bracken_kraken.report.tsv', sample=SAMPLE_INDEX),
    output:
        count_matrix = path + 'results/r-project/01_Taxa_Sample/{project}/bracken_taxonomic_profiles.count.tsv',
    threads: config["Params_r_sample"]["03_taxonomy"]["threads"]
    params:
        bracken_dir = path + 'results/r-sample/03_taxonomy/bracken/',
        suffix = ".bracken_kraken.report.tsv",
    conda:
        config["conda_envs"]["metalink-main"]
    log:
        path + 'logs/r-project/01_bracken_taxonomy_merge/{project}.log'
    benchmark:
        path + 'benchmarks/r-project/01_bracken_taxonomy_merge/{project}.benchmark'
    shell:
        """
        # python to extract seq num rather than relative abundance, generate the bracken_taxonomic_profiles.tsv
        python {Dir_script_metalink}/kraken2_merge.cmd.py \
            --input {params.bracken_dir} \
            --output {output.count_matrix} \
            --suffix {params.suffix} >> {log} 2>&1
        """

rule r_project_kraken2_taxonomic_merge:
    input:
        kraken = expand(path + 'results/r-sample/03_taxonomy/kraken2/{sample}.kraken2.report.tsv', sample=SAMPLE_INDEX),
    output:
        count_matrix = path + 'results/r-project/01_Taxa_Sample/{project}/kraken2_taxonomic_profiles.count.tsv',
    threads: config["Params_r_sample"]["03_taxonomy"]["threads"]
    params:
        kraken2_dir = path + 'results/r-sample/03_taxonomy/kraken2/',
        suffix = ".kraken2.report.tsv",
    conda:
        config["conda_envs"]["metalink-main"]
    log:
        path + 'logs/r-project/01_kraken2_taxonomy_merge/{project}.log'
    benchmark:
        path + 'benchmarks/r-project/01_kraken2_taxonomy_merge/{project}.benchmark'
    shell:
        """
        # python to extract seq num rather than relative abundance, generate the kraken2_taxonomic_profiles.tsv
        python {Dir_script_metalink}/kraken2_merge.cmd.py \
            --input {params.kraken2_dir} \
            --output {output.count_matrix} \
            --suffix {params.suffix} >> {log} 2>&1
        """

rule r_project_kraken2_functional_merge:
    input:
        kraken_report = expand(path + 'results/r-sample/04_function/kraken2/annotations/{sample}.kraken2.report.tsv', sample=SAMPLE_INDEX),
    output:
        count_matrix = path + 'results/r-project/02_KO_Sample/{project}/kraken2_functional_profiles.count.tsv',
    threads: config["Params_r_sample"]["04_function"]["threads"]
    params:
        kraken2_dir = path + 'results/r-sample/04_function/kraken2/annotations/',
        suffix = ".kraken2.report.tsv",
    conda:
        config["conda_envs"]["metalink-main"]
    log:
        path + 'logs/r-project/01_kraken2_function_merge/{project}.log'
    benchmark:
        path + 'benchmarks/r-project/01_kraken2_function_merge/{project}.benchmark'
    shell:
        """
        # python to extract seq num rather than relative abundance, generate the kraken2_functional_profiles.tsv
        python {Dir_script_metalink}/kraken2_merge.cmd.py \
            --input {params.kraken2_dir} \
            --output {output.count_matrix} \
            --suffix {params.suffix} >> {log} 2>&1
        """

rule r_project_00_FragGeneScanRs_stats:
    input:
        faa = expand(path + 'results/r-sample/02_FragGeneScanRs_stats/{sample}_read_s.faa.stats', sample=SAMPLE_INDEX),
        fna = expand(path + 'results/r-sample/02_FragGeneScanRs_stats/{sample}_read_s.fna.stats', sample=SAMPLE_INDEX),
        ffn = expand(path + 'results/r-sample/02_FragGeneScanRs_stats/{sample}_read_s.ffn.stats', sample=SAMPLE_INDEX),
    output:
        faa = temp(path + 'results/r-project/00_FragGeneScanRs_stats/{project}/02_FGS_faa_stats.tsv'),
        fna = temp(path + 'results/r-project/00_FragGeneScanRs_stats/{project}/02_FGS_fna_stats.tsv'),
        ffn = temp(path + 'results/r-project/00_FragGeneScanRs_stats/{project}/02_FGS_ffn_stats.tsv'),
        merged_stats = path + 'results/r-project/00_FragGeneScanRs_stats/{project}/02_FGS_stats.tsv',
    conda:
        config["conda_envs"]["metalink-main"]
    log:
        path + 'logs/r-project/00_FragGeneScanRs_stats/{project}/02_cat_FGS_stats.log',
    benchmark:
        path + 'benchmarks/r-project/00_FragGeneScanRs_stats/{project}/02_cat_FGS_stats.benchmark'
    shell:
        """
        # Process faa.stats files
        for file1 in {input.faa}; do
            if [ -s {output.faa} ]; then
                echo "pass"
            else
                echo -e "ID\tformat_faa\ttype_faa\tnum_seqs_faa\tsum_len_faa\tmin_len_faa\tavg_len_faa\tmax_len_faa" > {output.faa}
            fi
            sed 1d $file1 | awk -F '\\t' 'BEGIN {{OFS="\\t"}} {{gsub("_read_s.faa.gz", "", $1); print}}' >> {output.faa}
        done 2> {log}

        # Process fna.stats files
        for file1 in {input.fna}; do
            if [ -s {output.fna} ]; then
                echo "pass"
            else
                echo -e "ID\tformat_fna\ttype_fna\tnum_seqs_fna\tsum_len_fna\tmin_len_fna\tavg_len_fna\tmax_len_fna" > {output.fna}
            fi
            sed 1d $file1 | awk -F '\\t' 'BEGIN {{OFS="\\t"}} {{gsub("_read_s.fna.gz", "", $1); print}}' >> {output.fna}
        done 2>> {log}

        # Process ffn.stats files
        for file1 in {input.ffn}; do
            if [ -s {output.ffn} ]; then
                echo "pass"
            else
                echo -e "ID\tformat_ffn\ttype_ffn\tnum_seqs_ffn\tsum_len_ffn\tmin_len_ffn\tavg_len_ffn\tmax_len_ffn" > {output.ffn}
            fi
            sed 1d $file1 | awk -F '\\t' 'BEGIN {{OFS="\\t"}} {{gsub("_read_s.ffn.gz", "", $1); print}}' >> {output.ffn}
        done 2>> {log}

        # Merge all three stats files
        python -c '
import pandas as pd
df_fna = pd.read_csv("{output.fna}", sep="\\t", header=0)
df_faa = pd.read_csv("{output.faa}", sep="\\t", header=0)
df_ffn = pd.read_csv("{output.ffn}", sep="\\t", header=0)
merged_df = df_fna.merge(df_faa, on="ID", how="outer").merge(df_ffn, on="ID", how="outer")
merged_df.to_csv("{output.merged_stats}", index=False, sep="\\t")
' 2>> {log}
        """

rule r_project_01_cat_Taxa_abundance_dcast:
    input:
        Taxa_abundance = expand(path + 'results/r-sample/05_merge_tax_fun_depth/Taxa_abundance/{sample}.Taxa_abundance.tsv', sample=SAMPLE_INDEX),
    output:
        out_long = path + 'results/r-project/01_Taxa_Sample/{project}/Taxa_Sample.long.tsv',
    params:
        tmp_dir = path + 'results/r-project/tmp/{project}/',
        prefix_wide = path + 'results/r-project/01_Taxa_Sample/{project}/Taxa_Sample.wide.tsv',
    conda: config["conda_envs"]["metalink-main"]
    threads: config["Params_r_sample"]["05_merge_tax_fun_depth"]["threads"]
    log: path + 'logs/r-project/01_cat_Taxa_abundance_dcast/{project}.log'
    benchmark:
        path + 'benchmarks/r-project/01_cat_Taxa_abundance_dcast/{project}.benchmark'
    shell:
        """
        mkdir -p {params.tmp_dir}
        printf '%s\\n' {input.Taxa_abundance} > {params.tmp_dir}/Taxa_abundance.sample.list.raw 2> {log}

        > {params.tmp_dir}/Taxa_abundance.sample.list
        while IFS= read -r file; do
            if [ -s "$file" ]; then
                printf '%s\\n' "$file" >> {params.tmp_dir}/Taxa_abundance.sample.list
            else
                echo "[WARNING] File missing or empty, skipped: $file" >> {log}
            fi
        done < {params.tmp_dir}/Taxa_abundance.sample.list.raw
        Rscript {Dir_script_metalink}/cat_KO_depth_dcast.cmd.R \
            --sample_tsv_list {params.tmp_dir}/Taxa_abundance.sample.list \
            --id_var Taxonomy \
            --value_var Taxa_ReadCount,Taxa_RelativeAbundance \
            --out_long {output.out_long} \
            --out_wide {params.prefix_wide} \
            --threads {threads} \
            > {log} 2>&1
        """

rule r_project_01_Taxa_visualize:
    input:
        out_long = path + 'results/r-project/01_Taxa_Sample/{project}/Taxa_Sample.long.tsv',
    output:
        mk = path + 'markers/r-project/01_Taxa_visualization/{project}.mk',
    params:
        prefix = path + 'results/r-project/01_Taxa_visualization/{project}',
        prefix_wide = path + 'results/r-project/01_Taxa_Sample/{project}/Taxa_Sample.Taxa_ReadCount.wide.tsv',
        TopN = config["Params_r_project"]["01_Taxa_visualize"]["TopN"],
    threads: 1
    conda:
        config["conda_envs"]["metalink-main"]
    log:
        path + 'logs/r-project/01_Taxa_visualization/{project}.log'
    benchmark:
        path + 'benchmarks/r-project/01_Taxa_visualization/{project}.benchmark'
    shell:
        """
        mkdir -p $(dirname {params.prefix}) 2>>{log}
        missing=""
        for f in {input.out_long} {params.prefix_wide}; do
            [ -s "$f" ] || missing="$missing $f"
        done
        if [ -n "$missing" ]; then
            echo "Missing/empty inputs:$missing; skip reads project_01_Taxa_visualize for {wildcards.project}" >>{log}
        else
            echo "All input files non-empty, running reads project_01_Taxa_visualize." >>{log}
            Rscript {Dir_script_metalink}/plot_bar.cmd.R \
                --Type Taxa \
                --TopN {params.TopN} \
                --Input {params.prefix_wide} \
                --Output {params.prefix} >>{log} 2>&1
        fi
        touch {output.mk}
        """

rule r_project_02_cat_KO_depth_dcast:
    input:
        KO_depth = expand(path + 'results/r-sample/05_merge_tax_fun_depth/KO_depth/{sample}.KO_depth.tsv', sample=SAMPLE_INDEX),
    output:
        out_long = path + 'results/r-project/02_KO_Sample/{project}/KO_Sample.long.tsv',
    params:
        tmp_dir = path + 'results/r-project/tmp/{project}/',
        prefix_wide = path + 'results/r-project/02_KO_Sample/{project}/KO_Sample.wide.tsv',
    threads: config["Params_r_sample"]["05_merge_tax_fun_depth"]["threads"]
    conda:
        config["conda_envs"]["metalink-main"]
    log:
        path + 'logs/r-project/02_cat_KO_depth_dcast/{project}.log'
    benchmark:
        path + 'benchmarks/r-project/02_cat_KO_depth_dcast/{project}.benchmark'
    shell:
        """
        mkdir -p {params.tmp_dir}
        printf '%s\\n' {input.KO_depth} > {params.tmp_dir}/KO_depth.sample.list.raw 2> {log}

        > {params.tmp_dir}/KO_depth.sample.list
        while IFS= read -r file; do
            if [ -s "$file" ]; then
                printf '%s\\n' "$file" >> {params.tmp_dir}/KO_depth.sample.list
            else
                echo "[WARNING] File missing or empty, skipped: $file" >> {log}
            fi
        done < {params.tmp_dir}/KO_depth.sample.list.raw
        Rscript {Dir_script_metalink}/cat_KO_depth_dcast.cmd.R \
            --sample_tsv_list {params.tmp_dir}/KO_depth.sample.list \
            --id_var KO \
            --value_var KO_SumCoverage,KO_RelativeAbundance,KO_ORFNum \
            --out_long {output.out_long} \
            --out_wide {params.prefix_wide} \
            --threads  {threads} \
            > {log} 2>&1
        """

rule r_project_03_community_module:
    input:
        out_long = path + 'results/r-project/02_KO_Sample/{project}/KO_Sample.long.tsv',
    output:
        mk = path + 'markers/r-project/03_community_module/{project}.mk',
    params:
        wide_RelativeAbundance = path + 'results/r-project/02_KO_Sample/{project}/KO_Sample.KO_RelativeAbundance.wide.tsv',
        wide_SumCoverage = path + 'results/r-project/02_KO_Sample/{project}/KO_Sample.KO_SumCoverage.wide.tsv',
        pathway_split = "TRUE",
        community_output = path + 'results/r-project/03_Module_Sample/{project}/',
        min_sample_num_per_KO = config["Params_r_project"]["03_mclink_sample"]["min_sample_num_per_KO"],
        preva = config["Params_r_project"]["03_mclink_sample"]["preva"],
        abund = config["Params_r_project"]["03_mclink_sample"]["abund"],
        completeness_threshold = config["Params_r_project"]["03_mclink_sample"]["completeness_threshold"],
        plus_scale_method = config["Params_r_project"]["03_mclink_sample"]["plus_scale_method"],
        comma_scale_method = config["Params_r_project"]["03_mclink_sample"]["comma_scale_method"],
    threads: config["Params_r_sample"]["05_merge_tax_fun_depth"]["threads"]
    conda:
        config["conda_envs"]["metalink-main"]
    log:
        path + 'logs/r-project/03_community_module/{project}.log'
    benchmark:
        path + 'benchmarks/r-project/03_community_module/{project}.benchmark'
    shell:
        """
        mkdir -p {params.community_output} 2> {log}
        if [ -s {params.wide_RelativeAbundance} ] && [ -s {params.wide_SumCoverage} ]; then
            echo "KO_Sample.KO_RelativeAbundance.wide.tsv and KO_Sample.KO_SumCoverage.wide.tsv non-empty, running mclink with -f." >> {log}
            bash {Dir_script_metalink}/mclink-scripts/mclink.sample.sh \
                -c {threads} \
                -i {params.wide_RelativeAbundance} \
                -f {params.wide_SumCoverage} \
                -k {kegg_modules_ko_info_path} \
                -s {Dir_script_metalink}/mclink-scripts \
                -p {params.pathway_split} \
                -n {params.min_sample_num_per_KO} \
                -r {params.preva} \
                -a {params.abund} \
                -e {params.completeness_threshold} \
                -l {params.plus_scale_method} \
                -m {params.comma_scale_method} \
                -o {params.community_output} >>{log} 2>&1
        elif [ -s {params.wide_RelativeAbundance} ]; then
            echo "KO_Sample.KO_SumCoverage.wide.tsv empty/missing, running mclink without -f." >> {log}
            bash {Dir_script_metalink}/mclink-scripts/mclink.sample.sh \
                -c {threads} \
                -i {params.wide_RelativeAbundance} \
                -k {kegg_modules_ko_info_path} \
                -s {Dir_script_metalink}/mclink-scripts \
                -p {params.pathway_split} \
                -n {params.min_sample_num_per_KO} \
                -r {params.preva} \
                -a {params.abund} \
                -e {params.completeness_threshold} \
                -l {params.plus_scale_method} \
                -m {params.comma_scale_method} \
                -o {params.community_output} >>{log} 2>&1
        else
            echo "KO_Sample.KO_RelativeAbundance.wide.tsv empty/missing for {wildcards.project}; skip reads_project_03_community_module." >>{log}
        fi
        touch {output.mk}
        """

rule r_project_03_Fun_visualize:
    input:
        mk = path + 'markers/r-project/03_community_module/{project}.mk',
    output:
        mk = path + 'markers/r-project/03_Fun_visualization/{project}.mk',
    params:
        prefix = path + 'results/r-project/03_Fun_visualization/{project}',
        prefix_wide = path + 'results/r-project/03_Module_Sample/{project}/SumCoverage/mc_abundance_filtered_by_completeness.tsv',
        TopN = config["Params_r_project"]["03_Fun_visualize"]["TopN"],
    threads: 1
    conda:
        config["conda_envs"]["metalink-main"]
    log:
        path + 'logs/r-project/03_Fun_visualization/{project}.log'
    benchmark:
        path + 'benchmarks/r-project/03_Fun_visualization/{project}.benchmark'
    shell:
        """
        mkdir -p $(dirname {params.prefix}) 2>>{log}
        missing=""
        for f in {params.prefix_wide}; do
            [ -s "$f" ] || missing="$missing $f"
        done
        if [ -n "$missing" ]; then
            echo "Missing/empty inputs:$missing; skip reads project_03_Fun_visualize for {wildcards.project}" >>{log}
        else
            echo "All input files non-empty, running reads project_03_Fun_visualize." >>{log}
            Rscript {Dir_script_metalink}/plot_bar.cmd.R \
                --Type Module \
                --TopN {params.TopN} \
                --Input {params.prefix_wide} \
                --Output {params.prefix} >>{log} 2>&1
        fi
        touch {output.mk}
        """

rule r_project_05_merge_tax_KO_depth:
    input:
        ko_depth = expand(path + 'results/r-sample/05_merge_tax_fun_depth/KO_depth/{sample}.KO_depth.tsv', sample=SAMPLE_INDEX),
    output:
        out_long = path + 'results/r-project/05_merge_tax_KO_depth/{project}/{rank}.mean.long.tsv',
    params:
        detail = lambda wc: expand(path + 'results/r-sample/05_merge_tax_fun_depth/detail_for_each_seq/{sample}.detail_for_each_{rank}.tsv', sample=SAMPLE_INDEX, rank=wc.rank),  # 只展开当前 rank
        tmp_dir = path + 'results/r-project/tmp/{project}/',
        prefix_wide = path + 'results/r-project/05_merge_tax_KO_depth/{project}/{rank}.wide.tsv',
        preva = config["Params_r_project"]["05_merge_tax_KO_depth"]["preva"],
        abund = config["Params_r_project"]["05_merge_tax_KO_depth"]["abund"],
    threads: config["Params_r_sample"]["05_merge_tax_fun_depth"]["threads"]
    conda:
        config["conda_envs"]["metalink-main"]
    log:
        path + 'logs/r-project/05_merge_tax_KO_depth/{project}/{rank}.log'
    benchmark:
        path + 'benchmarks/r-project/05_merge_tax_KO_depth/{project}/{rank}.benchmark'
    shell:
        """
        mkdir -p {params.tmp_dir}
        printf '%s\\n' {params.detail} > {params.tmp_dir}/detail_for_each_{wildcards.rank}.sample.list.raw 2> {log}

        > {params.tmp_dir}/detail_for_each_{wildcards.rank}.sample.list
        while IFS= read -r file; do
            if [ -s "$file" ]; then
                printf '%s\\n' "$file" >> {params.tmp_dir}/detail_for_each_{wildcards.rank}.sample.list
            else
                echo "[WARNING] File missing or empty, skipped: $file" >> {log}
            fi
        done < {params.tmp_dir}/detail_for_each_{wildcards.rank}.sample.list.raw
        Rscript {Dir_script_metalink}/cat_KO_depth_mean.cmd.R\
            --sample_tsv_list {params.tmp_dir}/detail_for_each_{wildcards.rank}.sample.list \
            --id_var KO,{wildcards.rank} \
            --value_var ORF_Num,Taxa_KO_SumCoverage,Taxa_KO_RelativeAbundance \
            --out_long {output.out_long} \
            --out_wide {params.prefix_wide} \
            --min_prevalence {params.preva} \
            --min_abundance {params.abund} \
            --threads  {threads} \
            > {log} 2>&1
        """

rule r_project_05_merge_tax_Module_RelativeAbundance:
    input:
        mk = lambda wc: expand(path + 'markers/r-sample/06_module_completeness/{sample}.{rank}.mk', sample=SAMPLE_INDEX, rank=wc.rank),
    output:
        out_long = path + 'results/r-project/05_merge_tax_Module_RelativeAbundance/{project}/{rank}.mean.long.tsv',
    params:
        detail = lambda wc: expand(path + 'results/r-sample/06_module_completeness/{sample}/{rank}/RelativeAbundance/mc_abundance.tsv', sample=SAMPLE_INDEX, rank=wc.rank),
        prefix_wide = path + 'results/r-project/05_merge_tax_Module_RelativeAbundance/{project}/{rank}.wide.tsv',
        tmp_dir = path + 'results/r-project/tmp/{project}/',
        preva = config["Params_r_project"]["05_merge_tax_Module_RelativeAbundance"]["preva"],
        abund = config["Params_r_project"]["05_merge_tax_Module_RelativeAbundance"]["abund"],
    threads: config["Params_r_sample"]["05_merge_tax_fun_depth"]["threads"]
    conda:
        config["conda_envs"]["metalink-main"]
    log:
        path + 'logs/r-project/05_merge_tax_Module_RelativeAbundance/{project}/{rank}.log'
    benchmark:
        path + 'benchmarks/r-project/05_merge_tax_Module_RelativeAbundance/{project}/{rank}.benchmark'
    shell:
        """
        mkdir -p {params.tmp_dir}
        printf '%s\n' {params.detail} > {params.tmp_dir}/mc_abundance_RA_{wildcards.rank}.sample.list.raw 2> {log}

        > {params.tmp_dir}/mc_abundance_RA_{wildcards.rank}.sample.list
        while IFS= read -r file; do
            if [ -s "$file" ]; then
                printf '%s\n' "$file" >> {params.tmp_dir}/mc_abundance_RA_{wildcards.rank}.sample.list
            else
                echo "[WARNING] File missing or empty, skipped: $file" >> {log}
            fi
        done < {params.tmp_dir}/mc_abundance_RA_{wildcards.rank}.sample.list.raw

        Rscript {Dir_script_metalink}/cat_Module_depth_mean.cmd.R\
            --sample_tsv_list {params.tmp_dir}/mc_abundance_RA_{wildcards.rank}.sample.list \
            --id_var Module_Entry,{wildcards.rank} \
            --value_var Taxa_Module_RelativeAbundance \
            --out_long {output.out_long} \
            --out_wide {params.prefix_wide} \
            --min_prevalence {params.preva} \
            --min_abundance {params.abund} \
            --threads {threads} \
            >> {log} 2>&1
        """

rule r_project_05_merge_tax_Module_SumCoverage:
    input:
        mk = lambda wc: expand(path + 'markers/r-sample/06_module_completeness/{sample}.{rank}.mk', sample=SAMPLE_INDEX, rank=wc.rank),
    output:
        out_long = path + 'results/r-project/05_merge_tax_Module_SumCoverage/{project}/{rank}.mean.long.tsv',
    params:
        detail = lambda wc: expand(path + 'results/r-sample/06_module_completeness/{sample}/{rank}/SumCoverage/mc_abundance.tsv', sample=SAMPLE_INDEX, rank=wc.rank),
        prefix_wide = path + 'results/r-project/05_merge_tax_Module_SumCoverage/{project}/{rank}.wide.tsv',
        tmp_dir = path + 'results/r-project/tmp/{project}/',
        preva = config["Params_r_project"]["05_merge_tax_Module_SumCoverage"]["preva"],
        abund = config["Params_r_project"]["05_merge_tax_Module_SumCoverage"]["abund"],
    threads: config["Params_r_sample"]["05_merge_tax_fun_depth"]["threads"]
    conda:
        config["conda_envs"]["metalink-main"]
    log:
        path + 'logs/r-project/05_merge_tax_Module_SumCoverage/{project}/{rank}.log'
    benchmark:
        path + 'benchmarks/r-project/05_merge_tax_Module_SumCoverage/{project}/{rank}.benchmark'
    shell:
        """
        mkdir -p {params.tmp_dir}
        printf '%s\n' {params.detail} > {params.tmp_dir}/mc_abundance_SC_{wildcards.rank}.sample.list.raw 2> {log}

        > {params.tmp_dir}/mc_abundance_SC_{wildcards.rank}.sample.list
        while IFS= read -r file; do
            if [ -s "$file" ]; then
                printf '%s\n' "$file" >> {params.tmp_dir}/mc_abundance_SC_{wildcards.rank}.sample.list
            else
                echo "[WARNING] File missing or empty, skipped: $file" >> {log}
            fi
        done < {params.tmp_dir}/mc_abundance_SC_{wildcards.rank}.sample.list.raw

        Rscript {Dir_script_metalink}/cat_Module_depth_mean.cmd.R\
            --sample_tsv_list {params.tmp_dir}/mc_abundance_SC_{wildcards.rank}.sample.list \
            --id_var Module_Entry,{wildcards.rank} \
            --value_var Taxa_Module_SumCoverage \
            --out_long {output.out_long} \
            --out_wide {params.prefix_wide} \
            --min_prevalence {params.preva} \
            --min_abundance {params.abund} \
            --threads {threads} \
            >> {log} 2>&1
        """

rule r_project_06_module_completeness:
    input:
        KO_Taxa_RelativeAbundance = path + 'results/r-project/05_merge_tax_KO_depth/{project}/{rank}.mean.long.tsv',
    output:
        mk = path + 'markers/r-project/06_module_completeness/{project}/{rank}.mk',
    params:
        KO_Taxa_RelativeAbundance_ranked = path + 'results/r-project/05_merge_tax_KO_depth/{project}/{rank}.Taxa_KO_RelativeAbundance.mean.wide.tsv',
        KO_Taxa_SumCoverage_ranked = path + 'results/r-project/05_merge_tax_KO_depth/{project}/{rank}.Taxa_KO_SumCoverage.mean.wide.tsv',
        Module_Taxa_RA = path + 'results/r-project/06_module_completeness/{project}/{rank}/',
        pathway_split = "TRUE",
        min_KO_num_in_one_taxon = config["Params_r_project"]["06_mclink_taxon"]["min_KO_num_in_one_taxon"],
        min_mean_KO_RA_in_one_taxon = config["Params_r_project"]["06_mclink_taxon"]["min_mean_KO_RA_in_one_taxon"],
        completeness_threshold = config["Params_r_project"]["06_mclink_taxon"]["completeness_threshold"],
        plus_scale_method = config["Params_r_project"]["06_mclink_taxon"]["plus_scale_method"],
        comma_scale_method = config["Params_r_project"]["06_mclink_taxon"]["comma_scale_method"],
        TopN = config["Params_r_project"]["06_mclink_taxon"]["TopN"],
    threads: config["Params_r_sample"]["05_merge_tax_fun_depth"]["threads"]
    conda:
        config["conda_envs"]["metalink-main"]
    log:
        path + 'logs/r-project/06_module_completeness/{project}/{rank}.log'
    benchmark:
        path + 'benchmarks/r-project/06_module_completeness/{project}/{rank}.benchmark'
    shell:
        """
        mkdir -p {params.Module_Taxa_RA} 2> {log}
        if [ -s {input.KO_Taxa_RelativeAbundance} ]; then
            echo "KO_Taxa_RelativeAbundance non-empty, running mclink." >> {log}
            if [ -f "{params.KO_Taxa_RelativeAbundance_ranked}" ]; then
                inputfile="{params.KO_Taxa_RelativeAbundance_ranked}"
            else
                inputfile="{input.KO_Taxa_RelativeAbundance}"
            fi
            f_arg=""
            if [ -s "{params.KO_Taxa_SumCoverage_ranked}" ]; then
                f_arg="-f {params.KO_Taxa_SumCoverage_ranked}"
            fi
            bash {Dir_script_metalink}/mclink-scripts/mclink.taxa.sh \
                -c {threads} \
                -i $inputfile \
                $f_arg \
                -k {kegg_modules_ko_info_path} \
                -s {Dir_script_metalink}/mclink-scripts \
                -p {params.pathway_split} \
                -n {params.min_KO_num_in_one_taxon} \
                -a {params.min_mean_KO_RA_in_one_taxon} \
                -t {params.TopN} \
                -e {params.completeness_threshold} \
                -l {params.plus_scale_method} \
                -m {params.comma_scale_method} \
                -o {params.Module_Taxa_RA} >>{log} 2>&1
        else
            echo "KO_Taxa_RelativeAbundance empty/missing for {wildcards.project}; skip reads_project_module_completeness." >>{log}
        fi
        touch {output.mk}
        """
