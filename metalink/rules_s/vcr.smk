'''
/*
 * @Date: 2022-12-04 17:15:47
 * @LastEditors: LiuyangLi
 * @LastEditTime: 2024-07-20 13:03:29
 * @FilePath: vcr.smk
 * @Description:
 */
'''

import os
import sys

Dir_script_metalink = config["Directory"]["script"]["metalink"]
kegg_modules_ko_info_path = config["Directory"]["script"]["kegg_modules_ko_info_path"]
taxid_tax_path = config["Directory"]["database"]["taxid_tax_path"]
taxid_fun_path = config["Directory"]["database"]["taxid_fun_path"]
KO_length_ref_path = config["Directory"]["script"]["KO_length_ref_path"]
contigs_taxonomy_classify_method = config["Params_c_sample"]["01_taxonomy"]["contigs_taxonomy_classify_method"]
contigs_function_classify_method = config["Params_c_sample"]["02_function"]["contigs_function_classify_method"]


folder_checks = [
                path + 'results/',
                path + 'results/vcr-sample/',
                path + 'results/vcr-sample/04_function/',
                path + 'results/vcr-sample/04_function/emapper/',
                path + 'results/vcr-sample/04_function/emapper/annotations',
                path + 'results/vcr-project/',
                ]
for item in folder_checks:
    if not os.path.isdir(item): os.mkdir(item)


base_sample_mapping_method = config["Params_b_sample"]["03_mapping"]["method"]
#base_sample_mapping_method=bwamem,bwamem2,strobealign,minibwa
reads_function_classify_method = config["Params_vcr_sample"]["04_function"]["reads_function_classify_method"]
#reads_function_classify_method=kraken2,emapper
def get_vcr_sample_04_function_classify_rule_mk(wildcards):
    if reads_function_classify_method == "kraken2":
        return rules.vcr_sample_04_function_kraken2.output.mk
    elif reads_function_classify_method == "emapper":
        return rules.vcr_sample_04_function_emapper.output.mk
    else:
        raise ValueError(f"Unknown reads_function_classify_method: {reads_function_classify_method}")

_long_read_types = ['pacbio-raw', 'pacbio-corr', 'pacbio-hifi', 'nano-raw', 'nano-corr', 'nano-hq']

rule vcr_sample_01_unmapped_reads_fastq_to_fasta:
    input:
        stat_raw = path + 'results/b-sample/04_bam_summary/{sample}.contigs_bam_stats.tsv',
    output:
        unmapped_r_s = path + 'results/vcr-sample/01_unmapped_reads/{sample}_unmapped_s.fastq.gz',
        fna = temp(path + 'results/vcr-sample/02_FragGeneScanRs/{sample}_unmapped_s.fna'),
    threads: config["Params_vcr_sample"]["01_unmapped_reads_trim_merge"]["threads"]
    conda:
        config["conda_envs"]["metalink-main"]
    params:
        bam = path + 'results/b-sample/03_mapping/' + base_sample_mapping_method + '/{sample}.bam',
    priority: 1
    log:
        path + 'logs/vcr-sample/01_unmapped_reads/{sample}.log'
    benchmark:
        path + 'benchmarks/vcr-sample/01_unmapped_reads/{sample}.benchmark'
    shell:
        """
        samtools view -@ {threads} -f 4 -F 256 -b {params.bam} | \
            samtools sort -@ {threads} -n -o - | \
            samtools fastq -@ {threads} > {output.unmapped_r_s} 2>> {log}
        if [ -s {output.unmapped_r_s} ]; then
            seqkit fq2fa -j {threads} {output.unmapped_r_s} > {output.fna} 2>> {log}
        else
            echo "unmapped reads empty/missing for {wildcards.sample}; skip vcr_sample_01_unmapped_reads_fastq_to_fasta." >>{log}
            touch {output.fna}
        fi
        """

if READ_TYPE in _long_read_types:
    rule vcr_sample_02_FragGeneScanRs_3rd:
        input:
            fna = path + 'results/vcr-sample/02_FragGeneScanRs/{sample}_unmapped_s.fna',
        output:
            fna = path + 'results/vcr-sample/02_FragGeneScanRs/{sample}_unmapped_s.fna.gz',
            ffn = path + 'results/vcr-sample/02_FragGeneScanRs/{sample}_unmapped_s.ffn.gz',
            faa = path + 'results/vcr-sample/02_FragGeneScanRs/{sample}_unmapped_s.faa.gz',
            fna_length = path + 'results/vcr-sample/02_FragGeneScanRs/{sample}_unmapped_s.fna_length.gz',
            ffn_length = path + 'results/vcr-sample/02_FragGeneScanRs/{sample}_unmapped_s.ffn_length.gz',
        threads: config["Params_vcr_sample"]["02_FragGeneScanRs"]["threads"]
        priority: 2
        params:
            model = config["Params_vcr_sample"]["02_FragGeneScanRs"]["model"],
            ffn = path + 'results/vcr-sample/02_FragGeneScanRs/{sample}_unmapped_s.ffn',
            faa = path + 'results/vcr-sample/02_FragGeneScanRs/{sample}_unmapped_s.faa',
            fna_length = path + 'results/vcr-sample/02_FragGeneScanRs/{sample}_unmapped_s.fna_length',
            ffn_length = path + 'results/vcr-sample/02_FragGeneScanRs/{sample}_unmapped_s.ffn_length',
        conda:
            config["conda_envs"]["metalink-main"]
        log:
            path + 'logs/vcr-sample/02_FragGeneScanRs/{sample}.log'
        benchmark:
            path + 'benchmarks/vcr-sample/02_FragGeneScanRs/{sample}.benchmark'
        shell:
            """
            if [ -s {input.fna} ]; then
                FragGeneScanRs -p {threads} -w 0 -t {params.model} \
                    -s {input.fna} \
                    -n {params.ffn} \
                    -a {params.faa} 2> {log}
                seqkit fx2tab -l -n -i -H -j {threads} {input.fna} > {params.fna_length} 2>> {log}
                seqkit fx2tab -l -n -i -H -j {threads} {params.ffn} > {params.ffn_length} 2>> {log}
                pigz -p {threads} -f {input.fna} 2>> {log}
                pigz -p {threads} -f {params.faa} 2>> {log}
                pigz -p {threads} -f {params.ffn} 2>> {log}
                pigz -p {threads} -f {params.fna_length} 2>> {log}
                pigz -p {threads} -f {params.ffn_length} 2>> {log}
            else
                echo "unmapped_s.fna empty/missing for {wildcards.sample}; skip vcr_sample_02_FragGeneScanRs." >>{log}
                touch {output.fna} {output.ffn} {output.faa} {output.fna_length} {output.ffn_length}
            fi
            """
else:
    rule vcr_sample_02_FragGeneScanRs_2nd:
        input:
            fna = path + 'results/vcr-sample/02_FragGeneScanRs/{sample}_unmapped_s.fna',
        output:
            fna = path + 'results/vcr-sample/02_FragGeneScanRs/{sample}_unmapped_s.fna.gz',
            ffn = path + 'results/vcr-sample/02_FragGeneScanRs/{sample}_unmapped_s.ffn.gz',
            faa = path + 'results/vcr-sample/02_FragGeneScanRs/{sample}_unmapped_s.faa.gz',
            ffn_length = path + 'results/vcr-sample/02_FragGeneScanRs/{sample}_unmapped_s.ffn_length.gz',
        threads: config["Params_vcr_sample"]["02_FragGeneScanRs"]["threads"]
        priority: 2
        params:
            model = config["Params_vcr_sample"]["02_FragGeneScanRs"]["model"],
            ffn = path + 'results/vcr-sample/02_FragGeneScanRs/{sample}_unmapped_s.ffn',
            faa = path + 'results/vcr-sample/02_FragGeneScanRs/{sample}_unmapped_s.faa',
            ffn_length = path + 'results/vcr-sample/02_FragGeneScanRs/{sample}_unmapped_s.ffn_length',
        conda:
            config["conda_envs"]["metalink-main"]
        log:
            path + 'logs/vcr-sample/02_FragGeneScanRs/{sample}.log'
        benchmark:
            path + 'benchmarks/vcr-sample/02_FragGeneScanRs/{sample}.benchmark'
        shell:
            """
            if [ -s {input.fna} ]; then
                FragGeneScanRs -p {threads} -w 0 -t {params.model} \
                    -s {input.fna} \
                    -n {params.ffn} \
                    -a {params.faa} 2> {log}
                seqkit fx2tab -l -n -i -H -j {threads} {params.ffn} > {params.ffn_length} 2>> {log}
                pigz -p {threads} -f {input.fna} 2>> {log}
                pigz -p {threads} -f {params.faa} 2>> {log}
                pigz -p {threads} -f {params.ffn} 2>> {log}
                pigz -p {threads} -f {params.ffn_length} 2>> {log}
            else
                echo "unmapped_s.fna empty/missing for {wildcards.sample}; skip vcr_sample_02_FragGeneScanRs." >>{log}
                touch {output.fna} {output.ffn} {output.faa} {output.ffn_length}
            fi
            """

rule vcr_sample_02_FragGeneScanRs_stats:
    input:
        fna = path + 'results/vcr-sample/02_FragGeneScanRs/{sample}_unmapped_s.fna.gz',
        ffn = path + 'results/vcr-sample/02_FragGeneScanRs/{sample}_unmapped_s.ffn.gz',
        faa = path + 'results/vcr-sample/02_FragGeneScanRs/{sample}_unmapped_s.faa.gz',
    output:
        fna = path + 'results/vcr-sample/02_FragGeneScanRs_stats/{sample}_unmapped_s.fna.stats',
        ffn = path + 'results/vcr-sample/02_FragGeneScanRs_stats/{sample}_unmapped_s.ffn.stats',
        faa = path + 'results/vcr-sample/02_FragGeneScanRs_stats/{sample}_unmapped_s.faa.stats',
    conda:
        config["conda_envs"]["metalink-main"]
    threads: config["Params_vcr_sample"]["02_FragGeneScanRs_stats"]["threads"]
    log:
        path + 'logs/vcr-sample/02_FragGeneScanRs_stats/{sample}.compute_trimmed_stats.log'
    benchmark:
        path + 'benchmarks/vcr-sample/02_FragGeneScanRs_stats/{sample}.compute_trimmed_stats.benchmark'
    shell:
        """
        seqkit stat -T -b -j {threads} {input.fna} > {output.fna} 2> {log}
        seqkit stat -T -b -j {threads} {input.ffn} > {output.ffn} 2>> {log}
        seqkit stat -T -b -j {threads} {input.faa} > {output.faa} 2>> {log}
        """

rule vcr_sample_04_function_kraken2:
    priority: 20
    input:
        ffn = path + 'results/vcr-sample/02_FragGeneScanRs/{sample}_unmapped_s.ffn.gz',
    output:
        kraken_report = path + 'results/vcr-sample/04_function/kraken2/{sample}.kraken2.report.tsv',
        detail = path + 'results/vcr-sample/04_function/kraken2/{sample}.kraken2.annotations.gz',
        mk = path + 'markers/vcr-sample/04_function/kraken2/{sample}.mk',
    threads: config["Params_vcr_sample"]["04_function"]["threads"]
    params:
        databases_path = config["Directory"]["database"]["kraken2_function_db"],
        confidence = config["Params_vcr_sample"]["04_function_kraken2"]["confidence"],
        detail = path + 'results/vcr-sample/04_function/kraken2/{sample}.kraken2.read2classification.tsv',
        tmp = path + 'results/vcr-sample/04_function/kraken2/{sample}.kraken2.read2classification.tsv.gz',
        extra_args = config["Params_vcr_sample"]["04_function_kraken2"]["extra_args"],
    conda:
        config["conda_envs"]["metalink-main"]
    log:
        path + 'logs/vcr-sample/04_function/kraken2/{sample}.log'
    benchmark:
        path + 'benchmarks/vcr-sample/04_function/kraken2/{sample}.benchmark'
    shell:
        """
        if [ -s {input.ffn} ]; then
            kraken2 {params.extra_args} --threads {threads} \
                --db {params.databases_path} \
                --use-names \
                --confidence {params.confidence} \
                --report {output.kraken_report} \
                {input.ffn} 1> {params.detail} 2> {log}
            pigz -f -p {threads} {params.detail} > {params.tmp} 2>> {log}
            mv {params.tmp} {output.detail}
        else
            echo "unmapped_s.ffn.gz empty/missing for {wildcards.sample}; skip vcr_sample_04_function_kraken2." >>{log}
            touch {output.kraken_report}
            echo "" | pigz -c > {output.detail}
        fi
        touch {output.mk}
        """

rule vcr_sample_04_function_emapper:
    priority: 20
    input:
        faa = path + 'results/vcr-sample/02_FragGeneScanRs/{sample}_unmapped_s.faa.gz',
    output:
        faa = temp(path + 'results/vcr-sample/02_FragGeneScanRs/{sample}_unmapped_s.faa'),
        mk = path + 'markers/vcr-sample/04_function/emapper/{sample}.mk',
    threads: config["Params_vcr_sample"]["04_function"]["threads"]
    conda:
        config["conda_envs"]["metalink-annotation"]
    params:
        seed_ortholog_evalue = config["Params_vcr_sample"]["04_function_emapper"]["seed_ortholog_evalue"],
        method = config["Params_vcr_sample"]["04_function_emapper"]["method"],
        query_cover = config["Params_vcr_sample"]["04_function_emapper"]["query_cover"],
        data_dir = config["Directory"]["database"]["eggnog_db"],
        out = path + 'results/vcr-sample/04_function/emapper/annotations/{sample}',
        extra_args = config["Params_vcr_sample"]["04_function_emapper"]["extra_args"],
    log:
        path + 'logs/vcr-sample/04_function/emapper/{sample}.log'
    benchmark:
        path + 'benchmarks/vcr-sample/04_function/emapper/{sample}.benchmark'
    shell:
        """
        if [ -s {input.faa} ]; then
            pigz -d --keep -p {threads} -f {input.faa} > {output.faa} 2> {log}
            emapper.py \
                -i {output.faa} \
                --cpu {threads} \
                --seed_ortholog_evalue {params.seed_ortholog_evalue} \
                -m {params.method} \
                --query_cover {params.query_cover} \
                --data_dir {params.data_dir} \
                -o {params.out} --override {params.extra_args} >> {log} 2>&1
            pigz -p {threads} -r -f {params.out}.emapper.annotations {params.out}.emapper.hits {params.out}.emapper.seed_orthologs 2>> {log}
        else
            echo "unmapped_s.faa.gz empty/missing for {wildcards.sample}; skip vcr_sample_04_function_emapper." >>{log}
            touch {output.faa}
        fi
        touch {output.mk}
        """

rule vcr_sample_04_function_done:
    priority: 20
    input:
        mk = get_vcr_sample_04_function_classify_rule_mk,
    output:
        mk = path + 'markers/vcr-sample/04_function/classify_done/{sample}.mk'
    shell:
        "cp {input.mk} {output.mk}"

rule vcr_sample_03_taxonomy_kraken2_singleton:
    priority: 10
    resources: kraken2_taxonomy=1
    input:
        unmapped_s  = path + 'results/vcr-sample/01_unmapped_reads/{sample}_unmapped_s.fastq.gz',
    output:
        kraken_report = path + 'results/vcr-sample/03_taxonomy/kraken2/unmapped_reads_singleton/{sample}.kraken2.report.tsv',
        detail = path + 'results/vcr-sample/03_taxonomy/kraken2/unmapped_reads_singleton/{sample}.kraken2.read2classification.tsv.gz',
    threads: config["Params_vcr_sample"]["03_taxonomy"]["threads"]
    params:
        databases_path = config["Directory"]["database"]["kraken2_taxonomy_db"],
        confidence = config["Params_vcr_sample"]["03_taxonomy_kraken2"]["confidence"],
        detail = path + 'results/vcr-sample/03_taxonomy/kraken2/unmapped_reads_singleton/{sample}.kraken2.read2classification.tsv',
        extra_args = config["Params_vcr_sample"]["03_taxonomy_kraken2"]["extra_args"],
    conda:
        config["conda_envs"]["metalink-main"]
    log:
        path + 'logs/vcr-sample/03_taxonomy/kraken2/unmapped_reads_singleton/{sample}.log'
    benchmark:
        path + 'benchmarks/vcr-sample/03_taxonomy/kraken2/unmapped_reads_singleton/{sample}.benchmark'
    shell:
        """
        if [ -s {input.unmapped_s} ]; then
            kraken2 {params.extra_args} --threads {threads} \
                --db {params.databases_path} \
                --use-names \
                --confidence {params.confidence} \
                --report {output.kraken_report} \
                {input.unmapped_s} 1> {params.detail} 2> {log}
            pigz -f -p {threads} {params.detail} > {output.detail} 2>> {log}
        else
            echo "unmapped_s.fastq.gz empty/missing for {wildcards.sample}; skip vcr_sample_03_taxonomy_kraken2_singleton." >>{log}
            touch {output.kraken_report}
            echo "" | pigz -c > {output.detail}
        fi
        """

rule vcr_sample_05_merge_tax_fun_depth:
    priority: 30
    input:
        contigs_tax_mk = path + 'markers/c-sample/01_taxonomy/classify_done/{sample}.mk',
        contigs_fun_mk = path + 'markers/c-sample/02_function/classify_done/{sample}.mk',
        contigs_depth = path + 'results/b-sample/07_contigs/contigs_filter_depth/{sample}',
        contigs_ffn_length = path + 'results/b-sample/07_contigs/contigs_ffn_length/{sample}',
        reads_singleton_tax_detail = path + 'results/vcr-sample/03_taxonomy/kraken2/unmapped_reads_singleton/{sample}.kraken2.read2classification.tsv.gz',
        reads_ffn_length = path + 'results/vcr-sample/02_FragGeneScanRs/{sample}_unmapped_s.ffn_length.gz',
        reads_length = path + 'results/b-sample/01_trim_stats/{sample}_trimmed_s.fastq.stats',
        reads_fun_mk = path + 'markers/vcr-sample/04_function/classify_done/{sample}.mk',
        contigs_ORF_loci = path + 'results/b-sample/07_contigs/contigs_ORF_loci/{sample}',
        vMAGs_summar_mk = path + 'markers/v-sample/10_summary/{sample}.mk',
    output:
        summary_Tax = path + 'results/vcr-sample/05_merge_tax_fun_depth/summary/{sample}.summary_Tax.tsv',
        summary_ORF_KO = path + 'results/vcr-sample/05_merge_tax_fun_depth/summary/{sample}.summary_ORF_KO.tsv',
        summary_ORF_Tax_KO = path + 'results/vcr-sample/05_merge_tax_fun_depth/summary/{sample}.summary_ORF_Tax_KO.tsv',
        summary_ORF_Source = path + 'results/vcr-sample/05_merge_tax_fun_depth/summary/{sample}.summary_ORF_Source.tsv',
        summary_ORF_Tax_Source = path + 'results/vcr-sample/05_merge_tax_fun_depth/summary/{sample}.summary_ORF_Tax_Source.tsv',
        KO_depth = path + 'results/vcr-sample/05_merge_tax_fun_depth/KO_depth/{sample}.KO_depth.tsv',
        Taxa_abundance = path + 'results/vcr-sample/05_merge_tax_fun_depth/Taxa_abundance/{sample}.Taxa_abundance.tsv',
        KO_Taxa_SumCoverage = path + 'results/vcr-sample/05_merge_tax_fun_depth/KO_Taxa_SumCoverage/{sample}.KO_Taxa_SumCoverage.Taxa.tsv',
        KO_Taxa_RelativeAbundance = path + 'results/vcr-sample/05_merge_tax_fun_depth/KO_Taxa_RelativeAbundance/{sample}.KO_Taxa_RelativeAbundance.Taxa.tsv',
        mode_flag = path + 'results/vcr-sample/05_merge_tax_fun_depth/mode/{sample}.mode.tsv',
    params:
        taxid_tax = taxid_tax_path,
        taxid_fun = taxid_fun_path,
        in_KO_length_ref = KO_length_ref_path,
        contigs_taxonomy_classify_method = contigs_taxonomy_classify_method,
        contigs_function_classify_method = contigs_function_classify_method,
        reads_function_classify_method = reads_function_classify_method,
        contigs_tax_detail = path + 'results/c-sample/01_taxonomy/' + contigs_taxonomy_classify_method + '/{sample}.' + contigs_taxonomy_classify_method + '.contig2classification.tsv',
        contigs_fun_detail = path + 'results/c-sample/02_function/' + contigs_function_classify_method + '/annotations/{sample}.' + contigs_function_classify_method + '.annotations.gz',
        vMAGs_tax_detail = path + 'results/v-sample/10_summary/{sample}/vMAG_summary_info.tsv',
        vMAGs_tax_detail_filter = path + 'results/vcr-sample/05_merge_tax_fun_depth/{sample}.vMAG_summary_info.filter.tsv',
        vMAGs_ffn_length = path + 'results/v-sample/10_summary/{sample}/Virus_annotation_results.tsv',
        vMAGs_contig_length = path + 'results/v-sample/10_summary/{sample}/final_virus.contig_length.txt',
        reads_fun_detail = path + 'results/vcr-sample/04_function/' + reads_function_classify_method + '/{sample}.' + reads_function_classify_method + '.annotations.gz',
        merge_detail = path + 'results/vcr-sample/05_merge_tax_fun_depth/detail_for_each_seq/{sample}.detail_for_each_seq.tsv',
        read_type = READ_TYPE,
        minORFLen = config["Params_vcr_sample"]["05_merge_tax_fun_depth"]["minORFLen"],
        checkv_quality = config["Params_vcr_sample"]["05_merge_tax_fun_depth"]["checkv_quality"],
        KO_mode = config["Params_vcr_sample"]["05_merge_tax_fun_depth"]["KO_mode"],
        min_KO_num_in_one_taxon = config["Params_vcr_sample"]["05_merge_tax_fun_depth"]["min_KO_num_in_one_taxon"],
        rank = ",".join(RANKS),
        vc_consistency = path + 'results/vcr-sample/05_merge_tax_fun_depth/MAG_contig_consistency/{sample}.vc_consistency.tsv',
        out_CAMI_format = path + 'results/vcr-sample/05_merge_tax_fun_depth/Taxa_abundance/{sample}.Taxa_abundance.CAMI_format.tsv',
        CAMI_format_output = config["Params_vcr_sample"]["05_merge_tax_fun_depth"]["CAMI_format_output"],
        output_detail = config["Params_vcr_sample"]["05_merge_tax_fun_depth"]["output_detail"],
        reads_fna_length_arg = lambda wc: f'--in_read_fna_length {path}results/vcr-sample/02_FragGeneScanRs/{wc.sample}_unmapped_s.fna_length.gz' if READ_TYPE in _long_read_types else '',
    threads: config["Params_vcr_sample"]["05_merge_tax_fun_depth"]["threads"]
    conda:
        config["conda_envs"]["metalink-main"]
    log:
        path + 'logs/vcr-sample/05_merge_tax_fun_depth/{sample}.log'
    benchmark:
        path + 'benchmarks/vcr-sample/05_merge_tax_fun_depth/{sample}.benchmark'
    shell:
        """
        echo "Start vcr_sample_05_merge_tax_fun_depth for {wildcards.sample}" > {log}
        mkdir -p $(dirname {params.vc_consistency}) 2>{log}

        reads_miss=""
        for f in {input.reads_singleton_tax_detail} \
                {params.reads_fun_detail} {input.reads_ffn_length} {input.reads_length}; do
            [ -s "$f" ] || reads_miss="$reads_miss $f"
        done

        contigs_miss=""
        for f in {params.contigs_tax_detail} {params.contigs_fun_detail} \
                {input.contigs_ffn_length} {input.contigs_depth}; do
            [ -s "$f" ] || contigs_miss="$contigs_miss $f"
        done

        vMAGs_miss=""
        for f in {input.contigs_ORF_loci} \
                {params.vMAGs_tax_detail} {params.vMAGs_ffn_length} {params.vMAGs_contig_length}; do
            [ -s "$f" ] || vMAGs_miss="$vMAGs_miss $f"
        done

        if [ -n "$vMAGs_miss" ]; then
            echo "vMAGs inputs missing/empty; passing vMAGs-related processes for {wildcards.sample}" >> {log}
        else
            echo "Filtering vMAGs for {wildcards.sample}" >> {log}
            Rscript {Dir_script_metalink}/MAGs_v_quality_filter.cmd.R \
                --in_vMAGs_tax_detail {params.vMAGs_tax_detail} \
                --checkv_quality {params.checkv_quality} \
                --out_vMAGs_tax_detail {params.vMAGs_tax_detail_filter} >>{log} 2>&1

            for f in {params.vMAGs_tax_detail_filter}; do
                [ -s "$f" ] || vMAGs_miss="$vMAGs_miss $f"
            done
        fi

        if [ -n "$vMAGs_miss" ] && [ -n "$reads_miss" ] && [ -n "$contigs_miss" ]; then
            echo "ERROR: reads, contigs, and vMAGs inputs missing/empty for {wildcards.sample}; cannot merge taxonomy/function/depth." >> {log}
            echo "  missing/empty reads inputs:$reads_miss" >> {log}
            echo "  missing/empty contigs inputs:$contigs_miss" >> {log}
            echo "  missing/empty vMAGs inputs:$vMAGs_miss" >> {log}
            mkdir -p $(dirname {output.mode_flag})
            echo "ERROR [vcr_sample_05_merge_tax_fun_depth] sample {wildcards.sample}: reads, contigs, and vMAGs inputs missing/empty; check upstream results. See log: {log}" >&2
            exit 1
        elif [ -n "$vMAGs_miss" ] && [ -n "$reads_miss" ] && [ -z "$contigs_miss" ]; then
            echo "No vMAGs and reads, but contigs present; running contigs-only merge." >> {log}
            Rscript {Dir_script_metalink}/combine_taxonomy_function_depth.c.cmd.R \
                --in_taxid_tax {params.taxid_tax} \
                --in_taxid_fun {params.taxid_fun} \
                --contigs_taxonomy_classify_method {params.contigs_taxonomy_classify_method} \
                --contigs_function_classify_method {params.contigs_function_classify_method} \
                --in_contig_tax_detail {params.contigs_tax_detail} \
                --in_contig_fun_detail {params.contigs_fun_detail} \
                --in_contig_ffn_length {input.contigs_ffn_length} \
                --in_contig_depth {input.contigs_depth} \
                --in_read_length {input.reads_length} \
                --minORFLen {params.minORFLen} \
                --KO_mode {params.KO_mode} \
                --min_KO_num_in_one_taxon {params.min_KO_num_in_one_taxon} \
                --rank {params.rank} \
                --sample_name {wildcards.sample} \
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
                --detail {params.merge_detail} \
                --out_mode_flag {output.mode_flag} >>{log} 2>&1
            if [ "{params.output_detail}" = "True" ]; then
                pigz -p {threads} -r -f {params.merge_detail} > {params.merge_detail}.gz 2>> {log}
            fi

        elif [ -n "$vMAGs_miss" ] && [ -z "$reads_miss" ] && [ -n "$contigs_miss" ]; then
            echo "No vMAGs and contigs, but reads present; running reads-only merge." >> {log}
            Rscript {Dir_script_metalink}/combine_taxonomy_function_depth.r.cmd.R \
                --reads_function_classify_method {params.reads_function_classify_method} \
                --in_read_singleton_tax_detail {input.reads_singleton_tax_detail} \
                --in_taxid_tax {params.taxid_tax} \
                --in_taxid_fun {params.taxid_fun} \
                --in_KO_length_ref {params.in_KO_length_ref} \
                --in_read_fun_detail {params.reads_fun_detail} \
                --in_read_ffn_length {input.reads_ffn_length} \
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
                --detail {params.merge_detail} \
                --out_mode_flag {output.mode_flag} >>{log} 2>&1
            if [ "{params.output_detail}" = "True" ]; then
                pigz -p {threads} -r -f {params.merge_detail} > {params.merge_detail}.gz 2>> {log}
            fi

        elif [ -n "$vMAGs_miss" ] && [ -z "$reads_miss" ] && [ -z "$contigs_miss" ]; then
            echo "No vMAGs, but reads and contigs present; running contigs-reads merge." >> {log}
            Rscript {Dir_script_metalink}/combine_taxonomy_function_depth.cr.cmd.R \
                --in_taxid_tax {params.taxid_tax} \
                --in_taxid_fun {params.taxid_fun} \
                --in_KO_length_ref {params.in_KO_length_ref} \
                --contigs_taxonomy_classify_method {params.contigs_taxonomy_classify_method} \
                --contigs_function_classify_method {params.contigs_function_classify_method} \
                --in_contig_tax_detail {params.contigs_tax_detail} \
                --in_contig_fun_detail {params.contigs_fun_detail} \
                --reads_function_classify_method {params.reads_function_classify_method} \
                --in_read_singleton_tax_detail {input.reads_singleton_tax_detail} \
                --in_read_fun_detail {params.reads_fun_detail} \
                --in_contig_ffn_length {input.contigs_ffn_length} \
                --in_read_ffn_length {input.reads_ffn_length} \
                {params.reads_fna_length_arg} \
                --in_contig_depth {input.contigs_depth} \
                --in_read_length {input.reads_length} \
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
                --detail {params.merge_detail} \
                --out_mode_flag {output.mode_flag} >>{log} 2>&1
            if [ "{params.output_detail}" = "True" ]; then
                pigz -p {threads} -r -f {params.merge_detail} > {params.merge_detail}.gz 2>> {log}
            fi

        elif [ -z "$vMAGs_miss" ] && [ -n "$reads_miss" ] && [ -z "$contigs_miss" ]; then
            echo "No reads, but vMAGs and contigs present; running vMAGs-contigs merge." >> {log}
            Rscript {Dir_script_metalink}/combine_taxonomy_function_depth.vc.cmd.R \
                --in_taxid_tax {params.taxid_tax} \
                --in_taxid_fun {params.taxid_fun} \
                --contigs_taxonomy_classify_method {params.contigs_taxonomy_classify_method} \
                --contigs_function_classify_method {params.contigs_function_classify_method} \
                --in_contig_tax_detail {params.contigs_tax_detail} \
                --in_contig_fun_detail {params.contigs_fun_detail} \
                --in_contig_ffn_length {input.contigs_ffn_length} \
                --in_contig_depth {input.contigs_depth} \
                --in_contigs_ORF_loci {input.contigs_ORF_loci} \
                --in_vMAGs_tax_detail {params.vMAGs_tax_detail_filter} \
                --in_vMAGs_ffn_length {params.vMAGs_ffn_length} \
                --in_vMAGs_contig_length {params.vMAGs_contig_length} \
                --sample_name {wildcards.sample} \
                --in_read_length {input.reads_length} \
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
                --out_vc_consistency {params.vc_consistency} \
                --CAMI_format_output {params.CAMI_format_output} \
                --out_CAMI_format {params.out_CAMI_format} \
                --output_detail {params.output_detail} \
                --detail {params.merge_detail} \
                --out_mode_flag {output.mode_flag} >>{log} 2>&1
            if [ "{params.output_detail}" = "True" ]; then
                pigz -p {threads} -r -f {params.merge_detail} > {params.merge_detail}.gz 2>> {log}
            fi

        elif [ -z "$vMAGs_miss" ] && [ -z "$reads_miss" ] && [ -z "$contigs_miss" ]; then
            echo "vMAGs, reads and contigs present; running vMAGs-contigs-reads merge." >> {log}
            Rscript {Dir_script_metalink}/combine_taxonomy_function_depth.vcr.cmd.R \
                --in_taxid_tax {params.taxid_tax} \
                --in_taxid_fun {params.taxid_fun} \
                --in_KO_length_ref {params.in_KO_length_ref} \
                --contigs_taxonomy_classify_method {params.contigs_taxonomy_classify_method} \
                --contigs_function_classify_method {params.contigs_function_classify_method} \
                --in_contig_tax_detail {params.contigs_tax_detail} \
                --in_contig_fun_detail {params.contigs_fun_detail} \
                --reads_function_classify_method {params.reads_function_classify_method} \
                --in_read_singleton_tax_detail {input.reads_singleton_tax_detail} \
                --in_read_fun_detail {params.reads_fun_detail} \
                --in_contig_ffn_length {input.contigs_ffn_length} \
                --in_read_ffn_length {input.reads_ffn_length} \
                {params.reads_fna_length_arg} \
                --in_contig_depth {input.contigs_depth} \
                --in_contigs_ORF_loci {input.contigs_ORF_loci} \
                --in_vMAGs_tax_detail {params.vMAGs_tax_detail_filter} \
                --in_vMAGs_ffn_length {params.vMAGs_ffn_length} \
                --in_vMAGs_contig_length {params.vMAGs_contig_length} \
                --sample_name {wildcards.sample} \
                --in_read_length {input.reads_length} \
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
                --out_vc_consistency {params.vc_consistency} \
                --CAMI_format_output {params.CAMI_format_output} \
                --out_CAMI_format {params.out_CAMI_format} \
                --output_detail {params.output_detail} \
                --detail {params.merge_detail} \
                --out_mode_flag {output.mode_flag} >>{log} 2>&1
            if [ "{params.output_detail}" = "True" ]; then
                pigz -p {threads} -r -f {params.merge_detail} > {params.merge_detail}.gz 2>> {log}
            fi

        else
            echo "Error: Detect vMAGs but no available contigs for {wildcards.sample}" >> {log}
            exit 1
        fi
        """

rule vcr_sample_05_visualize:
    priority: 40
    input:
        summary_Tax = path + 'results/vcr-sample/05_merge_tax_fun_depth/summary/{sample}.summary_Tax.tsv',
        summary_ORF_KO = path + 'results/vcr-sample/05_merge_tax_fun_depth/summary/{sample}.summary_ORF_KO.tsv',
        summary_ORF_Tax_KO = path + 'results/vcr-sample/05_merge_tax_fun_depth/summary/{sample}.summary_ORF_Tax_KO.tsv',
        mode_flag = path + 'results/vcr-sample/05_merge_tax_fun_depth/mode/{sample}.mode.tsv',
    output:
        mk = path + 'markers/vcr-sample/05_visualize/{sample}.mk',
    params:
        prefix = path + 'results/vcr-sample/05_visualization/{sample}',
    threads: 1
    conda:
        config["conda_envs"]["metalink-main"]
    log:
        path + 'logs/vcr-sample/05_visualize/{sample}.log'
    benchmark:
        path + 'benchmarks/vcr-sample/05_visualize/{sample}.benchmark'
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
            tax_format=$(cat {input.mode_flag})
            echo "All summary files non-empty, running visualization with TaxFormat $tax_format." >>{log}
            Rscript {Dir_script_metalink}/plot_summary.cmd.R \
                --TaxFormat $tax_format \
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

rule vcr_sample_06_module_completeness:
    priority: 40
    input:
        KO_Taxa_RelativeAbundance = path + 'results/vcr-sample/05_merge_tax_fun_depth/KO_Taxa_RelativeAbundance/{sample}.KO_Taxa_RelativeAbundance.Taxa.tsv',
    output:
        mk = path + 'markers/vcr-sample/06_module_completeness/{sample}.{rank}.mk',
    params:
        KO_Taxa_RelativeAbundance_ranked = path + 'results/vcr-sample/05_merge_tax_fun_depth/KO_Taxa_RelativeAbundance/{sample}.KO_Taxa_RelativeAbundance.Taxa.{rank}.tsv',
        KO_Taxa_SumCoverage_ranked = path + 'results/vcr-sample/05_merge_tax_fun_depth/KO_Taxa_SumCoverage/{sample}.KO_Taxa_SumCoverage.Taxa.{rank}.tsv',
        KO_Taxa_SumCoverage = path + 'results/vcr-sample/05_merge_tax_fun_depth/KO_Taxa_SumCoverage/{sample}.KO_Taxa_SumCoverage.Taxa.tsv',
        Module_Taxa_RA = path + 'results/vcr-sample/06_module_completeness/{sample}/{rank}/',
        pathway_split = "TRUE",
        min_KO_num_in_one_taxon = config["Params_vcr_sample"]["06_mclink_taxon"]["min_KO_num_in_one_taxon"],
        min_mean_KO_RA_in_one_taxon = config["Params_vcr_sample"]["06_mclink_taxon"]["min_mean_KO_RA_in_one_taxon"],
        completeness_threshold = config["Params_vcr_sample"]["06_mclink_taxon"]["completeness_threshold"],
        plus_scale_method = config["Params_vcr_sample"]["06_mclink_taxon"]["plus_scale_method"],
        comma_scale_method = config["Params_vcr_sample"]["06_mclink_taxon"]["comma_scale_method"],
        TopN = config["Params_vcr_sample"]["06_mclink_taxon"]["TopN"],
    threads: config["Params_vcr_sample"]["05_merge_tax_fun_depth"]["threads"]
    conda:
        config["conda_envs"]["metalink-main"]
    log:
        path + 'logs/vcr-sample/06_module_completeness/{sample}.{rank}.log'
    benchmark:
        path + 'benchmarks/vcr-sample/06_module_completeness/{sample}.{rank}.benchmark'
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
            echo "KO_Taxa_RelativeAbundance empty/missing for {wildcards.sample}; skip vcr_sample_03_module_completeness." >>{log}
        fi
        touch {output.mk}
        """


rule vcr_project_00_FragGeneScanRs_stats:
    input:
        faa = expand(path + 'results/vcr-sample/02_FragGeneScanRs_stats/{sample}_unmapped_s.faa.stats', sample=SAMPLE_INDEX),
        fna = expand(path + 'results/vcr-sample/02_FragGeneScanRs_stats/{sample}_unmapped_s.fna.stats', sample=SAMPLE_INDEX),
        ffn = expand(path + 'results/vcr-sample/02_FragGeneScanRs_stats/{sample}_unmapped_s.ffn.stats', sample=SAMPLE_INDEX),
    output:
        faa = temp(path + 'results/vcr-project/00_FragGeneScanRs_stats/{project}/02_FGS_faa_stats.tsv'),
        fna = temp(path + 'results/vcr-project/00_FragGeneScanRs_stats/{project}/02_FGS_fna_stats.tsv'),
        ffn = temp(path + 'results/vcr-project/00_FragGeneScanRs_stats/{project}/02_FGS_ffn_stats.tsv'),
        merged_stats = path + 'results/vcr-project/00_FragGeneScanRs_stats/{project}/02_FGS_stats.tsv',
    conda:
        config["conda_envs"]["metalink-main"]
    log:
        path + 'logs/vcr-project/00_FragGeneScanRs_stats/{project}/02_cat_FGS_stats.log',
    benchmark:
        path + 'benchmarks/vcr-project/00_FragGeneScanRs_stats/{project}/02_cat_FGS_stats.benchmark'
    shell:
        """
        # Process faa.stats files
        for file1 in {input.faa}; do
            if [ -s {output.faa} ]; then
                echo "pass"
            else
                echo -e "ID\tformat_faa\ttype_faa\tnum_seqs_faa\tsum_len_faa\tmin_len_faa\tavg_len_faa\tmax_len_faa" > {output.faa}
            fi
            sed 1d $file1 | awk -F '\\t' 'BEGIN {{OFS="\\t"}} {{gsub("_unmapped_s.faa.gz", "", $1); print}}' >> {output.faa}
        done 2> {log}

        # Process fna.stats files
        for file1 in {input.fna}; do
            if [ -s {output.fna} ]; then
                echo "pass"
            else
                echo -e "ID\tformat_fna\ttype_fna\tnum_seqs_fna\tsum_len_fna\tmin_len_fna\tavg_len_fna\tmax_len_fna" > {output.fna}
            fi
            sed 1d $file1 | awk -F '\\t' 'BEGIN {{OFS="\\t"}} {{gsub("_unmapped_s.fna.gz", "", $1); print}}' >> {output.fna}
        done 2>> {log}

        # Process ffn.stats files
        for file1 in {input.ffn}; do
            if [ -s {output.ffn} ]; then
                echo "pass"
            else
                echo -e "ID\tformat_ffn\ttype_ffn\tnum_seqs_ffn\tsum_len_ffn\tmin_len_ffn\tavg_len_ffn\tmax_len_ffn" > {output.ffn}
            fi
            sed 1d $file1 | awk -F '\\t' 'BEGIN {{OFS="\\t"}} {{gsub("_unmapped_s.ffn.gz", "", $1); print}}' >> {output.ffn}
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

rule vcr_project_01_cat_Taxa_abundance_dcast:
    input:
        Taxa_abundance = expand(path + 'results/vcr-sample/05_merge_tax_fun_depth/Taxa_abundance/{sample}.Taxa_abundance.tsv', sample=SAMPLE_INDEX),
    output:
        out_long = path + 'results/vcr-project/01_Taxa_Sample/{project}/Taxa_Sample.long.tsv',
    params:
        tmp_dir = path + 'results/vcr-project/tmp/{project}/',
        prefix_wide = path + 'results/vcr-project/01_Taxa_Sample/{project}/Taxa_Sample.wide.tsv',
    conda: config["conda_envs"]["metalink-main"]
    threads: config["Params_vcr_sample"]["05_merge_tax_fun_depth"]["threads"]
    log: path + 'logs/vcr-project/01_cat_Taxa_abundance_dcast/{project}.log'
    benchmark:
        path + 'benchmarks/vcr-project/01_cat_Taxa_abundance_dcast/{project}.benchmark'
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

rule vcr_project_01_Taxa_visualize:
    input:
        out_long = path + 'results/vcr-project/01_Taxa_Sample/{project}/Taxa_Sample.long.tsv',
    output:
        mk = path + 'markers/vcr-project/01_Taxa_visualization/{project}.mk',
    params:
        prefix = path + 'results/vcr-project/01_Taxa_visualization/{project}',
        prefix_wide = path + 'results/vcr-project/01_Taxa_Sample/{project}/Taxa_Sample.Taxa_ReadCount.wide.tsv',
        TopN = config["Params_vcr_project"]["01_Taxa_visualize"]["TopN"],
    threads: 1
    conda:
        config["conda_envs"]["metalink-main"]
    log:
        path + 'logs/vcr-project/01_Taxa_visualization/{project}.log'
    benchmark:
        path + 'benchmarks/vcr-project/01_Taxa_visualization/{project}.benchmark'
    shell:
        """
        mkdir -p $(dirname {params.prefix}) 2>>{log}
        missing=""
        for f in {input.out_long} {params.prefix_wide}; do
            [ -s "$f" ] || missing="$missing $f"
        done
        if [ -n "$missing" ]; then
            echo "Missing/empty inputs:$missing; skip vcr project_01_Taxa_visualize for {wildcards.project}" >>{log}
        else
            echo "All input files non-empty, running vcr project_01_Taxa_visualize." >>{log}
            Rscript {Dir_script_metalink}/plot_bar.cmd.R \
                --Type Taxa \
                --TopN {params.TopN} \
                --Input {params.prefix_wide} \
                --Output {params.prefix} >>{log} 2>&1
        fi
        touch {output.mk}
        """

rule vcr_project_02_cat_KO_depth_dcast:
    input:
        KO_depth = expand(path + 'results/vcr-sample/05_merge_tax_fun_depth/KO_depth/{sample}.KO_depth.tsv', sample=SAMPLE_INDEX),
    output:
        out_long = path + 'results/vcr-project/02_KO_Sample/{project}/KO_Sample.long.tsv',
    params:
        tmp_dir = path + 'results/vcr-project/tmp/{project}/',
        prefix_wide = path + 'results/vcr-project/02_KO_Sample/{project}/KO_Sample.wide.tsv',
    threads: config["Params_vcr_sample"]["05_merge_tax_fun_depth"]["threads"]
    conda:
        config["conda_envs"]["metalink-main"]
    log:
        path + 'logs/vcr-project/02_cat_KO_depth_dcast/{project}.log'
    benchmark:
        path + 'benchmarks/vcr-project/02_cat_KO_depth_dcast/{project}.benchmark'
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

rule vcr_project_03_community_module:
    input:
        out_long = path + 'results/vcr-project/02_KO_Sample/{project}/KO_Sample.long.tsv',
    output:
        mk = path + 'markers/vcr-project/03_community_module/{project}.mk',
    params:
        wide_RelativeAbundance = path + 'results/vcr-project/02_KO_Sample/{project}/KO_Sample.KO_RelativeAbundance.wide.tsv',
        wide_SumCoverage = path + 'results/vcr-project/02_KO_Sample/{project}/KO_Sample.KO_SumCoverage.wide.tsv',
        pathway_split = "TRUE",
        community_output = path + 'results/vcr-project/03_Module_Sample/{project}/',
        min_sample_num_per_KO = config["Params_vcr_project"]["03_mclink_sample"]["min_sample_num_per_KO"],
        preva = config["Params_vcr_project"]["03_mclink_sample"]["preva"],
        abund = config["Params_vcr_project"]["03_mclink_sample"]["abund"],
        completeness_threshold = config["Params_vcr_project"]["03_mclink_sample"]["completeness_threshold"],
        plus_scale_method = config["Params_vcr_project"]["03_mclink_sample"]["plus_scale_method"],
        comma_scale_method = config["Params_vcr_project"]["03_mclink_sample"]["comma_scale_method"],
    threads: config["Params_vcr_sample"]["05_merge_tax_fun_depth"]["threads"]
    conda:
        config["conda_envs"]["metalink-main"]
    log:
        path + 'logs/vcr-project/03_community_module/{project}.log'
    benchmark:
        path + 'benchmarks/vcr-project/03_community_module/{project}.benchmark'
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
            echo "KO_Sample.KO_RelativeAbundance.wide.tsv empty/missing for {wildcards.project}; skip vcr_project_03_community_module." >>{log}
        fi
        touch {output.mk}
        """

rule vcr_project_03_Fun_visualize:
    input:
        mk = path + 'markers/vcr-project/03_community_module/{project}.mk',
    output:
        mk = path + 'markers/vcr-project/03_Fun_visualization/{project}.mk',
    params:
        prefix = path + 'results/vcr-project/03_Fun_visualization/{project}',
        prefix_wide = path + 'results/vcr-project/03_Module_Sample/{project}/SumCoverage/mc_abundance_filtered_by_completeness.tsv',
        TopN = config["Params_vcr_project"]["03_Fun_visualize"]["TopN"],
    threads: 1
    conda:
        config["conda_envs"]["metalink-main"]
    log:
        path + 'logs/vcr-project/03_Fun_visualization/{project}.log'
    benchmark:
        path + 'benchmarks/vcr-project/03_Fun_visualization/{project}.benchmark'
    shell:
        """
        mkdir -p $(dirname {params.prefix}) 2>>{log}
        missing=""
        for f in {params.prefix_wide}; do
            [ -s "$f" ] || missing="$missing $f"
        done
        if [ -n "$missing" ]; then
            echo "Missing/empty inputs:$missing; skip vcr project_03_Fun_visualize for {wildcards.project}" >>{log}
        else
            echo "All input files non-empty, running vcr project_03_Fun_visualize." >>{log}
            Rscript {Dir_script_metalink}/plot_bar.cmd.R \
                --Type Module \
                --TopN {params.TopN} \
                --Input {params.prefix_wide} \
                --Output {params.prefix} >>{log} 2>&1
        fi
        touch {output.mk}
        """

rule vcr_project_05_merge_tax_KO_depth:
    input:
        ko_depth = expand(path + 'results/vcr-sample/05_merge_tax_fun_depth/KO_depth/{sample}.KO_depth.tsv', sample=SAMPLE_INDEX),
    output:
        out_long = path + 'results/vcr-project/05_merge_tax_KO_depth/{project}/{rank}.mean.long.tsv',
    params:
        detail = lambda wc: expand(path + 'results/vcr-sample/05_merge_tax_fun_depth/detail_for_each_seq/{sample}.detail_for_each_{rank}.tsv', sample=SAMPLE_INDEX, rank=wc.rank),  # 只展开当前 rank
        tmp_dir = path + 'results/vcr-project/tmp/{project}/',
        prefix_wide = path + 'results/vcr-project/05_merge_tax_KO_depth/{project}/{rank}.wide.tsv',
        preva = config["Params_vcr_project"]["05_merge_tax_KO_depth"]["preva"],
        abund = config["Params_vcr_project"]["05_merge_tax_KO_depth"]["abund"],
    threads: config["Params_vcr_sample"]["05_merge_tax_fun_depth"]["threads"]
    conda:
        config["conda_envs"]["metalink-main"]
    log:
        path + 'logs/vcr-project/05_merge_tax_KO_depth/{project}/{rank}.log'
    benchmark:
        path + 'benchmarks/vcr-project/05_merge_tax_KO_depth/{project}/{rank}.benchmark'
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

rule vcr_project_05_merge_tax_Module_RelativeAbundance:
    input:
        mk = lambda wc: expand(path + 'markers/vcr-sample/06_module_completeness/{sample}.{rank}.mk', sample=SAMPLE_INDEX, rank=wc.rank),
    output:
        out_long = path + 'results/vcr-project/05_merge_tax_Module_RelativeAbundance/{project}/{rank}.mean.long.tsv',
    params:
        detail = lambda wc: expand(path + 'results/vcr-sample/06_module_completeness/{sample}/{rank}/RelativeAbundance/mc_abundance.tsv', sample=SAMPLE_INDEX, rank=wc.rank),
        prefix_wide = path + 'results/vcr-project/05_merge_tax_Module_RelativeAbundance/{project}/{rank}.wide.tsv',
        tmp_dir = path + 'results/vcr-project/tmp/{project}/',
        preva = config["Params_vcr_project"]["05_merge_tax_Module_RelativeAbundance"]["preva"],
        abund = config["Params_vcr_project"]["05_merge_tax_Module_RelativeAbundance"]["abund"],
    threads: config["Params_vcr_sample"]["05_merge_tax_fun_depth"]["threads"]
    conda:
        config["conda_envs"]["metalink-main"]
    log:
        path + 'logs/vcr-project/05_merge_tax_Module_RelativeAbundance/{project}/{rank}.log'
    benchmark:
        path + 'benchmarks/vcr-project/05_merge_tax_Module_RelativeAbundance/{project}/{rank}.benchmark'
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

rule vcr_project_05_merge_tax_Module_SumCoverage:
    input:
        mk = lambda wc: expand(path + 'markers/vcr-sample/06_module_completeness/{sample}.{rank}.mk', sample=SAMPLE_INDEX, rank=wc.rank),
    output:
        out_long = path + 'results/vcr-project/05_merge_tax_Module_SumCoverage/{project}/{rank}.mean.long.tsv',
    params:
        detail = lambda wc: expand(path + 'results/vcr-sample/06_module_completeness/{sample}/{rank}/SumCoverage/mc_abundance.tsv', sample=SAMPLE_INDEX, rank=wc.rank),
        prefix_wide = path + 'results/vcr-project/05_merge_tax_Module_SumCoverage/{project}/{rank}.wide.tsv',
        tmp_dir = path + 'results/vcr-project/tmp/{project}/',
        preva = config["Params_vcr_project"]["05_merge_tax_Module_SumCoverage"]["preva"],
        abund = config["Params_vcr_project"]["05_merge_tax_Module_SumCoverage"]["abund"],
    threads: config["Params_vcr_sample"]["05_merge_tax_fun_depth"]["threads"]
    conda:
        config["conda_envs"]["metalink-main"]
    log:
        path + 'logs/vcr-project/05_merge_tax_Module_SumCoverage/{project}/{rank}.log'
    benchmark:
        path + 'benchmarks/vcr-project/05_merge_tax_Module_SumCoverage/{project}/{rank}.benchmark'
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

rule vcr_project_06_module_completeness:
    input:
        KO_Taxa_RelativeAbundance = path + 'results/vcr-project/05_merge_tax_KO_depth/{project}/{rank}.mean.long.tsv',
    output:
        mk = path + 'markers/vcr-project/06_module_completeness/{project}/{rank}.mk',
    params:
        KO_Taxa_RelativeAbundance_ranked = path + 'results/vcr-project/05_merge_tax_KO_depth/{project}/{rank}.Taxa_KO_RelativeAbundance.mean.wide.tsv',
        KO_Taxa_SumCoverage_ranked = path + 'results/vcr-project/05_merge_tax_KO_depth/{project}/{rank}.Taxa_KO_SumCoverage.mean.wide.tsv',
        Module_Taxa_RA = path + 'results/vcr-project/06_module_completeness/{project}/{rank}/',
        pathway_split = "TRUE",
        min_KO_num_in_one_taxon = config["Params_vcr_project"]["06_mclink_taxon"]["min_KO_num_in_one_taxon"],
        min_mean_KO_RA_in_one_taxon = config["Params_vcr_project"]["06_mclink_taxon"]["min_mean_KO_RA_in_one_taxon"],
        completeness_threshold = config["Params_vcr_project"]["06_mclink_taxon"]["completeness_threshold"],
        plus_scale_method = config["Params_vcr_project"]["06_mclink_taxon"]["plus_scale_method"],
        comma_scale_method = config["Params_vcr_project"]["06_mclink_taxon"]["comma_scale_method"],
        TopN = config["Params_vcr_project"]["06_mclink_taxon"]["TopN"],
    threads: config["Params_vcr_sample"]["05_merge_tax_fun_depth"]["threads"]
    conda:
        config["conda_envs"]["metalink-main"]
    log:
        path + 'logs/vcr-project/06_module_completeness/{project}/{rank}.log'
    benchmark:
        path + 'benchmarks/vcr-project/06_module_completeness/{project}/{rank}.benchmark'
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
            echo "KO_Taxa_RelativeAbundance empty/missing for {wildcards.project}; skip vcr_project_module_completeness." >>{log}
        fi
        touch {output.mk}
        """

# The b-sample BAM is consumed by b_sample_04_bam_summary (upstream of the
# unmapped-read extraction), b_sample_07_contigs_prodigal (jgi depth) and the
# unmapped-read extraction itself (via params); the inputs below are their final
# downstream products, so the BAM is only deleted once all consumers are done.
rule vcr_sample_12_delete_bam:
    input:
        unmapped_r_s = path + 'results/vcr-sample/01_unmapped_reads/{sample}_unmapped_s.fastq.gz',
        contigs_fna = path + 'results/b-sample/07_contigs/contigs_fna/{sample}.fna',
    output:
        mk = path + 'markers/vcr-sample/12_delete_bam/{sample}.mk'
    threads: 1
    priority: 100
    params:
        bam = path + 'results/b-sample/03_mapping/' + base_sample_mapping_method + '/{sample}.bam',
    log:
        path + 'logs/vcr-sample/12_delete_bam/{sample}.log'
    benchmark:
        path + 'benchmarks/vcr-sample/12_delete_bam/{sample}.benchmark'
    shell:
        """
        if [ -e "{params.bam}" ];then rm -r {params.bam}*;else echo Folder trimmed.bam not exist.;fi > {log} 2>&1
        touch {output.mk}
        """
