'''
/*
 * @Date: 2022-12-04 17:15:47
 * @LastEditors: LiuyangLi
 * @LastEditTime: 2024-10-16 21:29:37
 * @FilePath: c.base.smk
 * @Description:
 */
'''
import os
import sys

Dir_script_metalink = config["Directory"]["script"]["metalink"]
kegg_modules_ko_info_path = config["Directory"]["script"]["kegg_modules_ko_info_path"]
kraken2_taxonomy_db = config["Directory"]["database"]["kraken2_taxonomy_db"]
CAT_taxonomy_db = config["Directory"]["database"]["CAT_taxonomy_db"]
CAT_taxonomy_tax = config["Directory"]["database"]["CAT_taxonomy_tax"]
taxid_tax_path = config["Directory"]["database"]["taxid_tax_path"]
taxid_fun_path = config["Directory"]["database"]["taxid_fun_path"]
contigs_taxonomy_classify_method = config["Params_c_sample"]["01_taxonomy"]["contigs_taxonomy_classify_method"]
contigs_function_classify_method = config["Params_c_sample"]["02_function"]["contigs_function_classify_method"]


folder_checks = [
                path + 'results/',
                path + 'results/c-sample/',
                path + 'results/c-sample/02_function/',
                path + 'results/c-sample/02_function/emapper/',
                path + 'results/c-sample/02_function/emapper/annotations/',
                path + 'results/c-project/',
                ]
for item in folder_checks:
    if not os.path.isdir(item): os.mkdir(item)

#contigs_taxonomy_classify_method=CAT_pack,kraken2
def get_contigs_sample_01_taxonomy_classify_rule_mk(wildcards):
    if contigs_taxonomy_classify_method == "CAT_pack":
        return rules.c_sample_01_taxonomy_CAT_pack.output.mk
    elif contigs_taxonomy_classify_method == "kraken2":
        return rules.c_sample_01_taxonomy_kraken2.output.mk
    else:
        raise ValueError(f"Unknown contigs_taxonomy_classify_method for c: {contigs_taxonomy_classify_method}")

##contigs_function_classify_method=emapper,kraken2
def get_contigs_sample_02_function_classify_rule_mk(wildcards):
    if contigs_function_classify_method == "emapper":
        return rules.c_sample_02_function_emapper.output.mk
    elif contigs_function_classify_method == "kraken2":
        return rules.c_sample_02_function_kraken2.output.mk
    else:
        raise ValueError(f"Unknown contigs_function_classify_method for c: {contigs_function_classify_method}")


rule c_sample_01_taxonomy_CAT_pack:
    input:
        contigs_fna = path + 'results/b-sample/07_contigs/contigs_fna/{sample}.fna',
        protein_faa = path + 'results/b-sample/07_contigs/contigs_protein_faa/{sample}.faa'
    output:
        mk = path + 'markers/c-sample/01_taxonomy/CAT_pack/{sample}.mk',
    params:
        db = CAT_taxonomy_db,
        tax = CAT_taxonomy_tax,
        out_dir= path + 'results/c-sample/01_taxonomy/CAT_pack/',
        CAT_range = config["Params_c_sample"]["01_taxonomy_CAT_pack"]["range"],
        CAT_fraction = config["Params_c_sample"]["01_taxonomy_CAT_pack"]["fraction"],
        extra_args = config["Params_c_sample"]["01_taxonomy_CAT_pack"]["extra_args"],
    conda: config["conda_envs"]["metalink-annotation"]
    priority: 1
    threads: config["Params_c_sample"]["01_taxonomy"]["threads"]
    log:
        path + 'logs/c-sample/01_taxonomy/CAT_pack/{sample}.log'
    benchmark:
        path + 'benchmarks/c-sample/01_taxonomy/CAT_pack/{sample}.benchmark'
    shell:
        """
        mkdir -p {params.out_dir}
        if [ -s {input.contigs_fna} ]; then
            echo "contigs_fna non-empty, running CAT_pack." > {log}
            CAT_pack contigs -n {threads} -r {params.CAT_range} -f {params.CAT_fraction} -c {input.contigs_fna} -p {input.protein_faa} -d {params.db} -t {params.tax} -o {params.out_dir}/{wildcards.sample}.CAT_pack --no_stars --force {params.extra_args} >> {log} 2>&1
            CAT_pack add_names -i {params.out_dir}/{wildcards.sample}.CAT_pack.ORF2LCA.txt -o {params.out_dir}/{wildcards.sample}.CAT_pack.ORF2LCA.name.txt -t {params.tax} --only_official --force >> {log} 2>&1
            CAT_pack add_names -i {params.out_dir}/{wildcards.sample}.CAT_pack.contig2classification.txt -o {params.out_dir}/{wildcards.sample}.CAT_pack.contig2classification.name.txt -t {params.tax} --only_official --force >> {log} 2>&1
            CAT_pack summarise -c {input.contigs_fna} -i {params.out_dir}/{wildcards.sample}.CAT_pack.contig2classification.name.txt -o {params.out_dir}/{wildcards.sample}.CAT_pack.summary.tsv --force >> {log} 2>&1
            mv {params.out_dir}/{wildcards.sample}.CAT_pack.ORF2LCA.name.txt {params.out_dir}/{wildcards.sample}.CAT_pack.ORF2LCA.tsv >> {log} 2>&1
            mv {params.out_dir}/{wildcards.sample}.CAT_pack.contig2classification.name.txt {params.out_dir}/{wildcards.sample}.CAT_pack.contig2classification.tsv >> {log} 2>&1
            rm -f {params.out_dir}/{wildcards.sample}.CAT_pack.ORF2LCA.txt >> {log} 2>&1
            rm -f {params.out_dir}/{wildcards.sample}.CAT_pack.contig2classification.txt >> {log} 2>&1
        else
            echo "contigs_fna empty/missing for {wildcards.sample}; skip contigs_sample_01_taxonomy_CAT_pack." >>{log}
            touch {params.out_dir}/{wildcards.sample}.CAT_pack.ORF2LCA.tsv \
                {params.out_dir}/{wildcards.sample}.CAT_pack.contig2classification.tsv \
                {params.out_dir}/{wildcards.sample}.CAT_pack.summary.tsv >> {log} 2>&1
        fi
        touch {output.mk}
        """

rule c_sample_01_taxonomy_kraken2:
    input:
        contigs_fna = path + 'results/b-sample/07_contigs/contigs_fna/{sample}.fna',
    output:
        mk = path + 'markers/c-sample/01_taxonomy/kraken2/{sample}.mk',
        kraken_report = path + 'results/c-sample/01_taxonomy/kraken2/{sample}.kraken2.report.tsv',
        detail = path + 'results/c-sample/01_taxonomy/kraken2/{sample}.kraken2.contig2classification.tsv',
    threads: config["Params_c_sample"]["01_taxonomy"]["threads"]
    params:
        databases_path = kraken2_taxonomy_db,
        confidence = config["Params_c_sample"]["01_taxonomy_kraken2"]["confidence"],
        unclassified_seq = path + 'results/c-sample/01_taxonomy/kraken2/{sample}.c.unclassified.fna',
        classified_seq = path + 'results/c-sample/01_taxonomy/kraken2/{sample}.c.classified.fna',
        out_dir= path + 'results/c-sample/01_taxonomy/kraken2/',
        extra_args = config["Params_c_sample"]["01_taxonomy_kraken2"]["extra_args"],
    conda: config["conda_envs"]["metalink-main"]
    resources: kraken2_taxonomy=1
    log:
        path + 'logs/c-sample/01_taxonomy/kraken2/{sample}.log'
    benchmark:
        path + 'benchmarks/c-sample/01_taxonomy/kraken2/{sample}.benchmark'
    shell:
        """
        mkdir -p {params.out_dir}
        if [ -s {input.contigs_fna} ]; then
            kraken2 {params.extra_args} --threads {threads} \
                --db {params.databases_path} \
                --use-names \
                --confidence {params.confidence} \
                --report {output.kraken_report} \
                {input.contigs_fna} 1> {output.detail} 2> {log}
        else
            echo "contigs_fna empty/missing for {wildcards.sample}; skip contigs_sample_01_taxonomy_kraken2." >>{log}
            touch {params.out_dir}/{wildcards.sample}.kraken2.report.tsv \
                {params.out_dir}/{wildcards.sample}.kraken2.contig2classification.tsv >> {log} 2>&1
        fi
        touch {output.mk}
        """

rule c_sample_01_taxonomy_done:
    input:
        mk = get_contigs_sample_01_taxonomy_classify_rule_mk,
    output:
        mk = path + 'markers/c-sample/01_taxonomy/classify_done/{sample}.mk',
    shell:
        "cp {input.mk} {output.mk}"

rule c_sample_02_function_emapper:
    input:
        protein_faa = path + 'results/b-sample/07_contigs/contigs_protein_faa/{sample}.faa'
    output:
        mk = path + 'markers/c-sample/02_function/emapper/{sample}.mk',
    threads: config["Params_c_sample"]["02_function"]["threads"]
    conda: config["conda_envs"]["metalink-annotation"]
    params:
        seed_ortholog_evalue = config["Params_c_sample"]["02_function_emapper"]["seed_ortholog_evalue"],
        method = config["Params_c_sample"]["02_function_emapper"]["method"],
        query_cover = config["Params_c_sample"]["02_function_emapper"]["query_cover"],
        data_dir = config["Directory"]["database"]["eggnog_db"],
        out = path + 'results/c-sample/02_function/emapper/annotations/{sample}',
        extra_args = config["Params_c_sample"]["02_function_emapper"]["extra_args"],
    log:
        path + 'logs/c-sample/02_function/emapper/{sample}.log'
    benchmark:
        path + 'benchmarks/c-sample/02_function/emapper/{sample}.benchmark'
    shell:
        """
        if [ -s {input.protein_faa} ]; then
            emapper.py \
                -i {input.protein_faa} \
                --cpu {threads} \
                --seed_ortholog_evalue {params.seed_ortholog_evalue} \
                -m {params.method} \
                --query_cover {params.query_cover} \
                --data_dir {params.data_dir} \
                -o {params.out} --override {params.extra_args} > {log} 2>&1
            pigz -p {threads} -r -f {params.out}.emapper.annotations {params.out}.emapper.hits {params.out}.emapper.seed_orthologs 2>> {log}
        else
            echo "protein_faa empty/missing for {wildcards.sample}; skip contigs_sample_02_function_emapper." >>{log}
        fi
        touch {output.mk}
        """

rule c_sample_02_function_kraken2:
    input:
        protein_ffn = path + 'results/b-sample/07_contigs/contigs_protein_ffn/{sample}.ffn'
    output:
        kraken_report = path + 'results/c-sample/02_function/kraken2/annotations/{sample}.kraken2.report.tsv',
        detail = path + 'results/c-sample/02_function/kraken2/annotations/{sample}.kraken2.annotations.gz',
        mk = path + 'markers/c-sample/02_function/kraken2/annotations/{sample}.mk',
    threads: config["Params_c_sample"]["02_function"]["threads"]
    params:
        databases_path = config["Directory"]["database"]["kraken2_function_db"],
        confidence = config["Params_c_sample"]["02_function_kraken2"]["confidence"],
        detail = path + 'results/c-sample/02_function/kraken2/annotations/{sample}.kraken2.mag2classification.tsv',
        tmp = path + 'results/c-sample/02_function/kraken2/annotations/{sample}.kraken2.mag2classification.tsv.gz',
        extra_args = config["Params_c_sample"]["02_function_kraken2"]["extra_args"],
    conda: config["conda_envs"]["metalink-main"]
    log:
        path + 'logs/c-sample/02_function/kraken2/{sample}.log'
    benchmark:
        path + 'benchmarks/c-sample/02_function/kraken2/{sample}.benchmark'
    shell:
        """
        if [ -s "{input.protein_ffn}" ]; then
            echo "Found protein ffn file: {input.protein_ffn}, running c_sample_02_function_kraken2..." >> {log}
            mkdir -p $(dirname {params.detail}) 2>> {log}
            kraken2 {params.extra_args} --threads {threads} \
                --db {params.databases_path} \
                --use-names \
                --confidence {params.confidence} \
                --report {output.kraken_report} \
                {input.protein_ffn} 1> {params.detail} 2>> {log}
            pigz -f -p {threads} {params.detail} 2>> {log}
            mv {params.tmp} {output.detail}
        else
            echo "Protein ffn file not found or empty: {input.protein_ffn}, skipping c_sample_02_function_kraken2." >> {log}
            touch {output.kraken_report}
            echo "" | pigz -c > {output.detail}
        fi
        touch {output.mk}
        """

rule c_sample_02_function_done:
    input:
        mk = get_contigs_sample_02_function_classify_rule_mk,
    output:
        mk = path + 'markers/c-sample/02_function/classify_done/{sample}.mk'
    shell:
        "cp {input.mk} {output.mk}"
