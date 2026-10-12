'''
/*
 * @Date: 2022-12-04 17:15:47
 * @LastEditors: LiuyangLi
 * @LastEditTime: 2024-07-20 13:03:29
 * @FilePath: m-s.smk
 * @Description:
 */
'''

import os
import sys

Dir_script_metalink = config["Directory"]["script"]["metalink"]
kegg_modules_ko_info_path = config["Directory"]["script"]["kegg_modules_ko_info_path"]
CAT_taxonomy_db = config["Directory"]["database"]["CAT_taxonomy_db"]
CAT_taxonomy_tax = config["Directory"]["database"]["CAT_taxonomy_tax"]
taxid_tax_path = config["Directory"]["database"]["taxid_tax_path"]
taxid_fun_path = config["Directory"]["database"]["taxid_fun_path"]
MAGs_taxonomy_classify_method = config["Params_m_sample"]["02_taxonomy"]["MAGs_taxonomy_classify_method"]
MAGs_function_classify_method = config["Params_m_sample"]["08_function"]["MAGs_function_classify_method"]
checkm2_db = config["Directory"]["database"]["checkm2_db"]
gtdbtk_db = config["Directory"]["database"]["gtdbtk_db"]

folder_checks = [
                path + 'results/',
                path + 'results/m-sample/',
                path + 'results/m-sample/tmp/',
                path + 'results/m-sample/01_checkm2/',
                path + 'results/m-sample/05_drep_eachsample_99/',
                path + 'results/m-sample/06_drep_eachsample_95/',
                path + 'results/m-sample/08_function/',
                path + 'results/m-sample/08_function/emapper/',
                path + 'results/m-project/',
                ]
for item in folder_checks:
    if not os.path.isdir(item): os.mkdir(item)

base_sample_mapping_method = config["Params_b_sample"]["03_mapping"]["method"]
#base_sample_mapping_method=bwamem,bwamem2,strobealign,minibwa
#MAGs_taxonomy_classify_method=BAT,gtdbtk
def get_m_sample_02_taxonomy_classify_rule_mk(wildcards):
    if MAGs_taxonomy_classify_method == "BAT":
        return rules.m_sample_02_taxonomy_BAT.output.mk
    elif MAGs_taxonomy_classify_method == "gtdbtk":
        return rules.m_sample_02_taxonomy_gtdbtk.output.mk
    else:
        raise ValueError(f"Unknown MAGs_taxonomy_classify_method for m: {MAGs_taxonomy_classify_method}")
##MAGs_function_classify_method=emapper,kraken2
def get_m_sample_08_function_classify_rule_mk(wildcards):
    if MAGs_function_classify_method == "emapper":
        return rules.m_sample_08_function_emapper.output.mk
    elif MAGs_function_classify_method == "kraken2":
        return rules.m_sample_08_function_kraken2.output.mk
    else:
        raise ValueError(f"Unknown MAGs_function_classify_method for m: {MAGs_function_classify_method}")


rule m_sample_01_checkm2:
    input:
        mk = path + 'markers/b-sample/08_bin_prodigal/{sample}.mk',
    output:
        mk = path + 'markers/m-sample/01_checkm2/{sample}.mk'
    params:
        in_path = path + 'results/b-sample/08_gather_bin_high_quality/bin_genome_fna_renamed/{sample}/',
        quality = path + 'results/m-sample/01_checkm2/{sample}/',
        genometype = "fna",
        checkm2_db = checkm2_db,
        extra_args = config["Params_m_sample"]["01_checkm2"]["extra_args"],
    threads: config["Params_m_sample"]["01_checkm2"]["threads"]
    conda: config["conda_envs"]["metalink-checkm2"]
    priority: 10
    log:
        path + 'logs/m-sample/01_checkm2/{sample}.log'
    benchmark:
        path + 'benchmarks/m-sample/01_checkm2/{sample}.benchmark'
    shell:
        """
        if [ -d "{params.in_path}" ] && [ -n "$(ls -A {params.in_path})" ]; then
            echo "Found bins in {params.in_path}, running checkm2..." > {log}
            echo "[CheckM2 CMD] checkm2 predict -t {threads} --input {params.in_path} -x {params.genometype} --force --output-directory {params.quality} {params.extra_args}" >> {log}
            mkdir -p {params.quality}
            checkm2 predict \
                -t {threads} \
                --input {params.in_path} \
                -x {params.genometype} \
                --database_path {params.checkm2_db} \
                --force \
                --output-directory {params.quality} {params.extra_args} >> {log} 2>&1
        else
            echo "No bins found in {params.in_path}, skipping checkm2." >> {log}
        fi
        touch {output.mk}
        """

rule m_sample_02_taxonomy_gtdbtk:
    priority: 10
    input:
        mk = path + 'markers/b-sample/08_bin_prodigal/{sample}.mk',
    output:
        mk = path + 'markers/m-sample/02_taxonomy/gtdbtk/{sample}.mk'
    params:
        in_path = path + 'results/b-sample/08_gather_bin_high_quality/bin_genome_fna_renamed/{sample}/',
        out_path = path + 'results/m-sample/02_taxonomy/gtdbtk/{sample}/',
        mash_db = temp(path + 'results/m-sample/tmp/{sample}/'),
        genometype = "fna",
        gtdbtk_db = gtdbtk_db,
        extra_args = config["Params_m_sample"]["02_taxonomy_gtdbtk"]["extra_args"],
    threads: config["Params_m_sample"]["02_taxonomy"]["threads"]
    conda: config["conda_envs"]["metalink-gtdbtk"]
    log:
        path + 'logs/m-sample/02_taxonomy/gtdbtk/{sample}.log'
    benchmark:
        path + 'benchmarks/m-sample/02_taxonomy/gtdbtk/{sample}.benchmark'
    shell:
        """
        if [ -d "{params.in_path}" ] && [ -n "$(ls -A {params.in_path})" ]; then
            echo "Found bins in {params.in_path}, running gtdbtk..." >> {log}
            mkdir -p {params.mash_db} 2> {log}
            mkdir -p {params.out_path} 2>> {log}
            echo "[GTDBTK CMD] gtdbtk classify_wf --mash_db {params.mash_db} --genome_dir {params.in_path} --out_dir {params.out_path} -x {params.genometype} --cpus {threads} --pplacer_cpus {threads} {params.extra_args}" >> {log}
            export GTDBTK_DATA_PATH={params.gtdbtk_db}
            gtdbtk classify_wf \
                --mash_db {params.mash_db} \
                --genome_dir {params.in_path} \
                --out_dir {params.out_path} \
                -x {params.genometype} \
                --cpus {threads} \
                --pplacer_cpus {threads} {params.extra_args} >> {log} 2>&1
        else
            echo "No bins found in {params.in_path}, skipping gtdbtk." >> {log}
        fi
        touch {output.mk}
        """

rule m_sample_02_taxonomy_BAT:
    input:
        mk1 = path + 'markers/b-sample/08_bin_prodigal/{sample}.mk',
    output:
        mk = path + 'markers/m-sample/02_taxonomy/BAT/{sample}.mk',
    params:
        out_dir= path + 'results/m-sample/02_taxonomy/BAT/{sample}/',
        bin_folder = path + 'results/b-sample/08_gather_bin_high_quality/bin_genome_fna_renamed/{sample}/',
        protein_faa = path + 'results/b-sample/08_gather_bin_high_quality/bin_protein_faa_renamed/{sample}/',
        bin_faa = temp(path + 'results/m-sample/tmp/{sample}/MAGs_protein_for_BAT.faa'),
        genometype = "fna",
        BAT_range = config["Params_m_sample"]["02_taxonomy_BAT"]["range"],
        BAT_fraction = config["Params_m_sample"]["02_taxonomy_BAT"]["fraction"],
        db = CAT_taxonomy_db,
        tax = CAT_taxonomy_tax,
        extra_args = config["Params_m_sample"]["02_taxonomy_BAT"]["extra_args"],
    conda: config["conda_envs"]["metalink-annotation"]
    priority: 10
    threads: config["Params_m_sample"]["02_taxonomy"]["threads"]
    log:
        path + 'logs/m-sample/02_taxonomy/{sample}.log'
    benchmark:
        path + 'benchmarks/m-sample/02_taxonomy/{sample}.benchmark'
    shell:
        """
        mkdir -p {params.out_dir} 2>{log}
        if [ -d "{params.bin_folder}" ] && [ -n "$(ls -A {params.bin_folder})" ] && [ -d "{params.protein_faa}" ] && [ -n "$(ls -A {params.protein_faa})" ]; then
            echo "Found bins in {params.bin_folder} and {params.protein_faa}, running BAT..." >> {log}
            mkdir -p $(dirname {params.bin_faa}) 2>> {log}
            cat {params.protein_faa}/*.faa > {params.bin_faa} 2>> {log}
            echo "Running BAT pipeline." >> {log} 2>&1
            CAT_pack bins -n {threads} -r {params.BAT_range} -f {params.BAT_fraction} -b {params.bin_folder} --bin_suffix {params.genometype} -p {params.bin_faa} -d {params.db} -t {params.tax} -o {params.out_dir}/BAT --no_stars --force {params.extra_args} >> {log} 2>&1
            CAT_pack add_names -i {params.out_dir}/BAT.ORF2LCA.txt -o {params.out_dir}/BAT.ORF2LCA.name.txt -t {params.tax} --only_official --force >> {log} 2>&1
            CAT_pack add_names -i {params.out_dir}/BAT.bin2classification.txt -o {params.out_dir}/BAT.bin2classification.name.txt -t {params.tax} --only_official --force >> {log} 2>&1
            CAT_pack summarise -i {params.out_dir}/BAT.bin2classification.name.txt -o {params.out_dir}/BAT.summary.txt --force >> {log} 2>&1
            mv {params.out_dir}/BAT.ORF2LCA.name.txt {params.out_dir}/BAT.ORF2LCA.tsv >> {log} 2>&1
            mv {params.out_dir}/BAT.bin2classification.name.txt {params.out_dir}/BAT.bin2classification.tsv >> {log} 2>&1
            rm -f {params.out_dir}/BAT.ORF2LCA.txt >> {log} 2>&1
            rm -f {params.out_dir}/BAT.bin2classification.txt >> {log} 2>&1
        else
            echo "contigs_fna empty/missing for {wildcards.sample}; skip m_sample_02_taxonomy_BAT." >>{log}
            touch {params.out_dir}/BAT.ORF2LCA.tsv \
                {params.out_dir}/BAT.bin2classification.tsv \
                {params.out_dir}/BAT.summary.tsv >> {log} 2>&1
        fi
        touch {output.mk}
        """

rule m_sample_02_taxonomy_done:
    priority: 10
    input:
        mk = get_m_sample_02_taxonomy_classify_rule_mk,
    output:
        mk = path + 'markers/m-sample/02_taxonomy/classify_done/{sample}.mk'
    shell:
        "cp {input.mk} {output.mk}"

rule m_sample_03_genome_stat:
    input:
        mk = path + 'markers/b-sample/08_bin_prodigal/{sample}.mk',
    output:
        mk = path + 'markers/m-sample/03_genome_stat/{sample}.mk'
    threads: 1
    conda: config["conda_envs"]["metalink-main"]
    params:
        in_path = path + 'results/b-sample/08_gather_bin_high_quality/bin_genome_fna_renamed/{sample}/',
        out_path = path + 'results/m-sample/03_genome_stat/',
        output = path + 'results/m-sample/03_genome_stat/{sample}.tsv',
        genometype = "fna",
    priority: 1
    log:
        path + 'logs/m-sample/03_genome_stat/{sample}.log'
    benchmark:
        path + 'benchmarks/m-sample/03_genome_stat/{sample}.benchmark'
    shell:
        """
        if [ -d "{params.in_path}" ] && [ -n "$(ls -A {params.in_path})" ]; then
            echo "Found bins in {params.in_path}, running fasta_stat.py..." >> {log}
            mkdir -p {params.out_path} 2> {log}
            python {Dir_script_metalink}/fasta_stat.py -i {params.in_path} -x {params.genometype} -o {params.output} >> {log} 2>&1
        else
            echo "No bins found in {params.in_path}, skipping fasta_stat.py." >> {log}
        fi
        touch {output.mk}
        """

rule m_sample_04_protein_stat:
    input:
        mk = path + 'markers/b-sample/08_bin_prodigal/{sample}.mk',
        mk2 = path + 'markers/m-sample/01_checkm2/{sample}.mk',
    output:
        mk = path + 'markers/m-sample/04_protein_stat/{sample}.mk'
    threads: 1
    conda: config["conda_envs"]["metalink-main"]
    params:
        in_path = path + 'results/b-sample/08_gather_bin_high_quality/bin_protein_faa_renamed/{sample}/',
        out_path = path + 'results/m-sample/04_protein_stat/',
        output = path + 'results/m-sample/04_protein_stat/{sample}.tsv',
        genometype = "faa",
    priority: 1
    log:
        path + 'logs/m-sample/04_protein_stat/{sample}.log'
    benchmark:
        path + 'benchmarks/m-sample/04_protein_stat/{sample}.benchmark'
    shell:
        """
        if [ -d "{params.in_path}" ] && [ -n "$(ls -A {params.in_path})" ]; then
            echo "Found bins in {params.in_path}, running fasta_stat_faa.py..." >> {log}
            mkdir -p {params.out_path} 2> {log}
            python {Dir_script_metalink}/fasta_stat_faa.py -i {params.in_path} -x {params.genometype} -o {params.output} >> {log} 2>&1
        else
            echo "No bins found in {params.in_path}, skipping fasta_stat_faa.py." >> {log}
        fi
        touch {output.mk}
        """

rule m_sample_05_drep_eachsample_99:
    input:
        mk = path + 'markers/b-sample/08_bin_prodigal/{sample}.mk',
    output:
        mk = path + 'markers/m-sample/05_drep_eachsample_99/{sample}.mk'
    threads: config["Params_m_sample"]["05_drep_eachsample_99"]["threads"]
    params:
        checkm = path + 'results/b-sample/08_gather_bin_high_quality/checkm_stats/{sample}.csv',
        g_path  = path + 'results/m-sample/tmp/05_drep_eachsample_99/{sample}.05_drep_eachsample_99.list',
        in_path = path + 'results/b-sample/08_gather_bin_high_quality/bin_genome_fa_renamed/{sample}/',
        out_path= path + 'results/m-sample/05_drep_eachsample_99/{sample}/',
        comp = config["Params_m_sample"]["05_drep_eachsample_99"]["comp"],
        con  = config["Params_m_sample"]["05_drep_eachsample_99"]["con"],
        sa   = config["Params_m_sample"]["05_drep_eachsample_99"]["sa"],
        extra_args = config["Params_m_sample"]["05_drep_eachsample_99"]["extra_args"],
    conda: config["conda_envs"]["metalink-annotation"]
    priority: 1
    log:
        path + 'logs/m-sample/05_drep_eachsample_99/{sample}.log'
    benchmark:
        path + 'benchmarks/m-sample/05_drep_eachsample_99/{sample}.benchmark'
    shell:
        r"""
        if [ -d "{params.out_path}" ]; then
            rm -r {params.out_path}
            echo "Removed old output directory {params.out_path}" > {log}
        fi
        if [ -d "{params.in_path}" ] && [ -n "$(ls -A {params.in_path})" ]; then
            echo "Found bins in {params.in_path}, running dRep 99..." >> {log}
            mkdir -p $(dirname {params.g_path}) {params.out_path}/dereplicated_genomes 2>> {log}
            find {params.in_path} -maxdepth 1 \( -type f -o -type l \) -name "*.fa" -print > {params.g_path} 2>> {log}
            n=$(wc -l < {params.g_path})
            DREP="dRep dereplicate {params.out_path} -g {params.g_path} -p {threads} \
                -comp {params.comp} -con {params.con} -sa {params.sa} {params.extra_args}"
            if [ -s "{params.checkm}" ]; then
                DREP="$DREP --genomeInfo {params.checkm}"
            fi
            if [ "$n" -eq 0 ]; then
                echo "No bins found in {params.g_path}, please check {params.g_path} and {params.in_path} for {wildcards.sample}." >> {log}
            elif [ "$n" -eq 1 ]; then
                echo "Only one genome found, copying to dereplicated_genomes without dRep 99." >> {log}
                cp $(cat {params.g_path}) {params.out_path}/dereplicated_genomes/
                genome=$(basename $(cat {params.g_path}))
                mkdir -p {params.out_path}/data_tables
                echo "genome,secondary_cluster,threshold,cluster_method,comparison_algorithm,primary_cluster" \
                    > {params.out_path}/data_tables/Cdb.csv
                echo "${{genome}},1_0,0,average,fastANI,1" \
                    >> {params.out_path}/data_tables/Cdb.csv
                echo "genome,cluster,score" \
                    > {params.out_path}/data_tables/Wdb.csv
                echo "${{genome}},1_0,100" \
                    >> {params.out_path}/data_tables/Wdb.csv
            else
                echo "$n bins found, running dRep." >> {log}
                echo "[dRep CMD] $DREP " >> {log}
                $DREP >> {log} 2>&1
            fi
        else
            echo "No bins found in {params.in_path}, skipping dRep 99." >> {log}
        fi
        touch {output.mk}
        """

rule m_sample_06_drep_eachsample_95:
    input:
        mk = path + "markers/m-sample/05_drep_eachsample_99/{sample}.mk"
    output:
        mk = path + 'markers/m-sample/06_drep_eachsample_95/{sample}.mk'
    threads: config["Params_m_sample"]["06_drep_eachsample_95"]["threads"]
    params:
        in_path  = path + 'results/m-sample/05_drep_eachsample_99/{sample}/dereplicated_genomes/',
        g_path   = path + 'results/m-sample/tmp/06_drep_eachsample_95/{sample}.MAGs_06_drep_eachsample_95.list',
        checkm   = path + 'results/m-sample/05_drep_eachsample_99/{sample}/data_tables/genomeInfo.csv',
        out_path = path + 'results/m-sample/06_drep_eachsample_95/{sample}/',
        comp = config["Params_m_sample"]["06_drep_eachsample_95"]["comp"],
        con  = config["Params_m_sample"]["06_drep_eachsample_95"]["con"],
        sa   = config["Params_m_sample"]["06_drep_eachsample_95"]["sa"],
        extra_args = config["Params_m_sample"]["06_drep_eachsample_95"]["extra_args"],
    conda: config["conda_envs"]["metalink-annotation"]
    priority: 1
    log:
        path + 'logs/m-sample/06_drep_eachsample_95/{sample}.log'
    benchmark:
        path + 'benchmarks/m-sample/06_drep_eachsample_95/{sample}.benchmark'
    shell:
        r"""
        if [ -d "{params.out_path}" ]; then
            rm -r {params.out_path}
            echo "Removed old output directory {params.out_path}" >> {log}
        fi
        if [ -d "{params.in_path}" ] && [ -n "$(ls -A {params.in_path})" ]; then
            echo "Found bins in {params.in_path}, running dRep..." >> {log}
            mkdir -p $(dirname {params.g_path}) {params.out_path}/dereplicated_genomes 2> {log}
            find {params.in_path} -maxdepth 1 \( -type f -o -type l \) -name "*.fa" -print > {params.g_path} 2>> {log}
            n=$(wc -l < {params.g_path})
            if [ "$n" -eq 0 ]; then
                echo "No bins found in {params.g_path}, please check {params.g_path} and {params.in_path} for {wildcards.sample}." >> {log}
            elif [ "$n" -eq 1 ]; then
                echo "Only one genome found, copying to dereplicated_genomes without dRep 95." >> {log}
                cp $(cat {params.g_path}) {params.out_path}/dereplicated_genomes/
                # Cdb.csv Wdb.csv
                genome=$(basename $(cat {params.g_path}))
                mkdir -p {params.out_path}/data_tables
                echo "genome,secondary_cluster,threshold,cluster_method,comparison_algorithm,primary_cluster" \
                    > {params.out_path}/data_tables/Cdb.csv
                echo "${{genome}},1_0,0,average,fastANI,1" \
                    >> {params.out_path}/data_tables/Cdb.csv
                echo "genome,cluster,score" \
                    > {params.out_path}/data_tables/Wdb.csv
                echo "${{genome}},1_0,100" \
                    >> {params.out_path}/data_tables/Wdb.csv
            else
                echo "$n bins found, running dRep." >> {log}
                echo "[dRep CMD] dRep dereplicate {params.out_path} -g {params.g_path} --genomeInfo {params.checkm} -p {threads} -comp {params.comp} -con {params.con} -sa {params.sa} {params.extra_args}" >> {log}
                dRep dereplicate {params.out_path} -g {params.g_path} --genomeInfo {params.checkm} -p {threads} -comp {params.comp} -con {params.con} -sa {params.sa} {params.extra_args} >> {log} 2>&1
            fi
        else
            echo "No bins found in {params.in_path}, skipping dRep 95." >> {log}
        fi
        touch {output.mk}
        """

rule m_sample_07_coverm_eachsample:
    input:
        mk = path + 'markers/b-sample/08_bin_prodigal/{sample}.mk',
    output:
        out_coverm = path + 'results/m-sample/07_coverm_eachsample/{sample}.tsv',
    threads: config["Params_m_sample"]["07_coverm_eachsample"]["threads"]
    conda: config["conda_envs"]["metalink-main"]
    priority: 80
    params:
        bam = path + 'results/b-sample/03_mapping/' + base_sample_mapping_method + '/{sample}.bam',
        MAGs_eachsample = path + 'results/b-sample/08_gather_bin_high_quality/bin_genome_fna_raw/{sample}/',
        genometype = "fa",
        method = config["Params_m_sample"]["07_coverm_eachsample"]["method"],
        min_read_percent_identity = config["Params_m_sample"]["07_coverm_eachsample"]["min_read_percent_identity"],
        min_read_aligned_length = config["Params_m_sample"]["07_coverm_eachsample"]["min_read_aligned_length"],
        contig_end_exclusion = config["Params_m_sample"]["07_coverm_eachsample"]["contig_end_exclusion"],
        extra_args = config["Params_m_sample"]["07_coverm_eachsample"]["extra_args"],
    log:
        path + 'logs/m-sample/07_coverm_eachsample/{sample}.log'
    benchmark:
        path + 'benchmarks/m-sample/07_coverm_eachsample/{sample}.benchmark'
    shell:
        """
        if [ -d "{params.MAGs_eachsample}" ]; then
            MAGs_eachsample={params.MAGs_eachsample}
        else
            echo "The directory MAGs_eachsample is not exist."
        fi > {log} 2>&1

        if [ -z "$(ls -A $MAGs_eachsample)" ]; then
            echo "The directory MAGs_eachsample is empty, pass..."
            touch {output.out_coverm}
        else
            echo "Calculating abundance..."
            coverm genome \
                -t {threads} \
                -m {params.method} \
                --min-read-percent-identity {params.min_read_percent_identity} \
                --min-read-aligned-length {params.min_read_aligned_length} \
                --contig-end-exclusion {params.contig_end_exclusion} \
                --bam-files {params.bam} \
                --genome-fasta-directory $MAGs_eachsample \
                -x {params.genometype} \
                -o {output.out_coverm} {params.extra_args}
        fi >> {log} 2>&1
        """

rule m_sample_08_function_emapper:
    priority: 20
    input:
        mk = path + 'markers/m-sample/01_checkm2/{sample}.mk',
    output:
        mk = path + 'markers/m-sample/08_function/emapper/{sample}.mk',
    threads: config["Params_m_sample"]["08_function"]["threads"]
    conda: config["conda_envs"]["metalink-annotation"]
    params:
        protein_faa = path + 'results/b-sample/08_gather_bin_high_quality/bin_protein_faa_renamed/{sample}/',
        all_faa = temp(path + 'results/m-sample/tmp/{sample}/MAGs_protein.faa'),
        seed_ortholog_evalue = config["Params_m_sample"]["08_function_emapper"]["seed_ortholog_evalue"],
        method = config["Params_m_sample"]["08_function_emapper"]["method"],
        query_cover = config["Params_m_sample"]["08_function_emapper"]["query_cover"],
        data_dir = config["Directory"]["database"]["eggnog_db"],
        out = path + 'results/m-sample/08_function/emapper/annotations/{sample}',
        extra_args = config["Params_m_sample"]["08_function_emapper"]["extra_args"],
    log:
        path + 'logs/m-sample/08_function/emapper/{sample}.log'
    benchmark:
        path + 'benchmarks/m-sample/08_function/emapper/{sample}.benchmark'
    shell:
        """
        if [ -d "{params.protein_faa}" ] && [ -n "$(ls -A {params.protein_faa})" ]; then
            echo "Found bins in {params.protein_faa}, running m_sample_08_function_emapper..." >> {log}
            mkdir -p $(dirname {params.all_faa}) 2> {log}
            cat {params.protein_faa}/*.faa > {params.all_faa} 2>> {log}
            mkdir -p $(dirname {params.out}) 2>> {log}
            emapper.py \
                -i {params.all_faa} \
                --cpu {threads} \
                --seed_ortholog_evalue {params.seed_ortholog_evalue} \
                -m {params.method} \
                --query_cover {params.query_cover} \
                --data_dir {params.data_dir} \
                -o {params.out} --override {params.extra_args} >> {log} 2>&1
            pigz -p {threads} -r -f {params.out}.emapper.annotations {params.out}.emapper.hits {params.out}.emapper.seed_orthologs 2>> {log}
        else
            echo "No bins found in {params.protein_faa}, skipping m_sample_08_function_emapper." >> {log}
        fi
        touch {output.mk}
        """

rule m_sample_08_function_kraken2:
    priority: 20
    input:
        mk1 = path + 'markers/b-sample/08_bin_prodigal/{sample}.mk',
    output:
        kraken_report = path + 'results/m-sample/08_function/kraken2/annotations/{sample}.kraken2.report.tsv',
        detail = path + 'results/m-sample/08_function/kraken2/annotations/{sample}.kraken2.annotations.gz',
        mk = path + 'markers/m-sample/08_function/kraken2/annotations/{sample}.mk',
    threads: config["Params_m_sample"]["08_function"]["threads"]
    params:
        protein_ffn = path + 'results/b-sample/08_gather_bin_high_quality/bin_protein_ffn_renamed/{sample}/',
        all_ffn = temp(path + 'results/m-sample/tmp/{sample}/MAGs_protein.ffn'),
        databases_path = config["Directory"]["database"]["kraken2_function_db"],
        confidence = config["Params_m_sample"]["08_function_kraken2"]["confidence"],
        detail = path + 'results/m-sample/08_function/kraken2/annotations/{sample}.kraken2.mag2classification.tsv',
        tmp = path + 'results/m-sample/08_function/kraken2/annotations/{sample}.kraken2.mag2classification.tsv.gz',
        extra_args = config["Params_m_sample"]["08_function_kraken2"]["extra_args"],
    conda: config["conda_envs"]["metalink-main"]
    log:
        path + 'logs/m-sample/08_function/kraken2/{sample}.log'
    benchmark:
        path + 'benchmarks/m-sample/08_function/kraken2/{sample}.benchmark'
    shell:
        """
        if [ -d "{params.protein_ffn}" ] && [ -n "$(ls {params.protein_ffn}/*.ffn 2>/dev/null)" ]; then
            echo "Found bins in {params.protein_ffn}, running m_sample_08_function_kraken2..." >> {log}
            mkdir -p $(dirname {params.all_ffn}) 2> {log}
            cat {params.protein_ffn}/*.ffn > {params.all_ffn} 2>> {log}
            mkdir -p $(dirname {params.detail}) 2>> {log}
            kraken2 {params.extra_args} --threads {threads} \
                --db {params.databases_path} \
                --use-names \
                --confidence {params.confidence} \
                --report {output.kraken_report} \
                {params.all_ffn} 1> {params.detail} 2> {log}
            pigz -f -p {threads} {params.detail} > {params.tmp} 2>> {log}
            mv {params.tmp} {output.detail}
        else
            echo "No bins found in {params.protein_ffn}, skipping m_sample_08_function_kraken2." >> {log}
            touch {output.kraken_report}
            echo "" | pigz -c > {output.detail}
        fi
        touch {output.mk}
        """


rule m_sample_08_function_done:
    priority: 20
    input:
        mk = get_m_sample_08_function_classify_rule_mk,
    output:
        mk = path + 'markers/m-sample/08_function/classify_done/{sample}.mk'
    shell:
        "cp {input.mk} {output.mk}"

rule m_sample_09_summary:
    input:
        mk1 = path + 'markers/m-sample/01_checkm2/{sample}.mk',
        mk2 = path + 'markers/m-sample/02_taxonomy/classify_done/{sample}.mk',
        mk3 = path + 'markers/m-sample/06_drep_eachsample_95/{sample}.mk',
        mk4 = path + 'markers/m-sample/03_genome_stat/{sample}.mk',
        mk5 = path + 'markers/m-sample/04_protein_stat/{sample}.mk',
        mk6 = path + 'markers/b-sample/08_bin_prodigal/{sample}.mk',
        coverm = path + 'results/m-sample/07_coverm_eachsample/{sample}.tsv',
    output:
        genome = path + 'results/m-sample/09_summary/{sample}/mMAG_summary_info.tsv',
        species = path + 'results/m-sample/09_summary/{sample}/mMAG_species_summary_info.tsv',
    threads: 1
    conda: config["conda_envs"]["metalink-main"]
    params:
        checkM = path + 'results/b-sample/08_gather_bin_high_quality/checkm_stats/{sample}.csv',
        checkM2 = path + 'results/m-sample/01_checkm2/{sample}/quality_report.tsv',
        gtdbtk_arc = path + 'results/m-sample/02_taxonomy/gtdbtk/{sample}/gtdbtk.ar53.summary.tsv',
        gtdbtk_bac = path + 'results/m-sample/02_taxonomy/gtdbtk/{sample}/gtdbtk.bac120.summary.tsv',
        taxid_tax = taxid_tax_path,
        MAGs_taxonomy_classify_method = MAGs_taxonomy_classify_method,
        MAGs_tax_detail = path + 'results/m-sample/02_taxonomy/' + MAGs_taxonomy_classify_method + '/{sample}/' + MAGs_taxonomy_classify_method + '.bin2classification.tsv',
        MAGs_03_genome_stat = path + 'results/m-sample/03_genome_stat/{sample}.tsv',
        MAGs_04_protein_stat = path + 'results/m-sample/04_protein_stat/{sample}.tsv',
        drep_99_Cdb = path + 'results/m-sample/05_drep_eachsample_99/{sample}/data_tables/Cdb.csv',
        drep_95_Cdb = path + 'results/m-sample/06_drep_eachsample_95/{sample}/data_tables/Cdb.csv',
        drep_95_Wdb = path + 'results/m-sample/06_drep_eachsample_95/{sample}/data_tables/Wdb.csv',
        drep_99_Wdb = path + 'results/m-sample/05_drep_eachsample_99/{sample}/data_tables/Wdb.csv',
    priority: 1
    log:
        path + 'logs/m-sample/09_summary/{sample}.log'
    benchmark:
        path + 'benchmarks/m-sample/09_summary/{sample}.benchmark'
    shell:
        """
        if [ "{params.MAGs_taxonomy_classify_method}" != "gtdbtk" ] && [ "{params.MAGs_taxonomy_classify_method}" != "BAT" ]; then
            echo "Error: Invalid MAGs_taxonomy_classify_method '{params.MAGs_taxonomy_classify_method}'. Only 'gtdbtk' or 'BAT' are allowed." >> {log}
            exit 1
        fi

        if [ "{params.MAGs_taxonomy_classify_method}" = "gtdbtk" ]; then
            echo "Processing GTDB-Tk taxonomy results for {wildcards.sample}..." >> {log}
            if [ -f {params.gtdbtk_arc} ]; then
                if [ -f {params.gtdbtk_bac} ]; then
                    cat {params.gtdbtk_arc} > {params.MAGs_tax_detail}
                    sed '1d' {params.gtdbtk_bac} >> {params.MAGs_tax_detail}
                else
                    echo "The directory of bacteria gtdbtk results is empty." >> {log}
                    cat {params.gtdbtk_arc} > {params.MAGs_tax_detail}
                fi
            else
                if [ -f {params.gtdbtk_bac} ]; then
                    echo "The directory of archaea gtdbtk results is empty." >> {log}
                    cat {params.gtdbtk_bac} > {params.MAGs_tax_detail}
                else
                    echo "The directory of archaea and bacteria gtdbtk results is empty." >> {log}
                    mkdir -p $(dirname {params.MAGs_tax_detail})
                    touch {params.MAGs_tax_detail}
                fi
            fi 2>> {log}
        else
            echo "Using BAT taxonomy results for {wildcards.sample}..." >> {log}
            if [ ! -f {params.MAGs_tax_detail} ]; then
                echo "Error: BAT taxonomy file not found: {params.MAGs_tax_detail}" >> {log}
                exit 1
            fi
        fi

        if [ ! -s {params.checkM2} ] || [ ! -s {params.MAGs_tax_detail} ] || [ ! -s {params.checkM} ] || \
        [ ! -s {params.MAGs_03_genome_stat} ] || [ ! -s {params.MAGs_04_protein_stat} ] || \
        [ ! -s {input.coverm} ] || [ ! -s {params.drep_99_Cdb} ] || [ ! -s {params.drep_95_Cdb} ] || \
        [ ! -s {params.drep_95_Wdb} ] || [ ! -s {params.drep_99_Wdb} ]; then
            echo "One or more input files empty/missing for {wildcards.sample}, skipping summary." >>{log}
            touch {output.genome} {output.species}
        else
            echo "All inputs non-empty, running combine_metadata_MAGs." >>{log}
            Rscript {Dir_script_metalink}/combine_metadata_MAGs.v6.cmd.R \
                --checkM2 {params.checkM2} \
                --in_taxid_tax {params.taxid_tax} \
                --taxonomy_classify_method {params.MAGs_taxonomy_classify_method} \
                --taxonomy {params.MAGs_tax_detail} \
                --checkM {params.checkM} \
                --genome_stat {params.MAGs_03_genome_stat} \
                --protein_stat {params.MAGs_04_protein_stat} \
                --in_coverm_trimmed_mean {input.coverm} \
                --drep_99_Cdb {params.drep_99_Cdb} \
                --drep_95_Cdb {params.drep_95_Cdb} \
                --drep_95_Wdb {params.drep_95_Wdb} \
                --drep_99_Wdb {params.drep_99_Wdb} \
                --out_genome {output.genome} \
                --out_species {output.species} >>{log} 2>&1
        fi
        """

# Extra inputs guaranteeing every consumer of the b-sample BAM in the DAG has
# finished before deletion:
# - mcr/vmcr: the unmapped-read extraction reads the BAM (via params); its
#   output and its own input (04_bam_summary stats) are added.
# - 07_contigs (jgi depth) reads the BAM and is in the DAG of every mode with
#   m-s.smk except pmr without project level (there nothing pulls it in, so
#   adding it would needlessly trigger the step).
m_sample_12_delete_metawrap_bam_extra_inputs = {}
if config.get('mode') in ('mcr', 'vmcr'):
    m_sample_12_delete_metawrap_bam_extra_inputs['stat_raw'] = path + 'results/b-sample/04_bam_summary/{sample}.contigs_bam_stats.tsv'
    m_sample_12_delete_metawrap_bam_extra_inputs['out_unmapped'] = path + 'results/' + config.get('mode') + '-sample/01_unmapped_reads/{sample}_unmapped_s.fastq.gz'
if config.get('mode') != 'pmr' or RUN_PROJECT_LEVEL:
    m_sample_12_delete_metawrap_bam_extra_inputs['contigs_fna'] = path + 'results/b-sample/07_contigs/contigs_fna/{sample}.fna'

rule m_sample_12_delete_metawrap_bam:
    input:
        out_coverm = path + 'results/m-sample/07_coverm_eachsample/{sample}.tsv',
        **m_sample_12_delete_metawrap_bam_extra_inputs,
    output:
        mk = path + 'markers/b-sample/10_delete_metawrap_bam/{sample}.mk'
    threads: 1
    priority: 100
    params:
        sn = '{sample}',
        bam = path + 'results/b-sample/03_mapping/' + base_sample_mapping_method + '/{sample}.bam',
    log:
        path + 'logs/b-sample/10_delete_metawrap_bam/{sample}.log'
    benchmark:
        path + 'benchmarks/b-sample/10_delete_metawrap_bam/{sample}.benchmark'
    shell:
        """
        if [ -e "{params.bam}" ];then rm -r {params.bam}*;else echo Folder trimmed.bam not exist.;fi > {log} 2>&1
        touch {output.mk}
        """

# The dRep data/ subdirectories (MASH/prodigal/checkM working files) are never
# consumed downstream; only dereplicated_genomes/ and data_tables/ are results.
# Delete data/ once the sample-level summary and, when the project level runs
# (modes with m-p.smk), the project-level dRep that gathers each sample's
# genomes have finished.
m_sample_12_delete_drep_data_inputs = [path + 'results/m-sample/09_summary/{sample}/mMAG_summary_info.tsv']
if RUN_PROJECT_LEVEL and config.get('mode') in ('m', 'pmr', 'vm'):
    m_sample_12_delete_drep_data_inputs += expand(path + 'markers/m-project/01_drep_project_99/{project}.mk', project=PROJECT_INDEX)
    m_sample_12_delete_drep_data_inputs += expand(path + 'markers/m-project/01_drep_project_95/{project}.mk', project=PROJECT_INDEX)

rule m_sample_12_delete_drep_data:
    input:
        m_sample_12_delete_drep_data_inputs,
    output:
        mk = path + 'markers/m-sample/12_delete_drep_data/{sample}.mk'
    threads: 1
    priority: 100
    params:
        drep_99_data = path + 'results/m-sample/05_drep_eachsample_99/{sample}/data',
        drep_95_data = path + 'results/m-sample/06_drep_eachsample_95/{sample}/data',
    log:
        path + 'logs/m-sample/12_delete_drep_data/{sample}.log'
    benchmark:
        path + 'benchmarks/m-sample/12_delete_drep_data/{sample}.benchmark'
    shell:
        """
        if [ -e "{params.drep_99_data}" ];then rm -r {params.drep_99_data};else echo Folder drep_eachsample_99.data not exist.;fi > {log} 2>&1
        if [ -e "{params.drep_95_data}" ];then rm -r {params.drep_95_data};else echo Folder drep_eachsample_95.data not exist.;fi >> {log} 2>&1
        touch {output.mk}
        """
