'''
/*
 * @Date: 2022-12-04 17:15:47
 * @LastEditors: LiuyangLi
 * @LastEditTime: 2024-07-20 13:03:29
 * @FilePath: m-p.smk
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

folder_checks = [
                path + 'results/',
                path + 'results/m-project/',
                ]
for item in folder_checks:
    if not os.path.isdir(item): os.mkdir(item)

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


rule m_project_01_drep_project_99:
    input:
        mk_06 = expand(path + 'markers/m-sample/05_drep_eachsample_99/{sample}.mk', sample=SAMPLE_INDEX),
    output:
        mk = path + 'markers/m-project/01_drep_project_99/{project}.mk',
    threads: config["Params_m_project"]["01_drep_project_99"]["threads"]
    params:
        in_path = expand(path + 'results/m-sample/05_drep_eachsample_99/{sample}/dereplicated_genomes/',sample=SAMPLE_INDEX),
        each_checkm = expand(path + 'results/b-sample/08_gather_bin_high_quality/checkm_stats/{sample}.csv',sample=SAMPLE_INDEX),
        g_path = path + 'results/m-project/tmp/{project}/genome_path_drep99.list',
        cat_checkm = path + 'results/m-project/tmp/{project}/genomeInfo_drep99.csv',
        out_path = path + 'results/m-project/01_drep_project_99/{project}/',
        comp = config["Params_m_project"]["01_drep_project_99"]["comp"],
        con = config["Params_m_project"]["01_drep_project_99"]["con"],
        sa = config["Params_m_project"]["01_drep_project_99"]["sa"],
        extra_args = config["Params_m_project"]["01_drep_project_99"]["extra_args"],
    conda: config["conda_envs"]["metalink-annotation"]
    priority: 1
    log:
        path + 'logs/m-project/01_drep_project_99/{project}.log'
    benchmark:
        path + 'benchmarks/m-project/01_drep_project_99/{project}.benchmark'
    shell:
        """
        if [ -d "{params.out_path}" ]; then
            rm -r {params.out_path}
            echo "Removed old output directory {params.out_path}" > {log}
        fi
        mkdir -p $(dirname {params.g_path}) {params.out_path}/dereplicated_genomes 2>> {log}
        rm -f "{params.g_path}"
        for d in {params.in_path}; do
            if [ -d "$d" ]; then
                find "$d" -maxdepth 1 \( -type f -o -type l \) -name "*.fa" -print >> "{params.g_path}"
            else
                echo "[WARN] Skip missing dereplicated_genomes dir: '$d'" >> {log}
            fi
        done 2>> {log}

        rm -f "{params.cat_checkm}" >>{log} 2>&1
        first=true
        for f in {params.each_checkm}; do
            if [ -s "$f" ]; then
                if $first; then
                    cat "$f" >> "{params.cat_checkm}"
                    first=false
                else
                    tail -n +2 "$f" >> "{params.cat_checkm}"
                fi
            else
                echo "[WARN] Skip invalid or empty each_checkm file: '$f'" >> {log}
            fi
        done 2>> {log}

        n=$(wc -l < {params.g_path})
        if [ "$n" -eq 0 ]; then
            echo "No genomes found for project {wildcards.project}, skipping dRep." >> {log}
        elif [ "$n" -eq 1 ]; then
            echo "Only one genome found, copying to dereplicated_genomes without dRep." >> {log}
            cp $(cat {params.g_path}) {params.out_path}/dereplicated_genomes/
            # Cdb.csv Wdb.csv
            genome=$(basename $(cat {params.g_path}))
            mkdir -p {params.out_path}/data_tables
            echo "genome,secondary_cluster,threshold,cluster_method,comparison_algorithm,primary_cluster" \\
                > {params.out_path}/data_tables/Cdb.csv
            echo "${{genome}},1_0,0,average,fastANI,1" \\
                >> {params.out_path}/data_tables/Cdb.csv
            echo "genome,cluster,score" \\
                > {params.out_path}/data_tables/Wdb.csv
            echo "${{genome}},1_0,100" \\
                >> {params.out_path}/data_tables/Wdb.csv
        else
            echo "$n genomes found, running dRep." >> {log}
            echo "[dRep CMD] dRep dereplicate {params.out_path} -g {params.g_path} --genomeInfo {params.cat_checkm} -p {threads} -comp {params.comp} -con {params.con} -sa {params.sa} {params.extra_args}" >> {log}
            dRep dereplicate {params.out_path} -g {params.g_path} --genomeInfo {params.cat_checkm} -p {threads} -comp {params.comp} -con {params.con} -sa {params.sa} {params.extra_args} >> {log} 2>&1
        fi
        touch {output.mk}
        """

rule m_project_01_drep_project_95:
    input:
        mk = path + 'markers/m-project/01_drep_project_99/{project}.mk',
    output:
        mk = path + 'markers/m-project/01_drep_project_95/{project}.mk',
    threads: config["Params_m_project"]["01_drep_project_95"]["threads"]
    conda: config["conda_envs"]["metalink-annotation"]
    params:
        in_path = path + 'results/m-project/01_drep_project_99/{project}/dereplicated_genomes/',
        g_path = path + 'results/m-project/tmp/{project}/genome_path_drep95.txt',
        checkm = path + 'results/m-project/01_drep_project_99/{project}/data_tables/genomeInfo.csv',
        out_path = path + 'results/m-project/01_drep_project_95/{project}/',
        comp = config["Params_m_project"]["01_drep_project_95"]["comp"],
        con = config["Params_m_project"]["01_drep_project_95"]["con"],
        sa = config["Params_m_project"]["01_drep_project_95"]["sa"],
        extra_args = config["Params_m_project"]["01_drep_project_95"]["extra_args"],
    priority: 1
    log:
        path + 'logs/m-project/01_drep_project_95/{project}.log'
    benchmark:
        path + 'benchmarks/m-project/01_drep_project_95/{project}.benchmark'
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
                echo "No bins found in {params.g_path}, please check {params.g_path} and {params.in_path} for {wildcards.project}." >> {log}
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

rule m_project_04_summary:
    input:
        mk = path + 'markers/m-project/01_drep_project_95/{project}.mk',
        sample_function_mk = expand(path + 'markers/m-sample/08_function/classify_done/{sample}.mk', sample=SAMPLE_INDEX),
        each_info = expand(path + 'results/m-sample/09_summary/{sample}/mMAG_summary_info.tsv', sample=SAMPLE_INDEX),
        each_MAGs_ffn_length = expand(path + 'results/b-sample/08_gather_bin_high_quality/bin_ffn_length/{sample}', sample=SAMPLE_INDEX),
    output:
        cat_info = path + 'results/m-project/04_summary/{project}/mMAG_summary_info.tsv',
    threads: 1
    conda: config["conda_envs"]["metalink-main"]
    params:
        cat_info = path + 'results/m-project/tmp/{project}/mMAG_summary_info_tmp.tsv',
        drep_project_99_Cdb = path + 'results/m-project/01_drep_project_99/{project}/data_tables/Cdb.csv',
        drep_project_95_Cdb = path + 'results/m-project/01_drep_project_95/{project}/data_tables/Cdb.csv',
        drep_project_95_Wdb = path + 'results/m-project/01_drep_project_95/{project}/data_tables/Wdb.csv',
        drep_project_99_Wdb = path + 'results/m-project/01_drep_project_99/{project}/data_tables/Wdb.csv',
        MAGs_function_classify_method = MAGs_function_classify_method,
        MAGs_fun_detail = expand(path + 'results/m-sample/08_function/' + MAGs_function_classify_method + '/annotations/{sample}.' + MAGs_function_classify_method + '.annotations.gz',sample=SAMPLE_INDEX),
        cat_fun_detail = path + 'results/m-project/04_summary/{project}/MAGs_protein.' + MAGs_function_classify_method + '.annotations',
        cat_info_bin_ffn_length = path + 'results/m-project/04_summary/{project}/bin_ffn_length',
    priority: 20
    log:
        path + 'logs/m-project/04_summary/{project}.log'
    benchmark:
        path + 'benchmarks/m-project/04_summary/{project}.benchmark'
    shell:
        """
        rm -f "{params.cat_fun_detail}" >{log} 2>&1
        if [ "{params.MAGs_function_classify_method}" = "emapper" ]; then
            #######  emapper
            first=true
            for f in {params.MAGs_fun_detail}; do
                if [ -s "$f" ]; then
                    num=3 # 要删除的行数
                    max=$(zcat "$f" | wc -l)
                    let sLine=max-num+1   # 删除的起始行
                    if $first; then
                        zcat "$f" | sed $sLine',$d' | sed '1,4d' >> "{params.cat_fun_detail}" #从起始行删除到最后行,再删除前5行
                        first=false
                    else
                        zcat "$f" | sed $sLine',$d' | sed '1,5d' >> "{params.cat_fun_detail}" #从起始行删除到最后行,再删除前5行
                    fi
                else
                    echo "[WARN] Skip invalid or empty emapper.annotations.gz file: '$f'" >> {log}
                fi
            done 2>> {log}
        elif [ "{params.MAGs_function_classify_method}" = "kraken2" ]; then
            #######  kraken2
            for f in {params.MAGs_fun_detail}; do
                if [ -s "$f" ]; then
                    zcat "$f" >> "{params.cat_fun_detail}" 2>> {log}
                else
                    echo "[WARN] Skip invalid or empty kraken2.annotations.gz file: '$f'" >> {log}
                fi
            done 2>> {log}
        else
            echo "[ERROR] Unknown MAGs_function_classify_method: '{params.MAGs_function_classify_method}'" >> {log}
        fi
        pigz -p {threads} -r -f {params.cat_fun_detail} 2>> {log}

        #######  bin_ffn_length
        rm -f "{params.cat_info_bin_ffn_length}" >>{log} 2>&1
        first=true
        for f in {input.each_MAGs_ffn_length}; do
            if [ -s "$f" ]; then
                if $first; then
                    cat "$f" >> "{params.cat_info_bin_ffn_length}"
                    first=false
                else
                    tail -n +2 "$f" >> "{params.cat_info_bin_ffn_length}"
                fi
            else
                echo "[WARN] Skip invalid or empty bin_ffn_length file: '$f'" >> {log}
            fi
        done 2>> {log}

        pigz -p {threads} -r -f {params.cat_info_bin_ffn_length} 2>> {log}

        #######  mMAG_summary_info
        rm -f "{params.cat_info}" >>{log} 2>&1
        first=true
        for f in {input.each_info}; do
            if [ -s "$f" ]; then
                if $first; then
                    cat "$f" >> "{params.cat_info}"
                    first=false
                else
                    tail -n +2 "$f" >> "{params.cat_info}"
                fi
            else
                echo "[WARN] Skip invalid or empty mMAG_summary_info.tsv file: '$f'" >> {log}
            fi
        done 2>> {log}

        if [ ! -s {params.cat_info} ] || [ ! -s {params.drep_project_99_Cdb} ] || \\
        [ ! -s {params.drep_project_95_Cdb} ] || [ ! -s {params.drep_project_95_Wdb} ]; then
            echo "One or more input files empty/missing for {wildcards.project}, skipping summary." >>{log}
            touch {output.cat_info}
        else
            Rscript {Dir_script_metalink}/combine_metadata_dRep99_95.v2.cmd.R \
                --metadata {params.cat_info} \
                --drep_99_Cdb {params.drep_project_99_Cdb} \
                --drep_95_Cdb {params.drep_project_95_Cdb} \
                --drep_95_Wdb {params.drep_project_95_Wdb} \
                --drep_99_Wdb {params.drep_project_99_Wdb} \
                -o {output.cat_info} >> {log} 2>&1
        fi
        """

# The project-level dRep data/ subdirectories (MASH/prodigal/checkM working
# files) are never consumed downstream; only dereplicated_genomes/ and
# data_tables/ are results. m_project_04_summary is the last step of the
# m-project chain that reads the dRep outputs (data_tables/, kept).
rule m_project_12_delete_drep_data:
    input:
        cat_info = path + 'results/m-project/04_summary/{project}/mMAG_summary_info.tsv',
    output:
        mk = path + 'markers/m-project/12_delete_drep_data/{project}.mk'
    threads: 1
    priority: 100
    params:
        drep_99_data = path + 'results/m-project/01_drep_project_99/{project}/data',
        drep_95_data = path + 'results/m-project/01_drep_project_95/{project}/data',
    log:
        path + 'logs/m-project/12_delete_drep_data/{project}.log'
    benchmark:
        path + 'benchmarks/m-project/12_delete_drep_data/{project}.benchmark'
    shell:
        """
        if [ -e "{params.drep_99_data}" ];then rm -r {params.drep_99_data};else echo Folder drep_project_99.data not exist.;fi > {log} 2>&1
        if [ -e "{params.drep_95_data}" ];then rm -r {params.drep_95_data};else echo Folder drep_project_95.data not exist.;fi >> {log} 2>&1
        touch {output.mk}
        """
