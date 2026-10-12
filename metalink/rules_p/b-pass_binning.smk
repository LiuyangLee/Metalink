'''
/*
 * @Date: 2022-12-04 17:15:47
 * @LastEditors: LiuyangLi
 * @LastEditTime: 2024-10-25 16:57:02
 * @FilePath: b-pass_binning.smk
 * @Description: base mode (pass metawrap-based binning)
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
                path + 'results/b-sample/06_bin_checkm/',
                path + 'results/b-sample/08_gather_bin_high_quality/',
                path + 'results/b-sample/08_gather_bin_high_quality/bin_genome_fna_renamed/',
                path + 'results/b-sample/08_gather_bin_high_quality/checkm_stats/',
                path + 'results/b-project/'
                ]

for item in folder_checks:
    if not os.path.isdir(item): os.mkdir(item)

rule b_sample_06_bin_checkm:
    input:
        in_path = lambda wildcards: SAMPLE_PATHS[wildcards.sample]["mMAG_dir"],
    output:
        mk = path + 'markers/b-sample/06_bin_checkm/{sample}.mk'
    threads: config["Params_b_sample"]["06_bin_checkm"]["threads"]
    conda: config["conda_envs"]["metalink-binning"]
    params:
        out_path = path + 'results/b-sample/06_bin_checkm/{sample}/',
        stats = path + 'results/b-sample/06_bin_checkm/{sample}/checkm_bins.Chdb.tsv',
        genometype = 'fa',
    priority: 60
    log:
        path + 'logs/b-sample/06_bin_checkm/{sample}.bin_checkm.log'
    benchmark:
        path + 'benchmarks/b-sample/06_bin_checkm/{sample}.bin_checkm.benchmark'
    shell:
        r"""
        if [ -d "{input.in_path}" ] && [ -n "$(find -L {input.in_path} -maxdepth 1 -type f \( -name '*.fasta' -o -name '*.fa' -o -name '*.fna' \) -print -quit)" ]; then
            echo "Found bins in {input.in_path}, running checkm..." > {log}
            mkdir -p {params.out_path} 2>> {log}
            checkm lineage_wf -x {params.genometype} {input.in_path} {params.out_path} -t {threads} --pplacer_threads {threads} >> {log} 2>&1
            checkm qa --tab_table {params.out_path}/lineage.ms {params.out_path} -f {params.stats} -t {threads} -o 2 >> {log} 2>&1
        else
            echo "No bins (.fasta/.fa/.fna) found in {input.in_path}, skipping checkm." >> {log}
        fi
        touch {output.mk}
        """

rule b_sample_08_gather_bin_checkm_stats:
    priority: 60
    input:
        mk_4 = path + 'markers/b-sample/06_bin_checkm/{sample}.mk',
    output:
        tsv = path + 'results/b-sample/08_gather_bin_high_quality/checkm_stats/{sample}.csv',
    params:
        sample = '{sample}',
        input_bin_dir = lambda wildcards: SAMPLE_PATHS[wildcards.sample]["mMAG_dir"],
        min_completeness = config["Params_b_sample"]["06_bin_checkm"].get("min_completeness", 50),
        max_contamination = config["Params_b_sample"]["06_bin_checkm"].get("max_contamination", 10),
        min_quality_score = config["Params_b_sample"]["06_bin_checkm"].get("min_quality_score", 50),
        stats = path + 'results/b-sample/06_bin_checkm/{sample}/checkm_bins.Chdb.tsv',
        bin_high_quality = path + 'results/b-sample/08_gather_bin_high_quality/bin_genome_fna_raw/{sample}/',
        bin_high_quality_renamed = path + 'results/b-sample/08_gather_bin_high_quality/bin_genome_fna_renamed/{sample}/',
    threads: config["Params_b_sample"]["08_bin_prodigal"]["threads"]
    conda: config["conda_envs"]["metalink-main"]
    log:
        path + 'logs/b-sample/08_gather_bin_high_quality/checkm_stats/{sample}.gather_bin_checkm_tsv.log'
    benchmark:
        path + 'benchmarks/b-sample/08_gather_bin_high_quality/checkm_stats/{sample}.gather_bin_checkm_tsv.benchmark'
    shell:
        """
        echo "[gather_bin_checkm_stats PARAMS] min_completeness={params.min_completeness}, max_contamination={params.max_contamination}, min_quality_score={params.min_quality_score}" > {log}
        if [ -d "{params.bin_high_quality}" ];then rm -r {params.bin_high_quality};else echo "Folder {params.bin_high_quality} not exist, gather_bin_checkm_stats good to go." >> {log};fi
        if [ -d "{params.bin_high_quality_renamed}" ];then rm -r {params.bin_high_quality_renamed};else echo "Folder {params.bin_high_quality_renamed} not exist, gather_bin_checkm_stats good to go." >> {log};fi
        STATS="{params.stats}"
        if [ -f "$STATS" ]; then
            mkdir -p {params.bin_high_quality}
            mkdir -p {params.bin_high_quality_renamed}
            python {Dir_script_metalink}/gather_bins_checkm.no_add_name.py {params.sample} "$STATS" {params.input_bin_dir} {params.bin_high_quality} {params.bin_high_quality_renamed} {output.tsv} {params.min_completeness} {params.max_contamination} {params.min_quality_score} >> {log} 2>&1
        else
            echo "File $STATS not exist" >&2
            touch {output.tsv}
        fi >> {log} 2>&1
        """

rule b_sample_08_bin_prodigal:
    priority: 80
    input:
        tsv = path + 'results/b-sample/08_gather_bin_high_quality/checkm_stats/{sample}.csv',
    output:
        bin_ffn_length = path + 'results/b-sample/08_gather_bin_high_quality/bin_ffn_length/{sample}',
        ORF_loci = path + 'results/b-sample/08_gather_bin_high_quality/bin_ORF_loci/{sample}',
        mk = path + 'markers/b-sample/08_bin_prodigal/{sample}.mk',
    conda: config["conda_envs"]["metalink-main"]
    threads: config["Params_b_sample"]["08_bin_prodigal"]["threads"]
    params:
        genome_fa = path + 'results/b-sample/08_gather_bin_high_quality/bin_genome_fa_renamed/{sample}',
        protein_ffn = path + 'results/b-sample/08_gather_bin_high_quality/bin_protein_ffn_renamed/{sample}',
        protein_faa = path + 'results/b-sample/08_gather_bin_high_quality/bin_protein_faa_renamed/{sample}',
        bin_high_quality_renamed = path + 'results/b-sample/08_gather_bin_high_quality/bin_genome_fna_renamed/{sample}',
        all_ffn = temp(path + 'results/b-sample/08_gather_bin_high_quality/tmp/{sample}.MAGs_protein.ffn'),
        extra_args = config["Params_b_sample"]["08_bin_prodigal"]["extra_args"],
    log:
        path + 'logs/b-sample/08_bin_prodigal/{sample}.gather_bin_checkm_tsv.log'
    benchmark:
        path + 'benchmarks/b-sample/08_bin_prodigal/{sample}.gather_bin_checkm_tsv.benchmark'
    shell:
        r"""
        echo "[bin_prodigal PARAMS] threads={threads}" > {log}
        if [ -d "{params.protein_ffn}" ];then rm -r {params.protein_ffn};else echo "Folder {params.protein_ffn} not exist, bin_prodigal good to go." >> {log};fi
        if [ -d "{params.protein_faa}" ];then rm -r {params.protein_faa};else echo "Folder {params.protein_faa} not exist, bin_prodigal good to go." >> {log};fi
        if [ -d "{params.genome_fa}" ];then rm -r {params.genome_fa};else echo "Folder {params.genome_fa} not exist, bin_prodigal good to go." >> {log};fi
        mkdir -p {params.protein_ffn} {params.protein_faa} {params.genome_fa}

        find -L {params.bin_high_quality_renamed} -maxdepth 1 -type f \( -name "*.fa" -o -name "*.fasta" -o -name "*.fna" \) | while IFS= read -r genome; do
            [ -e "$genome" ] || continue
            base=$(basename "${{genome%.*}}")
            ln -sf "$genome" {params.genome_fa}/"${{base}}".fa
            prodigal {params.extra_args} \
                -p single \
                -q \
                -i "$genome" \
                -d {params.protein_ffn}/"${{base}}".ffn \
                -a {params.protein_faa}/"${{base}}".faa -o /dev/null
        done >> {log} 2>&1

        mkdir -p $(dirname {params.all_ffn}) 2>> {log}
        if find {params.protein_ffn} -maxdepth 1 -name "*.ffn" -print -quit | grep -q .; then
            find {params.protein_ffn} -maxdepth 1 -name "*.ffn" -exec cat {{}} + > {params.all_ffn} 2>> {log}
            seqkit fx2tab -l -n -i -H -j {threads} {params.all_ffn} > {output.bin_ffn_length} 2>> {log}
        else
            echo "No .ffn files found for sample {wildcards.sample}" >> {log}
            touch {params.all_ffn}
            touch {output.bin_ffn_length}
        fi

        awk 'BEGIN {{OFS="\t"; print "ORF_name", "start", "end", "strand"}}
        />/ {{
            split($0, parts, " # ")
            orf_name = parts[1]
            sub(/^>/, "", orf_name)
            print orf_name, parts[2], parts[3], parts[4]
        }}' {params.all_ffn} > {output.ORF_loci}
        touch {output.mk}
        """

rule b_project_05_gather_project_bin_checkm_stats:
    input:
        tsv = expand(path + 'results/b-sample/08_gather_bin_high_quality/checkm_stats/{sample}.csv', sample=SAMPLE_INDEX),
    output:
        checkM = path + 'results/b-project/05_gather_project_bin_checkm_stats/{project}/05_checkm_high_quality_bins.csv',
        mk = path + 'markers/b-project/05_gather_project_bin_checkm_stats/{project}/05_gather_project_bin_checkm_stats.mk',
    threads: config["Params_b_project"]["05_gather_project_bin_checkm_stats"]["threads"]
    log:
        path + 'logs/b-project/05_gather_project_bin_checkm_stats/{project}/gather_project_bin_checkm_stats.log'
    benchmark:
        path + 'benchmarks/b-project/05_gather_project_bin_checkm_stats/{project}/gather_project_bin_checkm_stats.benchmark'
    shell:
        """
        touch {output.checkM}
        echo "genome,completeness,contamination,strain_heterogeneity,size,N50,QS" > {output.checkM} 2> {log}
        for f in {input.tsv}; do
            tail -n +2 "$f" >> {output.checkM}
        done 2>>{log}
        touch {output.mk}
        """
