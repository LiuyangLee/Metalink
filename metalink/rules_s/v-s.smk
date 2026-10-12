'''
/*
 * @Date: 2022-12-04 17:15:47
 * @LastEditors: LiuyangLi
 * @LastEditTime: 2024-10-25 16:57:02
 * @FilePath: v-s.smk
 * @Description: merged reads to assembly and short reads annotation
 */
'''
import os
import sys

Dir_script_metalink = config["Directory"]["script"]["metalink"]
kegg_modules_ko_info_path = config["Directory"]["script"]["kegg_modules_ko_info_path"]
conda_env_python = config["conda_envs"]["metalink-main"]
conda_env_seqkit = config["conda_envs"]["metalink-main"]
conda_env_checkv = config["conda_envs"]["metalink-annotation"]
conda_env_vrhyme = config["conda_envs"]["metalink-vrhyme"]
conda_env_virsorter2 = config["conda_envs"]["metalink-virsorter"]
conda_env_vibrant = config["conda_envs"]["metalink-vibrant"]
conda_env_genomad = config["conda_envs"]["metalink-genomad"]
conda_env_iphop = config["conda_envs"]["metalink-iphop"]
conda_env_gtdbtk = config["conda_envs"]["metalink-gtdbtk"]
gtdbtk_db = config["Directory"]["database"]["gtdbtk_db"]
checkv_db = config["Directory"]["database"]["checkv_db"]
genomad_db = config["Directory"]["database"]["genomad_db"]
virsorter2_db = config["Directory"]["database"]["virsorter2_db"]
vibrant_db = config["Directory"]["database"]["vibrant_db"]
iphop_db = config["Directory"]["database"]["iphop_db"]
kraken2_taxonomy_db = config["Directory"]["database"]["kraken2_taxonomy_db"]
CAT_taxonomy_virus_db = config["Directory"]["database"]["CAT_taxonomy_virus_db"]
CAT_taxonomy_virus_tax = config["Directory"]["database"]["CAT_taxonomy_virus_tax"]
identify_method = config["Params_v_sample"]["01_virus_identify"]["identify_method"]
virus_length_limit = config["Params_v_sample"]["01_virus_identify"]["virus_length_limit"]

path = config["Directory"]["work"]
path = os.path.join(path, '')
os.makedirs(path, exist_ok=True)
os.chdir(path)

folder_checks = [
                path + 'results/',
                path + 'results/v-sample/',
                path + 'results/v-sample/01_identification',
                path + 'results/v-sample/03_mapping',
                path + 'results/v-sample/03_taxonomy',
                path + 'results/v-sample/04_vRhyme',
                path + 'results/v-sample/05_CheckV',
                path + 'results/v-sample/06_vclust',
                path + 'results/v-sample/07_iPHoP',
                path + 'results/v-sample/08_taxonomy',
                path + 'results/v-sample/09_coverm',
                path + 'results/v-sample/10_summary',
                path + 'results/v-sample/11_visualization',
                ]
for item in folder_checks:
    if not os.path.isdir(item): os.mkdir(item)

#identify_method=vb,vs_vb,genomad,genomad_vs_vb
def get_v_sample_01_taxonomy_identify_rule_mk(wildcards):
    if identify_method == "genomad":
        return rules.v_sample_01_identify_genomad.output.mk
    elif identify_method == "vs_vb":
        return rules.v_sample_01_identify_vs_vb.output.mk
    elif identify_method == "genomad_vs_vb":
        return rules.v_sample_01_combine_genomad_vs_vb.output.mk
    else:
        raise ValueError(f"Unknown identify_method: {identify_method}")

#v_sample_mapping_method=bwamem,strobealign,minibwa
v_sample_mapping_method = config["Params_v_sample"]["03_mapping"]["method"]
rule v_sample_01_identify_done:
    input:
        mk = get_v_sample_01_taxonomy_identify_rule_mk,
    output:
        mk = path + 'markers/v-sample/01_identification/{sample}.identify_done.mk'
    shell:
        "cp {input.mk} {output.mk}"

rule v_sample_01_identify_genomad:
    input:
        metagenome = path + 'results/b-sample/07_contigs/contigs_fna/{sample}.fna',
    output:
        mk = path + 'markers/v-sample/01_identification/{sample}.genomad.mk',
    params:
        out_dir= path + 'results/v-sample/',
        identify_method = "genomad",
        virus_length_limit = virus_length_limit,
        conda_env_python = conda_env_python,
        conda_env_seqkit = conda_env_seqkit,
        conda_env_checkv = conda_env_checkv,
        conda_env_vibrant = conda_env_vibrant,
        conda_env_virsorter2 = conda_env_virsorter2,
        conda_env_genomad = conda_env_genomad,
        checkv_db = checkv_db,
        virsorter2_db = virsorter2_db,
        genomad_db = genomad_db,
        vibrant_db = vibrant_db,
    conda:
        config["conda_envs"]["metalink-main"]
    priority: 1
    threads: config["Params_v_sample"]["01_virus_identify"]["threads"]
    log:
        path + 'logs/v-sample/01_identification/{sample}.genomad.log'
    benchmark:
        path + 'benchmarks/v-sample/01_identification/{sample}.genomad.benchmark'
    shell:
        """
        mkdir -p {params.out_dir}
        python {Dir_script_metalink}/viwrap-scripts/Metalink.01_identify.py \
            --input_length_limit {params.virus_length_limit} \
            --input_metagenome {input.metagenome} \
            --out_dir {params.out_dir} \
            --identify_method {params.identify_method} \
            --threads {threads} \
            --conda_env_python {params.conda_env_python} \
            --conda_env_seqkit {params.conda_env_seqkit} \
            --conda_env_checkv {params.conda_env_checkv} \
            --conda_env_vibrant {params.conda_env_vibrant} \
            --conda_env_virsorter2 {params.conda_env_virsorter2} \
            --conda_env_genomad {params.conda_env_genomad} \
            --checkv_db {params.checkv_db} \
            --virsorter2_db {params.virsorter2_db} \
            --genomad_db {params.genomad_db} \
            --vibrant_db {params.vibrant_db} \
            --root_dir {Dir_script_metalink} > {log} 2>&1
        touch {output.mk}
        """

rule v_sample_01_identify_vs_vb:
    input:
        metagenome = path + 'results/b-sample/07_contigs/contigs_fna/{sample}.fna',
    output:
        mk = path + 'markers/v-sample/01_identification/{sample}.vs_vb.mk',
    params:
        out_dir= path + 'results/v-sample/',
        identify_method = "vs_vb",
        virus_length_limit = virus_length_limit,
        conda_env_python = conda_env_python,
        conda_env_seqkit = conda_env_seqkit,
        conda_env_checkv = conda_env_checkv,
        conda_env_vibrant = conda_env_vibrant,
        conda_env_virsorter2 = conda_env_virsorter2,
        conda_env_genomad = conda_env_genomad,
        checkv_db = checkv_db,
        virsorter2_db = virsorter2_db,
        genomad_db = genomad_db,
        vibrant_db = vibrant_db,
    conda:
        config["conda_envs"]["metalink-main"]
    priority: 1
    threads: config["Params_v_sample"]["01_virus_identify"]["threads"]
    log:
        path + 'logs/v-sample/01_identification/{sample}.vs_vb.log'
    benchmark:
        path + 'benchmarks/v-sample/01_identification/{sample}.vs_vb.benchmark'
    shell:
        """
        mkdir -p {params.out_dir}
        python {Dir_script_metalink}/viwrap-scripts/Metalink.01_identify.py \
            --input_length_limit {params.virus_length_limit} \
            --input_metagenome {input.metagenome} \
            --out_dir {params.out_dir} \
            --identify_method {params.identify_method} \
            --threads {threads} \
            --conda_env_python {params.conda_env_python} \
            --conda_env_seqkit {params.conda_env_seqkit} \
            --conda_env_checkv {params.conda_env_checkv} \
            --conda_env_vibrant {params.conda_env_vibrant} \
            --conda_env_virsorter2 {params.conda_env_virsorter2} \
            --conda_env_genomad {params.conda_env_genomad} \
            --checkv_db {params.checkv_db} \
            --virsorter2_db {params.virsorter2_db} \
            --genomad_db {params.genomad_db} \
            --vibrant_db {params.vibrant_db} \
            --root_dir {Dir_script_metalink} > {log} 2>&1
        touch {output.mk}
        """

rule v_sample_01_combine_genomad_vs_vb:
    input:
        mk_genomad = path + 'markers/v-sample/01_identification/{sample}.genomad.mk',
        mk_vs_vb   = path + 'markers/v-sample/01_identification/{sample}.vs_vb.mk',
        metagenome = path + 'results/b-sample/07_contigs/contigs_fna/{sample}.fna',
    output:
        mk = path + 'markers/v-sample/01_identification/{sample}.genomad_vs_vb.mk'
    params:
        out_dir= path + 'results/v-sample/',
        conda_env_python = conda_env_python,
        conda_env_seqkit = conda_env_seqkit,
    conda:
        config["conda_envs"]["metalink-main"]
    priority: 1
    threads: config["Params_v_sample"]["01_virus_identify"]["threads"]
    log:
        path + 'logs/v-sample/01_identification/{sample}.genomad_vs_vb.log'
    benchmark:
        path + 'benchmarks/v-sample/01_identification/{sample}.genomad_vs_vb.benchmark'
    shell:
        """
        mkdir -p {params.out_dir}
        python {Dir_script_metalink}/viwrap-scripts/Metalink.01_combine_genomad_vs_vb.py \
            --input_metagenome {input.metagenome} \
            --out_dir {params.out_dir} \
            --threads {threads} \
            --conda_env_python {params.conda_env_python} \
            --conda_env_seqkit {params.conda_env_seqkit} \
            --root_dir {Dir_script_metalink} > {log} 2>&1
        touch {output.mk}
        """

rule v_sample_03_index_bwamem:
    input:
        mk = get_v_sample_01_taxonomy_identify_rule_mk,
    output:
        mk = path + 'markers/v-sample/03_index/bwamem/{sample}.mk',
    params:
        source_viral_scaffold = path + 'results/v-sample/01_identification/{sample}/' + identify_method + '_result/final_virus.fna',
        viral_scaffold = path + 'results/v-sample/03_index/bwamem/{sample}.final_virus.fna',
        prefix       = path + 'results/v-sample/03_index/bwamem/{sample}.final_virus',
    conda:
        config["conda_envs"]["metalink-main"]
    priority: 60
    threads: config["Params_v_sample"]["03_mapping"]["threads"]
    log:
        path + 'logs/v-sample/03_index/bwamem/{sample}.03_index.log'
    benchmark:
        path + 'benchmarks/v-sample/03_index/bwamem/{sample}.03_index.benchmark'
    shell:
        """
        mkdir -p $(dirname {params.viral_scaffold}) 2>{log}
        if [ -s {params.source_viral_scaffold} ]; then
            echo "Running v_sample_03_index_bwamem." >>{log}
            ln -sf {params.source_viral_scaffold} {params.viral_scaffold} 2>> {log}
            bwa index -p {params.prefix} {params.viral_scaffold} 2>>{log}
        else
            echo "Empty viral scaffold for {wildcards.sample}; creating empty BAM index." >>{log}
        fi
        touch {output.mk}
        """

rule v_sample_03_mapping_bwamem:
    input:
        r_s = path + 'results/b-sample/01_trim/{sample}_trimmed_s.fastq.gz',
        mk = path + 'markers/v-sample/03_index/bwamem/{sample}.mk',
    output:
        bam = path + 'results/v-sample/03_mapping/bwamem/{sample}/{sample}.bam',
        mk = path + 'markers/v-sample/03_mapping/bwamem/{sample}/{sample}.mk',
    params:
        prefix = path + 'results/v-sample/03_index/bwamem/{sample}.final_virus',
        viral_scaffold = path + 'results/v-sample/03_index/bwamem/{sample}.final_virus.fna',
        sort_mem = config["Params_v_sample"]["03_mapping"].get("samtools_sort_mem", "768M"),
        mapping_extra_args = config["Params_v_sample"]["03_mapping"]["extra_args"],
    conda:
        config["conda_envs"]["metalink-main"]
    priority: 60
    threads: config["Params_v_sample"]["03_mapping"]["threads"]
    log:
        path + 'logs/v-sample/03_mapping/bwamem/{sample}/{sample}.virus_03_mapping.log'
    benchmark:
        path + 'benchmarks/v-sample/03_mapping/bwamem/{sample}/{sample}.virus_03_mapping.benchmark'
    shell:
        """
       if [ -s {params.viral_scaffold} ]; then
            echo "Running v_sample_03_mapping_bwamem." >>{log}
            (bwa mem -t {threads} {params.mapping_extra_args} \
                {params.prefix} \
                {input.r_s} \
                2> {log} \
                | samtools view -bS -h -@ {threads} \
                2>> {log} \
                | samtools sort -@ {threads} -m {params.sort_mem} \
                -o {output.bam}) >> {log} 2>&1
        else
            echo "Empty viral scaffold for {wildcards.sample}; creating empty BAM placeholder." >>{log}
            touch {output.bam}
        fi
        touch {output.mk}
        """

rule v_sample_03_index_bwamem2:
    input:
        mk = get_v_sample_01_taxonomy_identify_rule_mk,
    output:
        mk = path + 'markers/v-sample/03_index/bwamem2/{sample}.mk',
    params:
        source_viral_scaffold = path + 'results/v-sample/01_identification/{sample}/' + identify_method + '_result/final_virus.fna',
        viral_scaffold = path + 'results/v-sample/03_index/bwamem2/{sample}.final_virus.fna',
        prefix       = path + 'results/v-sample/03_index/bwamem2/{sample}.final_virus',
    conda:
        config["conda_envs"]["metalink-main"]
    priority: 60
    threads: config["Params_v_sample"]["03_mapping"]["threads"]
    log:
        path + 'logs/v-sample/03_index/bwamem2/{sample}.03_index.log'
    benchmark:
        path + 'benchmarks/v-sample/03_index/bwamem2/{sample}.03_index.benchmark'
    shell:
        """
        mkdir -p $(dirname {params.viral_scaffold}) 2>{log}
        if [ -s {params.source_viral_scaffold} ]; then
            echo "Running v_sample_03_index_bwamem2." >>{log}
            ln -sf {params.source_viral_scaffold} {params.viral_scaffold} 2>> {log}
            bwa-mem2 index -p {params.prefix} {params.viral_scaffold} >>{log} 2>&1
        else
            echo "Empty viral scaffold for {wildcards.sample}; creating empty BAM index." >>{log}
        fi
        touch {output.mk}
        """

rule v_sample_03_mapping_bwamem2:
    input:
        r_s = path + 'results/b-sample/01_trim/{sample}_trimmed_s.fastq.gz',
        mk = path + 'markers/v-sample/03_index/bwamem2/{sample}.mk',
    output:
        bam = path + 'results/v-sample/03_mapping/bwamem2/{sample}/{sample}.bam',
        mk = path + 'markers/v-sample/03_mapping/bwamem2/{sample}/{sample}.mk',
    params:
        prefix = path + 'results/v-sample/03_index/bwamem2/{sample}.final_virus',
        viral_scaffold = path + 'results/v-sample/03_index/bwamem2/{sample}.final_virus.fna',
        sort_mem = config["Params_v_sample"]["03_mapping"].get("samtools_sort_mem", "768M"),
        mapping_extra_args = config["Params_v_sample"]["03_mapping"]["extra_args"],
    conda:
        config["conda_envs"]["metalink-main"]
    priority: 60
    threads: config["Params_v_sample"]["03_mapping"]["threads"]
    log:
        path + 'logs/v-sample/03_mapping/bwamem2/{sample}/{sample}.virus_03_mapping.log'
    benchmark:
        path + 'benchmarks/v-sample/03_mapping/bwamem2/{sample}/{sample}.virus_03_mapping.benchmark'
    shell:
        """
       if [ -s {params.viral_scaffold} ]; then
            echo "Running v_sample_03_mapping_bwamem2." >>{log}
            (bwa-mem2 mem -t {threads} {params.mapping_extra_args} \
                {params.prefix} \
                {input.r_s} \
                2> {log} \
                | samtools view -bS -h -@ {threads} \
                2>> {log} \
                | samtools sort -@ {threads} -m {params.sort_mem} \
                -o {output.bam}) >> {log} 2>&1
        else
            echo "Empty viral scaffold for {wildcards.sample}; creating empty BAM placeholder." >>{log}
            touch {output.bam}
        fi
        touch {output.mk}
        """

rule v_sample_03_index_minibwa:
    input:
        mk = get_v_sample_01_taxonomy_identify_rule_mk,
    output:
        mk = path + 'markers/v-sample/03_index/minibwa/{sample}.mk',
    params:
        source_viral_scaffold = path + 'results/v-sample/01_identification/{sample}/' + identify_method + '_result/final_virus.fna',
        viral_scaffold = path + 'results/v-sample/03_index/minibwa/{sample}.final_virus.fna',
        prefix       = path + 'results/v-sample/03_index/minibwa/{sample}.final_virus',
    conda:
        config["conda_envs"]["metalink-main"]
    priority: 60
    threads: config["Params_v_sample"]["03_mapping"]["threads"]
    log:
        path + 'logs/v-sample/03_index/minibwa/{sample}.03_index.log'
    benchmark:
        path + 'benchmarks/v-sample/03_index/minibwa/{sample}.03_index.benchmark'
    shell:
        """
        mkdir -p $(dirname {params.viral_scaffold}) 2>{log}
        if [ -s {params.source_viral_scaffold} ]; then
            echo "Running v_sample_03_index_minibwa." >>{log}
            ln -sf {params.source_viral_scaffold} {params.viral_scaffold} 2>> {log}
            minibwa index {params.viral_scaffold} {params.prefix} -t {threads} 2>> {log}
        else
            echo "Empty viral scaffold for {wildcards.sample}; creating empty BAM index." >>{log}
        fi
        touch {output.mk}
        """

rule v_sample_03_mapping_minibwa:
    input:
        r_s = path + 'results/b-sample/01_trim/{sample}_trimmed_s.fastq.gz',
        mk = path + 'markers/v-sample/03_index/minibwa/{sample}.mk',
    output:
        bam = path + 'results/v-sample/03_mapping/minibwa/{sample}/{sample}.bam',
        mk = path + 'markers/v-sample/03_mapping/minibwa/{sample}/{sample}.mk',
    params:
        prefix = path + 'results/v-sample/03_index/minibwa/{sample}.final_virus',
        viral_scaffold = path + 'results/v-sample/03_index/minibwa/{sample}.final_virus.fna',
        sort_mem = config["Params_v_sample"]["03_mapping"].get("samtools_sort_mem", "768M"),
        mapping_extra_args = config["Params_v_sample"]["03_mapping"]["extra_args"],
    conda:
        config["conda_envs"]["metalink-main"]
    priority: 60
    threads: config["Params_v_sample"]["03_mapping"]["threads"]
    log:
        path + 'logs/v-sample/03_mapping/minibwa/{sample}/{sample}.virus_03_mapping.log'
    benchmark:
        path + 'benchmarks/v-sample/03_mapping/minibwa/{sample}/{sample}.virus_03_mapping.benchmark'
    shell:
        """
       if [ -s {params.viral_scaffold} ]; then
            echo "Running v_sample_03_mapping_minibwa." >>{log}
            (minibwa map -t {threads} {params.mapping_extra_args} \
                {params.prefix} \
                {input.r_s} \
                2> {log} \
                | samtools view -bS -h -@ {threads} \
                2>> {log} \
                | samtools sort -@ {threads} -m {params.sort_mem} \
                -o {output.bam}) >> {log} 2>&1
        else
            echo "Empty viral scaffold for {wildcards.sample}; creating empty BAM placeholder." >>{log}
            touch {output.bam}
        fi
        touch {output.mk}
        """

rule v_sample_03_index_strobealign:
    priority: 60
    input:
        mk = get_v_sample_01_taxonomy_identify_rule_mk,
        r_s = expand(path + 'results/b-sample/01_trim/{sample}_trimmed_s.fastq.gz', sample=SAMPLE_INDEX),
    output:
        mk = path + 'markers/v-sample/03_index/strobealign/{sample}.mk',
    params:
        source_viral_scaffold = path + 'results/v-sample/01_identification/{sample}/' + identify_method + '_result/final_virus.fna',
        viral_scaffold = path + 'results/v-sample/03_index/strobealign/{sample}.final_virus.fna',
        len_lst     = path + 'results/v-sample/03_index/strobealign/{sample}.unique_scaled_length.txt',
    conda:
        config["conda_envs"]["metalink-main"]
    threads: config["Params_v_sample"]["03_mapping"]["threads"]
    log:
        path + 'logs/v-sample/03_index/strobealign/{sample}.03_index.log'
    benchmark:
        path + 'benchmarks/v-sample/03_index/strobealign/{sample}.03_index.benchmark'
    shell:
        """
        mkdir -p $(dirname {params.viral_scaffold}) 2>{log}
        if [ -s {params.source_viral_scaffold} ]; then
            echo "Running v_sample_03_index_strobealign." >>{log}
            ln -sf {params.source_viral_scaffold} {params.viral_scaffold} 2>> {log}

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

            for fq in {input.r_s}
            do
                Length=`seqkit head -n 1000 $fq \
                    |seqkit fx2tab -l | cut -f 2 | awk '{{ if (length > max) max = length }} END {{ print max }}'`
                ScaledLength=$(process_number $Length)
                echo $ScaledLength
            done > {params.len_lst} 2>>{log}

            sort -nu {params.len_lst} -o {params.len_lst} 2>>{log}

            for L in $(cat {params.len_lst}); do
                strobealign --create-index -t {threads} -r $L {params.viral_scaffold}
            done >>{log} 2>&1
        else
            echo "Empty viral scaffold for {wildcards.sample}; skip indexing." >>{log}
        fi
        touch {output.mk}
        """

rule v_sample_03_mapping_strobealign:
    input:
        r_s = path + 'results/b-sample/01_trim/{sample}_trimmed_s.fastq.gz',
        mk = path + 'markers/v-sample/03_index/strobealign/{sample}.mk',
    output:
        bam = path + 'results/v-sample/03_mapping/strobealign/{sample}/{sample}.bam',
        mk = path + 'markers/v-sample/03_mapping/strobealign/{sample}/{sample}.mk',
    params:
        viral_scaffold = path + 'results/v-sample/03_index/strobealign/{sample}.final_virus.fna',
        sort_mem = config["Params_v_sample"]["03_mapping"].get("samtools_sort_mem", "768M"),
        mapping_extra_args = config["Params_v_sample"]["03_mapping"]["extra_args"],
    conda:
        config["conda_envs"]["metalink-main"]
    priority: 60
    threads: config["Params_v_sample"]["03_mapping"]["threads"]
    log:
        path + 'logs/v-sample/03_mapping/strobealign/{sample}/{sample}.virus_03_mapping.log'
    benchmark:
        path + 'benchmarks/v-sample/03_mapping/strobealign/{sample}/{sample}.virus_03_mapping.benchmark'
    shell:
        """
        if [ -s {params.viral_scaffold} ]; then
            echo "Running strobealign." >>{log}
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
        Length=`seqkit head -n 1000 {input.r_s} | seqkit fx2tab -l | cut -f 2 | awk '{{ if (length > max) max = length }} END {{ print max }}'`
        ScaledLength=$(process_number $Length)
            (strobealign --use-index -t {threads} -r $ScaledLength {params.mapping_extra_args} {params.viral_scaffold} {input.r_s} \
                2> {log} \
                | samtools sort -@ {threads} -m {params.sort_mem} -o {output.bam}) >> {log} 2>&1
        else
            echo "Empty viral scaffold for {wildcards.sample}; creating empty BAM placeholder." >>{log}
            touch {output.bam}
        fi
        touch {output.mk}
        """

rule v_sample_03_taxonomy:
    input:
        mk = get_v_sample_01_taxonomy_identify_rule_mk,
    output:
        mk = path + 'markers/v-sample/03_taxonomy/{sample}.mk',
    params:
        sample = "{sample}",
        viral_scaffold = path + 'results/v-sample/01_identification/{sample}/' + identify_method + '_result/final_virus.fna',
        viral_faa = path + 'results/v-sample/01_identification/{sample}/' + identify_method + '_result/final_virus.faa',
        CAT_range = config["Params_v_sample"]["03_taxonomy_CAT_pack"]["range"],
        CAT_fraction = config["Params_v_sample"]["03_taxonomy_CAT_pack"]["fraction"],
        db = CAT_taxonomy_virus_db,
        tax = CAT_taxonomy_virus_tax,
        out_dir= path + 'results/v-sample/03_taxonomy/{sample}/',
        extra_args = config["Params_v_sample"]["03_taxonomy_CAT_pack"]["extra_args"],
    conda:
        config["conda_envs"]["metalink-annotation"]
    priority: 10
    threads: config["Params_v_sample"]["03_taxonomy_CAT_pack"]["threads"]
    log:
        path + 'logs/v-sample/03_taxonomy/{sample}.virus_03_taxonomy.log'
    benchmark:
        path + 'benchmarks/v-sample/03_taxonomy/{sample}.virus_03_taxonomy.benchmark'
    shell:
        """
        mkdir -p {params.out_dir} 2>{log}
        if [ -s {params.viral_scaffold} ]; then
            echo "Running CAT_pack pipeline." >> {log} 2>&1
            CAT_pack contigs -n {threads} -r {params.CAT_range} -f {params.CAT_fraction} -c {params.viral_scaffold} -p {params.viral_faa} -d {params.db} -t {params.tax} -o {params.out_dir}/CAT_pack --no_stars --force {params.extra_args} >> {log} 2>&1
            CAT_pack add_names -i {params.out_dir}/CAT_pack.ORF2LCA.txt -o {params.out_dir}/CAT_pack.ORF2LCA.name.txt -t {params.tax} --only_official --force >> {log} 2>&1
            CAT_pack add_names -i {params.out_dir}/CAT_pack.contig2classification.txt -o {params.out_dir}/CAT_pack.contig2classification.name.txt -t {params.tax} --only_official --force >> {log} 2>&1
            CAT_pack summarise -c {params.viral_scaffold} -i {params.out_dir}/CAT_pack.contig2classification.name.txt -o {params.out_dir}/CAT_pack.summary.tsv --force >> {log} 2>&1
            mv {params.out_dir}/CAT_pack.ORF2LCA.name.txt {params.out_dir}/CAT_pack.ORF2LCA.tsv >> {log} 2>&1
            mv {params.out_dir}/CAT_pack.contig2classification.name.txt {params.out_dir}/CAT_pack.contig2classification.tsv >> {log} 2>&1
            rm -f {params.out_dir}/CAT_pack.ORF2LCA.txt >> {log} 2>&1
            rm -f {params.out_dir}/CAT_pack.contig2classification.txt >> {log} 2>&1
        else
            echo "Empty viral scaffold for {wildcards.sample}; skip CAT_pack taxonomy." >> {log} 2>&1
            touch {params.out_dir}/CAT_pack.ORF2LCA.tsv \
                {params.out_dir}/CAT_pack.contig2classification.tsv \
                {params.out_dir}/CAT_pack.summary.tsv >> {log} 2>&1
        fi
        touch {output.mk}
        """

rule v_sample_04_binning:
    input:
        metagenome = path + 'results/b-sample/07_contigs/contigs_fna/{sample}.fna',
        bam = path + 'results/v-sample/03_mapping/' + v_sample_mapping_method + '/{sample}/{sample}.bam',
    output:
        mk = path + 'markers/v-sample/04_vRhyme/{sample}.mk',
    params:
        viral_scaffold = path + 'results/v-sample/01_identification/{sample}/' + identify_method + '_result/final_virus.fna',
        scf2lytic_or_lyso_summary = path + 'results/v-sample/01_identification/{sample}/' + identify_method + '_result/scf2lytic_or_lyso.summary.txt',
        outdir_vrhyme = path + 'results/v-sample/04_vRhyme/{sample}/',
        identify_method = identify_method,
        conda_env_vrhyme = conda_env_vrhyme,
        conda_env_checkv = conda_env_checkv,
        checkv_db = checkv_db,
    conda:
        config["conda_envs"]["metalink-main"]
    priority: 1
    threads: config["Params_v_sample"]["04_binning"]["threads"]
    log:
        path + 'logs/v-sample/04_vRhyme/{sample}.virus_04_binning.log'
    benchmark:
        path + 'benchmarks/v-sample/04_vRhyme/{sample}.virus_04_binning.benchmark'
    shell:
        """
        python {Dir_script_metalink}/viwrap-scripts/Metalink.04_vRhyme.py \
            --viral_scaffold {params.viral_scaffold} \
            --scf2lytic_or_lyso_summary {params.scf2lytic_or_lyso_summary} \
            --input_metagenome {input.metagenome} \
            --input_bam {input.bam} \
            --outdir_vrhyme {params.outdir_vrhyme} \
            --threads {threads} \
            --conda_env_vrhyme {params.conda_env_vrhyme} \
            --conda_env_checkv {params.conda_env_checkv} \
            --checkv_db {params.checkv_db} \
            --root_dir {Dir_script_metalink} > {log} 2>&1
        touch {output.mk}
        """

rule v_sample_05_CheckV:
    input:
        mk = path + 'markers/v-sample/04_vRhyme/{sample}.mk',
    output:
        mk = path + 'markers/v-sample/05_CheckV/{sample}.mk',
    params:
        indir_vrhyme = path + 'results/v-sample/04_vRhyme/{sample}/',
        outdir_checkv = path + 'results/v-sample/05_CheckV/{sample}/',
        conda_env_checkv = conda_env_checkv,
        checkv_db = checkv_db,
    conda:
        config["conda_envs"]["metalink-main"]
    priority: 1
    threads: config["Params_v_sample"]["05_checkv"]["threads"]
    log:
        path + 'logs/v-sample/05_CheckV/{sample}.virus_05_CheckV.log'
    benchmark:
        path + 'benchmarks/v-sample/05_CheckV/{sample}.virus_05_CheckV.benchmark'
    shell:
        """
        python {Dir_script_metalink}/viwrap-scripts/Metalink.05_CheckV.py \
            --indir_vrhyme {params.indir_vrhyme} \
            --outdir_checkv {params.outdir_checkv} \
            --threads {threads} \
            --conda_env_checkv {params.conda_env_checkv} \
            --checkv_db {params.checkv_db} \
            --root_dir {Dir_script_metalink} > {log} 2>&1
        touch {output.mk}
        """

rule v_sample_06_vclust:
    input:
        mk = path + 'markers/v-sample/05_CheckV/{sample}.mk',
    output:
        fltr_txt = path + 'results/v-sample/06_vclust/{sample}/fltr.txt',
        ani_tsv = path + 'results/v-sample/06_vclust/{sample}/ani.tsv',
        ani_ids_tsv = path + 'results/v-sample/06_vclust/{sample}/ani.ids.tsv',
        ani_ids_scale_tsv = path + 'results/v-sample/06_vclust/{sample}/ani.ids.scale.tsv',
        clusters_tsv = path + 'results/v-sample/06_vclust/{sample}/clusters.tsv',
        unique_clusters = path + 'results/v-sample/06_vclust/{sample}/unique_clusters.tsv',
    params:
        Nlinked_viral_gn = path + 'results/v-sample/05_CheckV/{sample}/all_vRhyme_fasta.Nlinked_viral_gn.fna',
        ani = config["Params_v_sample"]["06_vclust"]["ani"],
        bin_folder = path + 'results/v-sample/04_vRhyme/{sample}/vRhyme_both_binned_unbinned/',
        vMAGs_folder = path + 'results/v-sample/06_vclust/{sample}/vMAGs_clusters_fna/',
        prefilter_extra_args = config["Params_v_sample"]["06_vclust"].get("prefilter_extra_args", ""),
        align_extra_args = config["Params_v_sample"]["06_vclust"].get("align_extra_args", ""),
        cluster_extra_args = config["Params_v_sample"]["06_vclust"].get("cluster_extra_args", ""),
    conda:
        config["conda_envs"]["metalink-vclust"]
    priority: 1
    threads: config["Params_v_sample"]["06_vclust"]["threads"]
    log:
        path + 'logs/v-sample/06_vclust/{sample}.06_vclust.log'
    benchmark:
        path + 'benchmarks/v-sample/06_vclust/{sample}.06_vclust.benchmark'
    shell:
        """
        if [ -s {params.Nlinked_viral_gn} ]; then
            n_seq=$(grep -c '^>' {params.Nlinked_viral_gn})
            mkdir -p {params.vMAGs_folder} 2>{log}
            if [ "$n_seq" -ge 2 ]; then
                echo "Running vclust pipeline." >> {log} 2>&1
                # Prefilter similar genome sequence pairs before conducting pairwise alignments.
                vclust prefilter -t {threads} -i {params.Nlinked_viral_gn} -o {output.fltr_txt} {params.prefilter_extra_args} > {log} 2>&1
                # Align similar genome sequence pairs and calculate pairwise ANI measures.
                vclust align -t {threads} -i {params.Nlinked_viral_gn} -o {output.ani_tsv} --filter {output.fltr_txt} {params.align_extra_args} >> {log} 2>&1

                python3 -c "
import csv,sys
in_f=sys.argv[1];out_f=sys.argv[2]
with open(in_f,newline='')as fin,open(out_f,'w',newline='')as fout:
    r=csv.DictReader(fin,delimiter='\t')
    w=csv.DictWriter(fout,fieldnames=r.fieldnames,delimiter='\t')
    w.writeheader()
    for row in r:
        if '_bin_' in row['id']:row['seq_len']=str(int(row['seq_len'])-1000)
        w.writerow(row)
" {output.ani_ids_tsv} {output.ani_ids_scale_tsv} >> {log} 2>&1

                # Cluster genome sequences based on given ANI measure and minimum threshold.
                vclust cluster -i {output.ani_tsv} -o {output.clusters_tsv} --ids {output.ani_ids_scale_tsv} --metric ani --ani {params.ani} -r {params.cluster_extra_args} >> {log} 2>&1

                awk -F'\\t' 'NR>1 && !seen[$2]++ {{print $2}}' {output.clusters_tsv} \
                    > {output.unique_clusters} 2>> {log}

                mkdir -p {params.vMAGs_folder} 2>> {log}
                awk '{{print $0".fna"}}' {output.unique_clusters} \
                    | xargs -r -I{{}} find {params.bin_folder} -type f -name {{}} -print0 \
                    | xargs -0 -r cp -t {params.vMAGs_folder} 2>> {log}
            else
                echo "Only one viral scaffold for {wildcards.sample}; skip vclust."
                vclust prefilter -t {threads} -i {params.Nlinked_viral_gn} -o {output.fltr_txt} {params.prefilter_extra_args} > {log} 2>&1
                seq_name=$(grep '^>' {params.Nlinked_viral_gn} | sed 's/^>//' | awk '{{print $1}}')
                printf "object\\tcluster\\n%s\\t%s\\n" "$seq_name" "$seq_name" > {output.clusters_tsv}
                printf "$seq_name" > {output.unique_clusters}
                touch {output.ani_tsv} {output.ani_ids_tsv} {output.ani_ids_scale_tsv}
            fi
        else
            echo "Empty viral scaffold for {wildcards.sample}; skip vclust." > {log} 2>&1
            touch {output.fltr_txt} {output.ani_tsv} {output.ani_ids_tsv} \
                {output.ani_ids_scale_tsv} {output.clusters_tsv} {output.unique_clusters}
        fi
        """

rule v_sample_07_iphop_none:
    input:
        mk = path + 'markers/v-sample/05_CheckV/{sample}.mk',
    output:
        mk = path + 'markers/v-sample/07_iPHoP/{sample}.mk',
    params:
        Nlinked_viral_gn = path + 'results/v-sample/05_CheckV/{sample}/all_vRhyme_fasta.Nlinked_viral_gn.fna',
        outdir_iphop= path + 'results/v-sample/07_iPHoP/{sample}/',
        conda_env_iphop = conda_env_iphop,
        conda_env_genomad = conda_env_genomad,
        conda_env_gtdbtk = conda_env_gtdbtk,
        gtdbtk_db = gtdbtk_db,
        iphop_db = iphop_db,
        genomad_db = genomad_db,
        custom_MAGs_dir = 'none',
        iPHoP_db_custom = 'none',
    conda:
        config["conda_envs"]["metalink-main"]
    priority: 1
    threads: config["Params_v_sample"]["07_iPHoP"]["threads"]
    log:
        path + 'logs/v-sample/07_iPHoP/{sample}.virus_07_iPHoP.log'
    benchmark:
        path + 'benchmarks/v-sample/07_iPHoP/{sample}.virus_07_iPHoP.benchmark'
    shell:
        """
        python {Dir_script_metalink}/viwrap-scripts/Metalink.07_iPHoP.py \
            --input_Nlinked_viral_gn {params.Nlinked_viral_gn} \
            --outdir_iphop {params.outdir_iphop} \
            --custom_MAGs_dir {params.custom_MAGs_dir} \
            --iPHoP_db_custom {params.iPHoP_db_custom} \
            --threads {threads} \
            --conda_env_iphop {params.conda_env_iphop} \
            --conda_env_genomad {params.conda_env_genomad} \
            --conda_env_gtdbtk {params.conda_env_gtdbtk} \
            --gtdbtk_db {params.gtdbtk_db} \
            --iphop_db {params.iphop_db} \
            --genomad_db {params.genomad_db} \
            --root_dir {Dir_script_metalink} > {log} 2>&1
        touch {output.mk}
        """

rule v_sample_08_taxonomy:
    input:
        mk = path + 'markers/v-sample/04_vRhyme/{sample}.mk',
    output:
        mk = path + 'markers/v-sample/08_taxonomy/{sample}.mk',
    params:
        sample = "{sample}",
        out_dir= path + 'results/v-sample/08_taxonomy/{sample}/',
        bin_folder = path + 'results/v-sample/04_vRhyme/{sample}/vRhyme_both_binned_unbinned/',
        bin_faa = path + 'results/v-sample/04_vRhyme/{sample}/all_vRhyme_faa.faa',
        BAT_range = config["Params_v_sample"]["08_taxonomy_BAT"]["range"],
        BAT_fraction = config["Params_v_sample"]["08_taxonomy_BAT"]["fraction"],
        db = CAT_taxonomy_virus_db,
        tax = CAT_taxonomy_virus_tax,
        extra_args = config["Params_v_sample"]["08_taxonomy_BAT"]["extra_args"],
    conda:
        config["conda_envs"]["metalink-annotation"]
    priority: 10
    threads: config["Params_v_sample"]["08_taxonomy_BAT"]["threads"]
    log:
        path + 'logs/v-sample/08_taxonomy/{sample}.virus_08_taxonomy.log'
    benchmark:
        path + 'benchmarks/v-sample/08_taxonomy/{sample}.virus_08_taxonomy.benchmark'
    shell:
        """
        mkdir -p {params.out_dir} 2>{log}
        if [ -s {params.bin_faa} ]; then
            echo "Running CAT_pack pipeline." >> {log} 2>&1
            CAT_pack bins -n {threads} -r {params.BAT_range} -f {params.BAT_fraction} -b {params.bin_folder} -p {params.bin_faa} -d {params.db} -t {params.tax} -o {params.out_dir}/BAT --no_stars --force {params.extra_args} >> {log} 2>&1
            CAT_pack add_names -i {params.out_dir}/BAT.ORF2LCA.txt -o {params.out_dir}/BAT.ORF2LCA.name.txt -t {params.tax} --only_official --force >> {log} 2>&1
            CAT_pack add_names -i {params.out_dir}/BAT.bin2classification.txt -o {params.out_dir}/BAT.bin2classification.name.txt -t {params.tax} --only_official --force >> {log} 2>&1
            CAT_pack summarise -i {params.out_dir}/BAT.bin2classification.name.txt -o {params.out_dir}/BAT.summary.txt --force >> {log} 2>&1
            mv {params.out_dir}/BAT.ORF2LCA.name.txt {params.out_dir}/BAT.ORF2LCA.tsv >> {log} 2>&1
            mv {params.out_dir}/BAT.bin2classification.name.txt {params.out_dir}/BAT.bin2classification.tsv >> {log} 2>&1
            rm -f {params.out_dir}/BAT.ORF2LCA.txt >> {log} 2>&1
            rm -f {params.out_dir}/BAT.bin2classification.txt >> {log} 2>&1
        else
            echo "Empty viral scaffold for {wildcards.sample}; skip CAT_pack taxonomy." >> {log} 2>&1
            touch {params.out_dir}/BAT.ORF2LCA.tsv \
                {params.out_dir}/BAT.bin2classification.tsv \
                {params.out_dir}/BAT.summary.tsv >> {log} 2>&1
        fi
        touch {output.mk}
        """

rule v_sample_09_coverm:
    input:
        bam = path + 'results/v-sample/03_mapping/' + v_sample_mapping_method + '/{sample}/{sample}.bam',
        mk = path + 'markers/v-sample/04_vRhyme/{sample}.mk',
    output:
        out_coverm = path + 'results/v-sample/09_coverm/{sample}/bin_coverage.tsv',
    threads: config["Params_v_sample"]["09_coverm_eachsample_vMAGs"]["threads"]
    conda:
        config["conda_envs"]["metalink-main"]
    priority: 80
    params:
        MAGs_eachsample = path + 'results/v-sample/04_vRhyme/{sample}/vRhyme_both_binned_unbinned_rawname/',
        genometype = "fna",
        method = config["Params_v_sample"]["09_coverm_eachsample_vMAGs"]["method"],
        min_read_percent_identity = config["Params_v_sample"]["09_coverm_eachsample_vMAGs"]["min_read_percent_identity"],
        min_read_aligned_length = config["Params_v_sample"]["09_coverm_eachsample_vMAGs"]["min_read_aligned_length"],
        contig_end_exclusion = config["Params_v_sample"]["09_coverm_eachsample_vMAGs"]["contig_end_exclusion"],
        extra_args = config["Params_v_sample"]["09_coverm_eachsample_vMAGs"]["extra_args"],
    log:
        path + 'logs/v-sample/09_coverm/{sample}.09_coverm.log'
    benchmark:
        path + 'benchmarks/v-sample/09_coverm/{sample}.09_coverm.benchmark'
    shell:
        """
        if [ -d "{params.MAGs_eachsample}" ]; then
            MAGs_eachsample={params.MAGs_eachsample}
        else
            echo "The directory vMAGs_folder is not exist."
        fi > {log} 2>&1

        if [ -z "$(ls -A $MAGs_eachsample)" ]; then
            echo "The directory vMAGs_folder is empty, pass."
            touch {output.out_coverm}
        else
            echo "Calculating abundance..."
            coverm genome \
                -t {threads} \
                -m {params.method} \
                --min-read-percent-identity {params.min_read_percent_identity} \
                --min-read-aligned-length {params.min_read_aligned_length} \
                --contig-end-exclusion {params.contig_end_exclusion} \
                --bam-files {input.bam} \
                --genome-fasta-directory $MAGs_eachsample \
                -x {params.genometype} \
                -o {output.out_coverm} {params.extra_args}
        fi >> {log} 2>&1
        """

rule v_sample_10_summary:
    input:
        metagenome = path + 'results/b-sample/07_contigs/contigs_fna/{sample}.fna',
        clusters_tsv = path + 'results/v-sample/06_vclust/{sample}/clusters.tsv',
        mk_07_iPHoP = path + 'markers/v-sample/07_iPHoP/{sample}.mk',
        mk_08_taxonomy = path + 'markers/v-sample/08_taxonomy/{sample}.mk',
        out_coverm = path + 'results/v-sample/09_coverm/{sample}/bin_coverage.tsv',
    output:
        mk = path + 'markers/v-sample/10_summary/{sample}.mk',
    params:
        out_dir= path + 'results/v-sample/',
        identify_method = identify_method,
        conda_env_seqkit = conda_env_seqkit,
        vibrant_db = vibrant_db,
        custom_MAGs_dir = 'none',
        iPHoP_db_custom = 'none',
    conda:
        config["conda_envs"]["metalink-main"]
    priority: 1
    threads: 8
    log:
        path + 'logs/v-sample/10_summary/{sample}.virus_10_summary.log'
    benchmark:
        path + 'benchmarks/v-sample/10_summary/{sample}.virus_10_summary.benchmark'
    shell:
        """
        python {Dir_script_metalink}/viwrap-scripts/Metalink.10_summary.py \
            --out_dir {params.out_dir} \
            --input_metagenome {input.metagenome} \
            --custom_MAGs_dir {params.custom_MAGs_dir} \
            --iPHoP_db_custom {params.iPHoP_db_custom} \
            --identify_method {params.identify_method} \
            --threads {threads} \
            --conda_env_seqkit {params.conda_env_seqkit} \
            --vibrant_db {params.vibrant_db} \
            --root_dir {Dir_script_metalink} >> {log} 2>&1
        touch {output.mk}
        """

rule v_sample_11_visualize:
    input:
        metagenome = path + 'results/b-sample/07_contigs/contigs_fna/{sample}.fna',
        mk = path + 'markers/v-sample/10_summary/{sample}.mk',
    output:
        mk = path + 'markers/v-sample/11_visualization/{sample}.mk',
    params:
        indir_summary = path + 'results/v-sample/10_summary/{sample}/',
        outdir_visualize = path + 'results/v-sample/11_visualization/{sample}/',
        vibrant_db = vibrant_db,
    conda:
        config["conda_envs"]["metalink-main"]
    priority: 1
    threads: 8
    log:
        path + 'logs/v-sample/11_visualization/{sample}.virus_11_visualize.log'
    benchmark:
        path + 'benchmarks/v-sample/11_visualization/{sample}.virus_11_visualize.benchmark'
    shell:
        """
        python {Dir_script_metalink}/viwrap-scripts/Metalink.11_visualize.py \
            --indir_summary {params.indir_summary} \
            --outdir_visualize {params.outdir_visualize} \
            --threads {threads} \
            --vibrant_db {params.vibrant_db} \
            --root_dir {Dir_script_metalink} > {log} 2>&1
        touch {output.mk}
        """

rule v_sample_12_delete_bam:
    input:
        out_coverm = path + 'results/v-sample/09_coverm/{sample}/bin_coverage.tsv',
    output:
        mk = path + 'markers/v-sample/12_delete_bam/{sample}.mk'
    threads: 1
    priority: 100
    params:
        bam = path + 'results/v-sample/03_mapping/' + v_sample_mapping_method + '/{sample}/{sample}.bam',
    log:
        path + 'logs/v-sample/12_delete_bam/{sample}.virus_12_delete_bam.log'
    benchmark:
        path + 'benchmarks/v-sample/12_delete_bam/{sample}.virus_12_delete_bam.benchmark'
    shell:
        """
        if [ -e "{params.bam}" ];then rm -r {params.bam}*;else echo Folder trimmed.bam not exist.;fi > {log} 2>&1
        touch {output.mk}
        """

# VirSorter2 writes its multi-GB iter-0 working directory inside the vs_vb
# result folder (used by the vs_vb and genomad_vs_vb identify methods); only
# the small final_* files above it are results. Delete iter-0 once
# identification is done; both the current nested and the flat layout are
# covered by existence checks.
rule v_sample_12_delete_vs2_iter0:
    input:
        mk = path + 'markers/v-sample/01_identification/{sample}.identify_done.mk',
    output:
        mk = path + 'markers/v-sample/12_delete_vs2_iter0/{sample}.mk'
    threads: 1
    priority: 100
    params:
        iter0_nested = path + 'results/v-sample/01_identification/{sample}/vs_vb_result/vs_vb_result/iter-0',
        iter0_flat = path + 'results/v-sample/01_identification/{sample}/vs_vb_result/iter-0',
    log:
        path + 'logs/v-sample/12_delete_vs2_iter0/{sample}.virus_12_delete_vs2_iter0.log'
    benchmark:
        path + 'benchmarks/v-sample/12_delete_vs2_iter0/{sample}.virus_12_delete_vs2_iter0.benchmark'
    shell:
        """
        if [ -e "{params.iter0_nested}" ];then rm -r {params.iter0_nested};else echo Folder vs2 iter-0 nested not exist.;fi > {log} 2>&1
        if [ -e "{params.iter0_flat}" ];then rm -r {params.iter0_flat};else echo Folder vs2 iter-0 flat not exist.;fi >> {log} 2>&1
        touch {output.mk}
        """
