'''
/*
 * @Date: 2022-12-04 17:15:47
 * @LastEditors: LiuyangLi
 * @LastEditTime: 2024-07-20 13:03:29
 * @FilePath: pmr_m.smk
 * @Description: project-mMAG-read chain: map sample reads to the project-level dereplicated mMAG reference, then integrate taxonomy/function/depth
 */
'''

import os
import sys

Dir_script_metalink = config["Directory"]["script"]["metalink"]
kegg_modules_ko_info_path = config["Directory"]["script"]["kegg_modules_ko_info_path"]
taxid_fun_path = config["Directory"]["database"]["taxid_fun_path"]
MAGs_taxonomy_classify_method = config["Params_m_sample"]["02_taxonomy"]["MAGs_taxonomy_classify_method"]
MAGs_function_classify_method = config["Params_m_sample"]["08_function"]["MAGs_function_classify_method"]

folder_checks = [
                path + 'results/',
                path + 'results/pmr_m-sample/',
                path + 'results/pmr_m-project/',
                ]
for item in folder_checks:
    if not os.path.isdir(item): os.mkdir(item)

pmr_m_sample_mapping_method = config["Params_pmr_m_sample"]["03_mapping"]["method"]
#pmr_m_sample_mapping_method=bwamem,bwamem2,strobealign,minibwa
def _pmr_m_get_mk(wildcards, rule_type):
    if pmr_m_sample_mapping_method == "strobealign":
        if rule_type == "index":
            return rules.pmr_m_project_03_index_strobealign.output.mk
        else:
            return rules.pmr_m_sample_03_mapping_strobealign.output.mk
    elif pmr_m_sample_mapping_method == "bwamem":
        if rule_type == "index":
            return rules.pmr_m_project_03_index_bwamem.output.mk
        else:
            return rules.pmr_m_sample_03_mapping_bwamem.output.mk
    elif pmr_m_sample_mapping_method == "bwamem2":
        if rule_type == "index":
            return rules.pmr_m_project_03_index_bwamem2.output.mk
        else:
            return rules.pmr_m_sample_03_mapping_bwamem2.output.mk
    elif pmr_m_sample_mapping_method == "minibwa":
        if rule_type == "index":
            return rules.pmr_m_project_03_index_minibwa.output.mk
        else:
            return rules.pmr_m_sample_03_mapping_minibwa.output.mk
    else:
        raise ValueError(f"Unknown mapping method: {pmr_m_sample_mapping_method}")

pmr_m_get_03_mapping_rule_mk = lambda wc: _pmr_m_get_mk(wc, "mapping")



def pmr_m_gated_input(wildcards, artifact):
    # The gated genome artifacts are content-managed side effects of the quality
    # gate. They are deliberately NOT declared outputs: Snakemake wipes declared
    # outputs before every rerun (prepare() -> remove_existing_output()), which
    # would defeat the gate's write-if-changed logic. Depend on the artifact once
    # it exists; on a fresh DAG fall back to the gate marker for ordering.
    f = path + 'results/pmr_m-project/02_quality_gate/{project}.{artifact}'.format(project=wildcards.project, artifact=artifact)
    if os.path.exists(f):
        return f
    return path + 'markers/pmr_m-project/02_quality_gate/{project}.mk'.format(project=wildcards.project)

rule pmr_m_project_02_quality_gate:
    input:
        m_tax = path + 'results/m-project/04_summary/{project}/mMAG_summary_info.tsv',
    output:
        m_filter = path + 'results/pmr_m-project/02_quality_gate/{project}.mMAG_summary_info.filter.tsv',
        mk = path + 'markers/pmr_m-project/02_quality_gate/{project}.mk',
    params:
        MAGs_taxonomy_classify_method = config["Params_m_sample"]["02_taxonomy"]["MAGs_taxonomy_classify_method"],
        MAGs_min_completeness1 = config["Params_pmr_m_project"]["02_quality_gate"]["MAGs_min_completeness1"],
        MAGs_max_contamination1 = config["Params_pmr_m_project"]["02_quality_gate"]["MAGs_max_contamination1"],
        MAGs_min_completeness2 = config["Params_pmr_m_project"]["02_quality_gate"]["MAGs_min_completeness2"],
        MAGs_max_contamination2 = config["Params_pmr_m_project"]["02_quality_gate"]["MAGs_max_contamination2"],
        drep_dir = path + 'results/m-project/01_drep_project_95/{project}/dereplicated_genomes/',
        gated_dir = path + 'results/pmr_m-project/02_quality_gate/{project}.mMAGs_gated',
        genome_list = path + 'results/pmr_m-project/02_quality_gate/{project}.mMAGs_gated.list',
    conda: config["conda_envs"]["metalink-main"]
    priority: 1
    threads: 1
    log:
        path + 'logs/pmr_m-project/02_quality_gate/{project}.02_quality_gate.log'
    benchmark:
        path + 'benchmarks/pmr_m-project/02_quality_gate/{project}.02_quality_gate.benchmark'
    shell:
        """
        mkdir -p $(dirname {log}) 2>/dev/null
        Rscript {Dir_script_metalink}/MAGs_m_quality_filter.cmd.R \
            --MAGs_taxonomy_classify_method {params.MAGs_taxonomy_classify_method} \
            --in_MAGs_tax_detail {input.m_tax} \
            --MAGs_min_completeness1 {params.MAGs_min_completeness1} \
            --MAGs_max_contamination1 {params.MAGs_max_contamination1} \
            --MAGs_min_completeness2 {params.MAGs_min_completeness2} \
            --MAGs_max_contamination2 {params.MAGs_max_contamination2} \
            --out_MAGs_tax_detail {output.m_filter} > {log} 2>&1

        # the gated genome set is sequence-level: rewrite it (and the gated
        # symlink dir) only when the set actually changed, so taxonomy-only
        # reruns do not retrigger the reference/index/mapping chain
        awk -F'\t' 'NR>1{{print $1}}' {output.m_filter} | sort -u > {params.genome_list}.tmp 2>> {log}
        if [ ! -d {params.gated_dir} ]; then _gate_reason="gated_dir missing";
        elif [ ! -f {params.genome_list} ]; then _gate_reason="genome_list missing";
        elif ! cmp -s {params.genome_list}.tmp {params.genome_list}; then _gate_reason="genome set changed";
        else _gate_reason=""; fi
        if [ -n "$_gate_reason" ]; then
            echo "[INFO] rewrite gated genomes: ${{_gate_reason}}" >> {log}
            mv {params.genome_list}.tmp {params.genome_list}
            mkdir -p {params.gated_dir} 2>> {log}
            rm -f {params.gated_dir}/*.fna 2>> {log}
            while read g; do
                for f in {params.drep_dir}/"$g".fa {params.drep_dir}/"$g".fna; do
                    [ -e "$f" ] || continue
                    ln -sf "$f" {params.gated_dir}/"$(basename "$f")"
                done
            done < {params.genome_list}
        else
            echo "[INFO] gated genome set unchanged; keep existing {params.genome_list} and {params.gated_dir}" >> {log}
            rm -f {params.genome_list}.tmp
        fi
        n_mMAG=$(wc -l < {params.genome_list})
        echo "[INFO] project {wildcards.project}: $n_mMAG mMAG(s) passed the quality gate." | tee -a {log} >&2
        if [ "$n_mMAG" -eq 0 ]; then
            echo "[WARNING] project {wildcards.project}: no mMAGs passed the quality gate; downstream reference/mapping will be empty." | tee -a {log} >&2
        fi
        touch {output.mk}
        """

rule pmr_m_project_03_index_bwamem:
    input:
        gated_list = lambda wildcards: pmr_m_gated_input(wildcards, 'mMAGs_gated.list'),
    output:
        mk = path + 'markers/pmr_m-project/03_index/bwamem/{project}.mk',
    params:
        MAGs_folder = path + 'results/pmr_m-project/02_quality_gate/{project}.mMAGs_gated/',
        MAGs_dRep   = path + 'results/pmr_m-project/03_index/bwamem/{project}.MAGs_dRep.fna',
        prefix      = path + 'results/pmr_m-project/03_index/bwamem/{project}.MAGs_dRep',
    conda: config["conda_envs"]["metalink-main"]
    priority: 1
    threads: config["Params_pmr_m_sample"]["03_mapping"]["threads"]
    log:
        path + 'logs/pmr_m-project/03_index/bwamem/{project}.03_index.log'
    benchmark:
        path + 'benchmarks/pmr_m-project/03_index/bwamem/{project}.03_index.benchmark'
    shell:
        """
        mkdir -p $(dirname {params.MAGs_dRep}) 2> {log}
        find -L {params.MAGs_folder} -type f -print0 | xargs -0 -r cat > {params.MAGs_dRep} 2>> {log}
        if [ ! -s {params.MAGs_dRep} ]; then
            echo "ERROR: project {wildcards.project}: no mMAGs passed the quality gate; reference {params.MAGs_dRep} is empty." >> {log}
            echo "ERROR [pmr_m_project_03_index] project {wildcards.project}: no mMAGs passed the quality gate; cannot build an empty reference. See log: {log}" >&2
            exit 1
        fi
        bwa index -p {params.prefix} {params.MAGs_dRep} 2>> {log}
        touch {output.mk}
        """

rule pmr_m_sample_03_mapping_bwamem:
    input:
        r1 = path + 'results/b-sample/01_trim/{sample}_trimmed_1.fastq.gz',
        r2 = path + 'results/b-sample/01_trim/{sample}_trimmed_2.fastq.gz',
        mk = path + 'markers/pmr_m-project/03_index/bwamem/{project}.mk',
    output:
        bam = path + 'results/pmr_m-sample/03_mapping/bwamem/{project}/{sample}.bam',
        mk = path + 'markers/pmr_m-sample/03_mapping/bwamem/{project}/{sample}.mk',
    params:
        prefix = path + 'results/pmr_m-project/03_index/bwamem/{project}.MAGs_dRep',
        sort_mem = config["Params_pmr_m_sample"]["03_mapping"].get("samtools_sort_mem", "2G"),
        mapping_extra_args = config["Params_pmr_m_sample"]["03_mapping"]["extra_args"],
    conda: config["conda_envs"]["metalink-main"]
    priority: 60
    threads: config["Params_pmr_m_sample"]["03_mapping"]["threads"]
    log:
        path + 'logs/pmr_m-sample/03_mapping/bwamem/{project}/{sample}.03_mapping.log'
    benchmark:
        path + 'benchmarks/pmr_m-sample/03_mapping/bwamem/{project}/{sample}.03_mapping.benchmark'
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

rule pmr_m_project_03_index_bwamem2:
    input:
        gated_list = lambda wildcards: pmr_m_gated_input(wildcards, 'mMAGs_gated.list'),
    output:
        mk = path + 'markers/pmr_m-project/03_index/bwamem2/{project}.mk',
    params:
        MAGs_folder = path + 'results/pmr_m-project/02_quality_gate/{project}.mMAGs_gated/',
        MAGs_dRep   = path + 'results/pmr_m-project/03_index/bwamem2/{project}.MAGs_dRep.fna',
        prefix      = path + 'results/pmr_m-project/03_index/bwamem2/{project}.MAGs_dRep',
    conda: config["conda_envs"]["metalink-main"]
    priority: 1
    threads: config["Params_pmr_m_sample"]["03_mapping"]["threads"]
    log:
        path + 'logs/pmr_m-project/03_index/bwamem2/{project}.03_index.log'
    benchmark:
        path + 'benchmarks/pmr_m-project/03_index/bwamem2/{project}.03_index.benchmark'
    shell:
        """
        mkdir -p $(dirname {params.MAGs_dRep}) 2> {log}
        find -L {params.MAGs_folder} -type f -print0 | xargs -0 -r cat > {params.MAGs_dRep} 2>> {log}
        if [ ! -s {params.MAGs_dRep} ]; then
            echo "ERROR: project {wildcards.project}: no mMAGs passed the quality gate; reference {params.MAGs_dRep} is empty." >> {log}
            echo "ERROR [pmr_m_project_03_index] project {wildcards.project}: no mMAGs passed the quality gate; cannot build an empty reference. See log: {log}" >&2
            exit 1
        fi
        bwa-mem2 index -p {params.prefix} {params.MAGs_dRep} >> {log} 2>&1
        touch {output.mk}
        """

rule pmr_m_sample_03_mapping_bwamem2:
    input:
        r1 = path + 'results/b-sample/01_trim/{sample}_trimmed_1.fastq.gz',
        r2 = path + 'results/b-sample/01_trim/{sample}_trimmed_2.fastq.gz',
        mk = path + 'markers/pmr_m-project/03_index/bwamem2/{project}.mk',
    output:
        bam = path + 'results/pmr_m-sample/03_mapping/bwamem2/{project}/{sample}.bam',
        mk = path + 'markers/pmr_m-sample/03_mapping/bwamem2/{project}/{sample}.mk',
    params:
        prefix = path + 'results/pmr_m-project/03_index/bwamem2/{project}.MAGs_dRep',
        sort_mem = config["Params_pmr_m_sample"]["03_mapping"].get("samtools_sort_mem", "2G"),
        mapping_extra_args = config["Params_pmr_m_sample"]["03_mapping"]["extra_args"],
    conda: config["conda_envs"]["metalink-main"]
    priority: 60
    threads: config["Params_pmr_m_sample"]["03_mapping"]["threads"]
    log:
        path + 'logs/pmr_m-sample/03_mapping/bwamem2/{project}/{sample}.03_mapping.log'
    benchmark:
        path + 'benchmarks/pmr_m-sample/03_mapping/bwamem2/{project}/{sample}.03_mapping.benchmark'
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

rule pmr_m_project_03_index_minibwa:
    input:
        gated_list = lambda wildcards: pmr_m_gated_input(wildcards, 'mMAGs_gated.list'),
    output:
        mk = path + 'markers/pmr_m-project/03_index/minibwa/{project}.mk',
    params:
        MAGs_folder = path + 'results/pmr_m-project/02_quality_gate/{project}.mMAGs_gated/',
        MAGs_dRep   = path + 'results/pmr_m-project/03_index/minibwa/{project}.MAGs_dRep.fna',
        prefix      = path + 'results/pmr_m-project/03_index/minibwa/{project}.MAGs_dRep',
    conda: config["conda_envs"]["metalink-main"]
    priority: 1
    threads: config["Params_pmr_m_sample"]["03_mapping"]["threads"]
    log:
        path + 'logs/pmr_m-project/03_index/minibwa/{project}.03_index.log'
    benchmark:
        path + 'benchmarks/pmr_m-project/03_index/minibwa/{project}.03_index.benchmark'
    shell:
        """
        mkdir -p $(dirname {params.MAGs_dRep}) 2> {log}
        find -L {params.MAGs_folder} -type f -print0 | xargs -0 -r cat > {params.MAGs_dRep} 2>> {log}
        if [ ! -s {params.MAGs_dRep} ]; then
            echo "ERROR: project {wildcards.project}: no mMAGs passed the quality gate; reference {params.MAGs_dRep} is empty." >> {log}
            echo "ERROR [pmr_m_project_03_index] project {wildcards.project}: no mMAGs passed the quality gate; cannot build an empty reference. See log: {log}" >&2
            exit 1
        fi
        minibwa index {params.MAGs_dRep} {params.prefix} -t {threads} 2>> {log}
        touch {output.mk}
        """

rule pmr_m_sample_03_mapping_minibwa:
    input:
        r1 = path + 'results/b-sample/01_trim/{sample}_trimmed_1.fastq.gz',
        r2 = path + 'results/b-sample/01_trim/{sample}_trimmed_2.fastq.gz',
        mk = path + 'markers/pmr_m-project/03_index/minibwa/{project}.mk',
    output:
        bam = path + 'results/pmr_m-sample/03_mapping/minibwa/{project}/{sample}.bam',
        mk = path + 'markers/pmr_m-sample/03_mapping/minibwa/{project}/{sample}.mk',
    params:
        prefix = path + 'results/pmr_m-project/03_index/minibwa/{project}.MAGs_dRep',
        sort_mem = config["Params_pmr_m_sample"]["03_mapping"].get("samtools_sort_mem", "2G"),
        mapping_extra_args = config["Params_pmr_m_sample"]["03_mapping"]["extra_args"],
    conda: config["conda_envs"]["metalink-main"]
    priority: 60
    threads: config["Params_pmr_m_sample"]["03_mapping"]["threads"]
    log:
        path + 'logs/pmr_m-sample/03_mapping/minibwa/{project}/{sample}.03_mapping.log'
    benchmark:
        path + 'benchmarks/pmr_m-sample/03_mapping/minibwa/{project}/{sample}.03_mapping.benchmark'
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

rule pmr_m_project_03_index_strobealign:
    input:
        gated_list = lambda wildcards: pmr_m_gated_input(wildcards, 'mMAGs_gated.list'),
        r1 = expand(path + 'results/b-sample/01_trim/{sample}_trimmed_1.fastq.gz', sample=SAMPLE_INDEX),
    output:
        mk = path + 'markers/pmr_m-project/03_index/strobealign/{project}.mk',
    params:
        MAGs_folder = path + 'results/pmr_m-project/02_quality_gate/{project}.mMAGs_gated/',
        MAGs_dRep   = path + 'results/pmr_m-project/03_index/strobealign/{project}.MAGs_dRep.fna',
        len_lst     = path + 'results/pmr_m-project/03_index/strobealign/{project}.unique_scaled_length.txt',
    conda: config["conda_envs"]["metalink-main"]
    threads: config["Params_pmr_m_sample"]["03_mapping"]["threads"]
    log:
        path + 'logs/pmr_m-project/03_index/strobealign/{project}.03_index.log'
    benchmark:
        path + 'benchmarks/pmr_m-project/03_index/strobealign/{project}.03_index.benchmark'
    shell:
        """
        mkdir -p $(dirname {params.MAGs_dRep}) 2> {log}
        find -L {params.MAGs_folder} -type f -print0 | xargs -0 -r cat > {params.MAGs_dRep} 2>> {log}
        if [ ! -s {params.MAGs_dRep} ]; then
            echo "ERROR: project {wildcards.project}: no mMAGs passed the quality gate; reference {params.MAGs_dRep} is empty." >> {log}
            echo "ERROR [pmr_m_project_03_index] project {wildcards.project}: no mMAGs passed the quality gate; cannot build an empty reference. See log: {log}" >&2
            exit 1
        fi

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
            strobealign --create-index -t {threads} -r $L {params.MAGs_dRep}
        done >> {log} 2>&1

        touch {output.mk}
        """

rule pmr_m_sample_03_mapping_strobealign:
    input:
        r1 = path + 'results/b-sample/01_trim/{sample}_trimmed_1.fastq.gz',
        r2 = path + 'results/b-sample/01_trim/{sample}_trimmed_2.fastq.gz',
        mk = path + 'markers/pmr_m-project/03_index/strobealign/{project}.mk',
    output:
        bam = path + 'results/pmr_m-sample/03_mapping/strobealign/{project}/{sample}.bam',
        mk = path + 'markers/pmr_m-sample/03_mapping/strobealign/{project}/{sample}.mk',
    params:
        MAGs_dRep = path + 'results/pmr_m-project/03_index/strobealign/{project}.MAGs_dRep.fna',
        sort_mem = config["Params_pmr_m_sample"]["03_mapping"].get("samtools_sort_mem", "2G"),
        mapping_extra_args = config["Params_pmr_m_sample"]["03_mapping"]["extra_args"],
    conda: config["conda_envs"]["metalink-main"]
    priority: 60
    threads: config["Params_pmr_m_sample"]["03_mapping"]["threads"]
    log:
        path + 'logs/pmr_m-sample/03_mapping/strobealign/{project}/{sample}.03_mapping.log'
    benchmark:
        path + 'benchmarks/pmr_m-sample/03_mapping/strobealign/{project}/{sample}.03_mapping.benchmark'
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
        (strobealign --use-index -t {threads} -r $ScaledLength {params.mapping_extra_args} {params.MAGs_dRep} {input.r1} {input.r2} \
            2> {log} \
            | samtools sort -@ {threads} -m {params.sort_mem} -o {output.bam}) >> {log} 2>&1
        touch {output.mk}
        """

rule pmr_m_sample_03_coverm_MAGs:
    input:
        mk = pmr_m_get_03_mapping_rule_mk,
        gated_list = lambda wildcards: pmr_m_gated_input(wildcards, 'mMAGs_gated.list'),
        bam = path + 'results/pmr_m-sample/03_mapping/' + pmr_m_sample_mapping_method + '/{project}/{sample}.bam',
    output:
        out_coverm = path + 'results/pmr_m-sample/03_coverm_MAGs/{project}/{sample}.tsv',
    threads: config["Params_pmr_m_sample"]["03_coverm_MAGs"]["threads"]
    conda: config["conda_envs"]["metalink-main"]
    priority: 80
    params:
        MAGs_folder = path + 'results/pmr_m-project/02_quality_gate/{project}.mMAGs_gated/',
        genometype = config["Params_pmr_m_sample"]["03_coverm_MAGs"]["type"],
        method = config["Params_pmr_m_sample"]["03_coverm_MAGs"]["method"],
        min_read_percent_identity = config["Params_pmr_m_sample"]["03_coverm_MAGs"]["min_read_percent_identity"],
        min_read_aligned_length = config["Params_pmr_m_sample"]["03_coverm_MAGs"]["min_read_aligned_length"],
        contig_end_exclusion = config["Params_pmr_m_sample"]["03_coverm_MAGs"]["contig_end_exclusion"],
        extra_args = config["Params_pmr_m_sample"]["03_coverm_MAGs"]["extra_args"],
    log:
        path + 'logs/pmr_m-sample/03_coverm_MAGs/{project}/{sample}.03_coverm_MAGs.log'
    benchmark:
        path + 'benchmarks/pmr_m-sample/03_coverm_MAGs/{project}/{sample}.03_coverm_MAGs.benchmark'
    shell:
        """
        if [ -d "{params.MAGs_folder}" ]; then
            MAGs_folder={params.MAGs_folder}
        else
            echo "The directory MAGs_folder is not exist."
        fi > {log} 2>&1

        if [ -z "$(ls -A $MAGs_folder)" ]; then
            echo "The directory MAGs_folder is empty, pass."
            touch {output.out_coverm}
        else
            echo "Calculating MAGs abundance..."
            coverm genome \
                -t {threads} \
                -m {params.method} \
                --min-read-percent-identity {params.min_read_percent_identity} \
                --min-read-aligned-length {params.min_read_aligned_length} \
                --contig-end-exclusion {params.contig_end_exclusion} \
                --bam-files {input.bam} \
                --genome-fasta-directory $MAGs_folder \
                -x {params.genometype} \
                -o {output.out_coverm} {params.extra_args}
        fi >> {log} 2>&1
        """

rule pmr_m_sample_03_coverm_contigs:
    priority: 80
    input:
        mk = pmr_m_get_03_mapping_rule_mk,
        gated_list = lambda wildcards: pmr_m_gated_input(wildcards, 'mMAGs_gated.list'),
        bam = path + 'results/pmr_m-sample/03_mapping/' + pmr_m_sample_mapping_method + '/{project}/{sample}.bam',
    output:
        out_coverm = path + 'results/pmr_m-sample/03_coverm_contigs/{project}/{sample}.tsv',
    conda: config["conda_envs"]["metalink-main"]
    threads: config["Params_pmr_m_sample"]["03_coverm_contigs"]["threads"]
    params:
        MAGs_folder = path + 'results/m-project/01_drep_project_95/{project}/dereplicated_genomes/',
        method = config["Params_pmr_m_sample"]["03_coverm_contigs"]["method"],
        min_read_percent_identity = config["Params_pmr_m_sample"]["03_coverm_contigs"]["min_read_percent_identity"],
        min_read_aligned_length = config["Params_pmr_m_sample"]["03_coverm_contigs"]["min_read_aligned_length"],
        contig_end_exclusion = config["Params_pmr_m_sample"]["03_coverm_contigs"]["contig_end_exclusion"],
        extra_args = config["Params_pmr_m_sample"]["03_coverm_contigs"]["extra_args"],
    log:
        path + 'logs/pmr_m-sample/03_coverm_contigs/{project}/{sample}.03_coverm_contigs.log'
    benchmark:
        path + 'benchmarks/pmr_m-sample/03_coverm_contigs/{project}/{sample}.03_coverm_contigs.benchmark'
    shell:
        """
        if [ -d "{params.MAGs_folder}" ]; then
            MAGs_folder={params.MAGs_folder}
        else
            echo "The directory MAGs_folder is not exist."
        fi > {log} 2>&1

        if [ -z "$(ls -A $MAGs_folder)" ]; then
            echo "The directory MAGs_folder is empty, pass."
            touch {output.out_coverm}
        else
            echo "Calculating genes abundance..."
            coverm contig \
                -t {threads} \
                -m {params.method} \
                --min-read-percent-identity {params.min_read_percent_identity} \
                --min-read-aligned-length {params.min_read_aligned_length} \
                --contig-end-exclusion {params.contig_end_exclusion} \
                --output-format dense \
                --bam-files {input.bam} \
                -o {output.out_coverm} {params.extra_args}
        fi >> {log} 2>&1
        """

rule pmr_m_project_04_merge:
    input:
        coverm = expand(path + 'results/pmr_m-sample/03_coverm_MAGs/{project}/{sample}.tsv', project=PROJECT_INDEX, sample=SAMPLE_INDEX),
        tax = path + 'results/pmr_m-project/02_quality_gate/{project}.mMAG_summary_info.filter.tsv',
    output:
        merge_file = path + 'results/pmr_m-project/04_summary/{project}/mMAG_species_abundance_matrix.tsv',
    threads: config["Params_pmr_m_sample"]["05_merge_tax_fun_depth"]["threads"]
    conda: config["conda_envs"]["metalink-main"]
    params:
        tmp_dir = path + 'results/pmr_m-project/04_summary/{project}/tmp_multi_MAGs_Sample/',
    log:
        path + 'logs/pmr_m-project/04_summary/{project}.04_merge.log'
    benchmark:
        path + 'benchmarks/pmr_m-project/04_summary/{project}.04_merge.benchmark'
    shell:
        """
        rm -rf {params.tmp_dir}
        mkdir -p {params.tmp_dir}
        for f in {input.coverm}; do
            if [ -s "$f" ]; then
                cp "$f" {params.tmp_dir}/
            else
                echo "[WARNING] File missing or empty, skipped: $f"
            fi
        done > {log} 2>&1
        Rscript {Dir_script_metalink}/merge_multi_MAGs_Sample_depth.cmd.R \
            -u {params.tmp_dir} \
            --tax_table {input.tax} \
            -o {output} >> {log} 2>&1
        rm -rf {params.tmp_dir}
        """

rule pmr_m_sample_05_merge_tax_fun_depth:
    priority: 30
    input:
        merge_file = path + 'results/pmr_m-project/04_summary/{project}/mMAG_species_abundance_matrix.tsv',
        MAGs_depth = path + 'results/pmr_m-sample/03_coverm_MAGs/{project}/{sample}.tsv',
        contigs_depth = path + 'results/pmr_m-sample/03_coverm_contigs/{project}/{sample}.tsv',
        mk_gate = path + 'markers/pmr_m-project/02_quality_gate/{project}.mk',
    output:
        summary_Tax = path + 'results/pmr_m-sample/05_merge_tax_fun_depth/{project}/summary/{sample}.summary_Tax.tsv',
        summary_ORF_KO = path + 'results/pmr_m-sample/05_merge_tax_fun_depth/{project}/summary/{sample}.summary_ORF_KO.tsv',
        summary_ORF_Tax_KO = path + 'results/pmr_m-sample/05_merge_tax_fun_depth/{project}/summary/{sample}.summary_ORF_Tax_KO.tsv',
        summary_ORF_Source = path + 'results/pmr_m-sample/05_merge_tax_fun_depth/{project}/summary/{sample}.summary_ORF_Source.tsv',
        summary_ORF_Tax_Source = path + 'results/pmr_m-sample/05_merge_tax_fun_depth/{project}/summary/{sample}.summary_ORF_Tax_Source.tsv',
        KO_depth = path + 'results/pmr_m-sample/05_merge_tax_fun_depth/{project}/KO_depth/{sample}.KO_depth.tsv',
        Taxa_abundance = path + 'results/pmr_m-sample/05_merge_tax_fun_depth/{project}/Taxa_abundance/{sample}.Taxa_abundance.tsv',
        KO_Taxa_SumCoverage = path + 'results/pmr_m-sample/05_merge_tax_fun_depth/{project}/KO_Taxa_SumCoverage/{sample}.KO_Taxa_SumCoverage.Taxa.tsv',
        KO_Taxa_RelativeAbundance = path + 'results/pmr_m-sample/05_merge_tax_fun_depth/{project}/KO_Taxa_RelativeAbundance/{sample}.KO_Taxa_RelativeAbundance.Taxa.tsv',
    params:
        taxid_fun = taxid_fun_path,
        MAGs_taxonomy_classify_method = MAGs_taxonomy_classify_method,
        MAGs_function_classify_method = MAGs_function_classify_method,
        MAGs_tax_detail = path + 'results/pmr_m-project/02_quality_gate/{project}.mMAG_summary_info.filter.tsv',
        MAGs_fun_detail = path + 'results/m-project/04_summary/{project}/MAGs_protein.' + MAGs_function_classify_method + '.annotations.gz',
        MAGs_ffn_length = path + 'results/m-project/04_summary/{project}/bin_ffn_length.gz',
        merge_detail = path + 'results/pmr_m-sample/05_merge_tax_fun_depth/{project}/detail_for_each_seq/{sample}.detail_for_each_seq.tsv',
        minORFLen = config["Params_pmr_m_sample"]["05_merge_tax_fun_depth"]["minORFLen"],
        KO_mode = config["Params_pmr_m_sample"]["05_merge_tax_fun_depth"]["KO_mode"],
        min_KO_num_in_one_taxon = config["Params_pmr_m_sample"]["05_merge_tax_fun_depth"]["min_KO_num_in_one_taxon"],
        rank = ",".join(RANKS),
        out_CAMI_format = path + 'results/pmr_m-sample/05_merge_tax_fun_depth/{project}/Taxa_abundance/{sample}.Taxa_abundance.CAMI_format.tsv',
        CAMI_format_output = config["Params_pmr_m_sample"]["05_merge_tax_fun_depth"]["CAMI_format_output"],
        output_detail = config["Params_pmr_m_sample"]["05_merge_tax_fun_depth"]["output_detail"],
    threads: config["Params_pmr_m_sample"]["05_merge_tax_fun_depth"]["threads"]
    conda: config["conda_envs"]["metalink-main"]
    log:
        path + 'logs/pmr_m-project/05_merge_tax_fun_depth/{project}/{sample}.log'
    benchmark:
        path + 'benchmarks/pmr_m-project/05_merge_tax_fun_depth/{project}/{sample}.benchmark'
    shell:
        """
        mMAGs_miss=""
        for f in {params.MAGs_tax_detail} {params.MAGs_fun_detail} \
                {params.MAGs_ffn_length} {input.MAGs_depth} {input.contigs_depth}; do
            [ -s "$f" ] || mMAGs_miss="$mMAGs_miss $f"
        done

        if [ -n "$mMAGs_miss" ]; then
            echo "mMAGs inputs missing/empty; passing mMAGs-related processes for {wildcards.sample}" >> {log}
        else
            for f in {params.MAGs_tax_detail}; do
                [ -s "$f" ] || mMAGs_miss="$mMAGs_miss $f"
            done
        fi

        if [ -n "$mMAGs_miss" ]; then
            echo "ERROR: mMAG_summary_info.filter.tsv empty or missing for {wildcards.sample}; cannot merge taxonomy/function/depth." >> {log}
            echo "  missing/empty mMAGs inputs:$mMAGs_miss" >> {log}
            echo "ERROR [pmr_m_sample_05_merge_tax_fun_depth] sample {wildcards.sample}: required mMAGs inputs (incl. mMAG_summary_info.filter.tsv) missing/empty; check upstream results. See log: {log}" >&2
            exit 1
        else
            echo "Found non-empty mMAG_summary_info.filter.tsv, running pmr-mode." >> {log}
            Rscript {Dir_script_metalink}/combine_taxonomy_function_depth.pmr_m.cmd.R \
                --in_taxid_fun {params.taxid_fun} \
                --MAGs_taxonomy_classify_method {params.MAGs_taxonomy_classify_method} \
                --MAGs_function_classify_method {params.MAGs_function_classify_method} \
                --in_MAGs_tax_detail {params.MAGs_tax_detail} \
                --in_MAGs_fun_detail {params.MAGs_fun_detail} \
                --in_MAGs_depth {input.MAGs_depth} \
                --in_MAGs_ffn_length {params.MAGs_ffn_length} \
                --in_contig_depth {input.contigs_depth} \
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
                --detail {params.merge_detail} >>{log} 2>&1
            if [ "{params.output_detail}" = "True" ]; then
                pigz -p {threads} -r -f {params.merge_detail} > {params.merge_detail}.gz 2>> {log}
            fi
        fi
        """

rule pmr_m_sample_05_visualize:
    priority: 40
    input:
        summary_Tax = path + 'results/pmr_m-sample/05_merge_tax_fun_depth/{project}/summary/{sample}.summary_Tax.tsv',
        summary_ORF_KO = path + 'results/pmr_m-sample/05_merge_tax_fun_depth/{project}/summary/{sample}.summary_ORF_KO.tsv',
        summary_ORF_Tax_KO = path + 'results/pmr_m-sample/05_merge_tax_fun_depth/{project}/summary/{sample}.summary_ORF_Tax_KO.tsv',
    output:
        mk = path + 'markers/pmr_m-sample/05_visualize/{project}/{sample}.mk',
    params:
        prefix = path + 'results/pmr_m-sample/05_visualization/{project}/{sample}',
    threads: 1
    conda: config["conda_envs"]["metalink-main"]
    log:
        path + 'logs/pmr_m-sample/05_visualize/{project}/{sample}.log'
    benchmark:
        path + 'benchmarks/pmr_m-sample/05_visualize/{project}/{sample}.benchmark'
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
            echo "All summary files non-empty, running pmr visualization." >>{log}
            Rscript {Dir_script_metalink}/plot_summary.cmd.R \
                --TaxFormat pmr \
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

rule pmr_m_sample_06_module_completeness:
    priority: 40
    input:
        KO_Taxa_RelativeAbundance = path + 'results/pmr_m-sample/05_merge_tax_fun_depth/{project}/KO_Taxa_RelativeAbundance/{sample}.KO_Taxa_RelativeAbundance.Taxa.tsv',
    output:
        mk = path + 'markers/pmr_m-sample/06_module_completeness/{project}/{sample}.{rank}.mk',
    params:
        KO_Taxa_RelativeAbundance_ranked = path + 'results/pmr_m-sample/05_merge_tax_fun_depth/{project}/KO_Taxa_RelativeAbundance/{sample}.KO_Taxa_RelativeAbundance.Taxa.{rank}.tsv',
        KO_Taxa_SumCoverage_ranked = path + 'results/pmr_m-sample/05_merge_tax_fun_depth/{project}/KO_Taxa_SumCoverage/{sample}.KO_Taxa_SumCoverage.Taxa.{rank}.tsv',
        KO_Taxa_SumCoverage = path + 'results/pmr_m-sample/05_merge_tax_fun_depth/{project}/KO_Taxa_SumCoverage/{sample}.KO_Taxa_SumCoverage.Taxa.tsv',
        Module_Taxa_RA = path + 'results/pmr_m-sample/06_module_completeness/{project}/{sample}/{rank}/',
        pathway_split = "TRUE",
        min_KO_num_in_one_taxon = config["Params_pmr_m_sample"]["06_mclink_taxon"]["min_KO_num_in_one_taxon"],
        min_mean_KO_RA_in_one_taxon = config["Params_pmr_m_sample"]["06_mclink_taxon"]["min_mean_KO_RA_in_one_taxon"],
        completeness_threshold = config["Params_pmr_m_sample"]["06_mclink_taxon"]["completeness_threshold"],
        plus_scale_method = config["Params_pmr_m_sample"]["06_mclink_taxon"]["plus_scale_method"],
        comma_scale_method = config["Params_pmr_m_sample"]["06_mclink_taxon"]["comma_scale_method"],
        TopN = config["Params_pmr_m_sample"]["06_mclink_taxon"]["TopN"],
    threads: config["Params_pmr_m_sample"]["05_merge_tax_fun_depth"]["threads"]
    conda: config["conda_envs"]["metalink-main"]
    log:
        path + 'logs/pmr_m-sample/06_module_completeness/{project}/{sample}.{rank}.log'
    benchmark:
        path + 'benchmarks/pmr_m-sample/06_module_completeness/{project}/{sample}.{rank}.benchmark'
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
            echo "KO_Taxa_RelativeAbundance empty/missing for {wildcards.sample}; skip pmr_m_sample_03_module_completeness." >>{log}
        fi
        touch {output.mk}
        """

rule pmr_m_project_01_cat_Taxa_abundance_dcast:
    input:
        Taxa_abundance = expand(path + 'results/pmr_m-sample/05_merge_tax_fun_depth/{project}/Taxa_abundance/{sample}.Taxa_abundance.tsv', project=PROJECT_INDEX, sample=SAMPLE_INDEX),
    output:
        out_long = path + 'results/pmr_m-project/01_Taxa_Sample/{project}/Taxa_Sample.long.tsv',
    params:
        tmp_dir = path + 'results/pmr_m-project/tmp/{project}/',
        prefix_wide = path + 'results/pmr_m-project/01_Taxa_Sample/{project}/Taxa_Sample.wide.tsv',
    conda: config["conda_envs"]["metalink-main"]
    threads: config["Params_pmr_m_sample"]["05_merge_tax_fun_depth"]["threads"]
    log: path + 'logs/pmr_m-project/01_cat_Taxa_abundance_dcast/{project}.log'
    benchmark:
        path + 'benchmarks/pmr_m-project/01_cat_Taxa_abundance_dcast/{project}.benchmark'
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
            --value_var Taxa_SumCoverage,Taxa_RelativeAbundance \
            --out_long {output.out_long} \
            --out_wide {params.prefix_wide} \
            --threads {threads} \
            > {log} 2>&1
        """

rule pmr_m_project_01_Taxa_visualize:
    input:
        out_long = path + 'results/pmr_m-project/01_Taxa_Sample/{project}/Taxa_Sample.long.tsv',
    output:
        mk = path + 'markers/pmr_m-project/01_Taxa_visualization/{project}.mk',
    params:
        prefix = path + 'results/pmr_m-project/01_Taxa_visualization/{project}',
        prefix_wide = path + 'results/pmr_m-project/01_Taxa_Sample/{project}/Taxa_Sample.Taxa_SumCoverage.wide.tsv',
        TopN = config["Params_pmr_m_project"]["01_Taxa_visualize"]["TopN"],
    threads: 1
    conda: config["conda_envs"]["metalink-main"]
    log:
        path + 'logs/pmr_m-project/01_Taxa_visualization/{project}.log'
    benchmark:
        path + 'benchmarks/pmr_m-project/01_Taxa_visualization/{project}.benchmark'
    shell:
        """
        mkdir -p $(dirname {params.prefix}) 2>>{log}
        missing=""
        for f in {input.out_long} {params.prefix_wide}; do
            [ -s "$f" ] || missing="$missing $f"
        done
        if [ -n "$missing" ]; then
            echo "Missing/empty inputs:$missing; skip pmr project_01_Taxa_visualize for {wildcards.project}" >>{log}
        else
            echo "All input files non-empty, running pmr project_01_Taxa_visualize." >>{log}
            Rscript {Dir_script_metalink}/plot_bar.cmd.R \
                --Type Taxa \
                --TopN {params.TopN} \
                --Input {params.prefix_wide} \
                --Output {params.prefix} >>{log} 2>&1
        fi
        touch {output.mk}
        """

rule pmr_m_project_02_cat_KO_depth_dcast:
    input:
        KO_depth = expand(path + 'results/pmr_m-sample/05_merge_tax_fun_depth/{project}/KO_depth/{sample}.KO_depth.tsv', project=PROJECT_INDEX, sample=SAMPLE_INDEX),
    output:
        out_long = path + 'results/pmr_m-project/02_KO_Sample/{project}/KO_Sample.long.tsv',
    params:
        tmp_dir = path + 'results/pmr_m-project/tmp/{project}/',
        prefix_wide = path + 'results/pmr_m-project/02_KO_Sample/{project}/KO_Sample.wide.tsv',
    threads: config["Params_pmr_m_sample"]["05_merge_tax_fun_depth"]["threads"]
    conda: config["conda_envs"]["metalink-main"]
    log:
        path + 'logs/pmr_m-project/02_cat_KO_depth_dcast/{project}.log'
    benchmark:
        path + 'benchmarks/pmr_m-project/02_cat_KO_depth_dcast/{project}.benchmark'
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

rule pmr_m_project_03_community_module:
    input:
        out_long = path + 'results/pmr_m-project/02_KO_Sample/{project}/KO_Sample.long.tsv',
    output:
        mk = path + 'markers/pmr_m-project/03_community_module/{project}.mk',
    params:
        wide_RelativeAbundance = path + 'results/pmr_m-project/02_KO_Sample/{project}/KO_Sample.KO_RelativeAbundance.wide.tsv',
        wide_SumCoverage = path + 'results/pmr_m-project/02_KO_Sample/{project}/KO_Sample.KO_SumCoverage.wide.tsv',
        pathway_split = "TRUE",
        community_output = path + 'results/pmr_m-project/03_Module_Sample/{project}/',
        min_sample_num_per_KO = config["Params_pmr_m_project"]["03_mclink_sample"]["min_sample_num_per_KO"],
        preva = config["Params_pmr_m_project"]["03_mclink_sample"]["preva"],
        abund = config["Params_pmr_m_project"]["03_mclink_sample"]["abund"],
        completeness_threshold = config["Params_pmr_m_project"]["03_mclink_sample"]["completeness_threshold"],
        plus_scale_method = config["Params_pmr_m_project"]["03_mclink_sample"]["plus_scale_method"],
        comma_scale_method = config["Params_pmr_m_project"]["03_mclink_sample"]["comma_scale_method"],
    threads: config["Params_pmr_m_sample"]["05_merge_tax_fun_depth"]["threads"]
    conda: config["conda_envs"]["metalink-main"]
    log:
        path + 'logs/pmr_m-project/03_community_module/{project}.log'
    benchmark:
        path + 'benchmarks/pmr_m-project/03_community_module/{project}.benchmark'
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
            echo "KO_Sample.KO_RelativeAbundance.wide.tsv empty/missing for {wildcards.project}; skip MAGs_project_03_community_module." >>{log}
        fi
        touch {output.mk}
        """

rule pmr_m_project_03_Fun_visualize:
    input:
        mk = path + 'markers/pmr_m-project/03_community_module/{project}.mk',
    output:
        mk = path + 'markers/pmr_m-project/03_Fun_visualization/{project}.mk',
    params:
        prefix = path + 'results/pmr_m-project/03_Fun_visualization/{project}',
        prefix_wide = path + 'results/pmr_m-project/03_Module_Sample/{project}/SumCoverage/mc_abundance_filtered_by_completeness.tsv',
        TopN = config["Params_pmr_m_project"]["03_Fun_visualize"]["TopN"],
    threads: 1
    conda: config["conda_envs"]["metalink-main"]
    log:
        path + 'logs/pmr_m-project/03_Fun_visualization/{project}.log'
    benchmark:
        path + 'benchmarks/pmr_m-project/03_Fun_visualization/{project}.benchmark'
    shell:
        """
        mkdir -p $(dirname {params.prefix}) 2>>{log}
        missing=""
        for f in {params.prefix_wide}; do
            [ -s "$f" ] || missing="$missing $f"
        done
        if [ -n "$missing" ]; then
            echo "Missing/empty inputs:$missing; skip pmr project_03_Fun_visualize for {wildcards.project}" >>{log}
        else
            echo "All input files non-empty, running pmr project_03_Fun_visualize." >>{log}
            Rscript {Dir_script_metalink}/plot_bar.cmd.R \
                --Type Module \
                --TopN {params.TopN} \
                --Input {params.prefix_wide} \
                --Output {params.prefix} >>{log} 2>&1
        fi
        touch {output.mk}
        """

rule pmr_m_project_05_merge_tax_KO_depth:
    input:
        ko_depth = expand(path + 'results/pmr_m-sample/05_merge_tax_fun_depth/{project}/KO_depth/{sample}.KO_depth.tsv', project=PROJECT_INDEX, sample=SAMPLE_INDEX),
    output:
        out_long = path + 'results/pmr_m-project/05_merge_tax_KO_depth/{project}/{rank}.mean.long.tsv',
    params:
        detail = lambda wc: expand(path + 'results/pmr_m-sample/05_merge_tax_fun_depth/{project}/detail_for_each_seq/{sample}.detail_for_each_{rank}.tsv', project=PROJECT_INDEX, sample=SAMPLE_INDEX, rank=wc.rank),  # 只展开当前 rank
        tmp_dir = path + 'results/pmr_m-project/tmp/{project}/',
        prefix_wide = path + 'results/pmr_m-project/05_merge_tax_KO_depth/{project}/{rank}.wide.tsv',
        preva = config["Params_pmr_m_project"]["05_merge_tax_KO_depth"]["preva"],
        abund = config["Params_pmr_m_project"]["05_merge_tax_KO_depth"]["abund"],
    threads: config["Params_pmr_m_sample"]["05_merge_tax_fun_depth"]["threads"]
    conda: config["conda_envs"]["metalink-main"]
    log:
        path + 'logs/pmr_m-project/05_merge_tax_KO_depth/{project}/{rank}.log'
    benchmark:
        path + 'benchmarks/pmr_m-project/05_merge_tax_KO_depth/{project}/{rank}.benchmark'
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

rule pmr_m_project_05_merge_tax_Module_RelativeAbundance:
    input:
        mk = lambda wc: expand(path + 'markers/pmr_m-sample/06_module_completeness/{project}/{sample}.{rank}.mk', project=PROJECT_INDEX, sample=SAMPLE_INDEX, rank=wc.rank),
    output:
        out_long = path + 'results/pmr_m-project/05_merge_tax_Module_RelativeAbundance/{project}/{rank}.mean.long.tsv',
    params:
        detail = lambda wc: expand(path + 'results/pmr_m-sample/06_module_completeness/{project}/{sample}/{rank}/RelativeAbundance/mc_abundance.tsv', project=PROJECT_INDEX, sample=SAMPLE_INDEX, rank=wc.rank),
        prefix_wide = path + 'results/pmr_m-project/05_merge_tax_Module_RelativeAbundance/{project}/{rank}.wide.tsv',
        tmp_dir = path + 'results/pmr_m-project/tmp/{project}/',
        preva = config["Params_pmr_m_project"]["05_merge_tax_Module_RelativeAbundance"]["preva"],
        abund = config["Params_pmr_m_project"]["05_merge_tax_Module_RelativeAbundance"]["abund"],
    threads: config["Params_pmr_m_sample"]["05_merge_tax_fun_depth"]["threads"]
    conda: config["conda_envs"]["metalink-main"]
    log:
        path + 'logs/pmr_m-project/05_merge_tax_Module_RelativeAbundance/{project}/{rank}.log'
    benchmark:
        path + 'benchmarks/pmr_m-project/05_merge_tax_Module_RelativeAbundance/{project}/{rank}.benchmark'
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

rule pmr_m_project_05_merge_tax_Module_SumCoverage:
    input:
        mk = lambda wc: expand(path + 'markers/pmr_m-sample/06_module_completeness/{project}/{sample}.{rank}.mk', project=PROJECT_INDEX, sample=SAMPLE_INDEX, rank=wc.rank),
    output:
        out_long = path + 'results/pmr_m-project/05_merge_tax_Module_SumCoverage/{project}/{rank}.mean.long.tsv',
    params:
        detail = lambda wc: expand(path + 'results/pmr_m-sample/06_module_completeness/{project}/{sample}/{rank}/SumCoverage/mc_abundance.tsv', project=PROJECT_INDEX, sample=SAMPLE_INDEX, rank=wc.rank),
        prefix_wide = path + 'results/pmr_m-project/05_merge_tax_Module_SumCoverage/{project}/{rank}.wide.tsv',
        tmp_dir = path + 'results/pmr_m-project/tmp/{project}/',
        preva = config["Params_pmr_m_project"]["05_merge_tax_Module_SumCoverage"]["preva"],
        abund = config["Params_pmr_m_project"]["05_merge_tax_Module_SumCoverage"]["abund"],
    threads: config["Params_pmr_m_sample"]["05_merge_tax_fun_depth"]["threads"]
    conda: config["conda_envs"]["metalink-main"]
    log:
        path + 'logs/pmr_m-project/05_merge_tax_Module_SumCoverage/{project}/{rank}.log'
    benchmark:
        path + 'benchmarks/pmr_m-project/05_merge_tax_Module_SumCoverage/{project}/{rank}.benchmark'
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

rule pmr_m_project_06_module_completeness:
    input:
        KO_Taxa_RelativeAbundance = path + 'results/pmr_m-project/05_merge_tax_KO_depth/{project}/{rank}.mean.long.tsv',
    output:
        mk = path + 'markers/pmr_m-project/06_module_completeness/{project}/{rank}.mk',
    params:
        KO_Taxa_RelativeAbundance_ranked = path + 'results/pmr_m-project/05_merge_tax_KO_depth/{project}/{rank}.Taxa_KO_RelativeAbundance.mean.wide.tsv',
        KO_Taxa_SumCoverage_ranked = path + 'results/pmr_m-project/05_merge_tax_KO_depth/{project}/{rank}.Taxa_KO_SumCoverage.mean.wide.tsv',
        Module_Taxa_RA = path + 'results/pmr_m-project/06_module_completeness/{project}/{rank}/',
        pathway_split = "TRUE",
        min_KO_num_in_one_taxon = config["Params_pmr_m_project"]["06_mclink_taxon"]["min_KO_num_in_one_taxon"],
        min_mean_KO_RA_in_one_taxon = config["Params_pmr_m_project"]["06_mclink_taxon"]["min_mean_KO_RA_in_one_taxon"],
        completeness_threshold = config["Params_pmr_m_project"]["06_mclink_taxon"]["completeness_threshold"],
        plus_scale_method = config["Params_pmr_m_project"]["06_mclink_taxon"]["plus_scale_method"],
        comma_scale_method = config["Params_pmr_m_project"]["06_mclink_taxon"]["comma_scale_method"],
        TopN = config["Params_pmr_m_project"]["06_mclink_taxon"]["TopN"],
    threads: config["Params_pmr_m_sample"]["05_merge_tax_fun_depth"]["threads"]
    conda: config["conda_envs"]["metalink-main"]
    log:
        path + 'logs/pmr_m-project/06_module_completeness/{project}/{rank}.log'
    benchmark:
        path + 'benchmarks/pmr_m-project/06_module_completeness/{project}/{rank}.benchmark'
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
            echo "KO_Taxa_RelativeAbundance empty/missing for {wildcards.project}; skip pmr_m_project_module_completeness." >>{log}
        fi
        touch {output.mk}
        """

rule pmr_m_sample_12_delete_bam:
    input:
        out_coverm = path + 'results/pmr_m-sample/03_coverm_MAGs/{project}/{sample}.tsv',
        out_coverm_contigs = path + 'results/pmr_m-sample/03_coverm_contigs/{project}/{sample}.tsv',
    output:
        mk = path + 'markers/pmr_m-sample/12_delete_bam/{project}/{sample}.mk'
    threads: 1
    priority: 100
    params:
        bam = path + 'results/pmr_m-sample/03_mapping/' + pmr_m_sample_mapping_method + '/{project}/{sample}.bam',
    log:
        path + 'logs/pmr_m-sample/12_delete_bam/{project}/{sample}.12_delete_bam.log'
    benchmark:
        path + 'benchmarks/pmr_m-sample/12_delete_bam/{project}/{sample}.12_delete_bam.benchmark'
    shell:
        """
        if [ -e "{params.bam}" ];then rm -r {params.bam}*;else echo Folder merged.bam not exist.;fi > {log} 2>&1
        touch {output.mk}
        """
