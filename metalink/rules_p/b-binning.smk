'''
/*
 * @Date: 2022-12-04 17:15:47
 * @LastEditors: LiuyangLi
 * @LastEditTime: 2024-10-25 16:57:02
 * @FilePath: b-binning.smk
 * @Description: base mode (include metawrap-based binning)
 */
'''
import os
import sys

Dir_script_metalink = config["Directory"]["script"]["metalink"]
kegg_modules_ko_info_path = config["Directory"]["script"]["kegg_modules_ko_info_path"]
kraken2_taxonomy_db = config["Directory"]["database"]["kraken2_taxonomy_db"]


def _fmt_num(x):
    """Format a number for file/directory names; drop trailing .0 from floats."""
    if isinstance(x, float) and x.is_integer():
        return str(int(x))
    return str(x)


folder_checks = [
                path + 'results/',
                path + 'results/b-sample/',
                path + 'results/b-sample/05_metawrap_binning/',
                path + 'results/b-sample/06_metawrap_bin_refinement/',
                path + 'results/b-sample/08_gather_bin_high_quality/',
                path + 'results/b-sample/08_gather_bin_high_quality/bin_genome_fna_renamed/',
                path + 'results/b-sample/08_gather_bin_high_quality/checkm_stats/',
                path + 'results/b-project/'
                ]

for item in folder_checks:
    if not os.path.isdir(item): os.mkdir(item)


rule b_sample_05_metawrap_binning:
    input:
        r1 = path + 'results/b-sample/01_trim/{sample}_trimmed_1.fastq.gz',
        r2 = path + 'results/b-sample/01_trim/{sample}_trimmed_2.fastq.gz',
        bam = path + 'results/b-sample/03_mapping/' + base_sample_mapping_method + '/{sample}.bam',
    output:
        mk = path + 'markers/b-sample/05_metawrap_binning/{sample}.mk'
    threads: config["Params_b_sample"]["05_metawrap_binning_depth"]["threads"]
    conda: config["conda_envs"]["metalink-binning"]
    params:
        min_read_percent_identity = config["Params_b_sample"]["05_metawrap_binning_depth"]["min_read_percent_identity"],
        contig_end_exclusion = config["Params_b_sample"]["05_metawrap_binning_depth"]["contig_end_exclusion"],
        min_contig_length = config["Params_b_sample"]["05_metawrap_binning_depth"]["min_contig_length"],
        ctg = path + 'results/b-sample/03_index/' + base_sample_mapping_method + '/{sample}.assembly_filter.fna',
        path = path + 'results/b-sample/05_metawrap_binning/{sample}/',
        binning_tools = config["Params_b_sample"]["05_metawrap_binning_depth"].get("binning_tools", "metabat2,maxbin2,concoct"),
        max_memory_gb = config["Params_b_sample"]["05_metawrap_binning_depth"].get("max_memory_gb", 512),
        extra_args = config["Params_b_sample"]["05_metawrap_binning_depth"]["extra_args"],
    priority: 60
    log:
        path + 'logs/b-sample/05_metawrap_binning/{sample}.metawrap_binning.log'
    benchmark:
        path + 'benchmarks/b-sample/05_metawrap_binning/{sample}.metawrap_binning.benchmark'
    shell:
        """
        # Validate at least 2 tools
        TOOL_COUNT=$(echo "{params.binning_tools}" | tr ',' '\n' | wc -l)
        if [ "$TOOL_COUNT" -lt 2 ]; then
            echo "ERROR: At least 2 binning tools required, got: {params.binning_tools}" > {log} 2>&1
            exit 1
        fi

        echo "[BINNING PARAMS] tools={params.binning_tools}, max_memory_gb={params.max_memory_gb}, min_contig_length={params.min_contig_length}, min_read_percent_identity={params.min_read_percent_identity}, contig_end_exclusion={params.contig_end_exclusion}" > {log}
        # Build tool flags from comma-separated list
        TOOL_FLAGS=""
        IFS=',' read -ra TOOLS <<< "{params.binning_tools}"
        for tool in "${{TOOLS[@]}}"; do
            case "$tool" in
                metabat2) TOOL_FLAGS="$TOOL_FLAGS --metabat2" ;;
                maxbin2)  TOOL_FLAGS="$TOOL_FLAGS --maxbin2" ;;
                concoct)  TOOL_FLAGS="$TOOL_FLAGS --concoct" ;;
            esac
        done

        if [ -f {input.bam} ]; then
            echo "Using previous bam file..."
            echo "[BINNING CMD] bash {Dir_script_metalink}/metawrap-scripts/binning.sh -a {params.ctg} -b {input.bam} -l {params.min_contig_length} --maxEdgeBases {params.contig_end_exclusion} --percentIdentity {params.min_read_percent_identity} -o {params.path} -t {threads} $TOOL_FLAGS -m {params.max_memory_gb} {input.r1} {input.r2} {params.extra_args}" >> {log}
            bash {Dir_script_metalink}/metawrap-scripts/binning.sh -a {params.ctg} -b {input.bam} \
                -l {params.min_contig_length} --maxEdgeBases {params.contig_end_exclusion} \
                --percentIdentity {params.min_read_percent_identity} \
                -o {params.path} -t {threads} $TOOL_FLAGS -m {params.max_memory_gb} \
                {input.r1} {input.r2} {params.extra_args} >> {log} 2>&1
        else
            echo "Bam file not exist, starting aligning reads to contigs..."
            echo "[BINNING CMD] bash {Dir_script_metalink}/metawrap-scripts/binning.sh -a {params.ctg} -l {params.min_contig_length} --maxEdgeBases {params.contig_end_exclusion} --percentIdentity {params.min_read_percent_identity} -o {params.path} -t {threads} $TOOL_FLAGS -m {params.max_memory_gb} {input.r1} {input.r2} {params.extra_args}" >> {log}
            bash {Dir_script_metalink}/metawrap-scripts/binning.sh -a {params.ctg} \
                -l {params.min_contig_length} --maxEdgeBases {params.contig_end_exclusion} \
                --percentIdentity {params.min_read_percent_identity} \
                -o {params.path} -t {threads} $TOOL_FLAGS -m {params.max_memory_gb} \
                {input.r1} {input.r2} {params.extra_args} >> {log} 2>&1
        fi >> {log} 2>&1
        touch {output.mk}
        """

rule b_sample_06_metawrap_bin_refinement:
    input:
        mk_1 = path + 'markers/b-sample/05_metawrap_binning/{sample}.mk',
    output:
        mk = path + 'markers/b-sample/06_metawrap_bin_refinement/{sample}.mk'
    threads: config["Params_b_sample"]["06_bin_checkm"]["threads"]
    conda: config["conda_envs"]["metalink-binning"]
    params:
        sn = '{sample}',
        path = path + 'results/b-sample/06_metawrap_bin_refinement/{sample}/',
        metabat2 = path + 'results/b-sample/05_metawrap_binning/{sample}/metabat2_bins/',
        maxbin2 = path + 'results/b-sample/05_metawrap_binning/{sample}/maxbin2_bins/',
        concoct = path + 'results/b-sample/05_metawrap_binning/{sample}/concoct_bins/',
        min_completeness = config["Params_b_sample"]["06_bin_checkm"].get("min_completeness", 50),
        max_contamination = config["Params_b_sample"]["06_bin_checkm"].get("max_contamination", 10),
        max_memory_gb = config["Params_b_sample"]["06_metawrap_bin_refinement"].get("max_memory_gb", 512),
        extra_args = config["Params_b_sample"]["06_metawrap_bin_refinement"]["extra_args"],
    priority: 60
    log:
        path + 'logs/b-sample/06_metawrap_bin_refinement/{sample}.metawrap_bin_refinement.log'
    benchmark:
        path + 'benchmarks/b-sample/06_metawrap_bin_refinement/{sample}.metawrap_bin_refinement.benchmark'
    shell:
        """
        echo "[MetaWrap bin_refinement PARAMS] min_completeness={params.min_completeness}, max_contamination={params.max_contamination}, max_memory_gb={params.max_memory_gb}" > {log}
        echo '[MetaWrap bin_refinement CMD] python {Dir_script_metalink}/run_metawrap_bin_refinement.py {params.metabat2} {params.maxbin2} {params.concoct} {params.path} {threads} {Dir_script_metalink} --memory-gb {params.max_memory_gb} --min-completeness {params.min_completeness} --max-contamination {params.max_contamination} --extra-args="{params.extra_args}"' >> {log}
        if [ -d "{params.path}" ];then rm -r {params.path};else echo "Folder not exist, MetaWrap bin_refinement good to go." >> {log};fi
        python {Dir_script_metalink}/run_metawrap_bin_refinement.py \
            {params.metabat2} \
            {params.maxbin2} \
            {params.concoct} \
            {params.path} \
            {threads} \
            {Dir_script_metalink}/metawrap-scripts/ \
            --memory-gb {params.max_memory_gb} \
            --min-completeness {params.min_completeness} \
            --max-contamination {params.max_contamination} \
            --extra-args="{params.extra_args}" >> {log} 2>&1
        touch {output.mk}
        """

rule b_sample_08_gather_bin_checkm_stats:
    priority: 60
    input:
        mk_4 = path + 'markers/b-sample/06_metawrap_bin_refinement/{sample}.mk',
    output:
        tsv = path + 'results/b-sample/08_gather_bin_high_quality/checkm_stats/{sample}.csv',
    params:
        sample = '{sample}',
        input_bin_dir = path + 'results/b-sample/06_metawrap_bin_refinement/{sample}/metawrap_' + _fmt_num(config["Params_b_sample"]["06_bin_checkm"].get("min_completeness", 50)) + '_' + _fmt_num(config["Params_b_sample"]["06_bin_checkm"].get("max_contamination", 10)) + '_bins',
        min_completeness = config["Params_b_sample"]["06_bin_checkm"].get("min_completeness", 50),
        max_contamination = config["Params_b_sample"]["06_bin_checkm"].get("max_contamination", 10),
        min_quality_score = config["Params_b_sample"]["06_bin_checkm"].get("min_quality_score", 50),
        stats = path + 'results/b-sample/06_metawrap_bin_refinement/{sample}/metawrap_' + _fmt_num(config["Params_b_sample"]["06_bin_checkm"].get("min_completeness", 50)) + '_' + _fmt_num(config["Params_b_sample"]["06_bin_checkm"].get("max_contamination", 10)) + '_bins.Chdb.tsv',
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
            python {Dir_script_metalink}/gather_bins_checkm.py {params.sample} "$STATS" {params.input_bin_dir} {params.bin_high_quality} {params.bin_high_quality_renamed} {output.tsv} {params.min_completeness} {params.max_contamination} {params.min_quality_score} >> {log} 2>&1
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

        find {params.bin_high_quality_renamed} -maxdepth 1 -type f \( -name "*.fa" -o -name "*.fasta" -o -name "*.fna" \) | while IFS= read -r genome; do
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

# The b-sample BAM is consumed by b_sample_04_bam_summary, b_sample_07_contigs_prodigal
# (jgi depth) and b_sample_05_metawrap_binning; the inputs below are their final
# downstream products, so the BAM is only deleted once all consumers are done.
rule b_sample_12_delete_bam:
    input:
        mk = path + 'markers/b-sample/08_bin_prodigal/{sample}.mk',
        stat_raw = path + 'results/b-sample/04_bam_summary/{sample}.contigs_bam_stats.tsv',
        contigs_fna = path + 'results/b-sample/07_contigs/contigs_fna/{sample}.fna',
    output:
        mk = path + 'markers/b-sample/12_delete_bam/{sample}.mk'
    threads: 1
    priority: 100
    params:
        bam = path + 'results/b-sample/03_mapping/' + base_sample_mapping_method + '/{sample}.bam',
    log:
        path + 'logs/b-sample/12_delete_bam/{sample}.log'
    benchmark:
        path + 'benchmarks/b-sample/12_delete_bam/{sample}.benchmark'
    shell:
        """
        if [ -e "{params.bam}" ];then rm -r {params.bam}*;else echo Folder trimmed.bam not exist.;fi > {log} 2>&1
        touch {output.mk}
        """

# The metaWRAP binning work_files/assembly copy is a binning-internal
# intermediate with no downstream consumer; delete it once binning and
# refinement have finished (08_bin_prodigal is the last step of that chain).
rule b_sample_12_delete_metawrap_work_files:
    input:
        mk = path + 'markers/b-sample/08_bin_prodigal/{sample}.mk',
    output:
        mk = path + 'markers/b-sample/12_delete_metawrap_work_files/{sample}.mk'
    threads: 1
    priority: 100
    params:
        asb = path + 'results/b-sample/05_metawrap_binning/{sample}/work_files/assembly'
    log:
        path + 'logs/b-sample/12_delete_metawrap_work_files/{sample}.log'
    benchmark:
        path + 'benchmarks/b-sample/12_delete_metawrap_work_files/{sample}.benchmark'
    shell:
        """
        if [ -e "{params.asb}.fa" ];then rm -r {params.asb}*;else echo Folder metabat2.assembly.fa not exist.;fi > {log} 2>&1
        touch {output.mk}
        """
