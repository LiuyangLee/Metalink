'''
/*
 * @Date: 2022-12-04 17:15:47
 * @LastEditors: LiuyangLi
 * @LastEditTime: 2024-07-20 13:03:29
 * @FilePath: vc.smk
 * @Description:
 */
'''

import os
import sys

Dir_script_metalink = config["Directory"]["script"]["metalink"]
kegg_modules_ko_info_path = config["Directory"]["script"]["kegg_modules_ko_info_path"]
taxid_tax_path = config["Directory"]["database"]["taxid_tax_path"]
taxid_fun_path = config["Directory"]["database"]["taxid_fun_path"]
contigs_taxonomy_classify_method = config["Params_c_sample"]["01_taxonomy"]["contigs_taxonomy_classify_method"]
contigs_function_classify_method = config["Params_c_sample"]["02_function"]["contigs_function_classify_method"]

folder_checks = [
                path + 'results/',
                path + 'results/vc-sample/',
                path + 'results/vc-sample/tmp/',
                ]
for item in folder_checks:
    if not os.path.isdir(item): os.mkdir(item)


rule vc_sample_05_merge_tax_fun_depth:
    priority: 30
    input:
        contigs_tax_mk = path + 'markers/c-sample/01_taxonomy/classify_done/{sample}.mk',
        contigs_fun_mk = path + 'markers/c-sample/02_function/classify_done/{sample}.mk',
        contigs_depth = path + 'results/b-sample/07_contigs/contigs_filter_depth/{sample}',
        contigs_ffn_length = path + 'results/b-sample/07_contigs/contigs_ffn_length/{sample}',
        reads_length = path + 'results/b-sample/01_trim_stats/{sample}_trimmed_1.fastq.stats',
        contigs_ORF_loci = path + 'results/b-sample/07_contigs/contigs_ORF_loci/{sample}',
        vMAGs_summar_mk = path + 'markers/v-sample/10_summary/{sample}.mk',
    output:
        summary_Tax = path + 'results/vc-sample/05_merge_tax_fun_depth/summary/{sample}.summary_Tax.tsv',
        summary_ORF_KO = path + 'results/vc-sample/05_merge_tax_fun_depth/summary/{sample}.summary_ORF_KO.tsv',
        summary_ORF_Tax_KO = path + 'results/vc-sample/05_merge_tax_fun_depth/summary/{sample}.summary_ORF_Tax_KO.tsv',
        summary_ORF_Source = path + 'results/vc-sample/05_merge_tax_fun_depth/summary/{sample}.summary_ORF_Source.tsv',
        summary_ORF_Tax_Source = path + 'results/vc-sample/05_merge_tax_fun_depth/summary/{sample}.summary_ORF_Tax_Source.tsv',
        KO_depth = path + 'results/vc-sample/05_merge_tax_fun_depth/KO_depth/{sample}.KO_depth.tsv',
        Taxa_abundance = path + 'results/vc-sample/05_merge_tax_fun_depth/Taxa_abundance/{sample}.Taxa_abundance.tsv',
        KO_Taxa_SumCoverage = path + 'results/vc-sample/05_merge_tax_fun_depth/KO_Taxa_SumCoverage/{sample}.KO_Taxa_SumCoverage.Taxa.tsv',
        KO_Taxa_RelativeAbundance = path + 'results/vc-sample/05_merge_tax_fun_depth/KO_Taxa_RelativeAbundance/{sample}.KO_Taxa_RelativeAbundance.Taxa.tsv',
        mode_flag = path + 'results/vc-sample/05_merge_tax_fun_depth/mode/{sample}.mode.tsv',
    params:
        taxid_tax = taxid_tax_path,
        taxid_fun = taxid_fun_path,
        contigs_taxonomy_classify_method = contigs_taxonomy_classify_method,
        contigs_function_classify_method = contigs_function_classify_method,
        contigs_tax_detail = path + 'results/c-sample/01_taxonomy/' + contigs_taxonomy_classify_method + '/{sample}.' + contigs_taxonomy_classify_method + '.contig2classification.tsv',
        contigs_fun_detail = path + 'results/c-sample/02_function/' + contigs_function_classify_method + '/annotations/{sample}.' + contigs_function_classify_method + '.annotations.gz',
        vMAGs_tax_detail = path + 'results/v-sample/10_summary/{sample}/vMAG_summary_info.tsv',
        vMAGs_tax_detail_filter = path + 'results/vc-sample/05_merge_tax_fun_depth/{sample}.vMAG_summary_info.filter.tsv',
        vMAGs_ffn_length = path + 'results/v-sample/10_summary/{sample}/Virus_annotation_results.tsv',
        vMAGs_contig_length = path + 'results/v-sample/10_summary/{sample}/final_virus.contig_length.txt',
        merge_detail = path + 'results/vc-sample/05_merge_tax_fun_depth/detail_for_each_seq/{sample}.detail_for_each_seq.tsv',
        minORFLen = config["Params_vc_sample"]["05_merge_tax_fun_depth"]["minORFLen"],
        checkv_quality = config["Params_vc_sample"]["05_merge_tax_fun_depth"]["checkv_quality"],
        KO_mode = config["Params_vc_sample"]["05_merge_tax_fun_depth"]["KO_mode"],
        min_KO_num_in_one_taxon = config["Params_vc_sample"]["05_merge_tax_fun_depth"]["min_KO_num_in_one_taxon"],
        rank = ",".join(RANKS),
        vc_consistency = path + 'results/vc-sample/05_merge_tax_fun_depth/MAG_contig_consistency/{sample}.vc_consistency.tsv',
        out_CAMI_format = path + 'results/vc-sample/05_merge_tax_fun_depth/Taxa_abundance/{sample}.Taxa_abundance.CAMI_format.tsv',
        CAMI_format_output = config["Params_vc_sample"]["05_merge_tax_fun_depth"]["CAMI_format_output"],
        output_detail = config["Params_vc_sample"]["05_merge_tax_fun_depth"]["output_detail"],
    threads: config["Params_vc_sample"]["05_merge_tax_fun_depth"]["threads"]
    conda: config["conda_envs"]["metalink-main"]
    log:
        path + 'logs/vc-sample/05_merge_tax_fun_depth/{sample}.log'
    benchmark:
        path + 'benchmarks/vc-sample/05_merge_tax_fun_depth/{sample}.benchmark'
    shell:
        """
        echo "Start vc_sample_05_merge_tax_fun_depth for {wildcards.sample}" > {log}
        mkdir -p $(dirname {params.vc_consistency}) 2>{log}

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

        if [ -n "$vMAGs_miss" ] && [ -n "$contigs_miss" ]; then
            echo "ERROR: both contigs and vMAGs inputs missing/empty for {wildcards.sample}; cannot merge taxonomy/function/depth." >> {log}
            echo "  missing/empty contigs inputs:$contigs_miss" >> {log}
            echo "  missing/empty vMAGs inputs:$vMAGs_miss" >> {log}
            echo "ERROR [vc_sample_05_merge_tax_fun_depth] sample {wildcards.sample}: both contigs and vMAGs inputs missing/empty; check upstream results. See log: {log}" >&2
            mkdir -p $(dirname {output.mode_flag})
            exit 1
        elif [ -n "$vMAGs_miss" ] && [ -z "$contigs_miss" ]; then
            echo "No vMAGs, but contigs present; running contigs-only merge." >> {log}
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

        elif [ -z "$vMAGs_miss" ] && [ -z "$contigs_miss" ]; then
            echo "vMAGs and contigs present; running vMAGs-contigs merge." >> {log}
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

        else
            echo "Error: Detect vMAGs but no available contigs for {wildcards.sample}" >> {log}
            exit 1
        fi
        """

rule vc_sample_05_visualize:
    priority: 40
    input:
        summary_Tax = path + 'results/vc-sample/05_merge_tax_fun_depth/summary/{sample}.summary_Tax.tsv',
        summary_ORF_KO = path + 'results/vc-sample/05_merge_tax_fun_depth/summary/{sample}.summary_ORF_KO.tsv',
        summary_ORF_Tax_KO = path + 'results/vc-sample/05_merge_tax_fun_depth/summary/{sample}.summary_ORF_Tax_KO.tsv',
        mode_flag = path + 'results/vc-sample/05_merge_tax_fun_depth/mode/{sample}.mode.tsv',
    output:
        mk = path + 'markers/vc-sample/05_visualize/{sample}.mk',
    params:
        prefix = path + 'results/vc-sample/05_visualization/{sample}',
    threads: 1
    conda: config["conda_envs"]["metalink-main"]
    log:
        path + 'logs/vc-sample/05_visualize/{sample}.log'
    benchmark:
        path + 'benchmarks/vc-sample/05_visualize/{sample}.benchmark'
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

rule vc_sample_06_module_completeness:
    priority: 40
    input:
        KO_Taxa_RelativeAbundance = path + 'results/vc-sample/05_merge_tax_fun_depth/KO_Taxa_RelativeAbundance/{sample}.KO_Taxa_RelativeAbundance.Taxa.tsv',
    output:
        mk = path + 'markers/vc-sample/06_module_completeness/{sample}.{rank}.mk',
    params:
        KO_Taxa_RelativeAbundance_ranked = path + 'results/vc-sample/05_merge_tax_fun_depth/KO_Taxa_RelativeAbundance/{sample}.KO_Taxa_RelativeAbundance.Taxa.{rank}.tsv',
        KO_Taxa_SumCoverage_ranked = path + 'results/vc-sample/05_merge_tax_fun_depth/KO_Taxa_SumCoverage/{sample}.KO_Taxa_SumCoverage.Taxa.{rank}.tsv',
        KO_Taxa_SumCoverage = path + 'results/vc-sample/05_merge_tax_fun_depth/KO_Taxa_SumCoverage/{sample}.KO_Taxa_SumCoverage.Taxa.tsv',
        Module_Taxa_RA = path + 'results/vc-sample/06_module_completeness/{sample}/{rank}/',
        pathway_split = "TRUE",
        min_KO_num_in_one_taxon = config["Params_vc_sample"]["06_mclink_taxon"]["min_KO_num_in_one_taxon"],
        min_mean_KO_RA_in_one_taxon = config["Params_vc_sample"]["06_mclink_taxon"]["min_mean_KO_RA_in_one_taxon"],
        completeness_threshold = config["Params_vc_sample"]["06_mclink_taxon"]["completeness_threshold"],
        plus_scale_method = config["Params_vc_sample"]["06_mclink_taxon"]["plus_scale_method"],
        comma_scale_method = config["Params_vc_sample"]["06_mclink_taxon"]["comma_scale_method"],
        TopN = config["Params_vc_sample"]["06_mclink_taxon"]["TopN"],
    threads: config["Params_vc_sample"]["05_merge_tax_fun_depth"]["threads"]
    conda: config["conda_envs"]["metalink-main"]
    log:
        path + 'logs/vc-sample/06_module_completeness/{sample}.{rank}.log'
    benchmark:
        path + 'benchmarks/vc-sample/06_module_completeness/{sample}.{rank}.benchmark'
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
            echo "KO_Taxa_RelativeAbundance empty/missing for {wildcards.sample}; skip vc_sample_03_module_completeness." >>{log}
        fi
        touch {output.mk}
        """

rule vc_project_01_cat_Taxa_abundance_dcast:
    input:
        Taxa_abundance = expand(path + 'results/vc-sample/05_merge_tax_fun_depth/Taxa_abundance/{sample}.Taxa_abundance.tsv', sample=SAMPLE_INDEX),
    output:
        out_long = path + 'results/vc-project/01_Taxa_Sample/{project}/Taxa_Sample.long.tsv',
    params:
        tmp_dir = path + 'results/vc-project/tmp/{project}/',
        prefix_wide = path + 'results/vc-project/01_Taxa_Sample/{project}/Taxa_Sample.wide.tsv',
    conda: config["conda_envs"]["metalink-main"]
    threads: config["Params_vc_sample"]["05_merge_tax_fun_depth"]["threads"]
    log: path + 'logs/vc-project/01_cat_Taxa_abundance_dcast/{project}.log'
    benchmark:
        path + 'benchmarks/vc-project/01_cat_Taxa_abundance_dcast/{project}.benchmark'
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

rule vc_project_01_Taxa_visualize:
    input:
        out_long = path + 'results/vc-project/01_Taxa_Sample/{project}/Taxa_Sample.long.tsv',
    output:
        mk = path + 'markers/vc-project/01_Taxa_visualization/{project}.mk',
    params:
        prefix = path + 'results/vc-project/01_Taxa_visualization/{project}',
        prefix_wide = path + 'results/vc-project/01_Taxa_Sample/{project}/Taxa_Sample.Taxa_ReadCount.wide.tsv',
        TopN = config["Params_vc_project"]["01_Taxa_visualize"]["TopN"],
    threads: 1
    conda: config["conda_envs"]["metalink-main"]
    log:
        path + 'logs/vc-project/01_Taxa_visualization/{project}.log'
    benchmark:
        path + 'benchmarks/vc-project/01_Taxa_visualization/{project}.benchmark'
    shell:
        """
        mkdir -p $(dirname {params.prefix}) 2>>{log}
        missing=""
        for f in {input.out_long} {params.prefix_wide}; do
            [ -s "$f" ] || missing="$missing $f"
        done
        if [ -n "$missing" ]; then
            echo "Missing/empty inputs:$missing; skip vc project_01_Taxa_visualize for {wildcards.project}" >>{log}
        else
            echo "All input files non-empty, running vc project_01_Taxa_visualize." >>{log}
            Rscript {Dir_script_metalink}/plot_bar.cmd.R \
                --Type Taxa \
                --TopN {params.TopN} \
                --Input {params.prefix_wide} \
                --Output {params.prefix} >>{log} 2>&1
        fi
        touch {output.mk}
        """

rule vc_project_02_cat_KO_depth_dcast:
    input:
        KO_depth = expand(path + 'results/vc-sample/05_merge_tax_fun_depth/KO_depth/{sample}.KO_depth.tsv', sample=SAMPLE_INDEX),
    output:
        out_long = path + 'results/vc-project/02_KO_Sample/{project}/KO_Sample.long.tsv',
    params:
        tmp_dir = path + 'results/vc-project/tmp/{project}/',
        prefix_wide = path + 'results/vc-project/02_KO_Sample/{project}/KO_Sample.wide.tsv',
    threads: config["Params_vc_sample"]["05_merge_tax_fun_depth"]["threads"]
    conda: config["conda_envs"]["metalink-main"]
    log:
        path + 'logs/vc-project/02_cat_KO_depth_dcast/{project}.log'
    benchmark:
        path + 'benchmarks/vc-project/02_cat_KO_depth_dcast/{project}.benchmark'
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

rule vc_project_03_community_module:
    input:
        out_long = path + 'results/vc-project/02_KO_Sample/{project}/KO_Sample.long.tsv',
    output:
        mk = path + 'markers/vc-project/03_community_module/{project}.mk',
    params:
        wide_RelativeAbundance = path + 'results/vc-project/02_KO_Sample/{project}/KO_Sample.KO_RelativeAbundance.wide.tsv',
        wide_SumCoverage = path + 'results/vc-project/02_KO_Sample/{project}/KO_Sample.KO_SumCoverage.wide.tsv',
        pathway_split = "TRUE",
        community_output = path + 'results/vc-project/03_Module_Sample/{project}/',
        min_sample_num_per_KO = config["Params_vc_project"]["03_mclink_sample"]["min_sample_num_per_KO"],
        preva = config["Params_vc_project"]["03_mclink_sample"]["preva"],
        abund = config["Params_vc_project"]["03_mclink_sample"]["abund"],
        completeness_threshold = config["Params_vc_project"]["03_mclink_sample"]["completeness_threshold"],
        plus_scale_method = config["Params_vc_project"]["03_mclink_sample"]["plus_scale_method"],
        comma_scale_method = config["Params_vc_project"]["03_mclink_sample"]["comma_scale_method"],
    threads: config["Params_vc_sample"]["05_merge_tax_fun_depth"]["threads"]
    conda: config["conda_envs"]["metalink-main"]
    log:
        path + 'logs/vc-project/03_community_module/{project}.log'
    benchmark:
        path + 'benchmarks/vc-project/03_community_module/{project}.benchmark'
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

rule vc_project_03_Fun_visualize:
    input:
        mk = path + 'markers/vc-project/03_community_module/{project}.mk',
    output:
        mk = path + 'markers/vc-project/03_Fun_visualization/{project}.mk',
    params:
        prefix = path + 'results/vc-project/03_Fun_visualization/{project}',
        prefix_wide = path + 'results/vc-project/03_Module_Sample/{project}/SumCoverage/mc_abundance_filtered_by_completeness.tsv',
        TopN = config["Params_vc_project"]["03_Fun_visualize"]["TopN"],
    threads: 1
    conda: config["conda_envs"]["metalink-main"]
    log:
        path + 'logs/vc-project/03_Fun_visualization/{project}.log'
    benchmark:
        path + 'benchmarks/vc-project/03_Fun_visualization/{project}.benchmark'
    shell:
        """
        mkdir -p $(dirname {params.prefix}) 2>>{log}
        missing=""
        for f in {params.prefix_wide}; do
            [ -s "$f" ] || missing="$missing $f"
        done
        if [ -n "$missing" ]; then
            echo "Missing/empty inputs:$missing; skip vc project_03_Fun_visualize for {wildcards.project}" >>{log}
        else
            echo "All input files non-empty, running vc project_03_Fun_visualize." >>{log}
            Rscript {Dir_script_metalink}/plot_bar.cmd.R \
                --Type Module \
                --TopN {params.TopN} \
                --Input {params.prefix_wide} \
                --Output {params.prefix} >>{log} 2>&1
        fi
        touch {output.mk}
        """

rule vc_project_05_merge_tax_KO_depth:
    input:
        ko_depth = expand(path + 'results/vc-sample/05_merge_tax_fun_depth/KO_depth/{sample}.KO_depth.tsv', sample=SAMPLE_INDEX),
    output:
        out_long = path + 'results/vc-project/05_merge_tax_KO_depth/{project}/{rank}.mean.long.tsv',
    params:
        detail = lambda wc: expand(path + 'results/vc-sample/05_merge_tax_fun_depth/detail_for_each_seq/{sample}.detail_for_each_{rank}.tsv', sample=SAMPLE_INDEX, rank=wc.rank),  # 只展开当前 rank
        tmp_dir = path + 'results/vc-project/tmp/{project}/',
        prefix_wide = path + 'results/vc-project/05_merge_tax_KO_depth/{project}/{rank}.wide.tsv',
        preva = config["Params_vc_project"]["05_merge_tax_KO_depth"]["preva"],
        abund = config["Params_vc_project"]["05_merge_tax_KO_depth"]["abund"],
    threads: config["Params_vc_sample"]["05_merge_tax_fun_depth"]["threads"]
    conda: config["conda_envs"]["metalink-main"]
    log:
        path + 'logs/vc-project/05_merge_tax_KO_depth/{project}/{rank}.log'
    benchmark:
        path + 'benchmarks/vc-project/05_merge_tax_KO_depth/{project}/{rank}.benchmark'
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

rule vc_project_05_merge_tax_Module_RelativeAbundance:
    input:
        mk = lambda wc: expand(path + 'markers/vc-sample/06_module_completeness/{sample}.{rank}.mk', sample=SAMPLE_INDEX, rank=wc.rank),
    output:
        out_long = path + 'results/vc-project/05_merge_tax_Module_RelativeAbundance/{project}/{rank}.mean.long.tsv',
    params:
        detail = lambda wc: expand(path + 'results/vc-sample/06_module_completeness/{sample}/{rank}/RelativeAbundance/mc_abundance.tsv', sample=SAMPLE_INDEX, rank=wc.rank),
        prefix_wide = path + 'results/vc-project/05_merge_tax_Module_RelativeAbundance/{project}/{rank}.wide.tsv',
        tmp_dir = path + 'results/vc-project/tmp/{project}/',
        preva = config["Params_vc_project"]["05_merge_tax_Module_RelativeAbundance"]["preva"],
        abund = config["Params_vc_project"]["05_merge_tax_Module_RelativeAbundance"]["abund"],
    threads: config["Params_vc_sample"]["05_merge_tax_fun_depth"]["threads"]
    conda: config["conda_envs"]["metalink-main"]
    log:
        path + 'logs/vc-project/05_merge_tax_Module_RelativeAbundance/{project}/{rank}.log'
    benchmark:
        path + 'benchmarks/vc-project/05_merge_tax_Module_RelativeAbundance/{project}/{rank}.benchmark'
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

rule vc_project_05_merge_tax_Module_SumCoverage:
    input:
        mk = lambda wc: expand(path + 'markers/vc-sample/06_module_completeness/{sample}.{rank}.mk', sample=SAMPLE_INDEX, rank=wc.rank),
    output:
        out_long = path + 'results/vc-project/05_merge_tax_Module_SumCoverage/{project}/{rank}.mean.long.tsv',
    params:
        detail = lambda wc: expand(path + 'results/vc-sample/06_module_completeness/{sample}/{rank}/SumCoverage/mc_abundance.tsv', sample=SAMPLE_INDEX, rank=wc.rank),
        prefix_wide = path + 'results/vc-project/05_merge_tax_Module_SumCoverage/{project}/{rank}.wide.tsv',
        tmp_dir = path + 'results/vc-project/tmp/{project}/',
        preva = config["Params_vc_project"]["05_merge_tax_Module_SumCoverage"]["preva"],
        abund = config["Params_vc_project"]["05_merge_tax_Module_SumCoverage"]["abund"],
    threads: config["Params_vc_sample"]["05_merge_tax_fun_depth"]["threads"]
    conda: config["conda_envs"]["metalink-main"]
    log:
        path + 'logs/vc-project/05_merge_tax_Module_SumCoverage/{project}/{rank}.log'
    benchmark:
        path + 'benchmarks/vc-project/05_merge_tax_Module_SumCoverage/{project}/{rank}.benchmark'
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

rule vc_project_06_module_completeness:
    input:
        KO_Taxa_RelativeAbundance = path + 'results/vc-project/05_merge_tax_KO_depth/{project}/{rank}.mean.long.tsv',
    output:
        mk = path + 'markers/vc-project/06_module_completeness/{project}/{rank}.mk',
    params:
        KO_Taxa_RelativeAbundance_ranked = path + 'results/vc-project/05_merge_tax_KO_depth/{project}/{rank}.Taxa_KO_RelativeAbundance.mean.wide.tsv',
        KO_Taxa_SumCoverage_ranked = path + 'results/vc-project/05_merge_tax_KO_depth/{project}/{rank}.Taxa_KO_SumCoverage.mean.wide.tsv',
        Module_Taxa_RA = path + 'results/vc-project/06_module_completeness/{project}/{rank}/',
        pathway_split = "TRUE",
        min_KO_num_in_one_taxon = config["Params_vc_project"]["06_mclink_taxon"]["min_KO_num_in_one_taxon"],
        min_mean_KO_RA_in_one_taxon = config["Params_vc_project"]["06_mclink_taxon"]["min_mean_KO_RA_in_one_taxon"],
        completeness_threshold = config["Params_vc_project"]["06_mclink_taxon"]["completeness_threshold"],
        plus_scale_method = config["Params_vc_project"]["06_mclink_taxon"]["plus_scale_method"],
        comma_scale_method = config["Params_vc_project"]["06_mclink_taxon"]["comma_scale_method"],
        TopN = config["Params_vc_project"]["06_mclink_taxon"]["TopN"],
    threads: config["Params_vc_sample"]["05_merge_tax_fun_depth"]["threads"]
    conda: config["conda_envs"]["metalink-main"]
    log:
        path + 'logs/vc-project/06_module_completeness/{project}/{rank}.log'
    benchmark:
        path + 'benchmarks/vc-project/06_module_completeness/{project}/{rank}.benchmark'
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
            echo "KO_Taxa_RelativeAbundance empty/missing for {wildcards.project}; skip vc_project_module_completeness." >>{log}
        fi
        touch {output.mk}
        """
