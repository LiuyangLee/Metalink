# control pipline endpoint
# should be one of the "taxonomy" | "function" | "taxonomy_and_function" | "all"
TARGET_STAGE = config.get("Params_c_sample", {}).get("00_target_stage", "all")

rule c_end_target:
    input:
        p_mk_11 =  expand(path + 'results/b-project/01_cat_raw_stats/{project}/01_raw_reads_stats.tsv', project=PROJECT_INDEX) if RUN_PROJECT_LEVEL else [],
        p_mk_12 =  expand(path + 'results/b-project/02_cat_trimmed_stats/{project}/02_clean_reads_stats.tsv', project=PROJECT_INDEX) if RUN_PROJECT_LEVEL else [],
        p_mk_13 =  expand(path + 'markers/b-project/03_contigs_stats/{project}/03_contigs_stats.mk', project=PROJECT_INDEX) if RUN_PROJECT_LEVEL else [],
        p_mk_14 =  expand(path + 'results/b-project/04_merge_reads_contigs_stats/{project}/04_reads_contigs_ffn_stats.tsv', project=PROJECT_INDEX) if RUN_PROJECT_LEVEL else [],
        p_bam_summary = expand(path + 'results/b-project/04_bam_summary/{project}.bam_mapping_summary.tsv', project=PROJECT_INDEX) if RUN_PROJECT_LEVEL else [],
        s_mk_02 = expand(path + 'markers/c-sample/01_taxonomy/classify_done/{sample}.mk', sample=SAMPLE_INDEX) if TARGET_STAGE in ["taxonomy", "taxonomy_and_function", "all"] else [],
        s_eg_03 = expand(path + 'markers/c-sample/02_function/classify_done/{sample}.mk', sample=SAMPLE_INDEX) if TARGET_STAGE in ["function", "taxonomy_and_function", "all"] else [],
        s_tax_03 = expand(path + 'results/c-sample/03_merge_tax_depth/Taxa_abundance/{sample}.Taxa_abundance.tsv', sample=SAMPLE_INDEX) if TARGET_STAGE in ["taxonomy", "taxonomy_and_function"] else [],
        s_fun_04 = expand(path + 'results/c-sample/04_merge_fun_depth/KO_depth/{sample}.KO_depth.tsv', sample=SAMPLE_INDEX) if TARGET_STAGE in ["function", "taxonomy_and_function"] else [],
        s_merge_05 = expand(path + 'results/c-sample/05_merge_tax_fun_depth/KO_Taxa_RelativeAbundance/{sample}.KO_Taxa_RelativeAbundance.Taxa.tsv', sample=SAMPLE_INDEX) if TARGET_STAGE == "all" else [],
        s_mc_05 = expand(path + 'markers/c-sample/05_visualize/{sample}.mk', sample=SAMPLE_INDEX) if TARGET_STAGE == "all" else [],
        s_mc_06 = expand(path + 'markers/c-sample/06_module_completeness/{sample}.{rank}.mk',sample=SAMPLE_INDEX, rank=RANKS) if TARGET_STAGE == "all" else [],
        s_Ta_out_long = expand(path + 'results/c-project/01_Taxa_Sample/{project}/Taxa_Sample.long.tsv', project=PROJECT_INDEX) if RUN_PROJECT_LEVEL and TARGET_STAGE in ["taxonomy", "taxonomy_and_function", "all"] else [],
        p_KO_out_long = expand(path + 'results/c-project/02_KO_Sample/{project}/KO_Sample.long.tsv', project=PROJECT_INDEX) if RUN_PROJECT_LEVEL and TARGET_STAGE in ["function", "taxonomy_and_function", "all"] else [],
        s_Ta_plot_mk = expand(path + 'markers/c-project/01_Taxa_visualization/{project}.mk', project=PROJECT_INDEX) if RUN_PROJECT_LEVEL and TARGET_STAGE in ["taxonomy", "taxonomy_and_function", "all"] else [],
        s_KO_plot_mk = expand(path + 'markers/c-project/03_Fun_visualization/{project}.mk', project=PROJECT_INDEX) if RUN_PROJECT_LEVEL and TARGET_STAGE == "all" else [],
        p_mc_mk = expand(path + 'markers/c-project/03_community_module/{project}.mk',project=PROJECT_INDEX) if RUN_PROJECT_LEVEL and TARGET_STAGE == "all" else [],
        p_mc_mk3 = expand(path + 'markers/c-project/06_module_completeness/{project}/{rank}.mk',project=PROJECT_INDEX, rank=RANKS) if RUN_PROJECT_LEVEL and TARGET_STAGE == "all" else [],
        p_out_long_mc = expand(path + 'results/c-project/05_merge_tax_KO_depth/{project}/{rank}.mean.long.tsv',rank=RANKS, project=PROJECT_INDEX) if RUN_PROJECT_LEVEL and TARGET_STAGE == "all" else [],
        p_out_long_Module_RelativeAbundance = expand(path + 'results/c-project/05_merge_tax_Module_RelativeAbundance/{project}/{rank}.mean.long.tsv',rank=RANKS, project=PROJECT_INDEX) if RUN_PROJECT_LEVEL and TARGET_STAGE == "all" else [],
        p_out_long_Module_SumCoverage = expand(path + 'results/c-project/05_merge_tax_Module_SumCoverage/{project}/{rank}.mean.long.tsv',rank=RANKS, project=PROJECT_INDEX) if RUN_PROJECT_LEVEL and TARGET_STAGE == "all" else [],
        s_b_mk_03 = expand(path + 'markers/b-sample/03_delete_assembly_intermediate/{sample}.mk', sample=SAMPLE_INDEX) if config.get('skip_level', 0) < 2 and not keep_intermediate('assembly') else [],
