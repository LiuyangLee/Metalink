# control pipline endpoint
# should be one of the "taxonomy" | "function" | "taxonomy_and_function" | "all"
TARGET_STAGE = config.get("Params_r_sample", {}).get("00_target_stage", "all")

rule r_end_target:
    input:
        p_mk_11 =  expand(path + 'results/b-project/01_cat_raw_stats/{project}/01_raw_reads_stats.tsv', project=PROJECT_INDEX) if RUN_PROJECT_LEVEL else [],
        p_mk_12 =  expand(path + 'results/b-project/02_cat_trimmed_stats/{project}/02_clean_reads_stats.tsv', project=PROJECT_INDEX) if RUN_PROJECT_LEVEL else [],
        s_trim_merge_1 = expand(path + 'results/r-sample/02_FragGeneScanRs/{sample}_both_merge_unmerge.fna.gz',sample=SAMPLE_INDEX) if TARGET_STAGE in ["function", "taxonomy_and_function", "all"] else [],
        s_trim_merge_2 = expand(path + 'results/r-sample/02_FragGeneScanRs/{sample}_both_merge_unmerge.ffn.gz',sample=SAMPLE_INDEX) if TARGET_STAGE in ["function", "taxonomy_and_function", "all"] else [],
        s_ffn_length = expand(path + 'results/r-sample/02_FragGeneScanRs/{sample}_both_merge_unmerge.ffn_length.gz',sample=SAMPLE_INDEX) if TARGET_STAGE in ["function", "taxonomy_and_function", "all"] else [],
        s_ffn_len  = expand(path + 'results/r-sample/02_FragGeneScanRs_stats/{sample}_both_merge_unmerge.fna.stats',sample=SAMPLE_INDEX) if TARGET_STAGE in ["function", "taxonomy_and_function", "all"] else [],
        s_bracken_kraken = expand(path + 'results/r-sample/03_taxonomy/bracken/{sample}.bracken_kraken.report.tsv',sample=SAMPLE_INDEX) if TARGET_STAGE in ["taxonomy", "function", "taxonomy_and_function", "all"] else [],
        s_reads_fun_mk = expand(path + 'markers/r-sample/04_function/classify_done/{sample}.mk',sample=SAMPLE_INDEX) if TARGET_STAGE in ["function", "taxonomy_and_function", "all"] else [],
        s_merge_05 = expand(path + 'results/r-sample/05_merge_tax_fun_depth/KO_Taxa_RelativeAbundance/{sample}.KO_Taxa_RelativeAbundance.Taxa.tsv', sample=SAMPLE_INDEX) if TARGET_STAGE == "all" else [],
        s_mc_05 = expand(path + 'markers/r-sample/05_visualize/{sample}.mk', sample=SAMPLE_INDEX) if TARGET_STAGE == "all" else [],
        s_mc_06 = expand(path + 'markers/r-sample/06_module_completeness/{sample}.{rank}.mk',rank=RANKS, sample=SAMPLE_INDEX) if TARGET_STAGE == "all" else [],
        merged_stats = expand(path + 'results/r-project/00_FragGeneScanRs_stats/{project}/02_FGS_stats.tsv',project=PROJECT_INDEX) if RUN_PROJECT_LEVEL and TARGET_STAGE in ["function", "taxonomy_and_function", "all"] else [],
        p_bracken_taxonomic_mk = expand(path + 'results/r-project/01_Taxa_Sample/{project}/bracken_taxonomic_profiles.count.tsv',project=PROJECT_INDEX) if RUN_PROJECT_LEVEL and TARGET_STAGE in ["taxonomy", "taxonomy_and_function", "all"] else [],
        p_kraken2_taxonomic_mk = expand(path + 'results/r-project/01_Taxa_Sample/{project}/kraken2_taxonomic_profiles.count.tsv',project=PROJECT_INDEX) if RUN_PROJECT_LEVEL and TARGET_STAGE in ["taxonomy", "taxonomy_and_function", "all"] else [],
        p_kraken2_functional_mk = expand(path + 'results/r-project/02_KO_Sample/{project}/kraken2_functional_profiles.count.tsv',project=PROJECT_INDEX) if RUN_PROJECT_LEVEL and TARGET_STAGE in ["function", "taxonomy_and_function", "all"] else [],
        p_Ta_out_long = expand(path + 'results/r-project/01_Taxa_Sample/{project}/Taxa_Sample.long.tsv',project=PROJECT_INDEX) if RUN_PROJECT_LEVEL and TARGET_STAGE == "all" else [],
        p_KO_out_long = expand(path + 'results/r-project/02_KO_Sample/{project}/KO_Sample.long.tsv',project=PROJECT_INDEX) if RUN_PROJECT_LEVEL and TARGET_STAGE == "all" else [],
        s_Ta_plot_mk = expand(path + 'markers/r-project/01_Taxa_visualization/{project}.mk', project=PROJECT_INDEX) if RUN_PROJECT_LEVEL and TARGET_STAGE == "all" else [],
        s_KO_plot_mk = expand(path + 'markers/r-project/03_Fun_visualization/{project}.mk', project=PROJECT_INDEX) if RUN_PROJECT_LEVEL and TARGET_STAGE == "all" else [],
        p_mc_mk = expand(path + 'markers/r-project/03_community_module/{project}.mk',project=PROJECT_INDEX) if RUN_PROJECT_LEVEL and TARGET_STAGE == "all" else [],
        p_mc_mk3 = expand(path + 'markers/r-project/06_module_completeness/{project}/{rank}.mk',project=PROJECT_INDEX, rank=RANKS) if RUN_PROJECT_LEVEL and TARGET_STAGE == "all" else [],
        p_out_long_mc = expand(path + 'results/r-project/05_merge_tax_KO_depth/{project}/{rank}.mean.long.tsv',rank=RANKS, project=PROJECT_INDEX) if RUN_PROJECT_LEVEL and TARGET_STAGE == "all" else [],
        p_out_long_Module_RelativeAbundance = expand(path + 'results/r-project/05_merge_tax_Module_RelativeAbundance/{project}/{rank}.mean.long.tsv',rank=RANKS, project=PROJECT_INDEX) if RUN_PROJECT_LEVEL and TARGET_STAGE == "all" else [],
        p_out_long_Module_SumCoverage = expand(path + 'results/r-project/05_merge_tax_Module_SumCoverage/{project}/{rank}.mean.long.tsv',rank=RANKS, project=PROJECT_INDEX) if RUN_PROJECT_LEVEL and TARGET_STAGE == "all" else [],
