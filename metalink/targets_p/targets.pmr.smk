# pmr targets are assembled conditionally on config['mags'] ('m' | 'v' | 'vm'):
#   m  : p+mMAG+r   -> b(binning) + m-s/m-p + pmr_m chain
#   v  : p+vMAG+r   -> b(no binning) + v-s + v-p + pmr_v chain
#   vm : p+vMAG+mMAG+r -> b(binning) + v-s + v-p + m-s/m-p + pmr_vm chain
_MAGS = config.get('mags', 'm')

rule pmr_end_target:
    input:
        # ---- b-project stats (all mags values) ----
        p_mk_11 =  expand(path + 'results/b-project/01_cat_raw_stats/{project}/01_raw_reads_stats.tsv', project=PROJECT_INDEX) if RUN_PROJECT_LEVEL else [],
        p_mk_12 =  expand(path + 'results/b-project/02_cat_trimmed_stats/{project}/02_clean_reads_stats.tsv', project=PROJECT_INDEX) if RUN_PROJECT_LEVEL else [],
        p_mk_13 =  expand(path + 'markers/b-project/03_contigs_stats/{project}/03_contigs_stats.mk', project=PROJECT_INDEX) if RUN_PROJECT_LEVEL else [],
        p_mk_14 =  expand(path + 'results/b-project/04_merge_reads_contigs_stats/{project}/04_reads_contigs_ffn_stats.tsv', project=PROJECT_INDEX) if RUN_PROJECT_LEVEL else [],
        # binning-derived b-project stats: only m/vm run the binning chain
        p_mk_15 =  expand(path + 'results/b-project/05_gather_project_bin_checkm_stats/{project}/05_checkm_high_quality_bins.csv', project=PROJECT_INDEX) if RUN_PROJECT_LEVEL and _MAGS in ('m', 'vm') else [],
        # ---- m chain (--mags m): pmr_m-sample / pmr_m-project products ----
        p_species_table = expand(path + 'results/pmr_m-project/04_summary/{project}/mMAG_species_abundance_matrix.tsv', project=PROJECT_INDEX) if RUN_PROJECT_LEVEL and _MAGS == 'm' else [],
        p_03_coverm_contigs = expand(path + 'results/pmr_m-sample/03_coverm_contigs/{project}/{sample}.tsv', project=PROJECT_INDEX, sample=SAMPLE_INDEX) if _MAGS == 'm' else [],
        s_merge_05 = expand(path + 'results/pmr_m-sample/05_merge_tax_fun_depth/{project}/KO_Taxa_RelativeAbundance/{sample}.KO_Taxa_RelativeAbundance.Taxa.tsv', project=PROJECT_INDEX, sample=SAMPLE_INDEX) if _MAGS == 'm' else [],
        s_mc_05 = expand(path + 'markers/pmr_m-sample/05_visualize/{project}/{sample}.mk', project=PROJECT_INDEX, sample=SAMPLE_INDEX) if _MAGS == 'm' else [],
        s_mc_06 = expand(path + 'markers/pmr_m-sample/06_module_completeness/{project}/{sample}.{rank}.mk',project=PROJECT_INDEX, rank=RANKS, sample=SAMPLE_INDEX) if _MAGS == 'm' else [],
        p_Ta_out_long = expand(path + 'results/pmr_m-project/01_Taxa_Sample/{project}/Taxa_Sample.long.tsv',project=PROJECT_INDEX) if RUN_PROJECT_LEVEL and _MAGS == 'm' else [],
        p_KO_out_long = expand(path + 'results/pmr_m-project/02_KO_Sample/{project}/KO_Sample.long.tsv',project=PROJECT_INDEX) if RUN_PROJECT_LEVEL and _MAGS == 'm' else [],
        s_Ta_plot_mk = expand(path + 'markers/pmr_m-project/01_Taxa_visualization/{project}.mk', project=PROJECT_INDEX) if RUN_PROJECT_LEVEL and _MAGS == 'm' else [],
        s_KO_plot_mk = expand(path + 'markers/pmr_m-project/03_Fun_visualization/{project}.mk', project=PROJECT_INDEX) if RUN_PROJECT_LEVEL and _MAGS == 'm' else [],
        p_mc_mk = expand(path + 'markers/pmr_m-project/03_community_module/{project}.mk',project=PROJECT_INDEX) if RUN_PROJECT_LEVEL and _MAGS == 'm' else [],
        p_mc_mk3 = expand(path + 'markers/pmr_m-project/06_module_completeness/{project}/{rank}.mk',project=PROJECT_INDEX, rank=RANKS) if RUN_PROJECT_LEVEL and _MAGS == 'm' else [],
        p_out_long_mc = expand(path + 'results/pmr_m-project/05_merge_tax_KO_depth/{project}/{rank}.mean.long.tsv',rank=RANKS, project=PROJECT_INDEX) if RUN_PROJECT_LEVEL and _MAGS == 'm' else [],
        p_out_long_Module_RelativeAbundance = expand(path + 'results/pmr_m-project/05_merge_tax_Module_RelativeAbundance/{project}/{rank}.mean.long.tsv',rank=RANKS, project=PROJECT_INDEX) if RUN_PROJECT_LEVEL and _MAGS == 'm' else [],
        p_out_long_Module_SumCoverage = expand(path + 'results/pmr_m-project/05_merge_tax_Module_SumCoverage/{project}/{rank}.mean.long.tsv',rank=RANKS, project=PROJECT_INDEX) if RUN_PROJECT_LEVEL and _MAGS == 'm' else [],
        # ---- v chain (--mags v): pmr_v-sample / pmr_v-project products ----
        v_cm_03 = expand(path + 'results/pmr_v-sample/03_coverm_MAGs/{project}/{sample}.tsv', project=PROJECT_INDEX, sample=SAMPLE_INDEX) if _MAGS == 'v' else [],
        v_03_coverm_contigs = expand(path + 'results/pmr_v-sample/03_coverm_contigs/{project}/{sample}.tsv', project=PROJECT_INDEX, sample=SAMPLE_INDEX) if _MAGS == 'v' else [],
        v_merge_05 = expand(path + 'results/pmr_v-sample/05_merge_tax_fun_depth/{project}/KO_Taxa_RelativeAbundance/{sample}.KO_Taxa_RelativeAbundance.Taxa.tsv', project=PROJECT_INDEX, sample=SAMPLE_INDEX) if _MAGS == 'v' else [],
        v_mc_05 = expand(path + 'markers/pmr_v-sample/05_visualize/{project}/{sample}.mk', project=PROJECT_INDEX, sample=SAMPLE_INDEX) if _MAGS == 'v' else [],
        v_mc_06 = expand(path + 'markers/pmr_v-sample/06_module_completeness/{project}/{sample}.{rank}.mk',project=PROJECT_INDEX, rank=RANKS, sample=SAMPLE_INDEX) if _MAGS == 'v' else [],
        v_species_table = expand(path + 'results/pmr_v-project/04_summary/{project}/vMAG_species_abundance_matrix.tsv', project=PROJECT_INDEX) if RUN_PROJECT_LEVEL and _MAGS == 'v' else [],
        v_Ta_out_long = expand(path + 'results/pmr_v-project/01_Taxa_Sample/{project}/Taxa_Sample.long.tsv',project=PROJECT_INDEX) if RUN_PROJECT_LEVEL and _MAGS == 'v' else [],
        v_KO_out_long = expand(path + 'results/pmr_v-project/02_KO_Sample/{project}/KO_Sample.long.tsv',project=PROJECT_INDEX) if RUN_PROJECT_LEVEL and _MAGS == 'v' else [],
        v_Ta_plot_mk = expand(path + 'markers/pmr_v-project/01_Taxa_visualization/{project}.mk', project=PROJECT_INDEX) if RUN_PROJECT_LEVEL and _MAGS == 'v' else [],
        v_KO_plot_mk = expand(path + 'markers/pmr_v-project/03_Fun_visualization/{project}.mk', project=PROJECT_INDEX) if RUN_PROJECT_LEVEL and _MAGS == 'v' else [],
        v_mc_mk = expand(path + 'markers/pmr_v-project/03_community_module/{project}.mk',project=PROJECT_INDEX) if RUN_PROJECT_LEVEL and _MAGS == 'v' else [],
        v_mc_mk3 = expand(path + 'markers/pmr_v-project/06_module_completeness/{project}/{rank}.mk',project=PROJECT_INDEX, rank=RANKS) if RUN_PROJECT_LEVEL and _MAGS == 'v' else [],
        v_out_long_mc = expand(path + 'results/pmr_v-project/05_merge_tax_KO_depth/{project}/{rank}.mean.long.tsv',rank=RANKS, project=PROJECT_INDEX) if RUN_PROJECT_LEVEL and _MAGS == 'v' else [],
        v_out_long_Module_RelativeAbundance = expand(path + 'results/pmr_v-project/05_merge_tax_Module_RelativeAbundance/{project}/{rank}.mean.long.tsv',rank=RANKS, project=PROJECT_INDEX) if RUN_PROJECT_LEVEL and _MAGS == 'v' else [],
        v_out_long_Module_SumCoverage = expand(path + 'results/pmr_v-project/05_merge_tax_Module_SumCoverage/{project}/{rank}.mean.long.tsv',rank=RANKS, project=PROJECT_INDEX) if RUN_PROJECT_LEVEL and _MAGS == 'v' else [],
        # ---- vm chain (--mags vm): pmr_vm-sample / pmr_vm-project products ----
        vm_species_table = expand(path + 'results/pmr_vm-project/04_summary/{project}/mMAG_vMAG_species_abundance_matrix.tsv', project=PROJECT_INDEX) if RUN_PROJECT_LEVEL and _MAGS == 'vm' else [],
        vm_03_coverm_contigs = expand(path + 'results/pmr_vm-sample/03_coverm_contigs/{project}/{sample}.tsv', project=PROJECT_INDEX, sample=SAMPLE_INDEX) if _MAGS == 'vm' else [],
        vm_merge_05 = expand(path + 'results/pmr_vm-sample/05_merge_tax_fun_depth/{project}/KO_Taxa_RelativeAbundance/{sample}.KO_Taxa_RelativeAbundance.Taxa.tsv', project=PROJECT_INDEX, sample=SAMPLE_INDEX) if _MAGS == 'vm' else [],
        vm_mc_05 = expand(path + 'markers/pmr_vm-sample/05_visualize/{project}/{sample}.mk', project=PROJECT_INDEX, sample=SAMPLE_INDEX) if _MAGS == 'vm' else [],
        vm_mc_06 = expand(path + 'markers/pmr_vm-sample/06_module_completeness/{project}/{sample}.{rank}.mk',project=PROJECT_INDEX, rank=RANKS, sample=SAMPLE_INDEX) if _MAGS == 'vm' else [],
        vm_Ta_out_long = expand(path + 'results/pmr_vm-project/01_Taxa_Sample/{project}/Taxa_Sample.long.tsv',project=PROJECT_INDEX) if RUN_PROJECT_LEVEL and _MAGS == 'vm' else [],
        vm_KO_out_long = expand(path + 'results/pmr_vm-project/02_KO_Sample/{project}/KO_Sample.long.tsv',project=PROJECT_INDEX) if RUN_PROJECT_LEVEL and _MAGS == 'vm' else [],
        vm_Ta_plot_mk = expand(path + 'markers/pmr_vm-project/01_Taxa_visualization/{project}.mk', project=PROJECT_INDEX) if RUN_PROJECT_LEVEL and _MAGS == 'vm' else [],
        vm_KO_plot_mk = expand(path + 'markers/pmr_vm-project/03_Fun_visualization/{project}.mk', project=PROJECT_INDEX) if RUN_PROJECT_LEVEL and _MAGS == 'vm' else [],
        vm_mc_mk = expand(path + 'markers/pmr_vm-project/03_community_module/{project}.mk',project=PROJECT_INDEX) if RUN_PROJECT_LEVEL and _MAGS == 'vm' else [],
        vm_mc_mk3 = expand(path + 'markers/pmr_vm-project/06_module_completeness/{project}/{rank}.mk',project=PROJECT_INDEX, rank=RANKS) if RUN_PROJECT_LEVEL and _MAGS == 'vm' else [],
        vm_out_long_mc = expand(path + 'results/pmr_vm-project/05_merge_tax_KO_depth/{project}/{rank}.mean.long.tsv',rank=RANKS, project=PROJECT_INDEX) if RUN_PROJECT_LEVEL and _MAGS == 'vm' else [],
        vm_out_long_Module_RelativeAbundance = expand(path + 'results/pmr_vm-project/05_merge_tax_Module_RelativeAbundance/{project}/{rank}.mean.long.tsv',rank=RANKS, project=PROJECT_INDEX) if RUN_PROJECT_LEVEL and _MAGS == 'vm' else [],
        vm_out_long_Module_SumCoverage = expand(path + 'results/pmr_vm-project/05_merge_tax_Module_SumCoverage/{project}/{rank}.mean.long.tsv',rank=RANKS, project=PROJECT_INDEX) if RUN_PROJECT_LEVEL and _MAGS == 'vm' else [],
        # ---- v-p products (--mags v/vm): pure project-level vMAG chain ----
        p_vc_01 = expand(path + 'results/v-project/06_vclust/{project}/clusters.tsv', project=PROJECT_INDEX) if RUN_PROJECT_LEVEL and _MAGS in ('v', 'vm') else [],
        p_vp_iphop = expand(path + 'markers/v-project/07_iPHoP/done/{project}.mk', project=PROJECT_INDEX) if RUN_PROJECT_LEVEL and _MAGS in ('v', 'vm') else [],
        p_vp_10 = expand(path + 'markers/v-project/10_summary/done/{project}.mk', project=PROJECT_INDEX) if RUN_PROJECT_LEVEL and _MAGS in ('v', 'vm') else [],
        p_vp_11 = expand(path + 'markers/v-project/11_visualization/done/{project}.mk', project=PROJECT_INDEX) if RUN_PROJECT_LEVEL and _MAGS in ('v', 'vm') else [],
        # ---- m-project summary (--mags vm): dRep refs + mMAG summary ----
        p_m_mk_drep99 = expand(path + 'markers/m-project/01_drep_project_99/{project}.mk', project=PROJECT_INDEX) if RUN_PROJECT_LEVEL and _MAGS == 'vm' else [],
        p_m_mk_drep95 = expand(path + 'markers/m-project/01_drep_project_95/{project}.mk', project=PROJECT_INDEX) if RUN_PROJECT_LEVEL and _MAGS == 'vm' else [],
        p_m_cat_info = expand(path + 'results/m-project/04_summary/{project}/mMAG_summary_info.tsv', project=PROJECT_INDEX) if RUN_PROJECT_LEVEL and _MAGS == 'vm' else [],
        # ---- keep_intermediate cleanup (conditioned by chain ownership) ----
        # metawrap binning intermediates: only m/vm run binning
        s_mk_10 = expand(path + 'markers/b-sample/10_delete_metawrap_bam/{sample}.mk', sample=SAMPLE_INDEX) if _MAGS in ('m', 'vm') and not keep_intermediate('bam') else [],
        s_mk_15 = expand(path + 'markers/b-sample/12_delete_metawrap_work_files/{sample}.mk', sample=SAMPLE_INDEX) if _MAGS in ('m', 'vm') and config.get('skip_level', 0) < 3 and not keep_intermediate('assembly') else [],
        # m dRep intermediates: only m/vm run the m chain
        s_mk_12 = expand(path + 'markers/m-sample/12_delete_drep_data/{sample}.mk', sample=SAMPLE_INDEX) if _MAGS in ('m', 'vm') and not keep_intermediate('drep_data') else [],
        p_mk_16 = expand(path + 'markers/m-project/12_delete_drep_data/{project}.mk', project=PROJECT_INDEX) if RUN_PROJECT_LEVEL and _MAGS in ('m', 'vm') and not keep_intermediate('drep_data') else [],
        # v-sample intermediates: only v/vm run the v chain
        s_mk_13 = expand(path + 'markers/v-sample/12_delete_bam/{sample}.mk', sample=SAMPLE_INDEX) if _MAGS in ('v', 'vm') and not keep_intermediate('bam') else [],
        s_mk_14 = expand(path + 'markers/v-sample/12_delete_vs2_iter0/{sample}.mk', sample=SAMPLE_INDEX) if _MAGS in ('v', 'vm') and not keep_intermediate('vs2_iter0') else [],
        # pmr_m project-reference mapping BAMs: only the m chain
        m_mk_12 = expand(path + 'markers/pmr_m-sample/12_delete_bam/{project}/{sample}.mk', project=PROJECT_INDEX, sample=SAMPLE_INDEX) if _MAGS == 'm' and not keep_intermediate('bam') else [],
        # pmr_v merged-mapping BAMs: only the v chain
        v_mk_12 = expand(path + 'markers/pmr_v-project/12_delete_bam/{project}/{sample}.mk', project=PROJECT_INDEX, sample=SAMPLE_INDEX) if RUN_PROJECT_LEVEL and _MAGS == 'v' and not keep_intermediate('bam') else [],
        # pmr_vm merged-mapping BAMs: only the vm chain
        vm_mk_12 = expand(path + 'markers/pmr_vm-sample/12_delete_bam/{project}/{sample}.mk', project=PROJECT_INDEX, sample=SAMPLE_INDEX) if _MAGS == 'vm' and not keep_intermediate('bam') else [],
        s_b_mk_03 = expand(path + 'markers/b-sample/03_delete_assembly_intermediate/{sample}.mk', sample=SAMPLE_INDEX) if config.get('skip_level', 0) < 2 and not keep_intermediate('assembly') else [],
