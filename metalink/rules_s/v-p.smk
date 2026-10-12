'''
/*
 * @Date: 2022-12-04 17:15:47
 * @LastEditors: LiuyangLi
 * @LastEditTime: 2024-10-25 16:57:02
 * @FilePath: v-p.smk
 * @Description: project-level vMAG aggregation: cross-sample vclust dereplication, iPHoP host prediction, summary and visualization
 */
'''
import os
import sys
import datetime

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

path = config["Directory"]["work"]
path = os.path.join(path, '')
os.makedirs(path, exist_ok=True)
os.chdir(path)

def _ml_log(level, msg):
    """Print a Metalink-style log message to stderr."""
    if os.environ.get('METALINK_SILENT'):
        return
    print(f"[{datetime.datetime.now().strftime('%Y-%m-%d %H:%M')} {level}] {msg}", file=sys.stderr)

_MODES_WITH_MMAG = ['vm']

def _mode_has_mmag():
    """True when the m-p dRep chain runs alongside v-p: vm mode, or pmr --mags vm."""
    mode = config.get('mode', '')
    return mode in _MODES_WITH_MMAG or (mode == 'pmr' and config.get('mags', 'm') == 'vm')

def _vp_chain_active():
    """True when the pure v-p (project-level vMAG) chain is part of the workflow:
    v/vm modes, or pmr with --mags v/vm."""
    mode = config.get('mode', '')
    return mode in ('v', 'vm') or (mode == 'pmr' and config.get('mags', 'm') in ('v', 'vm'))

def _get_iphop_warning_msg(param_key, level):
    mode = config.get('mode', '')
    if not _vp_chain_active():
        return ''
    custom = config[param_key]["07_iPHoP"]["custom_MAGs_dir"]
    if custom in ("auto_drep99", "auto_drep95") and not _mode_has_mmag():
        return f"[ERROR] custom_MAGs_dir='{custom}' requires the mMAG dRep chain (vm mode or pmr --mags vm); mode='{mode}' mags='{config.get('mags', 'm')}' does not provide it. Use 'none' or a user-defined path."
    if custom not in ("none", "auto_drep99", "auto_drep95"):
        if not os.path.isdir(custom):
            return f"[ERROR] custom_MAGs_dir '{custom}' does not exist or is not a directory."
        fna_files = [f for f in os.listdir(custom) if f.endswith('.fna')]
        if not fna_files:
            return f"[ERROR] custom_MAGs_dir '{custom}' contains no *.fna files."
    return ''

# Global config validation: strict check for modes that run the v-p chain
for _param_key, _level in [("Params_v_project", "project")]:
    _mode = config.get('mode', '')
    if not _vp_chain_active():
        continue
    _custom = config[_param_key]["07_iPHoP"]["custom_MAGs_dir"]
    if _custom == "none":
        _ml_log("INFO", f"Config validation: {_level}-level custom_MAGs_dir='none', iPHoP custom MAGs disabled.")
    elif _custom in ("auto_drep99", "auto_drep95") and not _mode_has_mmag():
        _ml_log("ERROR", f"Config validation: {_param_key}['07_iPHoP']['custom_MAGs_dir']='{_custom}' requires the mMAG dRep chain (vm mode or pmr --mags vm); mode='{_mode}' mags='{config.get('mags', 'm')}' does not provide it. Use 'none' or a user-defined path.")
        sys.exit(1)
    elif _custom not in ("none", "auto_drep99", "auto_drep95"):
        if not os.path.isdir(_custom):
            _ml_log("ERROR", f"Config validation: custom_MAGs_dir '{_custom}' does not exist or is not a directory.")
            sys.exit(1)
        _fna_files = [f for f in os.listdir(_custom) if f.endswith('.fna')]
        if not _fna_files:
            _ml_log("ERROR", f"Config validation: custom_MAGs_dir '{_custom}' contains no *.fna files.")
            sys.exit(1)
        _ml_log("INFO", f"Config validation: {_level}-level custom_MAGs_dir='{_custom}' is active for mode='{_mode}' ({len(_fna_files)} *.fna files found).")
    else:
        _ml_log("INFO", f"Config validation: {_level}-level custom_MAGs_dir='{_custom}' is active for mode='{_mode}'.")


folder_checks = [
                path + 'results/',
                path + 'results/v-project',
                path + 'results/v-project/05_CheckV',
                path + 'results/v-project/06_vclust',
                path + 'results/v-project/07_iPHoP',
                path + 'results/v-project/10_summary',
                path + 'results/v-project/11_visualization'
                ]
for item in folder_checks:
    if not os.path.isdir(item): os.mkdir(item)


def _v_project_07_iphop_get_suffix():
    mode = config.get('mode', '')
    if not _vp_chain_active():
        return "none"
    custom = config["Params_v_project"]["07_iPHoP"]["custom_MAGs_dir"]
    if custom == "none":
        return "none"
    if custom in ("auto_drep99", "auto_drep95") and not _mode_has_mmag():
        raise ValueError(f"[ERROR] custom_MAGs_dir='{custom}' requires the mMAG dRep chain (vm mode or pmr --mags vm); mode='{mode}' mags='{config.get('mags', 'm')}' does not provide it. Use 'none' or a user-defined path.")
    if custom not in ("none", "auto_drep99", "auto_drep95"):
        if not os.path.isdir(custom):
            raise ValueError(f"[ERROR] custom_MAGs_dir '{custom}' does not exist or is not a directory.")
        fna_files = [f for f in os.listdir(custom) if f.endswith('.fna')]
        if not fna_files:
            raise ValueError(f"[ERROR] custom_MAGs_dir '{custom}' contains no *.fna files.")
    if custom == "auto_drep99":
        return "auto_drep99"
    elif custom == "auto_drep95":
        return "auto_drep95"
    else:
        return "custom"

def _v_project_07_iphop_get_mk(wildcards):
    suffix = _v_project_07_iphop_get_suffix()
    if suffix == "auto_drep99":
        return rules.v_project_07_iphop_auto_drep99.output.mk
    elif suffix == "auto_drep95":
        return rules.v_project_07_iphop_auto_drep95.output.mk
    elif suffix == "none":
        return rules.v_project_07_iphop_none.output.mk
    else:
        return rules.v_project_07_iphop_custom.output.mk

v_project_get_07_iphop_rule_mk = lambda wc: _v_project_07_iphop_get_mk(wc)

def _v_project_07_iphop_get_outdir(wildcards):
    suffix = _v_project_07_iphop_get_suffix()
    return path + f'results/v-project/07_iPHoP/{suffix}/{wildcards.project}/'


rule v_project_01_vclust:
    input:
        mk = expand(path + 'markers/v-sample/10_summary/{sample}.mk', sample=SAMPLE_INDEX),
    output:
        Project_Nlinked_viral_gn = path + 'results/v-project/05_CheckV/{project}/all_vRhyme_fasta.Nlinked_viral_gn.fna',
        fltr_txt = path + 'results/v-project/06_vclust/{project}/fltr.txt',
        ani_tsv = path + 'results/v-project/06_vclust/{project}/ani.tsv',
        ani_ids_tsv = path + 'results/v-project/06_vclust/{project}/ani.ids.tsv',
        ani_ids_scale_tsv = path + 'results/v-project/06_vclust/{project}/ani.ids.scale.tsv',
        clusters_tsv = path + 'results/v-project/06_vclust/{project}/clusters.tsv',
        unique_clusters = path + 'results/v-project/06_vclust/{project}/unique_clusters.tsv',
        vMAGs_vclust = path + 'results/v-project/06_vclust/{project}/vMAGs_vclust.fna',
    params:
        Nlinked_viral_gn = expand(path + 'results/v-sample/05_CheckV/{sample}/all_vRhyme_fasta.Nlinked_viral_gn.fna', sample=SAMPLE_INDEX),
        ani = config["Params_v_project"]["06_vclust"]["ani"],
        bin_folder = expand(path + 'results/v-sample/04_vRhyme/{sample}/vRhyme_both_binned_unbinned/', sample=SAMPLE_INDEX),
        vMAGs_folder = path + 'results/v-project/06_vclust/{project}/vMAGs_clusters_fna/',
        prefilter_extra_args = config["Params_v_project"]["06_vclust"].get("prefilter_extra_args", ""),
        align_extra_args = config["Params_v_project"]["06_vclust"].get("align_extra_args", ""),
        cluster_extra_args = config["Params_v_project"]["06_vclust"].get("cluster_extra_args", ""),
    conda:
        config["conda_envs"]["metalink-vclust"]
    priority: 1
    threads: config["Params_v_project"]["06_vclust"]["threads"]
    log:
        path + 'logs/v-project/06_vclust/{project}.06_vclust.log'
    benchmark:
        path + 'benchmarks/v-project/06_vclust/{project}.06_vclust.benchmark'
    shell:
        """
        find {params.Nlinked_viral_gn} -type f -print0 | xargs -0 cat > {output.Project_Nlinked_viral_gn}

        if [ -s {output.Project_Nlinked_viral_gn} ]; then
            n_seq=$(grep -c '^>' {output.Project_Nlinked_viral_gn})
            mkdir -p {params.vMAGs_folder} 2>{log}
            if [ "$n_seq" -ge 2 ]; then
                echo "Running vclust pipeline." >> {log} 2>&1
                # Prefilter similar genome sequence pairs before conducting pairwise alignments.
                vclust prefilter -t {threads} -i {output.Project_Nlinked_viral_gn} -o {output.fltr_txt} {params.prefilter_extra_args} > {log} 2>&1

                # Align similar genome sequence pairs and calculate pairwise ANI measures.
                vclust align -t {threads} -i {output.Project_Nlinked_viral_gn} -o {output.ani_tsv} --filter {output.fltr_txt} {params.align_extra_args} >> {log} 2>&1

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

                while IFS= read -r c; do
                    hits=$(find {params.bin_folder} -type f -name "${{c}}.fna" 2>> {log})
                    nh=$(printf '%s\\n' "$hits" | grep -c .)
                    if [ "$nh" -eq 0 ]; then
                        echo "[WARNING] No genome file found for cluster: $c" >> {log}
                    elif [ "$nh" -gt 1 ]; then
                        echo "[WARNING] Multiple files ($nh) matched cluster: $c; all copied, later ones overwrite earlier in vMAGs_folder" >> {log}
                    fi
                    printf '%s\\n' "$hits" | while IFS= read -r f; do
                        [ -n "$f" ] && cp "$f" {params.vMAGs_folder}/
                    done
                done < {output.unique_clusters} 2>> {log}
            else
                echo "Only one viral scaffold for {wildcards.project}; skip vclust."
                vclust prefilter -t {threads} -i {output.Project_Nlinked_viral_gn} -o {output.fltr_txt} {params.prefilter_extra_args} > {log} 2>&1
                seq_name=$(grep '^>' {output.Project_Nlinked_viral_gn} | sed 's/^>//' | awk '{{print $1}}')
                printf "object\\tcluster\\n%s\\t%s\\n" "$seq_name" "$seq_name" > {output.clusters_tsv}
                printf "$seq_name" > {output.unique_clusters}
                touch {output.ani_tsv} {output.ani_ids_tsv} {output.ani_ids_scale_tsv}
                while IFS= read -r c; do
                    hits=$(find {params.bin_folder} -type f -name "${{c}}.fna" 2>> {log})
                    nh=$(printf '%s\\n' "$hits" | grep -c .)
                    if [ "$nh" -eq 0 ]; then
                        echo "[WARNING] No genome file found for cluster: $c" >> {log}
                    elif [ "$nh" -gt 1 ]; then
                        echo "[WARNING] Multiple files ($nh) matched cluster: $c; all copied, later ones overwrite earlier in vMAGs_folder" >> {log}
                    fi
                    printf '%s\\n' "$hits" | while IFS= read -r f; do
                        [ -n "$f" ] && cp "$f" {params.vMAGs_folder}/
                    done
                done < {output.unique_clusters} 2>> {log}
            fi
            find {params.vMAGs_folder} -type f -print0 | xargs -0 cat > {output.vMAGs_vclust} 2>> {log}
        else
            echo "Empty viral scaffold for {wildcards.project}; skip vclust." > {log} 2>&1
            mkdir -p {params.vMAGs_folder}
            touch {output.fltr_txt} {output.ani_tsv} {output.ani_ids_tsv} \
                {output.ani_ids_scale_tsv} {output.clusters_tsv} {output.unique_clusters} \
                {output.vMAGs_vclust}
        fi
        """

rule v_project_07_iphop_none:
    input:
        clusters_tsv = path + 'results/v-project/06_vclust/{project}/clusters.tsv',
    output:
        mk = path + 'markers/v-project/07_iPHoP/none/{project}.mk',
    params:
        Nlinked_viral_gn = path + 'results/v-project/05_CheckV/{project}/all_vRhyme_fasta.Nlinked_viral_gn.fna',
        outdir_iphop = path + 'results/v-project/07_iPHoP/none/{project}/results',
        out_dir = path + 'results/v-project/',
        conda_env_iphop = conda_env_iphop,
        conda_env_genomad = conda_env_genomad,
        conda_env_gtdbtk = conda_env_gtdbtk,
        gtdbtk_db = gtdbtk_db,
        iphop_db = iphop_db,
        genomad_db = genomad_db,
        custom_MAGs_dir = "none",
        iPHoP_db_custom = "none",
        Project = path + 'results/v-project/06_vclust/{project}',
        warning_msg = _get_iphop_warning_msg("Params_v_project", "project"),
    conda:
        config["conda_envs"]["metalink-main"]
    priority: 1
    threads: config["Params_v_project"]["07_iPHoP"]["threads"]
    log:
        path + 'logs/v-project/07_iPHoP/none/{project}.virus_07_iPHoP.log'
    benchmark:
        path + 'benchmarks/v-project/07_iPHoP/none/{project}.virus_07_iPHoP.benchmark'
    shell:
        """
        echo "[INFO] custom_MAGs_dir resolved to: {params.custom_MAGs_dir}" > {log} 2>&1
        if [ -n "{params.warning_msg}" ]; then
            echo "{params.warning_msg}" >> {log} 2>&1
        fi
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
            --root_dir {Dir_script_metalink} >> {log} 2>&1
        touch {output.mk}
        """

rule v_project_07_iphop_auto_drep99:
    input:
        clusters_tsv = path + 'results/v-project/06_vclust/{project}/clusters.tsv',
        mk_drep99 = path + 'markers/m-project/01_drep_project_99/{project}.mk',
    output:
        mk = path + 'markers/v-project/07_iPHoP/auto_drep99/{project}.mk',
    params:
        Nlinked_viral_gn = path + 'results/v-project/05_CheckV/{project}/all_vRhyme_fasta.Nlinked_viral_gn.fna',
        outdir_iphop = path + 'results/v-project/07_iPHoP/auto_drep99/{project}/results',
        out_dir = path + 'results/v-project/',
        conda_env_iphop = conda_env_iphop,
        conda_env_genomad = conda_env_genomad,
        conda_env_gtdbtk = conda_env_gtdbtk,
        gtdbtk_db = gtdbtk_db,
        iphop_db = iphop_db,
        genomad_db = genomad_db,
        custom_MAGs_dir = path + 'results/m-project/01_drep_project_99/{project}/dereplicated_genomes/',
        custom_MAGs_dir_fna = path + 'results/v-project/07_iPHoP/auto_drep99/{project}/dereplicated_genomes_fna_link',
        iPHoP_db_custom = path + 'results/v-project/07_iPHoP/auto_drep99/{project}/iphop_db_temp/',
        Project = path + 'results/v-project/06_vclust/{project}',
        warning_msg = _get_iphop_warning_msg("Params_v_project", "project"),
    conda:
        config["conda_envs"]["metalink-main"]
    priority: 1
    threads: config["Params_v_project"]["07_iPHoP"]["threads"]
    log:
        path + 'logs/v-project/07_iPHoP/auto_drep99/{project}.virus_07_iPHoP.log'
    benchmark:
        path + 'benchmarks/v-project/07_iPHoP/auto_drep99/{project}.virus_07_iPHoP.benchmark'
    shell:
        r"""
        echo "[INFO] custom_MAGs_dir resolved to: {params.custom_MAGs_dir_fna}" > {log} 2>&1
        if [ -n "{params.warning_msg}" ]; then
            echo "{params.warning_msg}" >> {log} 2>&1
        fi

        if [ -d "{params.iPHoP_db_custom}" ];then rm -r {params.iPHoP_db_custom};else echo "Folder {params.iPHoP_db_custom} not exist, v_project_07_iphop_auto_drep99 good to go." >> {log};fi
        if [ -d "{params.custom_MAGs_dir_fna}" ];then rm -r {params.custom_MAGs_dir_fna};else echo "Folder {params.custom_MAGs_dir_fna} not exist, v_project_07_iphop_auto_drep99 good to go." >> {log};fi

        mkdir -p {params.custom_MAGs_dir_fna} >> {log} 2>&1
        find {params.custom_MAGs_dir} -maxdepth 1 \( -type f -o -type l \) -name "*.fa" | while IFS= read -r genome; do
            [ -e $genome ] || continue
            base=$(basename $genome .fa)
            ln -sf $genome {params.custom_MAGs_dir_fna}/${{base}}.fna
        done >> {log} 2>&1

        python {Dir_script_metalink}/viwrap-scripts/Metalink.07_iPHoP.py \
            --input_Nlinked_viral_gn {params.Nlinked_viral_gn} \
            --outdir_iphop {params.outdir_iphop} \
            --custom_MAGs_dir {params.custom_MAGs_dir_fna} \
            --iPHoP_db_custom {params.iPHoP_db_custom} \
            --threads {threads} \
            --conda_env_iphop {params.conda_env_iphop} \
            --conda_env_genomad {params.conda_env_genomad} \
            --conda_env_gtdbtk {params.conda_env_gtdbtk} \
            --gtdbtk_db {params.gtdbtk_db} \
            --iphop_db {params.iphop_db} \
            --genomad_db {params.genomad_db} \
            --root_dir {Dir_script_metalink} >> {log} 2>&1
        touch {output.mk}
        """

rule v_project_07_iphop_auto_drep95:
    input:
        clusters_tsv = path + 'results/v-project/06_vclust/{project}/clusters.tsv',
        mk_drep95 = path + 'markers/m-project/01_drep_project_95/{project}.mk',
    output:
        mk = path + 'markers/v-project/07_iPHoP/auto_drep95/{project}.mk',
    params:
        Nlinked_viral_gn = path + 'results/v-project/05_CheckV/{project}/all_vRhyme_fasta.Nlinked_viral_gn.fna',
        outdir_iphop = path + 'results/v-project/07_iPHoP/auto_drep95/{project}/results',
        out_dir = path + 'results/v-project/',
        conda_env_iphop = conda_env_iphop,
        conda_env_genomad = conda_env_genomad,
        conda_env_gtdbtk = conda_env_gtdbtk,
        gtdbtk_db = gtdbtk_db,
        iphop_db = iphop_db,
        genomad_db = genomad_db,
        custom_MAGs_dir = path + 'results/m-project/01_drep_project_95/{project}/dereplicated_genomes',
        custom_MAGs_dir_fna = path + 'results/v-project/07_iPHoP/auto_drep95/{project}/dereplicated_genomes_fna_link',
        iPHoP_db_custom = path + 'results/v-project/07_iPHoP/auto_drep95/{project}/iphop_db_temp',
        Project = path + 'results/v-project/06_vclust/{project}',
        warning_msg = _get_iphop_warning_msg("Params_v_project", "project"),
    conda:
        config["conda_envs"]["metalink-main"]
    priority: 1
    threads: config["Params_v_project"]["07_iPHoP"]["threads"]
    log:
        path + 'logs/v-project/07_iPHoP/auto_drep95/{project}.virus_07_iPHoP.log'
    benchmark:
        path + 'benchmarks/v-project/07_iPHoP/auto_drep95/{project}.virus_07_iPHoP.benchmark'
    shell:
        r"""
        echo "[INFO] custom_MAGs_dir resolved to: {params.custom_MAGs_dir}" > {log} 2>&1
        if [ -n "{params.warning_msg}" ]; then
            echo "{params.warning_msg}" >> {log} 2>&1
        fi

        if [ -d "{params.iPHoP_db_custom}" ];then rm -r {params.iPHoP_db_custom};else echo "Folder {params.iPHoP_db_custom} not exist, v_project_07_iphop_auto_drep95 good to go." >> {log};fi
        if [ -d "{params.custom_MAGs_dir_fna}" ];then rm -r {params.custom_MAGs_dir_fna};else echo "Folder {params.custom_MAGs_dir_fna} not exist, v_project_07_iphop_auto_drep95 good to go." >> {log};fi

        mkdir -p {params.custom_MAGs_dir_fna} >> {log} 2>&1
        find {params.custom_MAGs_dir} -maxdepth 1 \( -type f -o -type l \) -name "*.fa" | while IFS= read -r genome; do
            [ -e $genome ] || continue
            base=$(basename $genome .fa)
            ln -sf $genome {params.custom_MAGs_dir_fna}/${{base}}.fna
        done >> {log} 2>&1

        python {Dir_script_metalink}/viwrap-scripts/Metalink.07_iPHoP.py \
            --input_Nlinked_viral_gn {params.Nlinked_viral_gn} \
            --outdir_iphop {params.outdir_iphop} \
            --custom_MAGs_dir {params.custom_MAGs_dir_fna} \
            --iPHoP_db_custom {params.iPHoP_db_custom} \
            --threads {threads} \
            --conda_env_iphop {params.conda_env_iphop} \
            --conda_env_genomad {params.conda_env_genomad} \
            --conda_env_gtdbtk {params.conda_env_gtdbtk} \
            --gtdbtk_db {params.gtdbtk_db} \
            --iphop_db {params.iphop_db} \
            --genomad_db {params.genomad_db} \
            --root_dir {Dir_script_metalink} >> {log} 2>&1
        touch {output.mk}
        """

rule v_project_07_iphop_custom:
    input:
        clusters_tsv = path + 'results/v-project/06_vclust/{project}/clusters.tsv',
    output:
        mk = path + 'markers/v-project/07_iPHoP/custom/{project}.mk',
    params:
        Nlinked_viral_gn = path + 'results/v-project/05_CheckV/{project}/all_vRhyme_fasta.Nlinked_viral_gn.fna',
        outdir_iphop = path + 'results/v-project/07_iPHoP/custom/{project}/results',
        out_dir = path + 'results/v-project/',
        conda_env_iphop = conda_env_iphop,
        conda_env_genomad = conda_env_genomad,
        conda_env_gtdbtk = conda_env_gtdbtk,
        gtdbtk_db = gtdbtk_db,
        iphop_db = iphop_db,
        genomad_db = genomad_db,
        custom_MAGs_dir_fna = config["Params_v_project"]["07_iPHoP"]["custom_MAGs_dir"],
        iPHoP_db_custom = path + 'results/v-project/07_iPHoP/custom/{project}/iphop_db_temp/',
        Project = path + 'results/v-project/06_vclust/{project}',
        warning_msg = _get_iphop_warning_msg("Params_v_project", "project"),
    conda:
        config["conda_envs"]["metalink-main"]
    priority: 1
    threads: config["Params_v_project"]["07_iPHoP"]["threads"]
    log:
        path + 'logs/v-project/07_iPHoP/custom/{project}.virus_07_iPHoP.log'
    benchmark:
        path + 'benchmarks/v-project/07_iPHoP/custom/{project}.virus_07_iPHoP.benchmark'
    shell:
        """
        echo "[INFO] custom_MAGs_dir resolved to: {params.custom_MAGs_dir_fna}" > {log} 2>&1
        if [ -n "{params.warning_msg}" ]; then
            echo "{params.warning_msg}" >> {log} 2>&1
        fi

        if [ -d "{params.iPHoP_db_custom}" ];then rm -r {params.iPHoP_db_custom};else echo "Folder {params.iPHoP_db_custom} not exist, v_project_07_iphop_custom good to go." >> {log};fi

        python {Dir_script_metalink}/viwrap-scripts/Metalink.07_iPHoP.py \
            --input_Nlinked_viral_gn {params.Nlinked_viral_gn} \
            --outdir_iphop {params.outdir_iphop} \
            --custom_MAGs_dir {params.custom_MAGs_dir_fna} \
            --iPHoP_db_custom {params.iPHoP_db_custom} \
            --threads {threads} \
            --conda_env_iphop {params.conda_env_iphop} \
            --conda_env_genomad {params.conda_env_genomad} \
            --conda_env_gtdbtk {params.conda_env_gtdbtk} \
            --gtdbtk_db {params.gtdbtk_db} \
            --iphop_db {params.iphop_db} \
            --genomad_db {params.genomad_db} \
            --root_dir {Dir_script_metalink} >> {log} 2>&1
        touch {output.mk}
        """

rule v_project_07_iphop_done:
    input:
        mk = v_project_get_07_iphop_rule_mk,
    output:
        mk = path + 'markers/v-project/07_iPHoP/done/{project}.mk',
    shell:
        "touch {output.mk}"

def _v_project_10_summary_get_mk(wildcards):
    suffix = _v_project_07_iphop_get_suffix()
    if suffix == "auto_drep99":
        return rules.v_project_10_summary_auto_drep99.output.mk
    elif suffix == "auto_drep95":
        return rules.v_project_10_summary_auto_drep95.output.mk
    elif suffix == "none":
        return rules.v_project_10_summary_none.output.mk
    else:
        return rules.v_project_10_summary_custom.output.mk

v_project_get_10_summary_rule_mk = lambda wc: _v_project_10_summary_get_mk(wc)

rule v_project_10_summary_none:
    input:
        clusters_tsv = path + 'results/v-project/06_vclust/{project}/clusters.tsv',
        mk_07_iPHoP = path + 'markers/v-project/07_iPHoP/done/{project}.mk',
        mk_04 = expand(path + 'markers/v-sample/04_vRhyme/{sample}.mk', sample=SAMPLE_INDEX),
        mk_10 = expand(path + 'markers/v-sample/10_summary/{sample}.mk', sample=SAMPLE_INDEX),
    output:
        mk = path + 'markers/v-project/10_summary/none/{project}.mk',
    params:
        indir_iphop = lambda wc: _v_project_07_iphop_get_outdir(wc) + 'results/',
        outdir_summary = path + 'results/v-project/10_summary/none/{project}/',
        custom_MAGs_dir_fna = lambda wc: (
            path + 'results/v-project/07_iPHoP/auto_drep99/' + wc.project + '/dereplicated_genomes_fna_link/'
            if config["Params_v_project"]["07_iPHoP"]["custom_MAGs_dir"] == "auto_drep99" and _mode_has_mmag()
            else (
                path + 'results/v-project/07_iPHoP/auto_drep95/' + wc.project + '/dereplicated_genomes_fna_link/'
                if config["Params_v_project"]["07_iPHoP"]["custom_MAGs_dir"] == "auto_drep95" and _mode_has_mmag()
                else (
                    "none"
                    if config["Params_v_project"]["07_iPHoP"]["custom_MAGs_dir"] == "none"
                    else config["Params_v_project"]["07_iPHoP"]["custom_MAGs_dir"]
                )
            )
        ),
        iPHoP_db_custom = lambda wc: _v_project_07_iphop_get_outdir(wc) + 'iphop_db_temp/',
        warning_msg = _get_iphop_warning_msg("Params_v_project", "project"),
        binned_folder             = expand(path + 'results/v-sample/04_vRhyme/{sample}/vRhyme_best_bins_fasta_modified/', sample=SAMPLE_INDEX),
        unbinned_folder           = expand(path + 'results/v-sample/04_vRhyme/{sample}/vRhyme_unbinned_viral_gn_fasta/', sample=SAMPLE_INDEX),
        vMAG_summary_info_sample = expand(path + 'results/v-sample/10_summary/{sample}/vMAG_summary_info.tsv', sample=SAMPLE_INDEX),
        Virus_annotation_results  = expand(path + 'results/v-sample/10_summary/{sample}/Virus_annotation_results.tsv', sample=SAMPLE_INDEX),
        AMG_pro2info              = expand(path + 'results/v-sample/10_summary/{sample}/AMG_results/AMG_pro2info.txt', sample=SAMPLE_INDEX),
        AMG_pros_faa              = expand(path + 'results/v-sample/10_summary/{sample}/AMG_results/AMG_pros.faa', sample=SAMPLE_INDEX),
        Gn2amg_statistics         = expand(path + 'results/v-sample/10_summary/{sample}/AMG_results/Gn2amg_statistics.txt', sample=SAMPLE_INDEX),
        Host_prediction_to_genome = expand(path + 'results/v-sample/10_summary/{sample}/Host_prediction_to_genome_m90.csv', sample=SAMPLE_INDEX),
        Host_prediction_to_genus  = expand(path + 'results/v-sample/10_summary/{sample}/Host_prediction_to_genus_m90.csv', sample=SAMPLE_INDEX),
        vMAGs_folder = path + 'results/v-project/10_summary/none/{project}/Virus_genomes_files/',
        AMGs_folder = path + 'results/v-project/10_summary/none/{project}/AMG_results/',
        vMAG_summary_info_project = path + 'results/v-project/10_summary/none/{project}/vMAG_summary_info.tsv',
        tmp_dir = path + 'results/v-project/10_summary/none/{project}/tmp/',
    conda:
        config["conda_envs"]["metalink-main"]
    priority: 1
    threads: 8
    log:
        path + 'logs/v-project/10_summary/none/{project}.virus_06_summary.log'
    benchmark:
        path + 'benchmarks/v-project/10_summary/none/{project}.virus_06_summary.benchmark'
    shell:
        """
        echo "[INFO] custom_MAGs_dir resolved to: {params.custom_MAGs_dir_fna}" >> {log} 2>&1
        if [ -n "{params.warning_msg}" ]; then
            echo "{params.warning_msg}" > {log} 2>&1
        fi
        if [ -d {params.tmp_dir} ]; then
            rm -rf {params.tmp_dir}
        fi >> {log} 2>&1
        mkdir -p {params.tmp_dir} >> {log} 2>&1

        if [ -d {params.vMAGs_folder} ]; then
            rm -rf {params.vMAGs_folder}
        fi >> {log} 2>&1
        mkdir -p {params.vMAGs_folder} >> {log} 2>&1

        if [ -d {params.AMGs_folder} ]; then
            rm -rf {params.AMGs_folder}
        fi >> {log} 2>&1
        mkdir -p {params.AMGs_folder} >> {log} 2>&1

        # 0. Write huge file lists into temp files (one per line)
        printf '%s\n' {params.binned_folder} > {params.tmp_dir}/Virus_genomes_files.lst 2>> {log}
        printf '%s\n' {params.unbinned_folder} >> {params.tmp_dir}/Virus_genomes_files.lst 2>> {log}
        printf '%s\n' {params.vMAG_summary_info_sample} > {params.tmp_dir}/vMAG_summary_info.lst 2>> {log}
        printf '%s\n' {params.Virus_annotation_results}  > {params.tmp_dir}/Virus_annotation_results.lst 2>> {log}
        printf '%s\n' {params.AMG_pro2info}              > {params.tmp_dir}/AMG_pro2info.lst 2>> {log}
        printf '%s\n' {params.AMG_pros_faa}              > {params.tmp_dir}/AMG_pros.faa.lst 2>> {log}
        printf '%s\n' {params.Gn2amg_statistics}         > {params.tmp_dir}/Gn2amg_statistics.lst 2>> {log}
        printf '%s\n' {params.Host_prediction_to_genome} > {params.tmp_dir}/Host_prediction_to_genome.lst 2>> {log}
        printf '%s\n' {params.Host_prediction_to_genus}  > {params.tmp_dir}/Host_prediction_to_genus.lst 2>> {log}

        # 1.1. Stream copy viral genomes (batch 1000 dirs)
        while IFS= read -r d; do
            if [ -d "$d" ]; then
                find "$d" -type f -print0
            else
                echo "[WARNING] Virus genomes dir missing, skipped: $d" >> {log}
            fi
        done < {params.tmp_dir}/Virus_genomes_files.lst \
            > {params.tmp_dir}/Virus_genomes_files.found 2>> {log}
        tr '\\0' '\\n' < {params.tmp_dir}/Virus_genomes_files.found | xargs -r -n1 basename | sort | uniq -d | \
            while IFS= read -r b; do
                echo "[WARNING] Duplicate genome basename across dirs, later copy overwrites earlier: $b" >> {log}
            done
        xargs -0 -r cp -t {params.vMAGs_folder} < {params.tmp_dir}/Virus_genomes_files.found 2>> {log}

        # 1.2. AMG_pros.faa
        cat {params.tmp_dir}/AMG_pros.faa.lst | xargs -r cat > {params.outdir_summary}/AMG_results/AMG_pros.faa 2>> {log}

        # 2. Stream merge + replace cluster column
        cat {params.tmp_dir}/vMAG_summary_info.lst | \
            xargs -r awk '
            BEGIN {{FS=OFS="\t"}}
            ARGIND==1 && FNR>1 {{map[$1]=$2; next}}
            FILENAME=="{input.clusters_tsv}" {{next}}
            FNR==1  {{if(++hdr==1) print}}
            FNR>1   {{if($1 in map) $2=map[$1]; print}}
            ' {input.clusters_tsv} - > {params.vMAG_summary_info_project} 2>> {log}

        # 3. Virus_annotation_results
        cat {params.tmp_dir}/Virus_annotation_results.lst | \
            xargs -r awk 'FNR==1 && !hdr++{{print}} FNR>1{{print}}' \
            > {params.outdir_summary}/Virus_annotation_results.tsv 2>> {log}

        # 4. AMG_pro2info
        cat {params.tmp_dir}/AMG_pro2info.lst | \
            xargs -r awk 'FNR==1 && !hdr++{{print}} FNR>1{{print}}' \
            > {params.outdir_summary}/AMG_results/AMG_pro2info.txt 2>> {log}

        # 5. Gn2amg_statistics
        cat {params.tmp_dir}/Gn2amg_statistics.lst | \
            xargs -r awk 'FNR==1 && !hdr++{{print}} FNR>1{{print}}' \
            > {params.outdir_summary}/AMG_results/Gn2amg_statistics.txt 2>> {log}

        # 6. Host_prediction_to_genome_m90
        cat {params.tmp_dir}/Host_prediction_to_genome.lst | \
            xargs -r awk 'FNR==1 && !hdr++{{print}} FNR>1{{print}}' \
            > {params.outdir_summary}/Host_prediction_to_genome_m90_without_customDB.csv 2>> {log}

        # 7. Host_prediction_to_genus_m90
        cat {params.tmp_dir}/Host_prediction_to_genus.lst | \
            xargs -r awk 'FNR==1 && !hdr++{{print}} FNR>1{{print}}' \
            > {params.outdir_summary}/Host_prediction_to_genus_m90_without_customDB.csv 2>> {log}

        # 8. Merge iPHoP host predictions with customDB (Host_prediction_to_genome/genus_m90)
        python {Dir_script_metalink}/viwrap-scripts/Metalink.10_summary_project.py \
            --indir_iphop {params.indir_iphop} \
            --outdir_summary {params.outdir_summary} \
            --custom_MAGs_dir {params.custom_MAGs_dir_fna} \
            --iPHoP_db_custom {params.iPHoP_db_custom} \
            --root_dir {Dir_script_metalink} >> {log} 2>&1
        touch {output.mk}
    """

rule v_project_10_summary_auto_drep99:
    input:
        clusters_tsv = path + 'results/v-project/06_vclust/{project}/clusters.tsv',
        mk_07_iPHoP = path + 'markers/v-project/07_iPHoP/done/{project}.mk',
        mk_04 = expand(path + 'markers/v-sample/04_vRhyme/{sample}.mk', sample=SAMPLE_INDEX),
        mk_10 = expand(path + 'markers/v-sample/10_summary/{sample}.mk', sample=SAMPLE_INDEX),
    output:
        mk = path + 'markers/v-project/10_summary/auto_drep99/{project}.mk',
    params:
        indir_iphop = lambda wc: _v_project_07_iphop_get_outdir(wc) + 'results/',
        outdir_summary = path + 'results/v-project/10_summary/auto_drep99/{project}/',
        custom_MAGs_dir_fna = lambda wc: (
            path + 'results/v-project/07_iPHoP/auto_drep99/' + wc.project + '/dereplicated_genomes_fna_link/'
            if config["Params_v_project"]["07_iPHoP"]["custom_MAGs_dir"] == "auto_drep99" and _mode_has_mmag()
            else (
                path + 'results/v-project/07_iPHoP/auto_drep95/' + wc.project + '/dereplicated_genomes_fna_link/'
                if config["Params_v_project"]["07_iPHoP"]["custom_MAGs_dir"] == "auto_drep95" and _mode_has_mmag()
                else (
                    "none"
                    if config["Params_v_project"]["07_iPHoP"]["custom_MAGs_dir"] == "none"
                    else config["Params_v_project"]["07_iPHoP"]["custom_MAGs_dir"]
                )
            )
        ),
        iPHoP_db_custom = lambda wc: _v_project_07_iphop_get_outdir(wc) + 'iphop_db_temp/',
        warning_msg = _get_iphop_warning_msg("Params_v_project", "project"),
        binned_folder             = expand(path + 'results/v-sample/04_vRhyme/{sample}/vRhyme_best_bins_fasta_modified/', sample=SAMPLE_INDEX),
        unbinned_folder           = expand(path + 'results/v-sample/04_vRhyme/{sample}/vRhyme_unbinned_viral_gn_fasta/', sample=SAMPLE_INDEX),
        vMAG_summary_info_sample = expand(path + 'results/v-sample/10_summary/{sample}/vMAG_summary_info.tsv', sample=SAMPLE_INDEX),
        Virus_annotation_results  = expand(path + 'results/v-sample/10_summary/{sample}/Virus_annotation_results.tsv', sample=SAMPLE_INDEX),
        AMG_pro2info              = expand(path + 'results/v-sample/10_summary/{sample}/AMG_results/AMG_pro2info.txt', sample=SAMPLE_INDEX),
        AMG_pros_faa              = expand(path + 'results/v-sample/10_summary/{sample}/AMG_results/AMG_pros.faa', sample=SAMPLE_INDEX),
        Gn2amg_statistics         = expand(path + 'results/v-sample/10_summary/{sample}/AMG_results/Gn2amg_statistics.txt', sample=SAMPLE_INDEX),
        Host_prediction_to_genome = expand(path + 'results/v-sample/10_summary/{sample}/Host_prediction_to_genome_m90.csv', sample=SAMPLE_INDEX),
        Host_prediction_to_genus  = expand(path + 'results/v-sample/10_summary/{sample}/Host_prediction_to_genus_m90.csv', sample=SAMPLE_INDEX),
        vMAGs_folder = path + 'results/v-project/10_summary/auto_drep99/{project}/Virus_genomes_files/',
        AMGs_folder = path + 'results/v-project/10_summary/auto_drep99/{project}/AMG_results/',
        vMAG_summary_info_project = path + 'results/v-project/10_summary/auto_drep99/{project}/vMAG_summary_info.tsv',
        tmp_dir = path + 'results/v-project/10_summary/auto_drep99/{project}/tmp/',
    conda:
        config["conda_envs"]["metalink-main"]
    priority: 1
    threads: 8
    log:
        path + 'logs/v-project/10_summary/auto_drep99/{project}.virus_06_summary.log'
    benchmark:
        path + 'benchmarks/v-project/10_summary/auto_drep99/{project}.virus_06_summary.benchmark'
    shell:
        """
        echo "[INFO] custom_MAGs_dir resolved to: {params.custom_MAGs_dir_fna}" >> {log} 2>&1
        if [ -n "{params.warning_msg}" ]; then
            echo "{params.warning_msg}" > {log} 2>&1
        fi
        if [ -d {params.tmp_dir} ]; then
            rm -rf {params.tmp_dir}
        fi >> {log} 2>&1
        mkdir -p {params.tmp_dir} >> {log} 2>&1

        if [ -d {params.vMAGs_folder} ]; then
            rm -rf {params.vMAGs_folder}
        fi >> {log} 2>&1
        mkdir -p {params.vMAGs_folder} >> {log} 2>&1

        if [ -d {params.AMGs_folder} ]; then
            rm -rf {params.AMGs_folder}
        fi >> {log} 2>&1
        mkdir -p {params.AMGs_folder} >> {log} 2>&1

        # 0. Write huge file lists into temp files (one per line)
        printf '%s\n' {params.binned_folder} > {params.tmp_dir}/Virus_genomes_files.lst 2>> {log}
        printf '%s\n' {params.unbinned_folder} >> {params.tmp_dir}/Virus_genomes_files.lst 2>> {log}
        printf '%s\n' {params.vMAG_summary_info_sample} > {params.tmp_dir}/vMAG_summary_info.lst 2>> {log}
        printf '%s\n' {params.Virus_annotation_results}  > {params.tmp_dir}/Virus_annotation_results.lst 2>> {log}
        printf '%s\n' {params.AMG_pro2info}              > {params.tmp_dir}/AMG_pro2info.lst 2>> {log}
        printf '%s\n' {params.AMG_pros_faa}              > {params.tmp_dir}/AMG_pros.faa.lst 2>> {log}
        printf '%s\n' {params.Gn2amg_statistics}         > {params.tmp_dir}/Gn2amg_statistics.lst 2>> {log}
        printf '%s\n' {params.Host_prediction_to_genome} > {params.tmp_dir}/Host_prediction_to_genome.lst 2>> {log}
        printf '%s\n' {params.Host_prediction_to_genus}  > {params.tmp_dir}/Host_prediction_to_genus.lst 2>> {log}

        # 1.1. Stream copy viral genomes (batch 1000 dirs)
        while IFS= read -r d; do
            if [ -d "$d" ]; then
                find "$d" -type f -print0
            else
                echo "[WARNING] Virus genomes dir missing, skipped: $d" >> {log}
            fi
        done < {params.tmp_dir}/Virus_genomes_files.lst \
            > {params.tmp_dir}/Virus_genomes_files.found 2>> {log}
        tr '\\0' '\\n' < {params.tmp_dir}/Virus_genomes_files.found | xargs -r -n1 basename | sort | uniq -d | \
            while IFS= read -r b; do
                echo "[WARNING] Duplicate genome basename across dirs, later copy overwrites earlier: $b" >> {log}
            done
        xargs -0 -r cp -t {params.vMAGs_folder} < {params.tmp_dir}/Virus_genomes_files.found 2>> {log}

        # 1.2. AMG_pros.faa
        cat {params.tmp_dir}/AMG_pros.faa.lst | xargs -r cat > {params.outdir_summary}/AMG_results/AMG_pros.faa 2>> {log}

        # 2. Stream merge + replace cluster column
        cat {params.tmp_dir}/vMAG_summary_info.lst | \
            xargs -r awk '
            BEGIN {{FS=OFS="\t"}}
            ARGIND==1 && FNR>1 {{map[$1]=$2; next}}
            FILENAME=="{input.clusters_tsv}" {{next}}
            FNR==1  {{if(++hdr==1) print}}
            FNR>1   {{if($1 in map) $2=map[$1]; print}}
            ' {input.clusters_tsv} - > {params.vMAG_summary_info_project} 2>> {log}

        # 3. Virus_annotation_results
        cat {params.tmp_dir}/Virus_annotation_results.lst | \
            xargs -r awk 'FNR==1 && !hdr++{{print}} FNR>1{{print}}' \
            > {params.outdir_summary}/Virus_annotation_results.tsv 2>> {log}

        # 4. AMG_pro2info
        cat {params.tmp_dir}/AMG_pro2info.lst | \
            xargs -r awk 'FNR==1 && !hdr++{{print}} FNR>1{{print}}' \
            > {params.outdir_summary}/AMG_results/AMG_pro2info.txt 2>> {log}

        # 5. Gn2amg_statistics
        cat {params.tmp_dir}/Gn2amg_statistics.lst | \
            xargs -r awk 'FNR==1 && !hdr++{{print}} FNR>1{{print}}' \
            > {params.outdir_summary}/AMG_results/Gn2amg_statistics.txt 2>> {log}

        # 6. Host_prediction_to_genome_m90
        cat {params.tmp_dir}/Host_prediction_to_genome.lst | \
            xargs -r awk 'FNR==1 && !hdr++{{print}} FNR>1{{print}}' \
            > {params.outdir_summary}/Host_prediction_to_genome_m90_without_customDB.csv 2>> {log}

        # 7. Host_prediction_to_genus_m90
        cat {params.tmp_dir}/Host_prediction_to_genus.lst | \
            xargs -r awk 'FNR==1 && !hdr++{{print}} FNR>1{{print}}' \
            > {params.outdir_summary}/Host_prediction_to_genus_m90_without_customDB.csv 2>> {log}

        # 8. Merge iPHoP host predictions with customDB (Host_prediction_to_genome/genus_m90)
        python {Dir_script_metalink}/viwrap-scripts/Metalink.10_summary_project.py \
            --indir_iphop {params.indir_iphop} \
            --outdir_summary {params.outdir_summary} \
            --custom_MAGs_dir {params.custom_MAGs_dir_fna} \
            --iPHoP_db_custom {params.iPHoP_db_custom} \
            --root_dir {Dir_script_metalink} >> {log} 2>&1
        touch {output.mk}
    """

rule v_project_10_summary_auto_drep95:
    input:
        clusters_tsv = path + 'results/v-project/06_vclust/{project}/clusters.tsv',
        mk_07_iPHoP = path + 'markers/v-project/07_iPHoP/done/{project}.mk',
        mk_04 = expand(path + 'markers/v-sample/04_vRhyme/{sample}.mk', sample=SAMPLE_INDEX),
        mk_10 = expand(path + 'markers/v-sample/10_summary/{sample}.mk', sample=SAMPLE_INDEX),
    output:
        mk = path + 'markers/v-project/10_summary/auto_drep95/{project}.mk',
    params:
        indir_iphop = lambda wc: _v_project_07_iphop_get_outdir(wc) + 'results/',
        outdir_summary = path + 'results/v-project/10_summary/auto_drep95/{project}/',
        custom_MAGs_dir_fna = lambda wc: (
            path + 'results/v-project/07_iPHoP/auto_drep99/' + wc.project + '/dereplicated_genomes_fna_link/'
            if config["Params_v_project"]["07_iPHoP"]["custom_MAGs_dir"] == "auto_drep99" and _mode_has_mmag()
            else (
                path + 'results/v-project/07_iPHoP/auto_drep95/' + wc.project + '/dereplicated_genomes_fna_link/'
                if config["Params_v_project"]["07_iPHoP"]["custom_MAGs_dir"] == "auto_drep95" and _mode_has_mmag()
                else (
                    "none"
                    if config["Params_v_project"]["07_iPHoP"]["custom_MAGs_dir"] == "none"
                    else config["Params_v_project"]["07_iPHoP"]["custom_MAGs_dir"]
                )
            )
        ),
        iPHoP_db_custom = lambda wc: _v_project_07_iphop_get_outdir(wc) + 'iphop_db_temp/',
        warning_msg = _get_iphop_warning_msg("Params_v_project", "project"),
        binned_folder             = expand(path + 'results/v-sample/04_vRhyme/{sample}/vRhyme_best_bins_fasta_modified/', sample=SAMPLE_INDEX),
        unbinned_folder           = expand(path + 'results/v-sample/04_vRhyme/{sample}/vRhyme_unbinned_viral_gn_fasta/', sample=SAMPLE_INDEX),
        vMAG_summary_info_sample = expand(path + 'results/v-sample/10_summary/{sample}/vMAG_summary_info.tsv', sample=SAMPLE_INDEX),
        Virus_annotation_results  = expand(path + 'results/v-sample/10_summary/{sample}/Virus_annotation_results.tsv', sample=SAMPLE_INDEX),
        AMG_pro2info              = expand(path + 'results/v-sample/10_summary/{sample}/AMG_results/AMG_pro2info.txt', sample=SAMPLE_INDEX),
        AMG_pros_faa              = expand(path + 'results/v-sample/10_summary/{sample}/AMG_results/AMG_pros.faa', sample=SAMPLE_INDEX),
        Gn2amg_statistics         = expand(path + 'results/v-sample/10_summary/{sample}/AMG_results/Gn2amg_statistics.txt', sample=SAMPLE_INDEX),
        Host_prediction_to_genome = expand(path + 'results/v-sample/10_summary/{sample}/Host_prediction_to_genome_m90.csv', sample=SAMPLE_INDEX),
        Host_prediction_to_genus  = expand(path + 'results/v-sample/10_summary/{sample}/Host_prediction_to_genus_m90.csv', sample=SAMPLE_INDEX),
        vMAGs_folder = path + 'results/v-project/10_summary/auto_drep95/{project}/Virus_genomes_files/',
        AMGs_folder = path + 'results/v-project/10_summary/auto_drep95/{project}/AMG_results/',
        vMAG_summary_info_project = path + 'results/v-project/10_summary/auto_drep95/{project}/vMAG_summary_info.tsv',
        tmp_dir = path + 'results/v-project/10_summary/auto_drep95/{project}/tmp/',
    conda:
        config["conda_envs"]["metalink-main"]
    priority: 1
    threads: 8
    log:
        path + 'logs/v-project/10_summary/auto_drep95/{project}.virus_06_summary.log'
    benchmark:
        path + 'benchmarks/v-project/10_summary/auto_drep95/{project}.virus_06_summary.benchmark'
    shell:
        """
        echo "[INFO] custom_MAGs_dir resolved to: {params.custom_MAGs_dir_fna}" >> {log} 2>&1
        if [ -n "{params.warning_msg}" ]; then
            echo "{params.warning_msg}" > {log} 2>&1
        fi
        if [ -d {params.tmp_dir} ]; then
            rm -rf {params.tmp_dir}
        fi >> {log} 2>&1
        mkdir -p {params.tmp_dir} >> {log} 2>&1

        if [ -d {params.vMAGs_folder} ]; then
            rm -rf {params.vMAGs_folder}
        fi >> {log} 2>&1
        mkdir -p {params.vMAGs_folder} >> {log} 2>&1

        if [ -d {params.AMGs_folder} ]; then
            rm -rf {params.AMGs_folder}
        fi >> {log} 2>&1
        mkdir -p {params.AMGs_folder} >> {log} 2>&1

        # 0. Write huge file lists into temp files (one per line)
        printf '%s\n' {params.binned_folder} > {params.tmp_dir}/Virus_genomes_files.lst 2>> {log}
        printf '%s\n' {params.unbinned_folder} >> {params.tmp_dir}/Virus_genomes_files.lst 2>> {log}
        printf '%s\n' {params.vMAG_summary_info_sample} > {params.tmp_dir}/vMAG_summary_info.lst 2>> {log}
        printf '%s\n' {params.Virus_annotation_results}  > {params.tmp_dir}/Virus_annotation_results.lst 2>> {log}
        printf '%s\n' {params.AMG_pro2info}              > {params.tmp_dir}/AMG_pro2info.lst 2>> {log}
        printf '%s\n' {params.AMG_pros_faa}              > {params.tmp_dir}/AMG_pros.faa.lst 2>> {log}
        printf '%s\n' {params.Gn2amg_statistics}         > {params.tmp_dir}/Gn2amg_statistics.lst 2>> {log}
        printf '%s\n' {params.Host_prediction_to_genome} > {params.tmp_dir}/Host_prediction_to_genome.lst 2>> {log}
        printf '%s\n' {params.Host_prediction_to_genus}  > {params.tmp_dir}/Host_prediction_to_genus.lst 2>> {log}

        # 1.1. Stream copy viral genomes (batch 1000 dirs)
        while IFS= read -r d; do
            if [ -d "$d" ]; then
                find "$d" -type f -print0
            else
                echo "[WARNING] Virus genomes dir missing, skipped: $d" >> {log}
            fi
        done < {params.tmp_dir}/Virus_genomes_files.lst \
            > {params.tmp_dir}/Virus_genomes_files.found 2>> {log}
        tr '\\0' '\\n' < {params.tmp_dir}/Virus_genomes_files.found | xargs -r -n1 basename | sort | uniq -d | \
            while IFS= read -r b; do
                echo "[WARNING] Duplicate genome basename across dirs, later copy overwrites earlier: $b" >> {log}
            done
        xargs -0 -r cp -t {params.vMAGs_folder} < {params.tmp_dir}/Virus_genomes_files.found 2>> {log}

        # 1.2. AMG_pros.faa
        cat {params.tmp_dir}/AMG_pros.faa.lst | xargs -r cat > {params.outdir_summary}/AMG_results/AMG_pros.faa 2>> {log}

        # 2. Stream merge + replace cluster column
        cat {params.tmp_dir}/vMAG_summary_info.lst | \
            xargs -r awk '
            BEGIN {{FS=OFS="\t"}}
            ARGIND==1 && FNR>1 {{map[$1]=$2; next}}
            FILENAME=="{input.clusters_tsv}" {{next}}
            FNR==1  {{if(++hdr==1) print}}
            FNR>1   {{if($1 in map) $2=map[$1]; print}}
            ' {input.clusters_tsv} - > {params.vMAG_summary_info_project} 2>> {log}

        # 3. Virus_annotation_results
        cat {params.tmp_dir}/Virus_annotation_results.lst | \
            xargs -r awk 'FNR==1 && !hdr++{{print}} FNR>1{{print}}' \
            > {params.outdir_summary}/Virus_annotation_results.tsv 2>> {log}

        # 4. AMG_pro2info
        cat {params.tmp_dir}/AMG_pro2info.lst | \
            xargs -r awk 'FNR==1 && !hdr++{{print}} FNR>1{{print}}' \
            > {params.outdir_summary}/AMG_results/AMG_pro2info.txt 2>> {log}

        # 5. Gn2amg_statistics
        cat {params.tmp_dir}/Gn2amg_statistics.lst | \
            xargs -r awk 'FNR==1 && !hdr++{{print}} FNR>1{{print}}' \
            > {params.outdir_summary}/AMG_results/Gn2amg_statistics.txt 2>> {log}

        # 6. Host_prediction_to_genome_m90
        cat {params.tmp_dir}/Host_prediction_to_genome.lst | \
            xargs -r awk 'FNR==1 && !hdr++{{print}} FNR>1{{print}}' \
            > {params.outdir_summary}/Host_prediction_to_genome_m90_without_customDB.csv 2>> {log}

        # 7. Host_prediction_to_genus_m90
        cat {params.tmp_dir}/Host_prediction_to_genus.lst | \
            xargs -r awk 'FNR==1 && !hdr++{{print}} FNR>1{{print}}' \
            > {params.outdir_summary}/Host_prediction_to_genus_m90_without_customDB.csv 2>> {log}

        # 8. Merge iPHoP host predictions with customDB (Host_prediction_to_genome/genus_m90)
        python {Dir_script_metalink}/viwrap-scripts/Metalink.10_summary_project.py \
            --indir_iphop {params.indir_iphop} \
            --outdir_summary {params.outdir_summary} \
            --custom_MAGs_dir {params.custom_MAGs_dir_fna} \
            --iPHoP_db_custom {params.iPHoP_db_custom} \
            --root_dir {Dir_script_metalink} >> {log} 2>&1
        touch {output.mk}
    """

rule v_project_10_summary_custom:
    input:
        clusters_tsv = path + 'results/v-project/06_vclust/{project}/clusters.tsv',
        mk_07_iPHoP = path + 'markers/v-project/07_iPHoP/done/{project}.mk',
        mk_04 = expand(path + 'markers/v-sample/04_vRhyme/{sample}.mk', sample=SAMPLE_INDEX),
        mk_10 = expand(path + 'markers/v-sample/10_summary/{sample}.mk', sample=SAMPLE_INDEX),
    output:
        mk = path + 'markers/v-project/10_summary/custom/{project}.mk',
    params:
        indir_iphop = lambda wc: _v_project_07_iphop_get_outdir(wc) + 'results/',
        outdir_summary = path + 'results/v-project/10_summary/custom/{project}/',
        custom_MAGs_dir_fna = lambda wc: (
            path + 'results/v-project/07_iPHoP/auto_drep99/' + wc.project + '/dereplicated_genomes_fna_link/'
            if config["Params_v_project"]["07_iPHoP"]["custom_MAGs_dir"] == "auto_drep99" and _mode_has_mmag()
            else (
                path + 'results/v-project/07_iPHoP/auto_drep95/' + wc.project + '/dereplicated_genomes_fna_link/'
                if config["Params_v_project"]["07_iPHoP"]["custom_MAGs_dir"] == "auto_drep95" and _mode_has_mmag()
                else (
                    "none"
                    if config["Params_v_project"]["07_iPHoP"]["custom_MAGs_dir"] == "none"
                    else config["Params_v_project"]["07_iPHoP"]["custom_MAGs_dir"]
                )
            )
        ),
        iPHoP_db_custom = lambda wc: _v_project_07_iphop_get_outdir(wc) + 'iphop_db_temp/',
        warning_msg = _get_iphop_warning_msg("Params_v_project", "project"),
        binned_folder             = expand(path + 'results/v-sample/04_vRhyme/{sample}/vRhyme_best_bins_fasta_modified/', sample=SAMPLE_INDEX),
        unbinned_folder           = expand(path + 'results/v-sample/04_vRhyme/{sample}/vRhyme_unbinned_viral_gn_fasta/', sample=SAMPLE_INDEX),
        vMAG_summary_info_sample = expand(path + 'results/v-sample/10_summary/{sample}/vMAG_summary_info.tsv', sample=SAMPLE_INDEX),
        Virus_annotation_results  = expand(path + 'results/v-sample/10_summary/{sample}/Virus_annotation_results.tsv', sample=SAMPLE_INDEX),
        AMG_pro2info              = expand(path + 'results/v-sample/10_summary/{sample}/AMG_results/AMG_pro2info.txt', sample=SAMPLE_INDEX),
        AMG_pros_faa              = expand(path + 'results/v-sample/10_summary/{sample}/AMG_results/AMG_pros.faa', sample=SAMPLE_INDEX),
        Gn2amg_statistics         = expand(path + 'results/v-sample/10_summary/{sample}/AMG_results/Gn2amg_statistics.txt', sample=SAMPLE_INDEX),
        Host_prediction_to_genome = expand(path + 'results/v-sample/10_summary/{sample}/Host_prediction_to_genome_m90.csv', sample=SAMPLE_INDEX),
        Host_prediction_to_genus  = expand(path + 'results/v-sample/10_summary/{sample}/Host_prediction_to_genus_m90.csv', sample=SAMPLE_INDEX),
        vMAGs_folder = path + 'results/v-project/10_summary/custom/{project}/Virus_genomes_files/',
        AMGs_folder = path + 'results/v-project/10_summary/custom/{project}/AMG_results/',
        vMAG_summary_info_project = path + 'results/v-project/10_summary/custom/{project}/vMAG_summary_info.tsv',
        tmp_dir = path + 'results/v-project/10_summary/custom/{project}/tmp/',
    conda:
        config["conda_envs"]["metalink-main"]
    priority: 1
    threads: 8
    log:
        path + 'logs/v-project/10_summary/custom/{project}.virus_06_summary.log'
    benchmark:
        path + 'benchmarks/v-project/10_summary/custom/{project}.virus_06_summary.benchmark'
    shell:
        """
        echo "[INFO] custom_MAGs_dir resolved to: {params.custom_MAGs_dir_fna}" >> {log} 2>&1
        if [ -n "{params.warning_msg}" ]; then
            echo "{params.warning_msg}" > {log} 2>&1
        fi
        if [ -d {params.tmp_dir} ]; then
            rm -rf {params.tmp_dir}
        fi >> {log} 2>&1
        mkdir -p {params.tmp_dir} >> {log} 2>&1

        if [ -d {params.vMAGs_folder} ]; then
            rm -rf {params.vMAGs_folder}
        fi >> {log} 2>&1
        mkdir -p {params.vMAGs_folder} >> {log} 2>&1

        if [ -d {params.AMGs_folder} ]; then
            rm -rf {params.AMGs_folder}
        fi >> {log} 2>&1
        mkdir -p {params.AMGs_folder} >> {log} 2>&1

        # 0. Write huge file lists into temp files (one per line)
        printf '%s\n' {params.binned_folder} > {params.tmp_dir}/Virus_genomes_files.lst 2>> {log}
        printf '%s\n' {params.unbinned_folder} >> {params.tmp_dir}/Virus_genomes_files.lst 2>> {log}
        printf '%s\n' {params.vMAG_summary_info_sample} > {params.tmp_dir}/vMAG_summary_info.lst 2>> {log}
        printf '%s\n' {params.Virus_annotation_results}  > {params.tmp_dir}/Virus_annotation_results.lst 2>> {log}
        printf '%s\n' {params.AMG_pro2info}              > {params.tmp_dir}/AMG_pro2info.lst 2>> {log}
        printf '%s\n' {params.AMG_pros_faa}              > {params.tmp_dir}/AMG_pros.faa.lst 2>> {log}
        printf '%s\n' {params.Gn2amg_statistics}         > {params.tmp_dir}/Gn2amg_statistics.lst 2>> {log}
        printf '%s\n' {params.Host_prediction_to_genome} > {params.tmp_dir}/Host_prediction_to_genome.lst 2>> {log}
        printf '%s\n' {params.Host_prediction_to_genus}  > {params.tmp_dir}/Host_prediction_to_genus.lst 2>> {log}

        # 1.1. Stream copy viral genomes (batch 1000 dirs)
        while IFS= read -r d; do
            if [ -d "$d" ]; then
                find "$d" -type f -print0
            else
                echo "[WARNING] Virus genomes dir missing, skipped: $d" >> {log}
            fi
        done < {params.tmp_dir}/Virus_genomes_files.lst \
            > {params.tmp_dir}/Virus_genomes_files.found 2>> {log}
        tr '\\0' '\\n' < {params.tmp_dir}/Virus_genomes_files.found | xargs -r -n1 basename | sort | uniq -d | \
            while IFS= read -r b; do
                echo "[WARNING] Duplicate genome basename across dirs, later copy overwrites earlier: $b" >> {log}
            done
        xargs -0 -r cp -t {params.vMAGs_folder} < {params.tmp_dir}/Virus_genomes_files.found 2>> {log}

        # 1.2. AMG_pros.faa
        cat {params.tmp_dir}/AMG_pros.faa.lst | xargs -r cat > {params.outdir_summary}/AMG_results/AMG_pros.faa 2>> {log}

        # 2. Stream merge + replace cluster column
        cat {params.tmp_dir}/vMAG_summary_info.lst | \
            xargs -r awk '
            BEGIN {{FS=OFS="\t"}}
            ARGIND==1 && FNR>1 {{map[$1]=$2; next}}
            FILENAME=="{input.clusters_tsv}" {{next}}
            FNR==1  {{if(++hdr==1) print}}
            FNR>1   {{if($1 in map) $2=map[$1]; print}}
            ' {input.clusters_tsv} - > {params.vMAG_summary_info_project} 2>> {log}

        # 3. Virus_annotation_results
        cat {params.tmp_dir}/Virus_annotation_results.lst | \
            xargs -r awk 'FNR==1 && !hdr++{{print}} FNR>1{{print}}' \
            > {params.outdir_summary}/Virus_annotation_results.tsv 2>> {log}

        # 4. AMG_pro2info
        cat {params.tmp_dir}/AMG_pro2info.lst | \
            xargs -r awk 'FNR==1 && !hdr++{{print}} FNR>1{{print}}' \
            > {params.outdir_summary}/AMG_results/AMG_pro2info.txt 2>> {log}

        # 5. Gn2amg_statistics
        cat {params.tmp_dir}/Gn2amg_statistics.lst | \
            xargs -r awk 'FNR==1 && !hdr++{{print}} FNR>1{{print}}' \
            > {params.outdir_summary}/AMG_results/Gn2amg_statistics.txt 2>> {log}

        # 6. Host_prediction_to_genome_m90
        cat {params.tmp_dir}/Host_prediction_to_genome.lst | \
            xargs -r awk 'FNR==1 && !hdr++{{print}} FNR>1{{print}}' \
            > {params.outdir_summary}/Host_prediction_to_genome_m90_without_customDB.csv 2>> {log}

        # 7. Host_prediction_to_genus_m90
        cat {params.tmp_dir}/Host_prediction_to_genus.lst | \
            xargs -r awk 'FNR==1 && !hdr++{{print}} FNR>1{{print}}' \
            > {params.outdir_summary}/Host_prediction_to_genus_m90_without_customDB.csv 2>> {log}

        # 8. Merge iPHoP host predictions with customDB (Host_prediction_to_genome/genus_m90)
        python {Dir_script_metalink}/viwrap-scripts/Metalink.10_summary_project.py \
            --indir_iphop {params.indir_iphop} \
            --outdir_summary {params.outdir_summary} \
            --custom_MAGs_dir {params.custom_MAGs_dir_fna} \
            --iPHoP_db_custom {params.iPHoP_db_custom} \
            --root_dir {Dir_script_metalink} >> {log} 2>&1
        touch {output.mk}
    """

rule v_project_10_summary_done:
    input:
        mk = v_project_get_10_summary_rule_mk,
    output:
        mk = path + 'markers/v-project/10_summary/done/{project}.mk',
    shell:
        "touch {output.mk}"

def _v_project_11_visualize_get_mk(wildcards):
    suffix = _v_project_07_iphop_get_suffix()
    if suffix == "auto_drep99":
        return rules.v_project_11_visualize_auto_drep99.output.mk
    elif suffix == "auto_drep95":
        return rules.v_project_11_visualize_auto_drep95.output.mk
    elif suffix == "none":
        return rules.v_project_11_visualize_none.output.mk
    else:
        return rules.v_project_11_visualize_custom.output.mk

v_project_get_11_visualize_rule_mk = lambda wc: _v_project_11_visualize_get_mk(wc)

rule v_project_11_visualize_none:
    input:
        mk = path + 'markers/v-project/10_summary/none/{project}.mk',
    output:
        mk = path + 'markers/v-project/11_visualization/none/{project}.mk',
    params:
        vibrant_db = vibrant_db,
        indir_summary = path + 'results/v-project/10_summary/none/{project}/',
        outdir_visualize = path + 'results/v-project/11_visualization/none/{project}/',
    conda:
        config["conda_envs"]["metalink-main"]
    priority: 1
    threads: 8
    log:
        path + 'logs/v-project/11_visualization/none/{project}.virus_11_visualize.log'
    benchmark:
        path + 'benchmarks/v-project/11_visualization/none/{project}.virus_11_visualize.benchmark'
    shell:
        """
        python {Dir_script_metalink}/viwrap-scripts/Metalink.11_visualize_project.py \
            --indir_summary {params.indir_summary} \
            --outdir_visualize {params.outdir_visualize} \
            --threads {threads} \
            --vibrant_db {params.vibrant_db} \
            --root_dir {Dir_script_metalink} > {log} 2>&1
        touch {output.mk}
        """

rule v_project_11_visualize_auto_drep99:
    input:
        mk = path + 'markers/v-project/10_summary/auto_drep99/{project}.mk',
    output:
        mk = path + 'markers/v-project/11_visualization/auto_drep99/{project}.mk',
    params:
        vibrant_db = vibrant_db,
        indir_summary = path + 'results/v-project/10_summary/auto_drep99/{project}/',
        outdir_visualize = path + 'results/v-project/11_visualization/auto_drep99/{project}/',
    conda:
        config["conda_envs"]["metalink-main"]
    priority: 1
    threads: 8
    log:
        path + 'logs/v-project/11_visualization/auto_drep99/{project}.virus_11_visualize.log'
    benchmark:
        path + 'benchmarks/v-project/11_visualization/auto_drep99/{project}.virus_11_visualize.benchmark'
    shell:
        """
        python {Dir_script_metalink}/viwrap-scripts/Metalink.11_visualize_project.py \
            --indir_summary {params.indir_summary} \
            --outdir_visualize {params.outdir_visualize} \
            --threads {threads} \
            --vibrant_db {params.vibrant_db} \
            --root_dir {Dir_script_metalink} > {log} 2>&1
        touch {output.mk}
        """

rule v_project_11_visualize_auto_drep95:
    input:
        mk = path + 'markers/v-project/10_summary/auto_drep95/{project}.mk',
    output:
        mk = path + 'markers/v-project/11_visualization/auto_drep95/{project}.mk',
    params:
        vibrant_db = vibrant_db,
        indir_summary = path + 'results/v-project/10_summary/auto_drep95/{project}/',
        outdir_visualize = path + 'results/v-project/11_visualization/auto_drep95/{project}/',
    conda:
        config["conda_envs"]["metalink-main"]
    priority: 1
    threads: 8
    log:
        path + 'logs/v-project/11_visualization/auto_drep95/{project}.virus_11_visualize.log'
    benchmark:
        path + 'benchmarks/v-project/11_visualization/auto_drep95/{project}.virus_11_visualize.benchmark'
    shell:
        """
        python {Dir_script_metalink}/viwrap-scripts/Metalink.11_visualize_project.py \
            --indir_summary {params.indir_summary} \
            --outdir_visualize {params.outdir_visualize} \
            --threads {threads} \
            --vibrant_db {params.vibrant_db} \
            --root_dir {Dir_script_metalink} > {log} 2>&1
        touch {output.mk}
        """

rule v_project_11_visualize_custom:
    input:
        mk = path + 'markers/v-project/10_summary/custom/{project}.mk',
    output:
        mk = path + 'markers/v-project/11_visualization/custom/{project}.mk',
    params:
        vibrant_db = vibrant_db,
        indir_summary = path + 'results/v-project/10_summary/custom/{project}/',
        outdir_visualize = path + 'results/v-project/11_visualization/custom/{project}/',
    conda:
        config["conda_envs"]["metalink-main"]
    priority: 1
    threads: 8
    log:
        path + 'logs/v-project/11_visualization/custom/{project}.virus_11_visualize.log'
    benchmark:
        path + 'benchmarks/v-project/11_visualization/custom/{project}.virus_11_visualize.benchmark'
    shell:
        """
        python {Dir_script_metalink}/viwrap-scripts/Metalink.11_visualize_project.py \
            --indir_summary {params.indir_summary} \
            --outdir_visualize {params.outdir_visualize} \
            --threads {threads} \
            --vibrant_db {params.vibrant_db} \
            --root_dir {Dir_script_metalink} > {log} 2>&1
        touch {output.mk}
        """

rule v_project_11_visualize_done:
    input:
        mk = v_project_get_11_visualize_rule_mk,
    output:
        mk = path + 'markers/v-project/11_visualization/done/{project}.mk',
    shell:
        "touch {output.mk}"
