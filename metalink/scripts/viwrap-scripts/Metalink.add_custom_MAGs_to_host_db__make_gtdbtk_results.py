#!/usr/bin/env python3

try:
    import warnings
    import sys
    import os
    from glob import glob
    warnings.filterwarnings("ignore")
except Exception as e:
    sys.stderr.write(str(e) + "\n\n")
    exit(1)

def add_custom_MAGs_to_host_db__make_gtdbtk_results(viwrap_outdir, custom_MAGs_dir, threads, gtdbtk_db):
    os.environ['GTDBTK_DATA_PATH'] = gtdbtk_db
    fasta_files = glob(os.path.join(custom_MAGs_dir, '*.fna'))
    print(fasta_files)
    if not fasta_files:
        # Step 1 Check the extension of input genomes
        sys.exit(f'Please make sure your there are input MAGs in {custom_MAGs_dir} and all of them end with ".fna"')
    else:
        bacteria_cmd = f'gtdbtk de_novo_wf --genome_dir {custom_MAGs_dir} --bacteria --outgroup_taxon p__Patescibacteria --out_dir {viwrap_outdir} --cpus {threads} --force --extension fna'
        print(bacteria_cmd)
        os.system(bacteria_cmd)

        archaea_cmd = f'gtdbtk de_novo_wf --genome_dir {custom_MAGs_dir} --archaea --outgroup_taxon p__Altiarchaeota --out_dir {viwrap_outdir} --cpus {threads} --force --extension fna'
        print(archaea_cmd)
        os.system(archaea_cmd)

viwrap_outdir, custom_MAGs_dir, threads, gtdbtk_db = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
add_custom_MAGs_to_host_db__make_gtdbtk_results(viwrap_outdir, custom_MAGs_dir, threads, gtdbtk_db)
