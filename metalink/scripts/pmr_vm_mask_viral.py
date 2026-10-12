#!/usr/bin/env python3
"""Mask viral regions on host MAG contigs for the pmr_vm merged reference.

For every final viral scaffold (per-sample vMAG_summary_info.tsv, scaffold
column), the parent-contig coordinates embedded in its ID
(|provirus_S_E from geNomad, |vs2region_S_E from the VS2 chain after
Metalink.embed_viral_coords.py) are used to:

  - interval scaffolds  -> mask [S, E] with N on the matching MAG contig
  - whole-contig viral  -> drop the entire MAG contig (MAG contamination)

MAG contigs are matched to viral parents by the (sample, contig_id) key:
viral parent "{sample}---{cid}"  <->  MAG contig "{sample}_bin.N---{cid}"
(coordinates are identical; only the contig name differs).

Host ORFs overlapping masked intervals (or sitting on dropped contigs) are
listed in excluded_ORFs.txt for combine_taxonomy_function_depth.pmr_vm.cmd.R
to drop from the host ORF tables (qaccver key).

Old-style names from runs predating the coordinate embedding
(||{i}_partial*, ||full_M, ||lt2gene_M) carry no coordinates and are skipped
with a warning; ||full / ||lt2gene without _M and plain names are treated as
whole-contig viral. Only stdlib is used.
"""

import argparse
import os
import re
import sys

PROV_RE = re.compile(r'\|(?:provirus|vs2region)_(\d+)_(\d+)$')
PARENT_RE = re.compile(r'^(.+)---(.+)$')
MAGCTG_RE = re.compile(r'^(.+)_bin\.\d+---(.+)$')


def classify_token(token, warns):
    """Classify one viral scaffold token.

    Returns ('interval', parent, S, E) or ('whole', parent) or None (skipped).
    """
    m = PROV_RE.search(token)
    if m:
        return ('interval', token[:m.start()], int(m.group(1)), int(m.group(2)))
    if '||' in token:
        base, suffix = token.split('||', 1)
        if suffix in ('full', 'lt2gene'):
            return ('whole', base)
        warns.append(f'[skip] old-style fragment without embedded coords: {token}')
        return None
    return ('whole', token)


def parse_parent(parent, warns):
    m = PARENT_RE.match(parent)
    if not m:
        warns.append(f'[skip] parent name has no "---" separator: {parent}')
        return None
    return m.group(1), m.group(2)


def index_mag_contigs(mags_folder, warns):
    """(sample, cid) -> (mag_path, contig_name); also mag_path -> [contig names]."""
    key2ctg = {}
    mag_files = []
    for f in sorted(os.listdir(mags_folder)):
        if not f.endswith(('.fa', '.fna', '.fasta')):
            continue
        path = os.path.join(mags_folder, f)
        mag_files.append(path)
        with open(path) as fh:
            for line in fh:
                if line.startswith('>'):
                    name = line[1:].split(None, 1)[0]
                    m = MAGCTG_RE.match(name)
                    if m:
                        key2ctg[(m.group(1), m.group(2))] = (path, name)
    return mag_files, key2ctg


def read_fasta(path):
    seqs, name = {}, None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('>'):
                name = line[1:].split(None, 1)[0]
                seqs[name] = []
            elif name is not None:
                seqs[name].append(line)
    return {k: ''.join(v) for k, v in seqs.items()}


def wrap(seq, width=60):
    return '\n'.join(seq[i:i + width] for i in range(0, len(seq), width))


def main():
    ap = argparse.ArgumentParser(
        description='Mask viral regions on host MAG contigs for the pmr_vm merged reference.')
    ap.add_argument('--mags_folder', required=True)
    ap.add_argument('--virus_summaries', required=True,
                    help='file listing per-sample vMAG_summary_info.tsv paths')
    ap.add_argument('--orf_loci', required=True,
                    help='file listing per-sample bin_ORF_loci paths')
    ap.add_argument('--out_genomes_dir', required=True)
    ap.add_argument('--out_intervals_bed', required=True)
    ap.add_argument('--out_excluded_orfs', required=True)
    ap.add_argument('--out_fully_viral_contigs', required=True)
    ap.add_argument('--mags_include_tsv', default=None,
                    help='filtered mMAG tax detail (mMAG_summary_info.filter.tsv); '
                         'if given, only listed genomes enter the masked set')
    ap.add_argument('--vmags_include_tsv', default=None,
                    help='filtered vMAG tax detail (vMAG_summary_info.filter.tsv); '
                         'if given, only scaffolds of listed genomes are masked/dropped')
    args = ap.parse_args()

    # optional quality-gate include lists (absent = keep everything)
    mags_include = None
    if args.mags_include_tsv and os.path.exists(args.mags_include_tsv):
        mags_include = set()
        with open(args.mags_include_tsv) as fh:
            for line in fh:
                line = line.rstrip('\n')
                if line and not line.startswith('genome'):
                    mags_include.add(line.split('\t')[0])
    vmag_scaffolds = None
    if args.vmags_include_tsv and os.path.exists(args.vmags_include_tsv):
        vmag_scaffolds = set()
        with open(args.vmags_include_tsv) as fh:
            header = fh.readline().rstrip('\n').split('\t')
            i_scf = header.index('scaffold') if 'scaffold' in header else None
            for line in fh:
                cols = line.rstrip('\n').split('\t')
                if i_scf is not None and len(cols) > i_scf:
                    for token in cols[i_scf].split(','):
                        token = token.strip()
                        if token:
                            vmag_scaffolds.add(token)

    warns = []

    # 1. collect viral intervals / whole-contig records keyed by (sample, cid)
    intervals = {}   # key -> list[(S, E)]
    whole = set()    # key
    n_scf = 0
    with open(args.virus_summaries) as fh:
        summary_files = [l.strip() for l in fh if l.strip()]
    for sf in summary_files:
        if not os.path.exists(sf):
            warns.append(f'[skip] virus summary missing: {sf}')
            continue
        with open(sf) as fh:
            header = fh.readline().rstrip('\n').split('\t')
            i_scf = header.index('scaffold')
            for line in fh:
                cols = line.rstrip('\n').split('\t')
                if len(cols) <= i_scf:
                    continue
                for token in cols[i_scf].split(','):
                    token = token.strip()
                    if not token:
                        continue
                    if vmag_scaffolds is not None and token not in vmag_scaffolds:
                        continue
                    n_scf += 1
                    rec = classify_token(token, warns)
                    if rec is None:
                        continue
                    key = parse_parent(rec[1], warns)
                    if key is None:
                        continue
                    if rec[0] == 'whole':
                        whole.add(key)
                    else:
                        intervals.setdefault(key, []).append((rec[2], rec[3]))

    # 2. index MAG contigs and assign hits
    mag_files, key2ctg = index_mag_contigs(args.mags_folder, warns)
    if mags_include is not None:
        before = len(mag_files)
        mag_files = [p for p in mag_files
                     if os.path.splitext(os.path.basename(p))[0] in mags_include]
        key2ctg = {k: v for k, v in key2ctg.items() if v[0] in set(mag_files)}
        print(f'[pmr_vm_mask_viral] quality gate: {len(mag_files)}/{before} MAGs pass the include list')
    ctg_intervals = {}   # contig_name -> list[(S, E)]
    ctg_drop = set()
    n_hit_int = n_hit_whole = 0
    for key, spans in intervals.items():
        if key in key2ctg and key not in whole:
            ctg_intervals.setdefault(key2ctg[key][1], []).extend(spans)
            n_hit_int += len(spans)
    for key in whole:
        if key in key2ctg:
            ctg_drop.add(key2ctg[key][1])
            n_hit_whole += 1
    for ctg in ctg_drop:
        ctg_intervals.pop(ctg, None)

    # 3. write masked genomes (touched) / symlinks (untouched)
    os.makedirs(args.out_genomes_dir, exist_ok=True)
    touched = {ctg for ctg in ctg_intervals} | ctg_drop
    n_masked = 0
    for path in mag_files:
        out = os.path.join(args.out_genomes_dir,
                           os.path.splitext(os.path.basename(path))[0] + '.fna')
        seqs = read_fasta(path)
        if not any(ctg in touched for ctg in seqs):
            if os.path.lexists(out):
                os.remove(out)
            # symlink so untouched MAGs track upstream updates (drift-following);
            # build_ref uses find -L so these are still picked up by `find -type f`
            os.symlink(os.path.abspath(path), out)
            continue
        n_masked += 1
        with open(out, 'w') as fh:
            for ctg, seq in seqs.items():
                if ctg in ctg_drop:
                    continue
                if ctg in ctg_intervals:
                    for s, e in sorted(ctg_intervals[ctg]):
                        if s < 1 or s > len(seq):
                            warns.append(f'[skip] interval out of range: {ctg} {s}-{e} (len {len(seq)})')
                            continue
                        if e > len(seq):
                            warns.append(f'[warn] interval end clamped: {ctg} {s}-{e} -> {e}={len(seq)}')
                            e = len(seq)
                        seq = seq[:s - 1] + 'N' * (e - s + 1) + seq[e:]
                fh.write(f'>{ctg}\n{wrap(seq)}\n')

    # 4. host ORF exclusion from per-sample ORF loci
    excluded = []
    n_loci = 0
    with open(args.orf_loci) as fh:
        loci_files = [l.strip() for l in fh if l.strip()]
    for lf in loci_files:
        if not os.path.exists(lf):
            warns.append(f'[skip] ORF loci missing: {lf}')
            continue
        with open(lf) as fh:
            for line in fh:
                line = line.rstrip('\n')
                if not line or line.startswith('ORF_name'):
                    continue
                n_loci += 1
                parts = line.split('\t')
                orf, start, end = parts[0], int(parts[1]), int(parts[2])
                ctg = orf.rsplit('_', 1)[0]
                if ctg in ctg_drop:
                    excluded.append(orf)
                elif ctg in ctg_intervals:
                    for s, e in ctg_intervals[ctg]:
                        if s <= end and e >= start:
                            excluded.append(orf)
                            break

    # 5. audit outputs
    with open(args.out_intervals_bed, 'w') as fh:
        for key, spans in sorted(intervals.items()):
            if key in key2ctg and key not in whole:
                ctg = key2ctg[key][1]
                for s, e in sorted(spans):
                    fh.write(f'{ctg}\t{s - 1}\t{e}\t{key[0]}---{key[1]}\tinterval\n')
        for key in sorted(whole):
            if key in key2ctg:
                fh.write(f'{key2ctg[key][1]}\t0\t.\t{key[0]}---{key[1]}\twhole\n')
    with open(args.out_excluded_orfs, 'w') as fh:
        for orf in excluded:
            fh.write(orf + '\n')
    with open(args.out_fully_viral_contigs, 'w') as fh:
        for ctg in sorted(ctg_drop):
            fh.write(ctg + '\n')

    print(f'[pmr_vm_mask_viral] viral scaffold tokens: {n_scf}; '
          f'intervals on MAG contigs: {n_hit_int}; whole-contig dropped: {n_hit_whole}; '
          f'MAGs rewritten: {n_masked} / {len(mag_files)}; '
          f'host ORFs excluded: {len(excluded)} / {n_loci}')
    for w in warns:
        print(f'[pmr_vm_mask_viral] {w}', file=sys.stderr)


if __name__ == '__main__':
    main()
