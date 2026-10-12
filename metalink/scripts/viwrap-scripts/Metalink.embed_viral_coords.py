#!/usr/bin/env python3
"""Embed parent-contig coordinates of VS2/CheckV viral fragments into sequence IDs.

Renames fragment scaffolds produced by the VirSorter2 (+VIBRANT) branch of
Metalink.01_identify.py:

    {parent}||0_partial     -> {parent}|vs2region_{S}_{E}
    {parent}||0_partial_1   -> {parent}|vs2region_{S}_{E}
    {parent}||full_1        -> {parent}|vs2region_{S}_{E}
    {parent}||lt2gene_1     -> {parent}|vs2region_{S}_{E}

S/E are 1-based inclusive coordinates on the parent contig, composed from
  - VS2 vs_vb_result/final-viral-boundary.tsv (trim_bp_start/end; the VS2-level
    cut, only meaningful for ||{i}_partial rows), and
  - CheckV CheckV_result/combined.fna header marks "start-end/length" (the
    CheckV-level cut, relative to its input sequence).
Every rename is asserted by exact substring matching of the fragment sequence
against the parent contig (exactly one hit at the composed position); on
assertion failure the old name is kept and a warning is logged.

Whole-contig sequences (||full / ||lt2gene without trailing _M) and geNomad
products (|provirus_S_E, plain names) are left untouched.

Rewritten in place inside --virsorter_outdir:
    final_vs2_virus.fna                 (record headers)
    final_vs2_virus.ffn / .faa          (protein-id prefixes in headers)
    final_vs2_virus.annotation.txt      (protein and scaffold columns)
    scf2lytic_or_lyso.summary.txt       (scaffold column)
Audit outputs:
    rename_map.tsv                      (old/new names, coords, category, status)
    stdout log                          (per-category counts, failures)

Only stdlib is used so any python interpreter can run this script.
"""

import argparse
import gzip
import os
import re
import sys
import tempfile

# {base} = name with the CheckV _M ordinal stripped; M = CheckV region ordinal
FRAG_RE = re.compile(r'^(?P<base>.+\|\|(?:full|lt2gene|\d+_partial))(?:_(?P<M>\d+))?$')
# CheckV combined.fna header mark: "start-end/input_length" (1-based)
HDR_MARK_RE = re.compile(r'^(\d+)-(\d+)/(\d+)$')
TAG = 'vs2region'


def open_maybe_gz(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path, 'rt')


def read_fasta(path):
    """Return {name: (desc, seq)}; name = first token of the header."""
    seqs = {}
    name = None
    with open_maybe_gz(path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('>'):
                parts = line[1:].split(None, 1)
                name = parts[0]
                seqs[name] = [parts[1] if len(parts) > 1 else '', []]
            elif name is not None:
                seqs[name][1].append(line)
    return {k: (v[0], ''.join(v[1])) for k, v in seqs.items()}


def read_boundary(path):
    """Return {seqname_new: (trim_bp_start, trim_bp_end)}."""
    table = {}
    with open_maybe_gz(path) as fh:
        header = fh.readline().rstrip('\n').split('\t')
        i_s, i_e = header.index('trim_bp_start'), header.index('trim_bp_end')
        i_new = header.index('seqname_new')
        for line in fh:
            c = line.rstrip('\n').split('\t')
            if len(c) > i_new and c[i_new]:
                table[c[i_new]] = (int(c[i_s]), int(c[i_e]))
    return table


def compose_coords(name, desc, boundary, warns):
    """Compose parent-contig (S, E) for a fragment scaffold name.

    Returns (parent, S, E, category) or None when the name is not a
    renameable fragment (whole-contig sequence or unrecognized pattern).
    """
    m = FRAG_RE.match(name)
    if not m:
        return None
    base, ord_m = m.group('base'), m.group('M')
    parent = base.split('||', 1)[0]

    if ord_m is None:
        if not base.endswith('_partial'):
            return None  # ||full / ||lt2gene without _M: whole contig
        if base not in boundary:
            warns.append(f'[skip] {name}: not found in final-viral-boundary.tsv')
            return None
        ts, te = boundary[base]
        return parent, ts, te, 'vs2_only'

    mark = HDR_MARK_RE.match(desc)
    if not mark:
        warns.append(f'[skip] {name}: no "start-end/len" mark in combined.fna header')
        return None
    cs, ce = int(mark.group(1)), int(mark.group(2))
    if base.endswith('_partial'):
        if base not in boundary:
            warns.append(f'[skip] {name}: {base} not found in final-viral-boundary.tsv')
            return None
        ts, _ = boundary[base]
        return parent, ts + cs - 1, ts + ce - 1, 'vs2+checkv'
    # ||full_M / ||lt2gene_M: the VS2-level sequence is the whole contig,
    # so CheckV coords are already parent-contig coords
    return parent, cs, ce, 'checkv_on_full'


def build_rename_map(metagenome_fna, boundary_tsv, combined_fna):
    """Return (rename_map, rows, warns).

    rename_map: {old_name: new_name}; rows: audit-table rows (list of str).
    """
    boundary = read_boundary(boundary_tsv)
    fragments = read_fasta(combined_fna)
    contigs = None  # lazy: only loaded if at least one candidate exists
    rename_map, rows, warns = {}, [], []

    for name, (desc, seq) in sorted(fragments.items()):
        composed = compose_coords(name, desc, boundary, warns)
        if composed is None:
            continue
        parent, s, e, category = composed
        if contigs is None:
            contigs = read_fasta(metagenome_fna)
        status = ''
        new_name = f'{parent}|{TAG}_{s}_{e}'
        if parent not in contigs:
            status = f'skipped: parent {parent} not in metagenome'
        else:
            pseq = contigs[parent][1]
            hits = pseq.count(seq)
            off = pseq.find(seq)
            if hits == 1 and off == s - 1 and off + len(seq) == e:
                status = 'renamed'
                rename_map[name] = new_name
            else:
                pos = 'miss' if off < 0 else f'{off + 1}-{off + len(seq)}'
                status = f'skipped: assertion failed (hits={hits}, match={pos}, composed={s}-{e})'
        rows.append((name, rename_map.get(name, ''), parent, str(s), str(e),
                     category, status))
        if status.startswith('skipped'):
            warns.append(f'[{status}] {name}')
    return rename_map, rows, warns


def _protein_rename(pid, rename_map):
    """Map a protein id ({scaffold}_{ordinal}) through the scaffold rename map."""
    for old in sorted(rename_map, key=len, reverse=True):
        if pid.startswith(old + '_'):
            return rename_map[old] + pid[len(old):]
    return pid


def _unique_tmp(path):
    """Unique temp file next to the target, so concurrent invocations never
    share a temp path (a shared fixed name led to torn interleaved writes)."""
    fd, tmp = tempfile.mkstemp(prefix=os.path.basename(path) + '.', suffix='.tmp',
                               dir=os.path.dirname(os.path.abspath(path)))
    os.close(fd)
    return tmp


def _rewrite_fasta_to_tmp(path, id_mapper):
    tmp = _unique_tmp(path)
    with open(path) as fin, open(tmp, 'w') as fout:
        for line in fin:
            if line.startswith('>'):
                head, sep, rest = line[1:].rstrip('\n').partition(' ')
                fout.write('>' + id_mapper(head) + (sep + rest if sep else '') + '\n')
            else:
                fout.write(line)
    return tmp


def _rewrite_annotation_to_tmp(path, rename_map):
    tmp = _unique_tmp(path)
    with open(path) as fin, open(tmp, 'w') as fout:
        for line in fin:
            parts = line.rstrip('\n').split('\t')
            if parts[0] != 'protein':
                if len(parts) >= 1:
                    parts[0] = _protein_rename(parts[0], rename_map)
                if len(parts) >= 2:
                    parts[1] = rename_map.get(parts[1], parts[1])
            fout.write('\t'.join(parts) + '\n')
    return tmp


def _rewrite_scf_summary_to_tmp(path, rename_map):
    tmp = _unique_tmp(path)
    with open(path) as fin, open(tmp, 'w') as fout:
        for line in fin:
            parts = line.rstrip('\n').split('\t')
            if parts[0] != 'scaffold':
                parts[0] = rename_map.get(parts[0], parts[0])
            fout.write('\t'.join(parts) + '\n')
    return tmp


def apply_rename(virsorter_outdir, rename_map):
    """Two-phase commit: write every rewritten file to its own unique temp file
    first, then atomically replace all of them. A kill during phase 1 leaves the
    originals untouched; during phase 2 it can only leave complete files (some
    old, some new), which a re-run reconciles because the rename is idempotent."""
    if not rename_map:
        return
    pending = []  # (tmp_path, final_path)
    try:
        fna = os.path.join(virsorter_outdir, 'final_vs2_virus.fna')
        if os.path.exists(fna) and os.path.getsize(fna) > 0:
            pending.append((_rewrite_fasta_to_tmp(fna, lambda h: rename_map.get(h, h)), fna))
        for ext in ('ffn', 'faa'):
            p = os.path.join(virsorter_outdir, f'final_vs2_virus.{ext}')
            if os.path.exists(p) and os.path.getsize(p) > 0:
                pending.append((_rewrite_fasta_to_tmp(p, lambda h: _protein_rename(h, rename_map)), p))
        ann = os.path.join(virsorter_outdir, 'final_vs2_virus.annotation.txt')
        if os.path.exists(ann) and os.path.getsize(ann) > 0:
            pending.append((_rewrite_annotation_to_tmp(ann, rename_map), ann))
        summ = os.path.join(virsorter_outdir, 'scf2lytic_or_lyso.summary.txt')
        if os.path.exists(summ) and os.path.getsize(summ) > 0:
            pending.append((_rewrite_scf_summary_to_tmp(summ, rename_map), summ))
    except BaseException:
        for tmp, _ in pending:
            if os.path.exists(tmp):
                os.unlink(tmp)
        raise
    for tmp, final in pending:
        os.replace(tmp, final)


def main():
    ap = argparse.ArgumentParser(
        description='Embed parent-contig coords of VS2/CheckV viral fragments into sequence IDs (|vs2region_S_E).')
    ap.add_argument('--input_metagenome', required=True,
                    help='metagenome contigs fasta used for substring-match assertion')
    ap.add_argument('--virsorter_outdir', required=True,
                    help='vs_vb output dir containing vs_vb_result/, CheckV_result/ and final_vs2_virus.*')
    ap.add_argument('--fail_loud', action='store_true',
                    help='exit non-zero if any fragment fails the substring-match assertion')
    args = ap.parse_args()

    boundary_tsv = os.path.join(args.virsorter_outdir, 'vs_vb_result', 'final-viral-boundary.tsv')
    combined_fna = os.path.join(args.virsorter_outdir, 'CheckV_result', 'combined.fna')
    for p in (args.input_metagenome, boundary_tsv, combined_fna):
        if not os.path.exists(p):
            sys.exit(f'[embed_viral_coords] required input missing: {p}')

    rename_map, rows, warns = build_rename_map(
        args.input_metagenome, boundary_tsv, combined_fna)
    apply_rename(args.virsorter_outdir, rename_map)

    map_file = os.path.join(args.virsorter_outdir, 'rename_map.tsv')
    tmp_map = _unique_tmp(map_file)
    with open(tmp_map, 'w') as fh:
        fh.write('old_name\tnew_name\tparent\tS\tE\tcategory\tstatus\n')
        for row in rows:
            fh.write('\t'.join(row) + '\n')
    os.replace(tmp_map, map_file)

    n_cat = {}
    for row in rows:
        n_cat[row[5]] = n_cat.get(row[5], 0) + 1
    print(f'[embed_viral_coords] fragments: {len(rows)} '
          f'({", ".join(f"{k}={v}" for k, v in sorted(n_cat.items()))}); '
          f'renamed: {len(rename_map)}; skipped: {len(rows) - len(rename_map)}')
    for w in warns:
        print(f'[embed_viral_coords] {w}', file=sys.stderr)
    print(f'[embed_viral_coords] rename map written to {map_file}')

    if args.fail_loud and len(rename_map) != len(rows):
        sys.exit('[embed_viral_coords] assertion failures encountered (--fail_loud)')


if __name__ == '__main__':
    main()
