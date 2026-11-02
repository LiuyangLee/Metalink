# Common initialization for all Metalink modes
# This file is included first by the main Snakefile
# It parses the input TSV and sets up SAMPLE_INDEX and SAMPLE_PATHS

PROJECT_INDEX = config["name_project"]
RUN_PROJECT_LEVEL = config["run_project"]
RANKS = config.get("taxonomy_rank", ["Phylum", "Class", "Order", "Family", "Genus", "Species"])
READ_TYPE = config.get("read_type", "paired")

ref_file = config["input_tsv"]
SAMPLE_INDEX = []
SAMPLE_PATHS = {}

# Parse input TSV file
try:
    with open(ref_file, "rt", encoding="utf-8-sig") as fh:
        first = next(fh, None)
        if first is None:
            raise ValueError("Input TSV file is empty: {}".format(ref_file))
        if "\t" not in first:
            raise ValueError(
                "Input file does not look like a TSV (no TAB found in the header line): {}\n"
                "  expected header columns like: sample, read_1, read_2, read_s, contig_file, mMAG_dir\n"
                "  first line: {}".format(ref_file, first.strip()[:200]))
        header = first.strip().split("\t")
        if "sample" not in header:
            raise ValueError(
                "Input TSV header must contain a 'sample' column: {}\n"
                "  header: {}".format(ref_file, ",".join(header)))
        col_sample = header.index("sample")
        # required columns depend on read_type and skip_level
        _req = ["read_1", "read_2"] if READ_TYPE in ("paired", "paired-rna") else ["read_s"]
        _skip = config.get("skip_level", 0)
        if _skip >= 2:
            _req.append("contig_file")
        if _skip >= 3:
            _req.append("mMAG_dir")
        _missing = [c for c in _req if c not in header]
        if _missing:
            raise ValueError(
                "Input TSV header is missing required column(s) {} for read_type={} and skip_level={}: {}\n"
                "  header: {}".format(",".join(_missing), READ_TYPE, _skip, ref_file, ",".join(header)))
        _req_idx = [(c, header.index(c)) for c in _req]
        for _ln, line in enumerate(fh, start=2):
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) != len(header):
                raise ValueError(
                    "Input TSV format error in {}: line {} has {} column(s) "
                    "but the header defines {} ({}). Please check the input "
                    "file format (e.g. missing/extra tab-separated fields).\n"
                    "  line {}: {}".format(
                        ref_file, _ln, len(parts), len(header),
                        ",".join(header), _ln, line))
            for _c, _ci in _req_idx:
                if parts[_ci] == "":
                    raise ValueError(
                        "Input TSV format error in {}: line {} has an empty value "
                        "in required column '{}'. Please check the input file.\n"
                        "  line {}: {}".format(ref_file, _ln, _c, _ln, line))
            sname = parts[col_sample]
            if "---" in sname:
                raise ValueError(
                    "Input TSV format error in {}: sample name '{}' (line {}) contains '---', "
                    "which is reserved as an internal delimiter (Metalink builds internal IDs "
                    "as sample---contig / genome---ORF). Please rename this sample.\n"
                    "  line {}: {}".format(ref_file, sname, _ln, _ln, line))
            SAMPLE_INDEX.append(sname)
            SAMPLE_PATHS[sname] = {
                header[i]: parts[i]
                for i in range(len(header))
                if i != col_sample
            }
except FileNotFoundError:
    print("Input TSV file is not found: {}".format(ref_file))
    raise

# ---------------------------------------------------------------------------
# Reserved delimiter check: '---' is used internally as {sample}---{contig}
# and {genome}---{ORF}. Sample names are checked above; here check user-provided
# FASTA headers (contig_file at skip>=2, mMAG_dir genomes at skip>=3) and mMAG
# file names. Fail fast with the exact location instead of cryptic split
# errors deep inside the merge steps.
# ---------------------------------------------------------------------------
import subprocess as _sp
import re as _re


def _triple_dash_headers(fpath):
    """Return up to 3 FASTA headers containing '---' (streaming, stops early)."""
    cmd = ['zgrep' if fpath.endswith('.gz') else 'grep', '-m', '3', '^>.*---', fpath]
    r = _sp.run(cmd, capture_output=True, text=True)
    return [l for l in r.stdout.splitlines() if l][:3]


def _check_fasta_headers(fpath):
    hits = _triple_dash_headers(fpath)
    if hits:
        raise ValueError(
            "FASTA file {} contains header(s) with '---' (reserved as an internal "
            "delimiter; Metalink builds internal IDs as sample---contig / genome---ORF):\n"
            "  {}".format(fpath, "\n  ".join(hits)) +
            "\nPlease rename these sequence headers (remove '---').")


_skip = config.get("skip_level", 0)
if _skip >= 2:
    for _s in SAMPLE_INDEX:
        _cf = SAMPLE_PATHS[_s].get("contig_file")
        if _cf:
            _check_fasta_headers(_cf)
if _skip >= 3:
    import glob as _glob
    # pmr pools all mMAG_dirs into ONE project-level mapping reference, so
    # contig IDs must also be unique ACROSS directories (except a genome
    # shared by several samples/directories with identical name and content).
    _mode = os.path.basename(str(config.get('workflow_target', ''))).replace('targets.', '').replace('.smk', '')
    _cid_owner = {}     # contig id -> genome basename (project-wide, pmr only)
    _genome_ids = {}    # genome basename -> frozenset of its contig ids (pmr only)
    for _s in SAMPLE_INDEX:
        _md = SAMPLE_PATHS[_s].get("mMAG_dir")
        if not _md:
            continue
        _seen_ids = {}
        for _f in sorted(_glob.glob(os.path.join(_md, "*"))):
            if not os.path.isfile(_f):
                continue
            _base = os.path.basename(_f)
            if "---" in _base:
                raise ValueError(
                    "mMAG file name '{}' in mMAG_dir '{}' (sample {}) contains '---', "
                    "which is reserved as an internal delimiter (internal IDs are built "
                    "as genome---ORF). Please rename this file.".format(_base, _md, _s))
            _check_fasta_headers(_f)
            # duplicate contig IDs across genomes in one mMAG_dir would collide
            # when the genomes are concatenated into a single mapping reference
            _genome = _re.sub(r'\.(fa|fna|fasta)(\.gz)?$', '', _base)
            _ids = set()
            _r = _sp.run(
                ['zgrep' if _f.endswith('.gz') else 'grep', '-h', '^>', _f],
                capture_output=True, text=True)
            for _h in _r.stdout.splitlines():
                _cid = _h[1:].split()[0] if len(_h) > 1 else _h
                _ids.add(_cid)
                if _cid in _seen_ids:
                    raise ValueError(
                        "Duplicate contig ID '{}' in mMAG_dir '{}' (sample {}): "
                        "found in both '{}' and '{}'. Contig IDs must be unique across "
                        "the genomes within one mMAG_dir, because the genomes are "
                        "concatenated into a single mapping reference. Please rename "
                        "the contig headers (e.g. prefix them with the genome name)."
                        .format(_cid, _md, _s, _seen_ids[_cid], _f))
                _seen_ids[_cid] = _f
                if _mode == 'pmr':
                    _owner = _cid_owner.get(_cid)
                    if _owner is not None and _owner != _genome:
                        raise ValueError(
                            "Duplicate contig ID '{}' across mMAG_dirs: found in genome "
                            "'{}' and genome '{}' (sample {}). pmr pools all mMAG_dirs "
                            "into one project-level mapping reference, so contig IDs "
                            "must be unique per genome across ALL directories. Please "
                            "rename the contig headers (e.g. prefix them with the "
                            "genome name). (Sharing the SAME genome file across samples "
                            "is fine.)"
                            .format(_cid, _owner, _genome, _s))
                    _cid_owner[_cid] = _genome
            if _mode == 'pmr':
                if _genome in _genome_ids and _genome_ids[_genome] != _ids:
                    raise ValueError(
                        "Genome name '{}' is used with different contig IDs in different "
                        "mMAG_dirs (e.g. sample {}). In pmr all directories are pooled "
                        "into one reference, so a genome name must always refer to the "
                        "same sequence set. Please rename one of them."
                        .format(_genome, _s))
                _genome_ids[_genome] = _ids

_KEEP_INTERMEDIATE_DEFAULTS = {'bam': True, 'drep_data': False, 'vs2_iter0': False, 'assembly': False}

def keep_intermediate(category):
    """Return True if the given intermediate-file category should be kept.

    Config value may be a dict (per-category), a legacy bool (applies to all
    categories), or absent (class defaults)."""
    ki = config.get('keep_intermediate', {})
    if isinstance(ki, bool):
        return ki
    return dict(ki).get(category, _KEEP_INTERMEDIATE_DEFAULTS[category])
