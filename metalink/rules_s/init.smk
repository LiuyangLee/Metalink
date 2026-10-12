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
            SAMPLE_INDEX.append(sname)
            SAMPLE_PATHS[sname] = {
                header[i]: parts[i]
                for i in range(len(header))
                if i != col_sample
            }
except FileNotFoundError:
    print("Input TSV file is not found: {}".format(ref_file))
    raise

_KEEP_INTERMEDIATE_DEFAULTS = {'bam': True, 'drep_data': False, 'vs2_iter0': False, 'assembly': False}

def keep_intermediate(category):
    """Return True if the given intermediate-file category should be kept.

    Config value may be a dict (per-category), a legacy bool (applies to all
    categories), or absent (class defaults)."""
    ki = config.get('keep_intermediate', {})
    if isinstance(ki, bool):
        return ki
    return dict(ki).get(category, _KEEP_INTERMEDIATE_DEFAULTS[category])
