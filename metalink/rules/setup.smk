import sys
import os
import hashlib
import shutil
import shlex
import yaml


# ---------------------------------------------------------------------------
# Metalink database setup workflow
#
# Downloads reference archives, validates checksums, extracts them, and
# optionally installs dependencies.  Database versions are taken from the
# bundle manifest (metalink/config/db_bundle.yaml) so that every Metalink
# release is paired with a well-defined, reproducible set of databases.
#
# Usage:
#   metalink download-database -d /path/to/DBDIR
#   metalink download-database -c /path/to/DBDIR/config/config.download.yaml
#
# Or directly with snakemake:
#   snakemake --snakefile metalink/rules/setup.smk --directory /path/to/DBDIR
# ---------------------------------------------------------------------------


def md5(fname):
    """Return the MD5 hex digest of a file, or None if it does not exist."""
    if not os.path.exists(fname):
        return None
    hash_md5 = hashlib.md5()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


# ---------------------------------------------------------------------------
# Load the database bundle manifest
# ---------------------------------------------------------------------------

def _resolve_bundle_path():
    """Return the path to the database bundle manifest.

    Priority:
        1. DB_bundle config value passed from the CLI
        2. The directory of the current Snakefile (more reliable than __file__
           inside Snakemake)
        3. The directory of this source file as a last resort
    """
    if 'DB_bundle' in config:
        return config['DB_bundle']
    try:
        # workflow.current_basedir is the directory of the current Snakefile.
        base = workflow.current_basedir
    except (NameError, AttributeError):
        try:
            base = os.path.dirname(os.path.abspath(workflow.snakefile))
        except (NameError, AttributeError):
            base = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(base, '..', 'config', 'db_bundle.yaml')


_BUNDLE_PATH = _resolve_bundle_path()
with open(_BUNDLE_PATH) as _bundle_f:
    DB_BUNDLE = yaml.safe_load(_bundle_f)

BUNDLE_VERSION = str(DB_BUNDLE.get('bundle_version', 'unknown'))


# Conda environment names come from the conda_envs mapping.  Priority:
#   1. the conda_envs section of the config file passed with --configfile
#      (generated configs carry the general-template mapping, possibly
#      renamed by the user)
#   2. the general config template next to the bundle manifest, so that a
#      renamed environment only needs to be updated in one place
#   3. the key itself
try:
    with open(os.path.join(os.path.dirname(_BUNDLE_PATH),
                           'template.config.general.yaml')) as _tpl_f:
        _CONDA_ENVS = (yaml.safe_load(_tpl_f) or {}).get('conda_envs', {}) or {}
except OSError:
    _CONDA_ENVS = {}

_CFG_CONDA_ENVS = config.get('conda_envs', {}) or {}


def _ce(key):
    """Resolve a conda env key (e.g. 'metalink-gtdbtk') to the env name."""
    name = _CFG_CONDA_ENVS.get(key) or _CONDA_ENVS.get(key)
    return str(name) if name else key


DB_DIR = os.path.abspath(os.getcwd())


# Size variant selection for Metalink-built kraken2 databases.
KRaken2_TAXONOMY_DB_SIZE = config.get('Kraken2_taxonomy_db_size', 'kraken2_raw')
KRaken2_FUNCTION_DB_SIZE = config.get('Kraken2_function_db_size', 'kraken2_raw')


# ---------------------------------------------------------------------------
# Database catalog derived from the bundle manifest
# ---------------------------------------------------------------------------

# URL config key prefixes for Metalink-built databases (must match the
# CLI/config keys used in earlier Metalink releases).
_METALINK_URL_PREFIX = {
    'kraken2_taxonomy': 'Kraken2_taxonomy',
    'kraken2_function': 'Kraken2_function',
    'cat_pack_db': 'Cat_pack',
}


def _resolve_metalink_url(db_name, archive, size=None):
    """Return the download URL for a Metalink-built database archive.

    Priority:
        1. Per-size config override: <Prefix>_db_url_<size>
        2. Generic config override: <Prefix>_db_url
        3. URL listed in the bundle manifest (per size variant when sizes
           exist, otherwise at the entry level)
    """
    prefix = _METALINK_URL_PREFIX.get(db_name)
    if prefix:
        url = config.get('{}_db_url_{}'.format(prefix, size), '') if size else ''
        if not url:
            url = config.get('{}_db_url'.format(prefix), '')
        if url:
            return url

    info = DB_BUNDLE['metalink_built_databases'][db_name]
    url = ''
    if size and size in info.get('size_variants', {}):
        url = info['size_variants'][size].get('url', '')
    if not url:
        url = info.get('url', '')
    if url:
        return url

    raise ValueError(
        "No download URL for database '{}'{}: add 'url' to its entry in {} "
        "or set '{}_db_url' in the config.".format(
            db_name,
            ' (size {})'.format(size) if size else '',
            _BUNDLE_PATH,
            prefix or db_name))


METALINK_DATABASES = {}
for _db_name, _db_info in DB_BUNDLE['metalink_built_databases'].items():
    _version = _db_info['version']
    _size = None
    if 'size_variants' in _db_info:
        _size = KRaken2_TAXONOMY_DB_SIZE if _db_name == 'kraken2_taxonomy' else KRaken2_FUNCTION_DB_SIZE
    _archive = _db_info['archive_prefix'].format(size=_size) if _size else _db_info['archive']
    _url = _resolve_metalink_url(_db_name, _archive, size=_size)
    _extract_dir = os.path.join(_db_info['dir'], _version)
    _markers = [m.format(size=_size) for m in _db_info['markers']]
    # Optional md5 from the manifest: per size variant when sizes exist,
    # otherwise at the entry level.  Empty/missing means no verification.
    # Coerce to str because an all-digit md5 parses as an int in YAML.
    if _size:
        _md5 = _db_info.get('size_variants', {}).get(_size, {}).get('md5')
    else:
        _md5 = _db_info.get('md5')
    _md5 = str(_md5).strip() if _md5 is not None else None
    _md5 = _md5 or None

    METALINK_DATABASES[_db_name] = {
        'archive': _archive,
        'dir': _extract_dir,
        'version': _version,
        'size': _size,
        'url_override': _url,
        'markers': _markers,
        'md5': _md5,
    }


THIRD_PARTY_DATABASES = {}
for _db_name, _db_info in DB_BUNDLE['third_party_databases'].items():
    THIRD_PARTY_DATABASES[_db_name] = {
        'target_dir': _db_info['target_dir'],
        'version': _db_info['version'],
        'markers': _db_info.get('markers', []),
    }


ALL_DB_NAMES = list(METALINK_DATABASES.keys()) + list(THIRD_PARTY_DATABASES.keys())


# ---------------------------------------------------------------------------
# Decide which databases to install
# ---------------------------------------------------------------------------

def parse_database_selection(selection):
    """Return a set of database names to install."""
    if selection is None or (isinstance(selection, str) and selection.strip() == ''):
        return set(ALL_DB_NAMES)

    # Config files carry a YAML list; --config Databases=... passes a
    # comma-separated string.
    tokens = selection.split(',') if isinstance(selection, str) else list(selection)

    selected = set()
    for token in tokens:
        token = str(token).strip()
        if not token:
            continue
        if token in ALL_DB_NAMES:
            selected.add(token)
        else:
            raise ValueError(
                f"Unknown database '{token}'. Choose from: "
                f"{', '.join(ALL_DB_NAMES)}"
            )
    return selected


SELECTED_DBS = parse_database_selection(config.get('Databases', ''))


def _setup_marker_name(db_name):
    """Return the versioned marker filename for a database."""
    info = METALINK_DATABASES.get(db_name)
    if info is not None:
        version = info['version']
        size = info.get('size')
        if size:
            return '{}_{}_{}.setup.mk'.format(db_name, version, size)
        return '{}_{}.setup.mk'.format(db_name, version)

    info = THIRD_PARTY_DATABASES.get(db_name, {})
    version = info.get('version')
    if version:
        return '{}_{}.setup.mk'.format(db_name, version)
    return '{}.setup.mk'.format(db_name)


# ---------------------------------------------------------------------------
# Top-level rule
# ---------------------------------------------------------------------------

setup_inputs = [os.path.join(DB_DIR, 'markers', _setup_marker_name(name)) for name in SELECTED_DBS]


rule all:
    input: setup_inputs


# ---------------------------------------------------------------------------
# Metalink-built database rules
# ---------------------------------------------------------------------------

for _db_name, _db_info in METALINK_DATABASES.items():
    _archive = _db_info['archive']
    _url = _db_info['url_override']
    _marker_name = _setup_marker_name(_db_name)

    rule:
        name: f"download_{_db_name}"
        output:
            os.path.join(DB_DIR, 'downloads', _archive)
        params:
            url=_url
        log:
            os.path.join(DB_DIR, 'logs', f"download_{_db_name}.log")
        shell:
            """
            (
            success=0
            mkdir -p $(dirname {output})
            # Progress is streamed into the log file (wget dot:giga = 1 MB per
            # dot); follow it live with: tail -f {log}
            for i in $(seq 1 5); do
                if command -v wget >/dev/null 2>&1; then
                    if wget --tries 2 --retry-connrefused --waitretry=60 --timeout=60 --progress=dot:giga -O {output} "{params.url}"; then
                        success=1
                        break
                    fi
                elif command -v curl >/dev/null 2>&1; then
                    if curl -L -C - -o {output} "{params.url}"; then
                        success=1
                        break
                    fi
                else
                    echo "Neither wget nor curl is available" >&2
                    exit 1
                fi
                echo "Download attempt $i/5 failed for {output}; retrying in 5s..." >&2
                sleep 5
            done
            if [ $success -ne 1 ]; then
                echo "ERROR: failed to download {params.url}" >&2
                exit 1
            fi
            # Fetch the sibling checksum file (<url>.md5); it is required so the
            # archive integrity can be verified before extraction.
            if command -v wget >/dev/null 2>&1; then
                wget --tries 2 --timeout=60 -q -O {output}.md5 "{params.url}.md5"
            else
                curl -fSL --max-time 60 -o {output}.md5 "{params.url}.md5"
            fi
            if [ ! -s {output}.md5 ]; then
                echo "ERROR: failed to download checksum file {params.url}.md5" >&2
                exit 1
            fi
            ) > {log} 2>&1
            """

    rule:
        name: f"setup_{_db_name}"
        input:
            os.path.join(DB_DIR, 'downloads', _archive)
        output:
            touch(os.path.join(DB_DIR, 'markers', _marker_name))
        params:
            db_name=_db_name
        log:
            os.path.join(DB_DIR, 'logs', f"setup_{_db_name}.log")
        run:
            _run_extract_and_verify(str(log), params.db_name)


def _run_extract_and_verify(log, db_name):
    """Run _extract_and_verify with stdout/stderr redirected to log."""
    os.makedirs(os.path.dirname(log), exist_ok=True)
    with open(log, 'w') as _log_f:
        old_stdout, old_stderr = sys.stdout, sys.stderr
        sys.stdout = sys.stderr = _log_f
        try:
            _extract_and_verify(db_name)
        finally:
            sys.stdout, sys.stderr = old_stdout, old_stderr


def _record_downloaded_size(db_name):
    """Record the downloaded size variant so later runs can auto-detect it."""
    info = METALINK_DATABASES.get(db_name, {})
    size = info.get('size')
    if not size:
        return
    size_key = '{}_db_size'.format(db_name)
    size_file = os.path.join(DB_DIR, '.metalink_db_sizes.yaml')
    sizes = {}
    if os.path.exists(size_file):
        try:
            with open(size_file) as f:
                for line in f:
                    line = line.strip()
                    if ':' in line and not line.startswith('#'):
                        k, v = line.split(':', 1)
                        sizes[k.strip()] = v.strip()
        except Exception:
            pass
    sizes[size_key] = size
    with open(size_file, 'w') as f:
        f.write('# Auto-generated by Metalink download-database.\n')
        for k, v in sorted(sizes.items()):
            f.write('{}: {}\n'.format(k, v))


def _record_installed_bundle():
    """Record the installed bundle version in DB_DIR."""
    bundle_file = os.path.join(DB_DIR, '.metalink_db_bundle.yaml')
    bundle_info = {'bundle_version': BUNDLE_VERSION}
    try:
        with open(bundle_file, 'w') as f:
            yaml.safe_dump(bundle_info, f, default_flow_style=False, sort_keys=True)
    except Exception:
        pass


def _extract_and_verify(db_name):
    """Extract the archive for db_name and verify markers/checksums."""
    info = METALINK_DATABASES[db_name]
    archive = os.path.join(DB_DIR, 'downloads', info['archive'])
    extract_dir = os.path.join(DB_DIR, info['dir'])
    expected_md5 = info.get('md5')

    # Checksum source: the manifest md5 field wins; otherwise the sibling
    # .md5 file fetched alongside the archive (standard md5sum output format).
    if expected_md5 is None:
        md5_file = archive + '.md5'
        if os.path.exists(md5_file):
            with open(md5_file) as f:
                parts = f.read().split()
            if parts:
                expected_md5 = parts[0].strip().lower() or None

    if expected_md5 is not None:
        actual_md5 = md5(archive)
        if actual_md5 != expected_md5:
            print(
                'ERROR: invalid checksum for {}: expected {}, got {}'.format(
                    info['archive'], expected_md5, actual_md5)
            )
            sys.exit(1)

    # Extract into a temporary directory first so we can handle all archive
    # layouts: (a) files directly at the archive root, (b) a single top-level
    # directory containing the payload, or (c) the archive carries the full
    # DBDIR-relative path (e.g. kraken2_taxonomy_db/v20250930/<size>/...).
    tmp_extract_dir = os.path.join(DB_DIR, 'downloads', db_name + '.tmp')
    if os.path.exists(tmp_extract_dir):
        shutil.rmtree(tmp_extract_dir)
    os.makedirs(tmp_extract_dir, exist_ok=True)
    shell("tar -xzf {} -C {}".format(
        shlex.quote(archive),
        shlex.quote(tmp_extract_dir),
    ))

    entries = os.listdir(tmp_extract_dir)
    if os.path.isdir(os.path.join(tmp_extract_dir, info['dir'])):
        # Layout (c): strip the DBDIR-relative prefix carried by the archive.
        src_dir = os.path.join(tmp_extract_dir, info['dir'])
    elif len(entries) == 1 and os.path.isdir(os.path.join(tmp_extract_dir, entries[0])):
        # Layout (b): single top-level directory; move its contents to the target.
        src_dir = os.path.join(tmp_extract_dir, entries[0])
    else:
        # Layout (a): files directly at the archive root.
        src_dir = tmp_extract_dir

    os.makedirs(extract_dir, exist_ok=True)
    for entry in os.listdir(src_dir):
        src = os.path.join(src_dir, entry)
        dst = os.path.join(extract_dir, entry)
        if os.path.exists(dst):
            shutil.rmtree(dst) if os.path.isdir(dst) else os.remove(dst)
        shutil.move(src, dst)

    shutil.rmtree(tmp_extract_dir)

    for marker in info['markers']:
        marker_file = os.path.join(extract_dir, marker)
        if not os.path.exists(marker_file):
            print(
                'ERROR: expected marker missing after extraction: {}'.format(marker_file)
            )
            sys.exit(1)

    # Remember which size variant was installed so the main workflow can use it
    # automatically without requiring the user to edit the config.
    _record_downloaded_size(db_name)

    # Record the bundle version that was just installed.
    _record_installed_bundle()

    # Remove the archive and its checksum file to save space.
    os.remove(archive)
    if os.path.exists(archive + '.md5'):
        os.remove(archive + '.md5')


# ---------------------------------------------------------------------------
# Third-party database rules
# ---------------------------------------------------------------------------

EGGNOG_CFG = DB_BUNDLE['third_party_databases']['eggnog']
EGGNOG_VERSION = EGGNOG_CFG['version']
EGGNOG_TARGET_DIR = os.path.join(DB_DIR, EGGNOG_CFG['target_dir'])
# Resolve {version} now so Snakemake does not treat it as a wildcard.
EGGNOG_BASE_URL = EGGNOG_CFG['base_url'].format(version=EGGNOG_VERSION)


rule download_eggnog:
    """Download eggNOG-mapper reference files."""
    output:
        db=os.path.join(EGGNOG_TARGET_DIR, 'eggnog.db.gz'),
        dmnd=os.path.join(EGGNOG_TARGET_DIR, 'eggnog_proteins.dmnd.gz'),
        taxa=os.path.join(EGGNOG_TARGET_DIR, 'eggnog.taxa.tar.gz'),
        mmseqs=os.path.join(EGGNOG_TARGET_DIR, 'mmseqs.tar.gz'),
    params:
        target_dir=EGGNOG_TARGET_DIR,
        base_url=EGGNOG_BASE_URL,
        mmseqs_url=EGGNOG_BASE_URL.replace('eggnog.embl.de', 'eggnog5.embl.de'),
    log:
        os.path.join(DB_DIR, 'logs', 'download_eggnog.log')
    shell:
        """
        (
        mkdir -p {params.target_dir}
        wget -c {params.base_url}eggnog.db.gz -O {output.db}
        wget -c {params.base_url}eggnog_proteins.dmnd.gz -O {output.dmnd}
        wget -c {params.base_url}eggnog.taxa.tar.gz -O {output.taxa}
        wget -c {params.mmseqs_url}mmseqs.tar.gz -O {output.mmseqs}
        ) > {log} 2>&1
        """


rule setup_eggnog:
    """Decompress and arrange eggNOG-mapper reference files."""
    input:
        rules.download_eggnog.output
    output:
        touch(os.path.join(DB_DIR, 'markers', _setup_marker_name('eggnog')))
    params:
        target_dir=EGGNOG_TARGET_DIR,
    threads: 20
    conda:
        _ce('metalink-annotation')
    log:
        os.path.join(DB_DIR, 'logs', 'setup_eggnog.log')
    shell:
        """
        (
        cd {params.target_dir}
        pigz -p {threads} -d eggnog.db.gz
        pigz -p {threads} -d eggnog_proteins.dmnd.gz
        tar -zxf eggnog.taxa.tar.gz
        tar -zxf mmseqs.tar.gz
        ) > {log} 2>&1
        """


GENOMAD_CFG = DB_BUNDLE['third_party_databases']['genomad']
GENOMAD_VERSION = GENOMAD_CFG['version']
GENOMAD_TARGET_DIR = os.path.join(DB_DIR, GENOMAD_CFG['target_dir'])
GENOMAD_SOURCE_DIR = os.path.join(DB_DIR, 'downloads', GENOMAD_CFG['download_subdir'])


rule download_and_setup_genomad:
    """Download the geNomad database."""
    output:
        touch(os.path.join(DB_DIR, 'markers', _setup_marker_name('genomad')))
    params:
        source_dir=GENOMAD_SOURCE_DIR,
        target_dir=GENOMAD_TARGET_DIR,
    conda:
        _ce('metalink-genomad')
    log:
        os.path.join(DB_DIR, 'logs', 'setup_genomad.log')
    shell:
        """
        (
        rm -rf {params.source_dir}
        genomad download-database {DB_DIR}/downloads
        rm -rf {params.target_dir}
        mkdir -p $(dirname {params.target_dir})
        mv {params.source_dir} {params.target_dir}
        ) > {log} 2>&1
        """


CHECKV_CFG = DB_BUNDLE['third_party_databases']['checkv']
CHECKV_VERSION = CHECKV_CFG['version']
CHECKV_TARGET_DIR = os.path.join(DB_DIR, CHECKV_CFG['target_dir'])
CHECKV_SOURCE_DIR = os.path.join(DB_DIR, 'downloads', CHECKV_CFG['download_subdir'])


rule download_and_setup_checkv:
    """Download the CheckV database and move it to a versioned directory."""
    output:
        touch(os.path.join(DB_DIR, 'markers', _setup_marker_name('checkv')))
    params:
        source_dir=CHECKV_SOURCE_DIR,
        target_dir=CHECKV_TARGET_DIR,
    conda:
        _ce('metalink-annotation')
    log:
        os.path.join(DB_DIR, 'logs', 'setup_checkv.log')
    shell:
        """
        (
        rm -rf {params.source_dir}
        rm -rf {params.target_dir}
        checkv download_database {DB_DIR}/downloads
        if [ -d {params.source_dir} ]; then
            mkdir -p $(dirname {params.target_dir})
            mv {params.source_dir} {params.target_dir}
        fi
        ) > {log} 2>&1
        """


VIBRANT_CFG = DB_BUNDLE['third_party_databases']['vibrant']
VIBRANT_VERSION = VIBRANT_CFG['version']
VIBRANT_TARGET_DIR = os.path.join(DB_DIR, VIBRANT_CFG['target_dir'])


rule download_and_setup_vibrant:
    """Download the VIBRANT database."""
    output:
        touch(os.path.join(DB_DIR, 'markers', _setup_marker_name('vibrant')))
    params:
        target_dir=VIBRANT_TARGET_DIR,
    conda:
        _ce('metalink-vibrant')
    log:
        os.path.join(DB_DIR, 'logs', 'setup_vibrant.log')
    shell:
        """
        (
        mkdir -p {params.target_dir}
        download-db.sh {params.target_dir}
        ) > {log} 2>&1
        """


IPHOP_CFG = DB_BUNDLE['third_party_databases']['iphop']
IPHOP_VERSION = IPHOP_CFG['version']
IPHOP_TARGET_DIR = os.path.join(DB_DIR, IPHOP_CFG['target_dir'])
IPHOP_ARCHIVE = os.path.basename(IPHOP_CFG['archive_url'])
IPHOP_EXTRACT_PATTERN = IPHOP_CFG['extracted_dir_pattern']


rule download_iphop:
    """Download the iPHoP host prediction database archive."""
    output:
        os.path.join(DB_DIR, 'downloads', IPHOP_ARCHIVE)
    params:
        archive_url=IPHOP_CFG['archive_url'],
    log:
        os.path.join(DB_DIR, 'logs', 'download_iphop.log')
    shell:
        """
        (
        mkdir -p $(dirname {output})
        wget -c {params.archive_url} \
        --no-check-certificate -O {output}
        ) > {log} 2>&1
        """


rule setup_iphop:
    """Extract and move the iPHoP database to a versioned directory."""
    input:
        rules.download_iphop.output
    output:
        touch(os.path.join(DB_DIR, 'markers', _setup_marker_name('iphop')))
    params:
        target_dir=IPHOP_TARGET_DIR,
        extract_pattern=IPHOP_EXTRACT_PATTERN,
        archive_name=IPHOP_ARCHIVE,
    conda:
        _ce('metalink-iphop')
    log:
        os.path.join(DB_DIR, 'logs', 'setup_iphop.log')
    shell:
        """
        (
        cd {DB_DIR}
        tar -zxvf downloads/{params.archive_name}
        rm -rf {params.target_dir}
        extracted=$(find . -maxdepth 1 -type d -name '{params.extract_pattern}' | head -n1)
        if [ -z "$extracted" ]; then
            echo "ERROR: could not find extracted iPHoP database directory" >&2
            exit 1
        fi
        mkdir -p $(dirname {params.target_dir})
        mv "$extracted" {params.target_dir}
        ) > {log} 2>&1
        """


GTDBTK_CFG = DB_BUNDLE['third_party_databases']['gtdbtk']
GTDBTK_VERSION = GTDBTK_CFG['version']
GTDBTK_TARGET_DIR = os.path.join(DB_DIR, GTDBTK_CFG['target_dir'])
GTDBTK_ARCHIVE_URL = GTDBTK_CFG['archive_url']
GTDBTK_ARCHIVE_NAME = os.path.basename(GTDBTK_ARCHIVE_URL)
GTDBTK_EXTRACTED_DIR = GTDBTK_CFG['extracted_dir']


rule download_gtdbtk:
    """Download the GTDB-Tk reference data archive."""
    output:
        os.path.join(DB_DIR, 'downloads', GTDBTK_ARCHIVE_NAME)
    params:
        archive_url=GTDBTK_CFG['archive_url'],
    log:
        os.path.join(DB_DIR, 'logs', 'download_gtdbtk.log')
    shell:
        """
        (
        mkdir -p $(dirname {output})
        wget -c {params.archive_url} -O {output}
        ) > {log} 2>&1
        """


rule setup_gtdbtk:
    """Extract and move the GTDB-Tk reference data to a versioned directory."""
    input:
        rules.download_gtdbtk.output
    output:
        touch(os.path.join(DB_DIR, 'markers', _setup_marker_name('gtdbtk')))
    params:
        target_dir=GTDBTK_TARGET_DIR,
        extracted_dir=GTDBTK_EXTRACTED_DIR,
        archive_name=GTDBTK_ARCHIVE_NAME,
    conda:
        _ce('metalink-gtdbtk')
    log:
        os.path.join(DB_DIR, 'logs', 'setup_gtdbtk.log')
    shell:
        """
        (
        cd {DB_DIR}
        tar -xvzf downloads/{params.archive_name}
        rm -rf {params.target_dir}
        mkdir -p $(dirname {params.target_dir})
        mv {params.extracted_dir} {params.target_dir}
        ) > {log} 2>&1
        """


CHECKM_CFG = DB_BUNDLE['third_party_databases']['checkm']
CHECKM_VERSION = CHECKM_CFG['version']
CHECKM_TARGET_DIR = os.path.join(DB_DIR, CHECKM_CFG['target_dir'])
CHECKM_ARCHIVE_URL = CHECKM_CFG['archive_url']
CHECKM_ARCHIVE_NAME = os.path.basename(CHECKM_ARCHIVE_URL)


rule download_checkm:
    """Download the CheckM reference data archive."""
    output:
        os.path.join(DB_DIR, 'downloads', CHECKM_ARCHIVE_NAME)
    params:
        archive_url=CHECKM_CFG['archive_url'],
    log:
        os.path.join(DB_DIR, 'logs', 'download_checkm.log')
    shell:
        """
        (
        mkdir -p $(dirname {output})
        wget -c {params.archive_url} -O {output}
        ) > {log} 2>&1
        """


rule setup_checkm:
    """Extract and configure the CheckM reference data."""
    input:
        rules.download_checkm.output
    output:
        touch(os.path.join(DB_DIR, 'markers', _setup_marker_name('checkm')))
    params:
        target_dir=CHECKM_TARGET_DIR,
        archive_name=CHECKM_ARCHIVE_NAME,
    conda:
        _ce('metalink-binning')
    log:
        os.path.join(DB_DIR, 'logs', 'setup_checkm.log')
    shell:
        """
        (
        cd {DB_DIR}
        rm -rf {params.target_dir}
        mkdir -p {params.target_dir}
        tar -zxvf downloads/{params.archive_name} -C {params.target_dir} --strip-components=1
        checkm data setRoot {params.target_dir}
        ) > {log} 2>&1
        """


VIRSORTER_CFG = DB_BUNDLE['third_party_databases']['virsorter']
VIRSORTER_VERSION = VIRSORTER_CFG['version']
VIRSORTER_TARGET_DIR = os.path.join(DB_DIR, VIRSORTER_CFG['target_dir'])


rule download_and_setup_virsorter:
    """Download and set up the VirSorter2 database."""
    output:
        touch(os.path.join(DB_DIR, 'markers', _setup_marker_name('virsorter')))
    params:
        target_dir=VIRSORTER_TARGET_DIR,
    threads: 4
    conda:
        _ce('metalink-virsorter')
    log:
        os.path.join(DB_DIR, 'logs', 'setup_virsorter.log')
    shell:
        """
        (
        virsorter setup -d {params.target_dir} -j {threads}
        ) > {log} 2>&1
        """


CHECKM2_CFG = DB_BUNDLE['third_party_databases']['checkm2']
CHECKM2_VERSION = CHECKM2_CFG['version']
CHECKM2_TARGET_DIR = os.path.join(DB_DIR, CHECKM2_CFG['target_dir'])
# Stage CheckM2 download in a tool-specific subdirectory so its md5sums.txt
# does not collide with other downloads in downloads/.
CHECKM2_STAGING_DIR = os.path.join(DB_DIR, 'downloads', 'checkm2')
CHECKM2_SOURCE_DIR = os.path.join(CHECKM2_STAGING_DIR, CHECKM2_CFG['download_subdir'])


rule download_and_setup_checkm2:
    """Download and arrange the CheckM2 database."""
    output:
        touch(os.path.join(DB_DIR, 'markers', _setup_marker_name('checkm2')))
    params:
        source_dir=CHECKM2_SOURCE_DIR,
        target_dir=CHECKM2_TARGET_DIR,
        staging_dir=CHECKM2_STAGING_DIR,
    conda:
        _ce('metalink-checkm2')
    log:
        os.path.join(DB_DIR, 'logs', 'setup_checkm2.log')
    shell:
        """
        (
        cd {DB_DIR}
        rm -rf {params.staging_dir}
        rm -rf {params.target_dir}
        mkdir -p {params.staging_dir}
        checkm2 database --download --path {params.staging_dir} --no_write_json_db
        mkdir -p $(dirname {params.target_dir})
        mv {params.source_dir} {params.target_dir}
        ) > {log} 2>&1
        """


# ---------------------------------------------------------------------------
# Error / success handlers
# ---------------------------------------------------------------------------

onerror:
    # Only clean up marker files that Snakemake has explicitly flagged as
    # incomplete.  Do not wipe the whole downloads directory or markers for
    # jobs that already succeeded -- that would force perfectly good database
    # setups to rerun after an unrelated failure.
    import base64

    incomplete_dir = os.path.join(DB_DIR, '.snakemake', 'incomplete')
    if os.path.isdir(incomplete_dir):
        for entry in os.listdir(incomplete_dir):
            try:
                path = base64.b64decode(entry.encode()).decode()
            except Exception:
                continue
            # Only remove files that live inside this database directory and
            # are marker outputs from the setup workflow.
            if (
                path.startswith(DB_DIR + os.sep)
                and path.endswith('.setup.mk')
                and os.path.exists(path)
            ):
                try:
                    os.remove(path)
                except OSError:
                    pass


onsuccess:
    _record_installed_bundle()
    print("All selected Metalink database setup finished")
