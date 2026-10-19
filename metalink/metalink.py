#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import os
import logging
import multiprocessing
import subprocess
import shlex
import tempfile
import time
import click
from pathlib import Path
from ruamel.yaml import YAML, CommentedMap
from ruamel.yaml.scalarstring import SingleQuotedScalarString
yaml = YAML(typ='rt')
yaml.default_flow_style = False
yaml.indent(mapping=2, sequence=4, offset=2)
yaml.preserve_quotes = True


def _load_db_bundle():
    """Load the Metalink database bundle manifest bundled with the code."""
    bundle_path = Path(__file__).with_name('config') / 'db_bundle.yaml'
    if not bundle_path.exists():
        return {}
    with open(bundle_path) as f:
        return yaml.load(f) or {}


_DB_BUNDLE = _load_db_bundle()


def set_logger(silent=False):
    level = logging.WARNING if silent else logging.INFO
    logging.basicConfig(
        level=level,
        datefmt="%Y-%m-%d %H:%M",
        format="[%(asctime)s %(levelname)s] %(message)s",
    )

def set_nested(d, path, val):
    """path: list of keys, e.g. ['Params','xx','yy']"""
    for k in path[:-1]:
        d = d.setdefault(k, {})
    d[path[-1]] = val

def _get_nested(d, path, default=None):
    """Safely retrieve a nested dict value by a list of keys."""
    for k in path:
        if not isinstance(d, dict) or k not in d:
            return default
        d = d[k]
    return d

def _warn_zero_kraken2_classification(work_dir, since, bracken_level=None,
                                      bracken_threshold=0):
    """Warn about problematic kraken2 reports from this run.

    Two cases are flagged (only reports written after `since`):
      1. 0% classified sequences at all (taxonomy and function reports).
      2. only when bracken_level is given (r mode runs bracken): taxonomy
         reports where reads were classified but none at the bracken level,
         so bracken abundance estimation was skipped.
    """
    import glob
    try:
        pattern = os.path.join(work_dir, 'results', '**',
                               '*.kraken2.report.tsv')
        zero = []
        no_species = []
        for rep in glob.glob(pattern, recursive=True):
            if os.path.getmtime(rep) < since:
                continue
            # The bracken-related check only applies to taxonomy reports of
            # runs that include a bracken step (r mode); LCFH function
            # reports reuse the D/P/.../S rank codes for functional levels.
            is_taxonomy = 'taxonomy' in rep and bracken_level is not None
            root_reads = None
            u_reads = None
            has_level_reads = False
            with open(rep) as f:
                for line in f:
                    fields = line.rstrip('\n').split('\t')
                    # kraken2 report: pct, clade reads, direct reads, rank code, ...
                    if len(fields) < 6:
                        continue
                    rank = fields[3].strip()
                    if rank == 'U':
                        try:
                            u_reads = float(fields[2])
                        except ValueError:
                            pass
                    elif rank == 'R':
                        # clade-level read count under root = total classified;
                        # more reliable than the rounded percentage column
                        try:
                            root_reads = float(fields[1])
                        except ValueError:
                            pass
                    elif is_taxonomy and rank == bracken_level:
                        # bracken keeps a node when its clade read count
                        # (field 1) is >= threshold; reads classified deeper
                        # are included there and mapped up during reestimation
                        try:
                            if float(fields[1]) >= max(int(bracken_threshold), 1):
                                has_level_reads = True
                        except (ValueError, TypeError):
                            pass
            # A report with zero classifications contains only the U line
            # (kraken2 omits the root row), so fall back to the U count.
            classified = root_reads if root_reads is not None else 0.0
            if root_reads is None and u_reads is None:
                continue
            if classified == 0:
                zero.append(rep)
            elif is_taxonomy and not has_level_reads:
                no_species.append(rep)
        def _label(rep):
            """sample name plus the classification stage the report belongs to"""
            sample = os.path.basename(rep).split('.kraken2.report.tsv')[0]
            if 'taxonomy' in rep:
                return '{} (taxonomy)'.format(sample)
            if 'function' in rep:
                return '{} (function)'.format(sample)
            return sample
        if zero:
            sample_names = ', '.join(_label(r) for r in zero)
            logging.warning(
                '%d kraken2 report(s) had 0%% classified sequences: %s\n'
                'This often happens when a reduced kraken2 database (e.g. '
                'kraken2_16GB) is used together with a high confidence threshold.\n'
                'Try lowering the kraken2 confidence in your config and rerun.\n'
                'See the FAQ "How do I choose the kraken2 confidence thresholds for '
                'taxonomic and functional classification?" in README_DETAIL.md '
                'for guidance.',
                len(zero), sample_names)
        if no_species:
            sample_names = ', '.join(_label(r) for r in no_species)
            logging.warning(
                '%d kraken2 taxonomy report(s) had classified reads but NONE '
                'classified at the bracken level (%s): %s\n'
                'Bracken abundance estimation (bracken -l %s) was skipped for them '
                '(empty Bracken outputs were produced).\n'
                'This often happens with a reduced kraken2 database and/or a high '
                'confidence threshold; consider lowering '
                '--reads-taxonomy-kraken2-confidence and rerun.\n'
                'See the FAQ "How do I choose the kraken2 confidence thresholds for '
                'taxonomic and functional classification?" in README_DETAIL.md '
                'for guidance.',
                len(no_species), bracken_level, sample_names, bracken_level)
    except Exception:
        pass

def _validate_cami_gtdbtk_compatibility(config, mode):
    """CAMI taxonomic profile output is disabled when mMAG taxonomy uses GTDB-Tk.

    When both are enabled in a mode that uses mMAGs, disable CAMI output and
    warn the user so the workflow does not produce malformed lineage columns.
    """
    modes_using_mmags = {'m', 'mc', 'mcr', 'pmr', 'vm', 'vmc', 'vmcr'}
    if mode not in modes_using_mmags:
        return False

    if mode == 'pmr':
        mags = config.get('mags', 'm')
        if mags == 'v':
            # vMAG taxonomy does not go through GTDB-Tk; no conflict.
            return False
        cami_key = ['Params_pmr_{}_sample'.format(mags), '05_merge_tax_fun_depth', 'CAMI_format_output']
    else:
        cami_key = ['Params_{}_sample'.format(mode), '05_merge_tax_fun_depth', 'CAMI_format_output']
    gtdbtk_key = ['Params_m_sample', '02_taxonomy', 'MAGs_taxonomy_classify_method']

    cami_enabled = _get_nested(config, cami_key, False)
    mag_tax_method = _get_nested(config, gtdbtk_key, 'gtdbtk')

    if cami_enabled and str(mag_tax_method).lower() == 'gtdbtk':
        logging.warning(
            'CAMI-format taxonomic profile output is disabled because mMAG taxonomy classify method'
            'is set to GTDB-Tk (%s). Use --mmags-taxonomy-classify-method BAT to enable it.',
            ' -> '.join(gtdbtk_key)
        )
        set_nested(config, cami_key, False)
        return True
    return False


def _migrate_pmr_quality_gate(config, mags):
    """Move quality-gate params from
    Params_pmr_{mags}_sample -> 05_merge_tax_fun_depth (<=0.3.5 configs) to
    Params_pmr_{mags}_project -> 02_quality_gate (quality gate now runs before
    coverm/reference building). Returns True when anything was migrated."""
    sample_key = 'Params_pmr_{}_sample'.format(mags)
    proj_key = 'Params_pmr_{}_project'.format(mags)
    merge_cfg = config.get(sample_key, {}).get('05_merge_tax_fun_depth')
    if not isinstance(merge_cfg, dict):
        return False
    keys = ['MAGs_min_completeness1', 'MAGs_max_contamination1',
            'MAGs_min_completeness2', 'MAGs_max_contamination2',
            'checkv_quality']
    present = [k for k in keys if k in merge_cfg]
    if not present:
        return False
    gate_cfg = config.setdefault(proj_key, {}).setdefault('02_quality_gate', {})
    for k in present:
        gate_cfg.setdefault(k, merge_cfg.pop(k))
    logging.warning(
        "Config migration: moved %s from %s -> 05_merge_tax_fun_depth to "
        "%s -> 02_quality_gate (quality gate moved before coverm/reference).",
        present, sample_key, proj_key)
    return True


def _validate_pmr_iphop_custom_mags_dir(config, mags):
    """pmr --mags v runs no m-p dereplication chain, so iPHoP custom_MAGs_dir
    cannot point at the auto-generated dereplicated mMAGs directories
    (auto_drep99/auto_drep95); only 'none' or a custom path is allowed.
    pmr --mags vm includes m-p, so all values are allowed there."""
    if mags != 'v':
        return
    iphop_key = ['Params_v_project', '07_iPHoP', 'custom_MAGs_dir']
    custom = _get_nested(config, iphop_key, 'none')
    if custom in ('auto_drep99', 'auto_drep95'):
        raise ValueError(
            "pmr --mags v does not run mMAG dereplication (m-p chain), so "
            "{} cannot be '{}' (no auto dRep mMAGs are produced). "
            "Use 'none' or an absolute path to a custom MAGs directory.".format(
                ' -> '.join(iphop_key), custom)
        )


def _validate_bracken_level(config, mode):
    """r mode only: bracken_level must be one of D/P/C/O/F/G/S; an invalid
    value would silently skip bracken for every sample."""
    if mode != 'r':
        return
    lvl = _get_nested(config, ['Params_r_sample', '03_taxonomy', 'bracken_level'])
    if lvl is None:
        return  # template default applies
    lvl = str(lvl).strip()
    if lvl not in set('DPCOFGS') or len(lvl) != 1:
        raise ValueError(
            "Invalid bracken_level {!r} in Params_r_sample -> 03_taxonomy; "
            "choose one of D/P/C/O/F/G/S.".format(lvl))
    config['Params_r_sample']['03_taxonomy']['bracken_level'] = lvl


def _validate_read_type_assembly_method(config, mode, skip_level=0):
    """Validate that read_type and assembly_method are compatible.

    read_type:
      - paired: must use megahit or spades
      - paired-rna: must use spades (runs rnaSPAdes)
      - single: must use megahit
      - single-rna: must use spades (runs rnaSPAdes)
      - pacbio-raw/corr/hifi, nano-raw/corr/hq: must use metaflye
    """
    if mode == 'r' or skip_level in (2, 3):
        return

    read_type = config.get('read_type')
    assembly_method = _get_nested(
        config, ['Params_b_sample', '02_assembly', 'method'], 'megahit'
    )
    long_read_types = {
        'pacbio-raw', 'pacbio-corr', 'pacbio-hifi',
        'nano-raw', 'nano-corr', 'nano-hq'
    }
    rna_read_types = {'paired-rna', 'single-rna'}
    if read_type in long_read_types and assembly_method != 'metaflye':
        raise ValueError(
            f"read_type '{read_type}' is a third-generation long read and requires "
            f"assembly method 'metaflye', got '{assembly_method}'"
        )
    if read_type in rna_read_types and assembly_method != 'spades':
        raise ValueError(
            f"read_type '{read_type}' requires assembly method 'spades' (rnaSPAdes), "
            f"got '{assembly_method}'"
        )
    if read_type == 'paired' and assembly_method not in {'megahit', 'spades'}:
        raise ValueError(
            f"read_type 'paired' supports assembly methods 'megahit' or 'spades', "
            f"got '{assembly_method}'"
        )
    if read_type == 'single' and assembly_method != 'megahit':
        raise ValueError(
            f"read_type 'single' supports assembly method 'megahit', got '{assembly_method}'"
        )


_SECTION_ORDER = {}

def get_section_name(prefix):
    rest = prefix.replace('Params_', '')
    parts = rest.split('_')
    mode_map = {
        'b': 'base',
        'c': 'contig',
        'cr': 'contig-read',
        'r': 'read',
        'm': 'mMAG',
        'mc': 'mMAG-contig',
        'mcr': 'mMAG-contig-read',
        'pmr': 'project-mMAG-read',
        'pmr_m': 'project-mMAG-read',
        'pmr_v': 'project-vMAG-read',
        'pmr_vm': 'project-vMAG-mMAG-read',
        'v': 'vMAG',
        'vc': 'vMAG-contig',
        'vcr': 'vMAG-contig-read',
        'vm': 'vMAG-mMAG',
        'vmc': 'vMAG-mMAG-contig',
        'vmcr': 'vMAG-mMAG-contig-read',
    }
    level_map = {
        'sample': 'sample-level',
        'project': 'project-level',
    }
    # Placeholder form shared by all pmr --mags chains, e.g.
    # Params_pmr_{mags}_sample: resolved to Params_pmr_{m,v,vm}_* at config time.
    if '{mags}' in parts:
        level_name = level_map.get(parts[-1], parts[-1])
        return f"pmr (mags-dependent) {level_name} options"
    # Longest-prefix match: two-segment mode keys (pmr_vm/pmr_m/pmr_v/...)
    # must be tried before single-segment ones, otherwise Params_pmr_m_sample
    # would be parsed as mode 'pmr' + level 'm'.
    if len(parts) >= 3 and '_'.join(parts[:2]) in mode_map:
        mode = '_'.join(parts[:2])
        level = parts[2]
        mode_name = mode_map.get(mode, mode.upper())
        level_name = level_map.get(level, level)
        return f"{mode_name} {level_name} options"
    if len(parts) >= 2:
        mode = parts[0]
        level = parts[1]
        mode_name = mode_map.get(mode, mode.upper())
        level_name = level_map.get(level, level)
        return f"{mode_name} {level_name} options"
    if rest == 'mags':
        return 'General options'
    return f"{rest} options"


def _resolve_mags_placeholder(yaml_path, mags):
    """Substitute the {mags} placeholder in a PARAMS_HELP yaml_path with the
    active --mags value (m/v/vm) for pmr mode."""
    return [p.replace('{mags}', mags) if isinstance(p, str) else p
            for p in yaml_path]


def _pmr_option_applicable(opt_name, yaml_path, cfg, mags):
    """Whether a pmr CLI option applies to the active --mags chain.

    Two checks: (1) type-specific threshold options only apply to the chains
    that use them; (2) the target section/subsection must exist in the
    per-mags config — options aimed at chains this mags does not run (e.g.
    virus_* under --mags m, binning/mmags_* under --mags v) would otherwise
    be written as dead keys nobody reads. Inapplicable options are skipped
    with a warning (returns False).
    """
    option_name = f"--{opt_name.replace('_', '-')}"
    if opt_name in _PMR_M_TYPE_OPTIONS and mags not in ('m', 'vm'):
        logging.warning(
            f"{option_name} only applies to --mags m/vm (mMAG quality thresholds); "
            f"ignored for --mags {mags}.")
        return False
    if opt_name in _PMR_V_TYPE_OPTIONS and mags not in ('v', 'vm'):
        logging.warning(
            f"{option_name} only applies to --mags v/vm (CheckV quality threshold); "
            f"ignored for --mags {mags}.")
        return False
    section = yaml_path[0]
    subsection = yaml_path[1] if len(yaml_path) > 1 else None
    if section not in cfg or (subsection is not None and subsection not in cfg.get(section, {})):
        logging.warning(
            f"{option_name} does not apply to --mags {mags} "
            f"(config has no {' -> '.join(yaml_path[:2])} section); ignored.")
        return False
    return True


class OrderedCommand(click.Command):
    """Command that groups options by their _section_name attribute."""

    def format_options(self, ctx, formatter):
        from collections import OrderedDict
        from click.formatting import measure_table, iter_rows, term_len, wrap_text

        sections = OrderedDict()
        for param in self.get_params(ctx):
            rv = param.get_help_record(ctx)
            if rv is not None:
                section_name = getattr(param, '_section_name', None)
                if section_name is None:
                    section_name = 'General options'
                if section_name not in sections:
                    sections[section_name] = []
                sections[section_name].append(rv)

        mode = ctx.invoked_subcommand or ctx.command.name
        ordered_sections = ['General options', 'Intermediate files & resources']
        if mode in _SECTION_ORDER:
            for sec in _SECTION_ORDER[mode]:
                if sec not in ordered_sections:
                    ordered_sections.append(sec)
        for sec in sections.keys():
            if sec not in ordered_sections:
                ordered_sections.append(sec)

        # Calculate global max option name length across ALL sections
        # for consistent alignment across the entire help page
        global_max_opt_len = 0
        for section_name in ordered_sections:
            if section_name in sections:
                for row in sections[section_name]:
                    if row[0]:
                        global_max_opt_len = max(global_max_opt_len, term_len(row[0]))
        # first_col = max option width + spacing (minimum 30 + 2)
        first_col = max(30, global_max_opt_len) + 2

        for section_name in ordered_sections:
            if section_name in sections:
                with formatter.section(section_name):
                    rows = sections[section_name]
                    for first, second in iter_rows(rows, 2):
                        formatter.write(f"{'':>{formatter.current_indent}}{first}")
                        if not second:
                            formatter.write(chr(10))
                            continue
                        spacing = first_col - formatter.current_indent - term_len(first)
                        if spacing > 0:
                            formatter.write(" " * spacing)
                        else:
                            formatter.write(chr(10))
                            formatter.write(" " * first_col)

                        text_width = max(formatter.width - first_col - 2, 10)
                        if second and second[0] == '\b':
                            # Click convention: a help string starting with \b
                            # should not be re-formatted. Strip the marker and
                            # wrap each original line independently so that
                            # hand-written lists/indentation are preserved while
                            # still ensuring no individual line exceeds the
                            # terminal width.
                            content = second[1:]
                            if content.startswith('\n'):
                                content = content[1:]
                            wrapped_text = '\n'.join(
                                wrap_text(line, text_width, preserve_paragraphs=True)
                                for line in content.split('\n')
                            )
                        else:
                            wrapped_text = wrap_text(second, text_width, preserve_paragraphs=True)
                        lines = wrapped_text.splitlines()

                        if lines:
                            formatter.write(lines[0] + chr(10))
                            for line in lines[1:]:
                                formatter.write(" " * first_col + line + chr(10))
                        else:
                            formatter.write(chr(10))

BASE_KEYS = [
    'mode', 'mags', 'threads', 'jobs',
    'workflow_target', 'workflow_snakefile',
    'workflow_includes', 'workflow_common_includes',
    'input_tsv', 'Directory', 'conda_envs',
    'run_project', 'name_project', 'read_type',
    'taxonomy_rank', 'skip_level', 'kraken2_taxonomy_jobs', 'keep_intermediate'
]

MODE_SPECIFIC = {
    'b':    [],
    'c':    ['Params_c_sample', 'Params_c_project'],
    'cr':   ['Params_c_sample', 'Params_c_project', 'Params_cr_sample', 'Params_cr_project'],
    'r':    ['Params_r_sample', 'Params_r_project'],
    'm':    ['Params_m_sample', 'Params_m_project'],
    'mc':   ['Params_m_sample', 'Params_m_project', 'Params_c_sample', 'Params_c_project', 'Params_mc_sample', 'Params_mc_project'],
    'mcr':  ['Params_m_sample', 'Params_m_project', 'Params_c_sample', 'Params_c_project', 'Params_mcr_sample', 'Params_mcr_project'],
    'pmr_m':  ['Params_m_sample', 'Params_m_project', 'Params_pmr_m_sample', 'Params_pmr_m_project'],
    'pmr_v':  ['Params_v_sample', 'Params_v_project', 'Params_pmr_v_sample', 'Params_pmr_v_project'],
    'pmr_vm': ['Params_m_sample', 'Params_m_project', 'Params_v_sample', 'Params_v_project', 'Params_pmr_vm_sample', 'Params_pmr_vm_project'],
    'v':    ['Params_v_sample', 'Params_v_project'],
    'vc':   ['Params_v_sample', 'Params_v_project', 'Params_c_sample', 'Params_c_project', 'Params_vc_sample', 'Params_vc_project'],
    'vcr':  ['Params_v_sample', 'Params_v_project', 'Params_c_sample', 'Params_c_project', 'Params_vcr_sample', 'Params_vcr_project'],
    'vm':   ['Params_v_sample', 'Params_v_project', 'Params_m_sample', 'Params_m_project', 'Params_vm_sample', 'Params_vm_project'],
    'vmc':  ['Params_v_sample', 'Params_v_project', 'Params_m_sample', 'Params_m_project', 'Params_c_sample', 'Params_c_project', 'Params_vmc_sample', 'Params_vmc_project'],
    'vmcr': ['Params_v_sample', 'Params_v_project', 'Params_m_sample', 'Params_m_project', 'Params_c_sample', 'Params_c_project', 'Params_vmcr_sample', 'Params_vmcr_project'],
}

def _merge_configs(general_cfg, mode_cfg):
    """
    Merge general config with mode-specific config.
    Mode-specific values override general values for the same keys.
    """
    merged = CommentedMap()
    # First copy all from general
    for key in general_cfg:
        merged[key] = general_cfg[key]
    # Then override/add from mode-specific
    for key in mode_cfg:
        merged[key] = mode_cfg[key]
    return merged


def _compact_yaml(filepath):
    """Add logical blank lines between config sections while preserving
    blank lines inside YAML multiline string scalars (| / >)."""
    with open(filepath, 'r') as f:
        content = f.read()
    lines = content.split('\n')
    result = []
    base_keys_set = set(BASE_KEYS)
    prev_was_base = False
    prev_was_params = False
    for line in lines:
        stripped = line.strip()
        is_top_level = (
            stripped != ''
            and not line.startswith(' ')
            and not line.startswith('#')
            and ':' in line
        )
        if is_top_level:
            key = stripped.split(':')[0]
            is_base = key in base_keys_set
            is_params = key.startswith('Params_')
            # Insert a section separator only if one isn't already present.
            needs_blank = (prev_was_base and is_params) or (prev_was_params and is_params)
            if needs_blank and result and result[-1].strip() != '':
                result.append('')
            prev_was_base = is_base
            prev_was_params = is_params
        # Preserve every original line, including blank lines (they may be
        # part of multiline string scalars).
        result.append(line)
    with open(filepath, 'w') as f:
        f.write('\n'.join(result) + '\n')


# Workflow include lists per skip_level:
#   0: input raw reads (trim + assembly + binning)
#   1: input clean reads (pass trim)
#   2: input clean reads and contigs (pass trim and assembly)
#   3: input clean reads, contigs and mMAGs (pass trim, assembly, and binning)
WORKFLOW_MAP = {
    0: {
        'b':    ['b-trim.smk', 'b-assembly.smk', 'b-binning.smk'],
        'c':    ['b-trim.smk', 'b-assembly.smk', 'c.base.smk', 'c.mtf.smk'],
        'cr':   ['b-trim.smk', 'b-assembly.smk', 'c.base.smk', 'cr.smk'],
        'r':    ['b-trim.smk', 'r.smk'],
        'm':    ['b-trim.smk', 'b-assembly.smk', 'b-binning.smk', 'm-s.smk', 'm-p.smk', 'm.mtf.smk'],
        'mc':   ['b-trim.smk', 'b-assembly.smk', 'b-binning.smk', 'm-s.smk', 'c.base.smk', 'mc.smk'],
        'mcr':  ['b-trim.smk', 'b-assembly.smk', 'b-binning.smk', 'm-s.smk', 'c.base.smk', 'mcr.smk'],
        'pmr':  ['b-trim.smk', 'b-assembly.smk', 'b-binning.smk', 'm-s.smk', 'm-p.smk', 'pmr_m.smk'],
        'v':    ['b-trim.smk', 'b-assembly.smk', 'v-s.smk', 'v-p.smk'],
        'vc':   ['b-trim.smk', 'b-assembly.smk', 'v-s.smk', 'c.base.smk', 'vc.smk'],
        'vcr':  ['b-trim.smk', 'b-assembly.smk', 'v-s.smk', 'c.base.smk', 'vcr.smk'],
        'vm':   ['b-trim.smk', 'b-assembly.smk', 'b-binning.smk', 'v-s.smk', 'v-p.smk', 'm-s.smk', 'm-p.smk', 'vm.smk'],
        'vmc':  ['b-trim.smk', 'b-assembly.smk', 'b-binning.smk', 'v-s.smk', 'm-s.smk', 'c.base.smk', 'vmc.smk'],
        'vmcr': ['b-trim.smk', 'b-assembly.smk', 'b-binning.smk', 'v-s.smk', 'm-s.smk', 'c.base.smk', 'vmcr.smk'],
    },
    1: {
        'b':    ['b-pass_trim.smk', 'b-assembly.smk', 'b-binning.smk'],
        'c':    ['b-pass_trim.smk', 'b-assembly.smk', 'c.base.smk', 'c.mtf.smk'],
        'cr':   ['b-pass_trim.smk', 'b-assembly.smk', 'c.base.smk', 'cr.smk'],
        'r':    ['b-pass_trim.smk', 'r.smk'],
        'm':    ['b-pass_trim.smk', 'b-assembly.smk', 'b-binning.smk', 'm-s.smk', 'm-p.smk', 'm.mtf.smk'],
        'mc':   ['b-pass_trim.smk', 'b-assembly.smk', 'b-binning.smk', 'm-s.smk', 'c.base.smk', 'mc.smk'],
        'mcr':  ['b-pass_trim.smk', 'b-assembly.smk', 'b-binning.smk', 'm-s.smk', 'c.base.smk', 'mcr.smk'],
        'pmr':  ['b-pass_trim.smk', 'b-assembly.smk', 'b-binning.smk', 'm-s.smk', 'm-p.smk', 'pmr_m.smk'],
        'v':    ['b-pass_trim.smk', 'b-assembly.smk', 'v-s.smk', 'v-p.smk'],
        'vc':   ['b-pass_trim.smk', 'b-assembly.smk', 'v-s.smk', 'c.base.smk', 'vc.smk'],
        'vcr':  ['b-pass_trim.smk', 'b-assembly.smk', 'v-s.smk', 'c.base.smk', 'vcr.smk'],
        'vm':   ['b-pass_trim.smk', 'b-assembly.smk', 'b-binning.smk', 'v-s.smk', 'v-p.smk', 'm-s.smk', 'm-p.smk', 'vm.smk'],
        'vmc':  ['b-pass_trim.smk', 'b-assembly.smk', 'b-binning.smk', 'v-s.smk', 'm-s.smk', 'c.base.smk', 'vmc.smk'],
        'vmcr': ['b-pass_trim.smk', 'b-assembly.smk', 'b-binning.smk', 'v-s.smk', 'm-s.smk', 'c.base.smk', 'vmcr.smk'],
    },
    2: {
        'b':    ['b-pass_trim.smk', 'b-pass_assembly.smk', 'b-binning.smk'],
        'c':    ['b-pass_trim.smk', 'b-pass_assembly.smk', 'c.base.smk', 'c.mtf.smk'],
        'cr':   ['b-pass_trim.smk', 'b-pass_assembly.smk', 'c.base.smk', 'cr.smk'],
        'r':    ['b-pass_trim.smk', 'r.smk'],
        'm':    ['b-pass_trim.smk', 'b-pass_assembly.smk', 'b-binning.smk', 'm-s.smk', 'm-p.smk', 'm.mtf.smk'],
        'mc':   ['b-pass_trim.smk', 'b-pass_assembly.smk', 'b-binning.smk', 'm-s.smk', 'c.base.smk', 'mc.smk'],
        'mcr':  ['b-pass_trim.smk', 'b-pass_assembly.smk', 'b-binning.smk', 'm-s.smk', 'c.base.smk', 'mcr.smk'],
        'pmr':  ['b-pass_trim.smk', 'b-pass_assembly.smk', 'b-binning.smk', 'm-s.smk', 'm-p.smk', 'pmr_m.smk'],
        'v':    ['b-pass_trim.smk', 'b-pass_assembly.smk', 'v-s.smk', 'v-p.smk'],
        'vc':   ['b-pass_trim.smk', 'b-pass_assembly.smk', 'v-s.smk', 'c.base.smk', 'vc.smk'],
        'vcr':  ['b-pass_trim.smk', 'b-pass_assembly.smk', 'v-s.smk', 'c.base.smk', 'vcr.smk'],
        'vm':   ['b-pass_trim.smk', 'b-pass_assembly.smk', 'b-binning.smk', 'v-s.smk', 'v-p.smk', 'm-s.smk', 'm-p.smk', 'vm.smk'],
        'vmc':  ['b-pass_trim.smk', 'b-pass_assembly.smk', 'b-binning.smk', 'v-s.smk', 'm-s.smk', 'c.base.smk', 'vmc.smk'],
        'vmcr': ['b-pass_trim.smk', 'b-pass_assembly.smk', 'b-binning.smk', 'v-s.smk', 'm-s.smk', 'c.base.smk', 'vmcr.smk'],
    },
    3: {
        'b':    ['b-pass_trim.smk', 'b-pass_assembly.smk', 'b-pass_binning.smk'],
        'c':    ['b-pass_trim.smk', 'b-pass_assembly.smk', 'c.base.smk', 'c.mtf.smk'],
        'cr':   ['b-pass_trim.smk', 'b-pass_assembly.smk', 'c.base.smk', 'cr.smk'],
        'r':    ['b-pass_trim.smk', 'r.smk'],
        'm':    ['b-pass_trim.smk', 'b-pass_assembly.smk', 'b-pass_binning.smk', 'm-s.smk', 'm-p.smk', 'm.mtf.smk'],
        'mc':   ['b-pass_trim.smk', 'b-pass_assembly.smk', 'b-pass_binning.smk', 'm-s.smk', 'c.base.smk', 'mc.smk'],
        'mcr':  ['b-pass_trim.smk', 'b-pass_assembly.smk', 'b-pass_binning.smk', 'm-s.smk', 'c.base.smk', 'mcr.smk'],
        'pmr':  ['b-pass_trim.smk', 'b-pass_assembly.smk', 'b-pass_binning.smk', 'm-s.smk', 'm-p.smk', 'pmr_m.smk'],
        'v':    ['b-pass_trim.smk', 'b-pass_assembly.smk', 'v-s.smk', 'v-p.smk'],
        'vc':   ['b-pass_trim.smk', 'b-pass_assembly.smk', 'v-s.smk', 'c.base.smk', 'vc.smk'],
        'vcr':  ['b-pass_trim.smk', 'b-pass_assembly.smk', 'v-s.smk', 'c.base.smk', 'vcr.smk'],
        'vm':   ['b-pass_trim.smk', 'b-pass_assembly.smk', 'b-pass_binning.smk', 'v-s.smk', 'v-p.smk', 'm-s.smk', 'm-p.smk', 'vm.smk'],
        'vmc':  ['b-pass_trim.smk', 'b-pass_assembly.smk', 'b-pass_binning.smk', 'v-s.smk', 'm-s.smk', 'c.base.smk', 'vmc.smk'],
        'vmcr': ['b-pass_trim.smk', 'b-pass_assembly.smk', 'b-pass_binning.smk', 'v-s.smk', 'm-s.smk', 'c.base.smk', 'vmcr.smk'],
    },
}

_SINGLE_LIKE_READ_TYPES = {'single', 'single-rna', 'pacbio-raw', 'pacbio-corr', 'pacbio-hifi',
                           'nano-raw', 'nano-corr', 'nano-hq'}


# pmr --mags v: vMAG chain has no metawrap binning, so skip_level 3 == 2.
_PMR_V_SMKS = {
    0: ['b-trim.smk', 'b-assembly.smk', 'v-s.smk', 'v-p.smk', 'pmr_v.smk'],
    1: ['b-pass_trim.smk', 'b-assembly.smk', 'v-s.smk', 'v-p.smk', 'pmr_v.smk'],
    2: ['b-pass_trim.smk', 'b-pass_assembly.smk', 'v-s.smk', 'v-p.smk', 'pmr_v.smk'],
    3: ['b-pass_trim.smk', 'b-pass_assembly.smk', 'v-s.smk', 'v-p.smk', 'pmr_v.smk'],
}

# pmr --mags vm: binning + v-s/v-p + m-s/m-p + combined-reference pmr_vm chain.
_PMR_VM_SMKS = {
    0: ['b-trim.smk', 'b-assembly.smk', 'b-binning.smk', 'v-s.smk', 'v-p.smk', 'm-s.smk', 'm-p.smk', 'pmr_vm.smk'],
    1: ['b-pass_trim.smk', 'b-assembly.smk', 'b-binning.smk', 'v-s.smk', 'v-p.smk', 'm-s.smk', 'm-p.smk', 'pmr_vm.smk'],
    2: ['b-pass_trim.smk', 'b-pass_assembly.smk', 'b-binning.smk', 'v-s.smk', 'v-p.smk', 'm-s.smk', 'm-p.smk', 'pmr_vm.smk'],
    3: ['b-pass_trim.smk', 'b-pass_assembly.smk', 'b-pass_binning.smk', 'v-s.smk', 'v-p.smk', 'm-s.smk', 'm-p.smk', 'pmr_vm.smk'],
}

_PMR_MAGS_VALUES = ('m', 'v', 'vm')

# pmr mode loads a per-mags config template (template.config.pmr_{m,v,vm}.yaml)
# instead of the usual template.config.{mode}.yaml.
_PMR_MAGS_TEMPLATE = {'m': 'pmr_m', 'v': 'pmr_v', 'vm': 'pmr_vm'}

# pmr type-specific options: mMAG CheckM thresholds apply to --mags m/vm,
# the CheckV quality threshold applies to --mags v/vm.
_PMR_M_TYPE_OPTIONS = (
    'mmags_min_completeness1', 'mmags_max_contamination1',
    'mmags_min_completeness2', 'mmags_max_contamination2',
)
_PMR_V_TYPE_OPTIONS = ('vmags_checkv_quality',)


def _generated_config_name(mode, cli_overrides):
    """Filename of the config written into <workdir>/config/.

    pmr gets a per-mags name (config.pmr_{m,v,vm}.yaml) so configs generated
    for different --mags values in the same working directory do not
    overwrite each other. Invalid mags values are rejected later in
    make_config before anything is written.
    """
    if mode == 'pmr':
        mags = (cli_overrides or {}).get('mags') or 'm'
        return 'config.pmr_{}.yaml'.format(mags)
    return 'config.{}.yaml'.format(mode)


def _workflow_smk_files(mode, skip_level, mags=None):
    if skip_level not in WORKFLOW_MAP:
        raise ValueError(f"skip_level must be one of [0/1/2/3], got: {skip_level}")
    if mode == 'pmr':
        mags = mags or 'm'
        if mags == 'v':
            return list(_PMR_V_SMKS[skip_level])
        if mags == 'vm':
            return list(_PMR_VM_SMKS[skip_level])
        if mags != 'm':
            raise ValueError(f"mags must be one of {list(_PMR_MAGS_VALUES)}, got: {mags}")
    return WORKFLOW_MAP[skip_level][mode]


def _rules_dir_for(read_type):
    return 'rules_s' if read_type in _SINGLE_LIKE_READ_TYPES else 'rules_p'


def _targets_dir_for(read_type):
    return 'targets_s' if read_type in _SINGLE_LIKE_READ_TYPES else 'targets_p'


def _infer_skip_level(config):
    """skip_level recorded in the config wins; fall back to inferring it from
    workflow_includes for configs generated by older versions."""
    skip_level = config.get('skip_level')
    if skip_level is not None:
        return skip_level
    names = {os.path.basename(str(x)) for x in (config.get('workflow_includes') or [])}
    if 'b-pass_binning.smk' in names:
        return 3
    if 'b-pass_assembly.smk' in names:
        return 2
    if 'b-pass_trim.smk' in names:
        return 1
    return 0


def make_config(mode, input_tsv, threads, cli_overrides, config_f, db_dir=None,
                working_dir=None, project_name=None, cli_provided=None):

    # --mags must be settled before the mode template is loaded: pmr uses a
    # per-mags template (pmr_m/pmr_v/pmr_vm). Non-pmr modes ignore 'mags'.
    mags = None
    if mode == 'pmr':
        mags = cli_overrides.get('mags') or 'm'
        if mags not in _PMR_MAGS_VALUES:
            raise ValueError(
                "--mags must be one of {}, got: {}".format(list(_PMR_MAGS_VALUES), mags)
            )

    config_dir = Path(__file__).with_name('config')
    general_template = config_dir / 'template.config.general.yaml'
    _template_name = _PMR_MAGS_TEMPLATE[mags] if mode == 'pmr' else mode
    mode_template = config_dir / 'template.config.{}.yaml'.format(_template_name)

    if not general_template.exists():
        raise FileNotFoundError('missing general template {}'.format(general_template))
    if not mode_template.exists():
        raise FileNotFoundError('missing mode template {}'.format(mode_template))

    # Load general config and mode-specific config, then merge (mode overrides general)
    with open(general_template, 'r') as f:
        general_cfg = yaml.load(f)
    with open(mode_template, 'r') as f:
        mode_cfg = yaml.load(f)

    cfg = _merge_configs(general_cfg, mode_cfg)

    # Track which config values came from explicit CLI options for logging
    _cli_param_log = {}
    if cli_provided is None:
        cli_provided = set()

    def _log_config_source(name, option_name, config_key, value):
        if name in cli_provided:
            _cli_param_log[config_key] = (option_name, value)

    # DBDIR is stored under Directory.database so all database settings live
    # in one place. It can still be overridden from the command line via -d.
    _directory = cfg.setdefault('Directory', {})
    _db_section = _directory.setdefault('database', {})

    if db_dir is not None:
        _db_section['DBDIR'] = str(Path(db_dir).resolve())
        _log_config_source('db_dir', '-d/--db-dir',
                           'Directory.database.DBDIR', _db_section['DBDIR'])

    # Resolve kraken2 database size selection to concrete paths.
    # Priority: 1) explicit config value, 2) size recorded by download-database,
    # 3) default 'kraken2_raw'. Users can still override by providing the
    # *_DB_path keys explicitly.
    def _load_downloaded_sizes(db_dir):
        """Return size choices previously written by download-database."""
        size_file = os.path.join(db_dir, '.metalink_db_sizes.yaml')
        if not os.path.exists(size_file):
            return {}
        try:
            with open(size_file) as f:
                data = yaml.load(f)
            if isinstance(data, dict):
                return data
        except Exception:
            pass
        return {}

    _dbdir_for_sizes = str(_db_section.get('DBDIR', 'db'))
    if not os.path.isabs(_dbdir_for_sizes):
        _dbdir_for_sizes = str((Path(__file__).parent / _dbdir_for_sizes).resolve())
    _downloaded_sizes = _load_downloaded_sizes(_dbdir_for_sizes)

    _kraken2_tax_size = _db_section.get('kraken2_taxonomy_db_size', '') or \
                        _downloaded_sizes.get('kraken2_taxonomy_db_size', '') or \
                        'kraken2_raw'
    _kraken2_fun_size = _db_section.get('kraken2_function_db_size', '') or \
                        _downloaded_sizes.get('kraken2_function_db_size', '') or \
                        'kraken2_raw'

    # With a reduced kraken2 database (anything other than the full
    # kraken2_raw), default the kraken2 classification confidence to 0:
    # sparser k-mer coverage means the full-database thresholds (0.6/0.7)
    # would leave most sequences unclassified. Explicit CLI values are
    # applied afterwards and take precedence over this default.
    for _opt_name, (_yaml_path, _default, *_) in PARAMS_HELP.get(mode, {}).items():
        if _yaml_path[-1] != 'confidence':
            continue
        _section = _yaml_path[-2]
        if 'taxonomy_kraken2' in _section and _kraken2_tax_size != 'kraken2_raw':
            pass
        elif 'function_kraken2' in _section and _kraken2_fun_size != 'kraken2_raw':
            pass
        else:
            continue
        if mode == 'pmr':
            _yaml_path = _resolve_mags_placeholder(_yaml_path, mags)
        if _get_nested(cfg, _yaml_path) is not None:
            set_nested(cfg, _yaml_path, 0.0)
            logging.info(
                'Setting config %s to 0.0 (reduced kraken2 database: %s)',
                ' -> '.join(_yaml_path),
                _kraken2_tax_size if 'taxonomy' in _section else _kraken2_fun_size)

    # Resolve kraken2 database build versions from the bundle manifest.
    _kraken2_tax_version = _db_section.get('kraken2_taxonomy_db_version', '') or \
                           _get_nested(_DB_BUNDLE, ['metalink_built_databases', 'kraken2_taxonomy', 'version'], '') or \
                           'v20250930'
    _kraken2_fun_version = _db_section.get('kraken2_function_db_version', '') or \
                           _get_nested(_DB_BUNDLE, ['metalink_built_databases', 'kraken2_function', 'version'], '') or \
                           'v20250930'

    # Collect derived database paths so they can be inserted at a clean position
    # in the config (after cat_pack_db) instead of appended at the end.
    _derived_paths = []

    if not _db_section.get('kraken2_taxonomy_db'):
        _derived_paths.append(('kraken2_taxonomy_db',
                               'kraken2_taxonomy_db/{}/{}/'.format(_kraken2_tax_version, _kraken2_tax_size)))
    if not _db_section.get('kraken2_function_db'):
        _derived_paths.append(('kraken2_function_db',
                               'kraken2_function_db/{}/{}/'.format(_kraken2_fun_version, _kraken2_fun_size)))

    # Temporarily inject kraken2 DB paths so dependent paths can use them.
    for _key, _val in _derived_paths:
        _db_section[_key] = _val

    def _derive_if_missing(key, base_key, *rel_parts, trailing_slash=False):
        if _db_section.get(key):
            return
        base = _db_section.get(base_key, '')
        if not base:
            return
        derived = os.path.join(base, *rel_parts)
        if trailing_slash and not derived.endswith('/'):
            derived += '/'
        _derived_paths.append((key, derived))
        # Make it available for further derivations immediately.
        _db_section[key] = derived

    # Remove the temporarily injected keys so we can re-insert them in order.
    for _key, _val in _derived_paths:
        if _key in _db_section:
            del _db_section[_key]

    # Insert derived paths right after cat_pack_db to keep the config tidy.
    # Use single-quoted scalars so they match the surrounding path entries.
    _insert_pos = list(_db_section.keys()).index('cat_pack_db') + 1
    for _key, _val in _derived_paths:
        _db_section.insert(_insert_pos, _key, SingleQuotedScalarString(_val))
        _insert_pos += 1

    # Resolve database paths that are relative to DBDIR.
    # Absolute paths are preserved so users can override individual tool databases.
    _db_section.setdefault('DBDIR', 'db')
    _dbdir = str(_db_section['DBDIR']).rstrip('/')
    for _db_key, _db_val in list(_db_section.items()):
        if _db_key == 'DBDIR':
            continue
        # Size selectors are plain strings, not filesystem paths.
        if _db_key.endswith('_db_size'):
            continue
        if not isinstance(_db_val, str):
            continue
        if os.path.isabs(_db_val):
            continue
        # Resolve relative paths (including empty string) under DBDIR.
        _db_section[_db_key] = os.path.join(_dbdir, _db_val) if _db_val else _dbdir

    # Resolve the Metalink scripts directory if it is relative.
    # Absolute paths are preserved so users can point to a custom installation.
    _script_section = _directory.setdefault('script', {})
    _script_path = _script_section.get('metalink', 'scripts/')
    if isinstance(_script_path, str) and not os.path.isabs(_script_path):
        _script_section['metalink'] = str((Path(__file__).parent / _script_path).resolve())

    # Resolve the KEGG module KO info table path if it is relative.
    # Default location is under the Metalink installation scripts directory.
    _kegg_path = _script_section.get(
        'kegg_modules_ko_info_path',
        'scripts/mclink-scripts/modules_ko_info.tsv'
    )
    if isinstance(_kegg_path, str) and not os.path.isabs(_kegg_path):
        _script_section['kegg_modules_ko_info_path'] = str((Path(__file__).parent / _kegg_path).resolve())

    # Resolve the KO reference length table path if it is relative.
    # The table ships with the Metalink scripts directory, so it is available
    # for every reads function classification method (kraken2 or emapper),
    # independent of the kraken2 function database.
    _ko_len_path = _script_section.get(
        'KO_length_ref_path',
        'scripts/KO_faa_seq_length_summary.tsv.gz'
    )
    if isinstance(_ko_len_path, str) and not os.path.isabs(_ko_len_path):
        _script_section['KO_length_ref_path'] = str((Path(__file__).parent / _ko_len_path).resolve())

    cfg['mode'] = mode
    if mode == 'pmr':
        cfg['mags'] = mags
        _log_config_source('mags', '--mags', 'mags', mags)
    cfg['input_tsv'] = str(Path(input_tsv).resolve())
    _log_config_source('input_tsv', '-i/--input-tsv', 'input_tsv', cfg['input_tsv'])
    cfg['threads'] = threads

    cfg['run_project'] = cli_overrides.get('run_project', False)
    _log_config_source('run_project', '--run-project', 'run_project', cfg['run_project'])

    # keep_intermediate is a per-category dict in the template; fall back to
    # the built-in defaults if a (legacy) config made it a bool or scalar.
    _KI_DEFAULTS = {'bam': True, 'drep_data': False, 'vs2_iter0': False, 'assembly': False}
    _ki = cfg.get('keep_intermediate')
    cfg['keep_intermediate'] = dict(_ki) if isinstance(_ki, dict) else dict(_KI_DEFAULTS)
    for _flag, _cat in (('delete_bam', 'bam'),
                        ('keep_drep_data', 'drep_data'),
                        ('keep_vs2_iter0', 'vs2_iter0'),
                        ('keep_assembly_intermediate', 'assembly')):
        if cli_overrides.get(_flag, False):
            cfg['keep_intermediate'][_cat] = _flag.startswith('keep_')
        _log_config_source(_flag, '--' + _flag.replace('_', '-'),
                           'keep_intermediate.' + _cat, cfg['keep_intermediate'][_cat])

    _rt = cli_overrides.get('read_type', 'paired')
    _valid_read_types = (
        'paired', 'single', 'paired-rna', 'single-rna',
        'pacbio-raw', 'pacbio-corr', 'pacbio-hifi',
        'nano-raw', 'nano-corr', 'nano-hq'
    )
    if _rt not in _valid_read_types:
        raise ValueError(
            "--read-type must be one of {}, got: {}".format(_valid_read_types, _rt)
        )
    cfg['read_type'] = _rt
    _log_config_source('read_type', '--read-type', 'read_type', cfg['read_type'])

    _taxonomy_rank = cli_overrides.get('taxonomy_rank') or cfg.get('taxonomy_rank') or 'Phylum,Class,Order,Family,Genus,Species'
    cfg['taxonomy_rank'] = [r.strip() for r in _taxonomy_rank.split(',')] if isinstance(_taxonomy_rank, str) else _taxonomy_rank
    _log_config_source('taxonomy_rank', '--taxonomy-rank', 'taxonomy_rank', cfg['taxonomy_rank'])

    if working_dir is not None:
        cfg.setdefault('Directory', {})['work'] = str(Path(working_dir).resolve())

    if project_name is not None:
        cfg['name_project'] = [project_name]
        _log_config_source('project_name', '--project-name', 'name_project', cfg['name_project'])

    # skip_level affects which smk files are included; record it in the config
    _skip_level = cli_overrides.get('skip_level', 0)
    if _skip_level not in WORKFLOW_MAP:
        raise ValueError(f"--skip-level must be one of [0/1/2/3], got: {_skip_level}")
    cfg['skip_level'] = _skip_level
    _log_config_source('skip_level', '--skip-level', 'skip_level', _skip_level)

    # concurrent cap for kraken2 taxonomy jobs; applied by the Snakefile as a
    # Snakemake global resource limit (kraken2_taxonomy). Absent = unlimited.
    _ktj = cli_overrides.get('kraken2_taxonomy_jobs')
    if _ktj is not None:
        cfg['kraken2_taxonomy_jobs'] = int(_ktj)
        _log_config_source('kraken2_taxonomy_jobs', '--kraken2-taxonomy-jobs', 'kraken2_taxonomy_jobs', _ktj)

    for opt_name, (yaml_path, default, *_ ) in PARAMS_HELP[mode].items():
        user_val = cli_overrides.get(opt_name)
        if user_val is not None:
            option_name = f"--{opt_name.replace('_', '-')}"
            if mode == 'pmr':
                yaml_path = _resolve_mags_placeholder(yaml_path, mags)
                if not _pmr_option_applicable(opt_name, yaml_path, cfg, mags):
                    continue
            logging.info(f'Setting config {" -> ".join(yaml_path)} to {user_val!r} (from CLI {option_name})')
            set_nested(cfg, yaml_path, user_val)

    for config_key, (option_name, value) in sorted(_cli_param_log.items()):
        logging.info(f'Setting config {config_key} to {value!r} (from CLI {option_name})')

    # skip_level was validated above; it determines which smk files are included
    skip_level = _skip_level
    smk_files = _workflow_smk_files(mode, skip_level, mags)

    # Validate read_type and assembly_method compatibility for modes that actually run assembly
    _validate_read_type_assembly_method(cfg, mode, skip_level)
    _validate_bracken_level(cfg, mode)

    _rules_dir = _rules_dir_for(cfg['read_type'])
    _targets_dir = _targets_dir_for(cfg['read_type'])
    cfg['workflow_common_includes'] = ["{}/init.smk".format(_rules_dir)]
    cfg['workflow_includes'] = ["{}/{}".format(_rules_dir, f) for f in smk_files]
    cfg['workflow_target'] = "{}/targets.{}.smk".format(_targets_dir, mode)
    cfg['workflow_snakefile'] = str(Path(__file__).parent / 'Snakefile')

    ordered_cfg = CommentedMap()
    for key in BASE_KEYS:
        if key in cfg:
            ordered_cfg[key] = cfg.pop(key)
    for p in ['Params_b_sample', 'Params_b_project']:
        if p in cfg:
            ordered_cfg[p] = cfg.pop(p)
    _sections_key = 'pmr_{}'.format(mags) if mode == 'pmr' else mode
    for p in MODE_SPECIFIC.get(_sections_key, []):
        if p in cfg:
            ordered_cfg[p] = cfg.pop(p)
    for key in sorted(cfg.keys()):
        ordered_cfg[key] = cfg.pop(key)

    cfg = ordered_cfg

    # Enforce parameter compatibility before writing the config file.
    if mode == 'pmr':
        _validate_pmr_iphop_custom_mags_dir(cfg, mags)
    _validate_cami_gtdbtk_compatibility(cfg, mode)

    with open(config_f, 'w') as f:
        f.write('# Metalink v{}\n'.format(__version__))
        f.write('# Note: when using this config with -c/--configfile, config-source options\n')
        f.write('#       (-i, -d, -w, --read-type, --skip-level)\n')
        f.write('#       must be edited here; they cannot be overridden from the command line.\n')
        f.write('#       Mode-specific execution parameters (e.g., --trim-sample-threads) can still be\n')
        f.write('#       overridden from the command line.\n')
        # Keep long database paths on a single line instead of folding them.
        _old_width = getattr(yaml, 'width', None)
        yaml.width = 4096
        yaml.dump(cfg, f)
        yaml.width = _old_width
    _compact_yaml(config_f)

def load_configfile(path):
    with open(path, 'r') as f:
        cfg = yaml.load(f)
    if cfg is None:
        raise ValueError(f"config file {path} is empty or with error format. Please check if template.config.general.yaml and template.config.<mode>.yaml are present and non-empty in the config directory.")
    return cfg

from metalink import __version__

# Make the package version visible to workflow scripts (R/Python helpers run
# via Snakemake shell jobs) so output headers can record it without
# hardcoding.  Scripts read METALINK_VERSION and fall back to "unknown" when
# invoked outside the CLI (e.g. direct snakemake runs).
os.environ.setdefault('METALINK_VERSION', __version__)

PARAMS_HELP = {
    'b': {
        'raw_stats_sample_threads': (
            ['Params_b_sample', '00_raw_stats', 'threads'], 4, int,
            '00_raw_stats: Threads for seqkit.'),
        'trim_sample_threads': (
            ['Params_b_sample', '01_trim', 'threads'], 16, int,
            '01_trim: Threads for fastp (next-generation short reads) or fastplong (third-generation long reads).'),
        'trim_length_required': (
            ['Params_b_sample', '01_trim', 'length_required'], 15, int,
            '01_trim: Reads shorter than length_required will be discarded'),
        'trim_n_base_limit': (
            ['Params_b_sample', '01_trim', 'n_base_limit'], 5, int,
            "01_trim: if one read's number of N base is >n_base_limit, then this read/pair is discarded."),
        'trim_compression': (
            ['Params_b_sample', '01_trim', 'compression'], 4, int,
            '01_trim: Compression level for gzip output [1-9]. 1 is fastest, 9 is smallest.'),
        'trim_extra_args': (
            ['Params_b_sample', '01_trim', 'extra_args'], '', str,
            '01_trim: Extra arguments forwarded to fastp or fastplong for raw reads trimming (e.g., "--reads_to_process 0 -Q").'),
        'trim_stats_sample_threads': (
            ['Params_b_sample', '01_trim_stats', 'threads'], 16, int,
            '01_trim_stats: Threads.'),
        'assembly_method': (
            ['Params_b_sample', '02_assembly', 'method'], 'megahit', str,
            '\b\n02_assembly: Underlying assembly software using one of [megahit/spades/metaflye].\n'
            'Must be compatible with --read-type (config key read_type):\n'
            'megahit: MEGAHIT (second-generation sequencing; for paired or single DNA reads).\n'
            'spades: SPAdes (second-generation sequencing; metaSPAdes for paired DNA reads, rnaSPAdes for paired-rna or single-rna RNA reads).\n'
            'metaflye: metaFlye (third-generation sequencing; required for pacbio-raw/corr/hifi or nano-raw/corr/hq long reads).'),
        'assembly_threads': (
            ['Params_b_sample', '02_assembly', 'threads'], 40, int,
            '02_assembly: Threads for [megahit/spades/metaflye].'),
        'megahit_presets': (
            ['Params_b_sample', '02_megahit', 'presets'], 'meta-sensitive', str,
            '\b\n02_megahit: MEGAHIT presets using one of [meta-sensitive/meta-large].\n'
            "meta-sensitive: '--min-count 1 --k-list 21,29,39,49,...,129,141'.\n"
            "meta-large: '--k-min 27 --k-max 127 --k-step 10'.\n"
            '(large & complex metagenomes, like soil).'),
        'megahit_k_list': (
            ['Params_b_sample', '02_megahit', 'k_list'], '', str,
            '02_megahit: MEGAHIT --k-list (optional). Explicit comma-separated k-mer sizes (e.g. "21,29,39,59,79,99,119,141"). Overrides --megahit-presets. If set, --megahit-k-min, --megahit-k-max, --megahit-k-step must NOT be set. Optional --megahit-min-count can be used together.'),
        'megahit_min_count': (
            ['Params_b_sample', '02_megahit', 'min_count'], 2, int,
            '02_megahit: MEGAHIT --min-count (optional). Minimum multiplicity for filtering (k_min+1)-mers. Can be used with --megahit-k-list or --megahit-k-min/--k-max/--k-step.'),
        'megahit_k_min': (
            ['Params_b_sample', '02_megahit', 'k_min'], '', str,
            '02_megahit: MEGAHIT --k-min (optional). Minimum kmer size (<= 255), must be odd number. If set with --k-max and --k-step, overrides --megahit-presets and uses custom k-mer mode.'),
        'megahit_k_max': (
            ['Params_b_sample', '02_megahit', 'k_max'], '', str,
            '02_megahit: MEGAHIT --k-max (optional). Maximum kmer size (<= 255), must be odd number.'),
        'megahit_k_step': (
            ['Params_b_sample', '02_megahit', 'k_step'], '', str,
            '02_megahit: MEGAHIT --k-step (optional). Increment of kmer size of each iteration (<= 28), must be even number.'),
        'megahit_extra_args': (
            ['Params_b_sample', '02_megahit', 'extra_args'], '', str,
            '02_megahit: Extra arguments forwarded to MEGAHIT (e.g. "--no-mercy --tmp-dir /path/to/your/tmp/dir"). Do NOT use this for k-mer settings (--k-list, --k-min, --presets, etc.).'),
        'spades_extra_args': (
            ['Params_b_sample', '02_spades', 'extra_args'], '', str,
            '02_spades: Extra arguments forwarded to SPAdes. Do NOT use this for input/output or thread settings.'),
        'metaflye_extra_args': (
            ['Params_b_sample', '02_metaflye', 'extra_args'], '', str,
            '02_metaflye: Extra arguments forwarded to Flye for third-generation long reads assembly (e.g., "--iterations 2"). Do NOT use this for read type or input/output settings.'),
        'mapping_threads': (
            ['Params_b_sample', '03_mapping', 'threads'], 20, int,
            '03_mapping: Threads for [minibwa/bwamem/bwamem2/strobealign].'),
        'mapping_method': (
            ['Params_b_sample', '03_mapping', 'method'], 'minibwa', str,
            '\b\n03_mapping: Underlying mapping software used [minibwa/bwamem/bwamem2/strobealign].\n'
            'minibwa: minibwa using default parameters.\n'
            'bwamem: bwa mem using default parameters.\n'
            'bwamem2: bwa-mem2 using default parameters.\n'
            'strobealign: strobealign using default parameters.'),
        'mapping_samtools_sort_mem': (
            ['Params_b_sample', '03_mapping', 'samtools_sort_mem'], '768M', str,
            '03_mapping: Set maximum memory per thread; suffix K/M/G recognized, e.g. 2G, 768M.'),
        'mapping_extra_args': (
            ['Params_b_sample', '03_mapping', 'extra_args'], '', str,
            '03_mapping: Extra arguments forwarded to the active mapping aligner (minibwa map, bwa mem, bwa-mem2 mem, or strobealign --use-index), (e.g. "-B 8" for minibwa map).'),
        'binning_metawrap_depth_threads': (
            ['Params_b_sample', '05_metawrap_binning_depth', 'threads'], 8, int,
            '05_metawrap_binning_depth: Threads.'),
        'binning_min_read_percent_identity': (
            ['Params_b_sample', '05_metawrap_binning_depth', 'min_read_percent_identity'], 97, float,
            '05_metawrap_binning_depth: Exclude reads by overall percent identity [0-100] during the binning process, e.g. 97 for 97%.'),
        'binning_contig_end_exclusion': (
            ['Params_b_sample', '05_metawrap_binning_depth', 'contig_end_exclusion'], 75, int,
            '05_metawrap_binning_depth: Exclude bases at the ends of contig sequences used for depth calculation and bining.'),
        'binning_min_contig_length': (
            ['Params_b_sample', '05_metawrap_binning_depth', 'min_contig_length'], 1500, int,
            '05_metawrap_binning_depth: Minimum contig length for binning.'),
        'binning_tools': (
            ['Params_b_sample', '05_metawrap_binning_depth', 'binning_tools'], 'metabat2,maxbin2,concoct', str,
            '\b\n05_metawrap_binning_depth: Comma-separated list of binning tools, choose 2 or 3 from [metabat2,maxbin2,concoct].\n'
            'At least 2 tools recommended for bin_refinement.'),
        'binning_max_memory_gb': (
            ['Params_b_sample', '05_metawrap_binning_depth', 'max_memory_gb'], 512, int,
            '05_metawrap_binning_depth: Maximum memory (GB) for MetaWRAP binning.'),
        'binning_metawrap_extra_args': (
            ['Params_b_sample', '05_metawrap_binning_depth', 'extra_args'], '', str,
            '05_metawrap_binning_depth: Extra arguments forwarded to metaWRAP binning (e.g., "-m 5").'),
        'metawrap_refinement_max_memory_gb': (
            ['Params_b_sample', '06_metawrap_bin_refinement', 'max_memory_gb'], 512, int,
            '06_metawrap_bin_refinement: Maximum memory in GB for MetaWRAP bin_refinement.'),
        'metawrap_refinement_extra_args': (
            ['Params_b_sample', '06_metawrap_bin_refinement', 'extra_args'], '', str,
            '06_metawrap_bin_refinement: Extra arguments forwarded to metaWRAP bin_refinement (e.g., "--remove-ambiguous").'),
        'bin_checkm_threads': (
            ['Params_b_sample', '06_bin_checkm', 'threads'], 8, int,
            '06_bin_checkm or 06_metawrap_bin_refinement: Threads.'),
        'bin_min_completeness': (
            ['Params_b_sample', '06_bin_checkm', 'min_completeness'], 50, float,
            '06_bin_checkm or 06_metawrap_bin_refinement: Minimum completeness (%) bin filtering, e.g. 50 for 50%.'),
        'bin_max_contamination': (
            ['Params_b_sample', '06_bin_checkm', 'max_contamination'], 10, float,
            '06_bin_checkm or 06_metawrap_bin_refinement: Maximum contamination (%) bin filtering, e.g. 10 for 10%.'),
        'bin_min_quality_score': (
            ['Params_b_sample', '06_bin_checkm', 'min_quality_score'], 50, float,
            '06_bin_checkm: Minimum quality score for bin filtering (GTDB-like), QS = completeness - 5*contamination, e.g. 50 for QS>=50.'),
        'contigs_depth_sample_threads': (
            ['Params_b_sample', '07_contigs_depth', 'threads'], 8, int,
            '07_contigs_depth: Threads.'),
        'depth_min_read_percent_identity': (
            ['Params_b_sample', '07_contigs_depth', 'min_read_percent_identity'], 95, float,
            '\b\n07_contigs_depth: Exclude reads by overall percent identity [0-100] for contig depth estimation\n'
            'using jgi_summarize_bam_contig_depths --percentIdentity, e.g. 95 for 95%.\n'
            'Reads with percent identity below this threshold are excluded from depth calculation.\n'
            'Higher values increase stringency, retaining only high-quality alignments.'),
        'depth_contig_end_exclusion': (
            ['Params_b_sample', '07_contigs_depth', 'contig_end_exclusion'], 75, int,
            '\b\n07_contigs_depth: Exclude bases at the ends of contig sequences for depth estimation using jgi_summarize_bam_contig_depths --maxEdgeBases.\n'
            'Bases within this distance from either end of a contig are ignored when calculating depth,\n'
            'preventing edge effects from partial alignments.'),
        'depth_min_contig_length': (
            ['Params_b_sample', '07_contigs_depth', 'min_contig_length'], 300, int,
            '\b\n07_contigs_depth: Minimum contig length for depth estimation [bp], e.g. 300 for 300 bp.\n'
            'Contigs shorter than this threshold are excluded from depth calculation.'),
        'contigs_prodigal_extra_args': (
            ['Params_b_sample', '07_contigs_depth', 'extra_args'], '', str,
            '07_contigs_depth: Extra arguments forwarded to prodigal for contig gene prediction (e.g., "-m -n").'),
        'bin_prodigal_threads': (
            ['Params_b_sample', '08_bin_prodigal', 'threads'], 8, int,
            '08_bin_prodigal: Threads.'),
        'bin_prodigal_extra_args': (
            ['Params_b_sample', '08_bin_prodigal', 'extra_args'], '', str,
            '08_bin_prodigal: Extra arguments forwarded to prodigal for MAG gene prediction (e.g., "-m -n").'),
        'raw_stats_project_threads': (
            ['Params_b_project', '01_cat_raw_stats', 'threads'], 4, int,
            '01_cat_raw_stats: Threads.'),
        'trimmed_stats_project_threads': (
            ['Params_b_project', '02_cat_trimmed_stats', 'threads'], 4, int,
            '02_cat_trimmed_stats: Threads.'),
        'contigs_stats_project_threads': (
            ['Params_b_project', '03_contigs_stats', 'threads'], 4, int,
            '03_contigs_stats: Threads.'),
        'merge_stats_project_threads': (
            ['Params_b_project', '04_merge_reads_contigs_stats', 'threads'], 4, int,
            '04_merge_reads_contigs_stats: Threads.'),
        'gather_checkm_project_threads': (
            ['Params_b_project', '05_gather_project_bin_checkm_stats', 'threads'], 4, int,
            '05_gather_project_bin_checkm_stats: Threads.'),
    },
    'c': {
        'raw_stats_sample_threads': (
            ['Params_b_sample', '00_raw_stats', 'threads'], 4, int,
            '00_raw_stats: Threads for seqkit.'),
        'trim_sample_threads': (
            ['Params_b_sample', '01_trim', 'threads'], 16, int,
            '01_trim: Threads for fastp (next-generation short reads) or fastplong (third-generation long reads).'),
        'trim_length_required': (
            ['Params_b_sample', '01_trim', 'length_required'], 15, int,
            '01_trim: Reads shorter than length_required will be discarded'),
        'trim_n_base_limit': (
            ['Params_b_sample', '01_trim', 'n_base_limit'], 5, int,
            "01_trim: if one read's number of N base is >n_base_limit, then this read/pair is discarded."),
        'trim_compression': (
            ['Params_b_sample', '01_trim', 'compression'], 4, int,
            '01_trim: Compression level for gzip output [1-9]. 1 is fastest, 9 is smallest.'),
        'trim_extra_args': (
            ['Params_b_sample', '01_trim', 'extra_args'], '', str,
            '01_trim: Extra arguments forwarded to fastp or fastplong for raw reads trimming (e.g., "--reads_to_process 0 -Q").'),
        'trim_stats_sample_threads': (
            ['Params_b_sample', '01_trim_stats', 'threads'], 16, int,
            '01_trim_stats: Threads.'),
        'assembly_method': (
            ['Params_b_sample', '02_assembly', 'method'], 'megahit', str,
            '\b\n02_assembly: Underlying assembly software using one of [megahit/spades/metaflye].\n'
            'Must be compatible with --read-type (config key read_type):\n'
            'megahit: MEGAHIT (second-generation sequencing; for paired or single DNA reads).\n'
            'spades: SPAdes (second-generation sequencing; metaSPAdes for paired DNA reads, rnaSPAdes for paired-rna or single-rna RNA reads).\n'
            'metaflye: metaFlye (third-generation sequencing; required for pacbio-raw/corr/hifi or nano-raw/corr/hq long reads).'),
        'assembly_threads': (
            ['Params_b_sample', '02_assembly', 'threads'], 40, int,
            '02_assembly: Threads for [megahit/spades/metaflye].'),
        'megahit_presets': (
            ['Params_b_sample', '02_megahit', 'presets'], 'meta-sensitive', str,
            '\b\n02_megahit: MEGAHIT presets using one of [meta-sensitive/meta-large].\n'
            "meta-sensitive: '--min-count 1 --k-list 21,29,39,49,...,129,141'.\n"
            "meta-large: '--k-min 27 --k-max 127 --k-step 10'.\n"
            '(large & complex metagenomes, like soil).'),
        'megahit_k_list': (
            ['Params_b_sample', '02_megahit', 'k_list'], '', str,
            '02_megahit: MEGAHIT --k-list (optional). Explicit comma-separated k-mer sizes (e.g. "21,29,39,59,79,99,119,141"). Overrides --megahit-presets. If set, --megahit-k-min, --megahit-k-max, --megahit-k-step must NOT be set. Optional --megahit-min-count can be used together.'),
        'megahit_min_count': (
            ['Params_b_sample', '02_megahit', 'min_count'], 2, int,
            '02_megahit: MEGAHIT --min-count (optional). Minimum multiplicity for filtering (k_min+1)-mers. Can be used with --megahit-k-list or --megahit-k-min/--k-max/--k-step.'),
        'megahit_k_min': (
            ['Params_b_sample', '02_megahit', 'k_min'], '', str,
            '02_megahit: MEGAHIT --k-min (optional). Minimum kmer size (<= 255), must be odd number. If set with --k-max and --k-step, overrides --megahit-presets and uses custom k-mer mode.'),
        'megahit_k_max': (
            ['Params_b_sample', '02_megahit', 'k_max'], '', str,
            '02_megahit: MEGAHIT --k-max (optional). Maximum kmer size (<= 255), must be odd number.'),
        'megahit_k_step': (
            ['Params_b_sample', '02_megahit', 'k_step'], '', str,
            '02_megahit: MEGAHIT --k-step (optional). Increment of kmer size of each iteration (<= 28), must be even number.'),
        'megahit_extra_args': (
            ['Params_b_sample', '02_megahit', 'extra_args'], '', str,
            '02_megahit: Extra arguments forwarded to MEGAHIT (e.g. "--no-mercy --tmp-dir /path/to/your/tmp/dir"). Do NOT use this for k-mer settings (--k-list, --k-min, --presets, etc.).'),
        'spades_extra_args': (
            ['Params_b_sample', '02_spades', 'extra_args'], '', str,
            '02_spades: Extra arguments forwarded to SPAdes. Do NOT use this for input/output or thread settings.'),
        'metaflye_extra_args': (
            ['Params_b_sample', '02_metaflye', 'extra_args'], '', str,
            '02_metaflye: Extra arguments forwarded to Flye for third-generation long reads assembly (e.g., "--iterations 2"). Do NOT use this for read type or input/output settings.'),
        'mapping_threads': (
            ['Params_b_sample', '03_mapping', 'threads'], 20, int,
            '03_mapping: Threads for [minibwa/bwamem/bwamem2/strobealign].'),
        'mapping_method': (
            ['Params_b_sample', '03_mapping', 'method'], 'minibwa', str,
            '\b\n03_mapping: Underlying mapping software used [minibwa/bwamem/bwamem2/strobealign].\n'
            'minibwa: minibwa using default parameters.\n'
            'bwamem: bwa mem using default parameters.\n'
            'bwamem2: bwa-mem2 using default parameters.\n'
            'strobealign: strobealign using default parameters.'),
        'mapping_samtools_sort_mem': (
            ['Params_b_sample', '03_mapping', 'samtools_sort_mem'], '768M', str,
            '03_mapping: Set maximum memory per thread; suffix K/M/G recognized, e.g. 2G, 768M.'),
        'mapping_extra_args': (
            ['Params_b_sample', '03_mapping', 'extra_args'], '', str,
            '03_mapping: Extra arguments forwarded to the active mapping aligner (minibwa map, bwa mem, bwa-mem2 mem, or strobealign --use-index), (e.g. "-B 8" for minibwa map).'),
        'contigs_depth_sample_threads': (
            ['Params_b_sample', '07_contigs_depth', 'threads'], 8, int,
            '07_contigs_depth: Threads.'),
        'depth_min_read_percent_identity': (
            ['Params_b_sample', '07_contigs_depth', 'min_read_percent_identity'], 95, float,
            '\b\n07_contigs_depth: Exclude reads by overall percent identity [0-100] for contig depth estimation\n'
            'using jgi_summarize_bam_contig_depths --percentIdentity, e.g. 95 for 95%.\n'
            'Reads with percent identity below this threshold are excluded from depth calculation.\n'
            'Higher values increase stringency, retaining only high-quality alignments.'),
        'depth_contig_end_exclusion': (
            ['Params_b_sample', '07_contigs_depth', 'contig_end_exclusion'], 75, int,
            '\b\n07_contigs_depth: Exclude bases at the ends of contig sequences for depth estimation using jgi_summarize_bam_contig_depths --maxEdgeBases.\n'
            'Bases within this distance from either end of a contig are ignored when calculating depth,\n'
            'preventing edge effects from partial alignments.'),
        'depth_min_contig_length': (
            ['Params_b_sample', '07_contigs_depth', 'min_contig_length'], 300, int,
            '\b\n07_contigs_depth: Minimum contig length for depth estimation [bp], e.g. 300 for 300 bp.\n'
            'Contigs shorter than this threshold are excluded from depth calculation.'),
        'contigs_prodigal_extra_args': (
            ['Params_b_sample', '07_contigs_depth', 'extra_args'], '', str,
            '07_contigs_depth: Extra arguments forwarded to prodigal for contig gene prediction (e.g., "-m -n").'),
        'raw_stats_project_threads': (
            ['Params_b_project', '01_cat_raw_stats', 'threads'], 4, int,
            '01_cat_raw_stats: Threads.'),
        'trimmed_stats_project_threads': (
            ['Params_b_project', '02_cat_trimmed_stats', 'threads'], 4, int,
            '02_cat_trimmed_stats: Threads.'),
        'contigs_stats_project_threads': (
            ['Params_b_project', '03_contigs_stats', 'threads'], 4, int,
            '03_contigs_stats: Threads.'),
        'merge_stats_project_threads': (
            ['Params_b_project', '04_merge_reads_contigs_stats', 'threads'], 4, int,
            '04_merge_reads_contigs_stats: Threads.'),
        'contigs_taxonomy_classify_method': (
            ['Params_c_sample', '01_taxonomy', 'contigs_taxonomy_classify_method'], 'kraken2', str,
            '01_taxonomy: Contigs taxonomy classifier using one of [CAT_pack/kraken2].'
        ),
        'contigs_taxonomy_threads': (
            ['Params_c_sample', '01_taxonomy', 'threads'], 60, int,
            '01_taxonomy: Threads for [CAT_pack/kraken2].'
        ),
        'contigs_taxonomy_cat_range': (
            ['Params_c_sample', '01_taxonomy_CAT_pack', 'range'], 10, float,
            '\b\n01_taxonomy_CAT_pack: Percentage range [0-100] of DIAMOND hits relative to top bit-score to include in LCA calculation.\n'
            'Higher values yield more conservative classifications by incorporating divergent hits, while lower values enable specific assignments.\n'
            'Applicable when --contigs-taxonomy-classify-method is set to CAT_pack.'
        ),
        'contigs_taxonomy_cat_fraction': (
            ['Params_c_sample', '01_taxonomy_CAT_pack', 'fraction'], 0.50, float,
            '\b\n01_taxonomy_CAT_pack: Minimum bit-score support fraction [0-0.99].\n'
            'Higher values enforce stricter consensus (one classification per sequence), while lower values allow multiple lineage detections.\n'
            'Applicable when --contigs-taxonomy-classify-method is set to CAT_pack.'
        ),
        'contigs_taxonomy_cat_pack_extra_args': (
            ['Params_c_sample', '01_taxonomy_CAT_pack', 'extra_args'], '', str,
            '01_taxonomy_CAT_pack: Extra arguments forwarded to CAT_pack contigs (e.g., "--block_size 12.0").'),
        'contigs_taxonomy_kraken2_confidence': (
            ['Params_c_sample', '01_taxonomy_kraken2', 'confidence'], 0.7, float,
            '\b\n01_taxonomy_kraken2: Confidence threshold [0-1] for taxonomic classification, applicable when --contigs-taxonomy-classify-method is set to kraken2.\n'
            'Default value was tuned on the CAMI2 marine dataset. For samples with higher complexity, a lower confidence threshold is usually recommended.'
        ),
        'contigs_taxonomy_kraken2_extra_args': (
            ['Params_c_sample', '01_taxonomy_kraken2', 'extra_args'], '', str,
            '01_taxonomy_kraken2: Extra arguments forwarded to kraken2 (e.g., "--memory-mapping").'),
        'contigs_function_classify_method': (
            ['Params_c_sample', '02_function', 'contigs_function_classify_method'], 'kraken2', str,
            '02_function: Contigs function classifier using one of [kraken2/emapper].'
        ),
        'contigs_function_threads': (
            ['Params_c_sample', '02_function', 'threads'], 60, int,
            '02_function: Threads for [kraken2/emapper].'),
        'contigs_function_emapper_evalue': (
            ['Params_c_sample', '02_function_emapper', 'seed_ortholog_evalue'], 1e-10, float,
            '02_function_emapper: Report only alignments with e-value below or equal to the threshold, applicable when --contigs-function-classify-method is set to emapper.'
        ),
        'contigs_function_emapper_method': (
            ['Params_c_sample', '02_function_emapper', 'method'], 'diamond', str,
            '\b\n02_function_emapper: search seed orthologs using one of [diamond/mmseqs], applicable when --contigs-function-classify-method is set to emapper.\n'
            'diamond: using default parameters of diamond mode of emapper.\n'
            'mmseqs: using default parameters of mmseqs mode of emapper'
        ),
        'contigs_function_emapper_query_cover': (
            ['Params_c_sample', '02_function_emapper', 'query_cover'], 30, float,
            '\b\n02_function_emapper: Report only alignments with query cover above or equal to the given percentage [0-100],.\n'
            'applicable when --contigs-function-classify-method is set to emapper'
        ),
        'contigs_function_emapper_extra_args': (
            ['Params_c_sample', '02_function_emapper', 'extra_args'], '', str,
            '02_function_emapper: Extra arguments forwarded to eggNOG-mapper (e.g., "--pident 30").'),
        'contigs_function_kraken2_confidence': (
            ['Params_c_sample', '02_function_kraken2', 'confidence'], 0.05, float,
            '02_function_kraken2: Confidence threshold [0-1] for functional classification, applicable when --contigs-function-classify-method is set to kraken2.'
        ),
        'contigs_function_kraken2_extra_args': (
            ['Params_c_sample', '02_function_kraken2', 'extra_args'], '', str,
            '02_function_kraken2: Extra arguments forwarded to kraken2 (e.g., "--memory-mapping").'),
        'merge_min_orf_length': (
            ['Params_c_sample', '05_merge_tax_fun_depth', 'minORFLen'], 50, int,
            '05_merge_tax_fun_depth: Minimum ORF length for finally used ORFs of contigs.'
        ),
        'merge_multiple_ko_handling_mode': (
            ['Params_c_sample', '05_merge_tax_fun_depth', 'KO_mode'], 'ignore_multiple_KOs', str,
            '\b\n05_merge_tax_fun_depth: Method using one of [split_multiple_KOs/keep_multiple_KOs/ignore_multiple_KOs] when an ORF has multiple KO annotations:.\n'
            'split_multiple_KOs: evenly distributing the total abundance across all assigned KOs;\n'
            'keep_multiple_KOs: retain original multiple KO annotations;\n'
            'ignore_multiple_KOs: ignore the ORFs with multiple KO annotations and only use the ORFs with single KO annotation.'
        ),
        'merge_min_ko_num_per_taxon': (
            ['Params_c_sample', '05_merge_tax_fun_depth', 'min_KO_num_in_one_taxon'], 100, int,
            '\b\n05_merge_tax_fun_depth: Minimum number of KOs per taxon, used to filter outputs in KO_Taxa_RelativeAbundance and KO_Taxa_SumCoverage.\n'
            'This threshold filters out taxa with insufficient KOs,\n'
            'retaining only KO-rich taxa to ensure robust metabolic reconstruction\n'
            'while preventing excessive output file sizes caused by sparse matrices.'
        ),
        'merge_tax_fun_depth_threads': (
            ['Params_c_sample', '05_merge_tax_fun_depth', 'threads'], 20, int,
            '05_merge_tax_fun_depth: Threads.'),
        'merge_cami_format_output': (
            ['Params_c_sample', '05_merge_tax_fun_depth', 'CAMI_format_output'], False, bool,
            '05_merge_tax_fun_depth: Report a taxonomic profile using the CAMI output format (flag only).'
        ),
        'merge_output_detail': (
            ['Params_c_sample', '05_merge_tax_fun_depth', 'output_detail'], False, bool,
            '05_merge_tax_fun_depth: Whether to output detail tables (flag only).'),
        'mclink_sample_min_ko_num_per_taxon': (
            ['Params_c_sample', '06_mclink_taxon', 'min_KO_num_in_one_taxon'], 100, int,
            '06_mclink_taxon: Minimum number of KOs per taxon for calculating taxon-level metabolic pathway/module completeness or abundance per sample. Taxa with fewer than this number of detected KOs in a given sample will be filtered out. Calculated independently for each sample.'
        ),
        'mclink_sample_min_ko_abundance_per_taxon': (
            ['Params_c_sample', '06_mclink_taxon', 'min_mean_KO_RA_in_one_taxon'], 0, float,
            '06_mclink_taxon: Minimum mean KO abundance [0-1e6] per taxon for calculating taxon-level metabolic pathway/module completeness or abundance per sample. Taxa with mean KO abundance below this threshold in a given sample will be filtered out. Note: GPM/TPM sum of all taxa in each sample is 1e6.'
        ),
        'mclink_sample_topn_taxon': (
            ['Params_c_sample', '06_mclink_taxon', 'TopN'], 100, int,
            '06_mclink_taxon: Number of top taxa (ranked by mean KO abundance) retained for calculating taxon-level metabolic pathway/module completeness or abundance per sample. Set to 0 to disable filtering.'
        ),
        'mclink_sample_completeness_per_taxon': (
            ['Params_c_sample', '06_mclink_taxon', 'completeness_threshold'], 0.5, float,
            '06_mclink_taxon: Completeness threshold [0-1] of modules for calculating taxon-level metabolic pathway/module abundance. Higher values retain only more complete profiles.'),
        'mclink_sample_plus_scale_method_per_taxon': (
            ['Params_c_sample', '06_mclink_taxon', 'plus_scale_method'], 'mean', str,
            "06_mclink_taxon: Aggregation method for '+'-separated KOs in a module [mean/max/min]. 'mean': average abundance across all required KOs; 'max': highest KO abundance; 'min': lowest KO abundance."),
        'mclink_sample_comma_scale_method_per_taxon': (
            ['Params_c_sample', '06_mclink_taxon', 'comma_scale_method'], 'sum', str,
            "06_mclink_taxon: Aggregation method for ','-separated alternative KOs in a module [sum/max]. 'sum': cumulative abundance of all alternative KOs; 'max': abundance of the highest single KO among alternatives. Note: this applies to abundance calculation only; module completeness is always computed with 'max'."),
        'barplot_topn_taxon': (
            ['Params_c_project', '01_Taxa_visualize', 'TopN'], 20, int,
            '01_Taxa_visualize: Number of top taxa to display in stacked barplots. Remaining categories are merged into Others. Maximum 29.'
        ),
        'barplot_topn_fun': (
            ['Params_c_project', '03_Fun_visualize', 'TopN'], 20, int,
            '03_Fun_visualize: Number of top functional categories to display in stacked barplots. Remaining categories are merged into Others. Maximum 29.'
        ),
        'mclink_min_ko_num_per_sample': (
            ['Params_c_project', '03_mclink_sample', 'min_sample_num_per_KO'], 1, int,
            '03_mclink_sample: Minimum number of samples in which a KO must be detected for calculating community-level metabolic pathway/module completeness or abundance.'
        ),
        'mclink_min_ko_prevalence_per_sample': (
            ['Params_c_project', '03_mclink_sample', 'preva'], 0, float,
            '03_mclink_sample: Minimum detection rate [0-1] of a KO across all samples for calculating community-level metabolic pathway/module completeness or abundance.'
        ),
        'mclink_min_ko_abundance_per_sample': (
            ['Params_c_project', '03_mclink_sample', 'abund'], 0, float,
            '03_mclink_sample: Minimum mean relative abundance [0-1e6] of a KO across all samples for calculating community-level metabolic pathway/module completeness or abundance. KOs with mean abundance below this threshold are excluded. Note: GPM/TPM sum of a sample is 1e6.'
        ),
        'mclink_completeness_per_sample': (
            ['Params_c_project', '03_mclink_sample', 'completeness_threshold'], 0.5, float,
            '03_mclink_sample: Completeness threshold [0-1] of modules for calculating community-level metabolic pathway/module abundance. Higher values retain only more complete profiles.'
        ),
        'mclink_plus_scale_method_per_sample': (
            ['Params_c_project', '03_mclink_sample', 'plus_scale_method'], 'mean', str,
            "03_mclink_sample: Aggregation method for '+'-separated KOs in a module [mean/max/min]. 'mean': average abundance across all required KOs; 'max': highest KO abundance; 'min': lowest KO abundance."
        ),
        'mclink_comma_scale_method_per_sample': (
            ['Params_c_project', '03_mclink_sample', 'comma_scale_method'], 'sum', str,
            "03_mclink_sample: Aggregation method for ','-separated alternative KOs in a module [sum/max]. 'sum': cumulative abundance of all alternative KOs; 'max': abundance of the highest single KO among alternatives. Note: this applies to abundance calculation only; module completeness is always computed with 'max'."
        ),
        'merge_min_taxa_ko_prevalence_across_samples': (
            ['Params_c_project', '05_merge_tax_KO_depth', 'preva'], 0.10, float,
            '05_merge_tax_KO_depth: Minimum prevalence [0-1] of a Taxon-KO across samples. Taxa-KOs detected in fewer samples than this proportion are excluded.'
        ),
        'merge_min_taxa_ko_abundance_across_samples': (
            ['Params_c_project', '05_merge_tax_KO_depth', 'abund'], 0.0, float,
            '05_merge_tax_KO_depth: Minimum mean relative abundance [0-1e6] of a Taxon-KO across samples. Taxa-KOs with mean abundance below this threshold are excluded. Note: GPM/TPM sum of a sample is 1e6.'
        ),
        'merge_min_taxa_module_prevalence_across_samples': (
            ['Params_c_project', '05_merge_tax_Module_RelativeAbundance', 'preva'], 0.10, float,
            '05_merge_tax_Module_RelativeAbundance: Minimum prevalence [0-1] of a Taxon-Module across samples. Taxa-Modules detected in fewer samples than this proportion are excluded.'
        ),
        'merge_min_taxa_module_abundance_across_samples': (
            ['Params_c_project', '05_merge_tax_Module_RelativeAbundance', 'abund'], 0.0, float,
            '05_merge_tax_Module_RelativeAbundance: Minimum mean relative abundance [0-1e6] of a Taxon-Module across samples. Taxa-Modules with mean abundance below this threshold are excluded. Note: GPM/TPM sum of a sample is 1e6.'
        ),
        'merge_min_taxa_module_sc_prevalence_across_samples': (
            ['Params_c_project', '05_merge_tax_Module_SumCoverage', 'preva'], 0.10, float,
            '05_merge_tax_Module_SumCoverage: Minimum prevalence [0-1] of a Taxon-Module across samples. Taxa-Modules detected in fewer samples than this proportion are excluded.'
        ),
        'merge_min_taxa_module_sc_abundance_across_samples': (
            ['Params_c_project', '05_merge_tax_Module_SumCoverage', 'abund'], 0.0, float,
            '05_merge_tax_Module_SumCoverage: Minimum mean sum coverage [0-1e6] of a Taxon-Module across samples. Taxa-Modules with mean sum coverage below this threshold are excluded. Note: GPM/TPM sum of a sample is 1e6.'
        ),
        'mclink_project_min_ko_num_per_taxon': (
            ['Params_c_project', '06_mclink_taxon', 'min_KO_num_in_one_taxon'], 100, int,
            '06_mclink_taxon: Minimum number of KOs per taxon for calculating taxon-level metabolic pathway/module completeness on the project-averaged taxon-function matrix. Taxa with fewer than this number of detected KOs in the project-averaged matrix will be filtered out. KO abundances are first averaged across all samples to generate a single project-level matrix before calculation.'
        ),
        'mclink_project_min_ko_abundance_per_taxon': (
            ['Params_c_project', '06_mclink_taxon', 'min_mean_KO_RA_in_one_taxon'], 0, float,
            '06_mclink_taxon: Minimum mean KO abundance [0-1e6] per taxon for calculating taxon-level metabolic pathway/module completeness on the project-averaged taxon-function matrix. Taxa with mean KO abundance below this threshold in the project-averaged matrix will be filtered out. Note: GPM/TPM sum of all taxa is 1e6 per sample; values are averaged across samples prior to filtering and calculation.'
        ),
        'mclink_project_topn_taxon': (
            ['Params_c_project', '06_mclink_taxon', 'TopN'], 100, int,
            '06_mclink_taxon: Number of top taxa (ranked by mean KO abundance in the project-averaged taxon-function matrix) retained for calculating taxon-level metabolic pathway/module completeness. Set to 0 to disable filtering.'
        ),
        'mclink_project_completeness_per_taxon': (
            ['Params_c_project', '06_mclink_taxon', 'completeness_threshold'], 0.5, float,
            '06_mclink_taxon: Completeness threshold [0-1] of modules for calculating taxon-level metabolic pathway/module abundance on the project-averaged taxon-function matrix. Higher values retain only more complete profiles.'
        ),
        'mclink_project_plus_scale_method_per_taxon': (
            ['Params_c_project', '06_mclink_taxon', 'plus_scale_method'], 'mean', str,
            "06_mclink_taxon: Aggregation method for '+'-separated KOs in a module [mean/max/min] on the project-averaged taxon-function matrix. 'mean': average abundance across all required KOs; 'max': highest KO abundance; 'min': lowest KO abundance."
        ),
        'mclink_project_comma_scale_method_per_taxon': (
            ['Params_c_project', '06_mclink_taxon', 'comma_scale_method'], 'sum', str,
            "06_mclink_taxon: Aggregation method for ','-separated alternative KOs in a module [sum/max] on the project-averaged taxon-function matrix. 'sum': cumulative abundance of all alternative KOs; 'max': abundance of the highest single KO among alternatives. Note: this applies to abundance calculation only; module completeness is always computed with 'max'."
        ),
    },
    'cr': {
        'raw_stats_sample_threads': (
            ['Params_b_sample', '00_raw_stats', 'threads'], 4, int,
            '00_raw_stats: Threads for seqkit.'),
        'trim_sample_threads': (
            ['Params_b_sample', '01_trim', 'threads'], 16, int,
            '01_trim: Threads for fastp (next-generation short reads) or fastplong (third-generation long reads).'),
        'trim_length_required': (
            ['Params_b_sample', '01_trim', 'length_required'], 15, int,
            '01_trim: Reads shorter than length_required will be discarded'),
        'trim_n_base_limit': (
            ['Params_b_sample', '01_trim', 'n_base_limit'], 5, int,
            "01_trim: if one read's number of N base is >n_base_limit, then this read/pair is discarded."),
        'trim_compression': (
            ['Params_b_sample', '01_trim', 'compression'], 4, int,
            '01_trim: Compression level for gzip output [1-9]. 1 is fastest, 9 is smallest.'),
        'trim_extra_args': (
            ['Params_b_sample', '01_trim', 'extra_args'], '', str,
            '01_trim: Extra arguments forwarded to fastp or fastplong for raw reads trimming (e.g., "--reads_to_process 0 -Q").'),
        'trim_stats_sample_threads': (
            ['Params_b_sample', '01_trim_stats', 'threads'], 16, int,
            '01_trim_stats: Threads.'),
        'assembly_method': (
            ['Params_b_sample', '02_assembly', 'method'], 'megahit', str,
            '\b\n02_assembly: Underlying assembly software using one of [megahit/spades/metaflye].\n'
            'Must be compatible with --read-type (config key read_type):\n'
            'megahit: MEGAHIT (second-generation sequencing; for paired or single DNA reads).\n'
            'spades: SPAdes (second-generation sequencing; metaSPAdes for paired DNA reads, rnaSPAdes for paired-rna or single-rna RNA reads).\n'
            'metaflye: metaFlye (third-generation sequencing; required for pacbio-raw/corr/hifi or nano-raw/corr/hq long reads).'),
        'assembly_threads': (
            ['Params_b_sample', '02_assembly', 'threads'], 40, int,
            '02_assembly: Threads for [megahit/spades/metaflye].'),
        'megahit_presets': (
            ['Params_b_sample', '02_megahit', 'presets'], 'meta-sensitive', str,
            '\b\n02_megahit: MEGAHIT presets using one of [meta-sensitive/meta-large].\n'
            "meta-sensitive: '--min-count 1 --k-list 21,29,39,49,...,129,141'.\n"
            "meta-large: '--k-min 27 --k-max 127 --k-step 10'.\n"
            '(large & complex metagenomes, like soil).'),
        'megahit_k_list': (
            ['Params_b_sample', '02_megahit', 'k_list'], '', str,
            '02_megahit: MEGAHIT --k-list (optional). Explicit comma-separated k-mer sizes (e.g. "21,29,39,59,79,99,119,141"). Overrides --megahit-presets. If set, --megahit-k-min, --megahit-k-max, --megahit-k-step must NOT be set. Optional --megahit-min-count can be used together.'),
        'megahit_min_count': (
            ['Params_b_sample', '02_megahit', 'min_count'], 2, int,
            '02_megahit: MEGAHIT --min-count (optional). Minimum multiplicity for filtering (k_min+1)-mers. Can be used with --megahit-k-list or --megahit-k-min/--k-max/--k-step.'),
        'megahit_k_min': (
            ['Params_b_sample', '02_megahit', 'k_min'], '', str,
            '02_megahit: MEGAHIT --k-min (optional). Minimum kmer size (<= 255), must be odd number. If set with --k-max and --k-step, overrides --megahit-presets and uses custom k-mer mode.'),
        'megahit_k_max': (
            ['Params_b_sample', '02_megahit', 'k_max'], '', str,
            '02_megahit: MEGAHIT --k-max (optional). Maximum kmer size (<= 255), must be odd number.'),
        'megahit_k_step': (
            ['Params_b_sample', '02_megahit', 'k_step'], '', str,
            '02_megahit: MEGAHIT --k-step (optional). Increment of kmer size of each iteration (<= 28), must be even number.'),
        'megahit_extra_args': (
            ['Params_b_sample', '02_megahit', 'extra_args'], '', str,
            '02_megahit: Extra arguments forwarded to MEGAHIT (e.g. "--no-mercy --tmp-dir /path/to/your/tmp/dir"). Do NOT use this for k-mer settings (--k-list, --k-min, --presets, etc.).'),
        'spades_extra_args': (
            ['Params_b_sample', '02_spades', 'extra_args'], '', str,
            '02_spades: Extra arguments forwarded to SPAdes. Do NOT use this for input/output or thread settings.'),
        'metaflye_extra_args': (
            ['Params_b_sample', '02_metaflye', 'extra_args'], '', str,
            '02_metaflye: Extra arguments forwarded to Flye for third-generation long reads assembly (e.g., "--iterations 2"). Do NOT use this for read type or input/output settings.'),
        'mapping_threads': (
            ['Params_b_sample', '03_mapping', 'threads'], 20, int,
            '03_mapping: Threads for [minibwa/bwamem/bwamem2/strobealign].'),
        'mapping_method': (
            ['Params_b_sample', '03_mapping', 'method'], 'minibwa', str,
            '\b\n03_mapping: Underlying mapping software used [minibwa/bwamem/bwamem2/strobealign].\n'
            'minibwa: minibwa using default parameters.\n'
            'bwamem: bwa mem using default parameters.\n'
            'bwamem2: bwa-mem2 using default parameters.\n'
            'strobealign: strobealign using default parameters.'),
        'mapping_samtools_sort_mem': (
            ['Params_b_sample', '03_mapping', 'samtools_sort_mem'], '768M', str,
            '03_mapping: Set maximum memory per thread; suffix K/M/G recognized, e.g. 2G, 768M.'),
        'mapping_extra_args': (
            ['Params_b_sample', '03_mapping', 'extra_args'], '', str,
            '03_mapping: Extra arguments forwarded to the active mapping aligner (minibwa map, bwa mem, bwa-mem2 mem, or strobealign --use-index), (e.g. "-B 8" for minibwa map).'),
        'contigs_depth_sample_threads': (
            ['Params_b_sample', '07_contigs_depth', 'threads'], 8, int,
            '07_contigs_depth: Threads.'),
        'depth_min_read_percent_identity': (
            ['Params_b_sample', '07_contigs_depth', 'min_read_percent_identity'], 95, float,
            '\b\n07_contigs_depth: Exclude reads by overall percent identity [0-100] for contig depth estimation\n'
            'using jgi_summarize_bam_contig_depths --percentIdentity, e.g. 95 for 95%.\n'
            'Reads with percent identity below this threshold are excluded from depth calculation.\n'
            'Higher values increase stringency, retaining only high-quality alignments.'),
        'depth_contig_end_exclusion': (
            ['Params_b_sample', '07_contigs_depth', 'contig_end_exclusion'], 75, int,
            '\b\n07_contigs_depth: Exclude bases at the ends of contig sequences for depth estimation using jgi_summarize_bam_contig_depths --maxEdgeBases.\n'
            'Bases within this distance from either end of a contig are ignored when calculating depth,\n'
            'preventing edge effects from partial alignments.'),
        'depth_min_contig_length': (
            ['Params_b_sample', '07_contigs_depth', 'min_contig_length'], 300, int,
            '\b\n07_contigs_depth: Minimum contig length for depth estimation [bp], e.g. 300 for 300 bp.\n'
            'Contigs shorter than this threshold are excluded from depth calculation.'),
        'contigs_prodigal_extra_args': (
            ['Params_b_sample', '07_contigs_depth', 'extra_args'], '', str,
            '07_contigs_depth: Extra arguments forwarded to prodigal for contig gene prediction (e.g., "-m -n").'),
        'raw_stats_project_threads': (
            ['Params_b_project', '01_cat_raw_stats', 'threads'], 4, int,
            '01_cat_raw_stats: Threads.'),
        'trimmed_stats_project_threads': (
            ['Params_b_project', '02_cat_trimmed_stats', 'threads'], 4, int,
            '02_cat_trimmed_stats: Threads.'),
        'contigs_stats_project_threads': (
            ['Params_b_project', '03_contigs_stats', 'threads'], 4, int,
            '03_contigs_stats: Threads.'),
        'merge_stats_project_threads': (
            ['Params_b_project', '04_merge_reads_contigs_stats', 'threads'], 4, int,
            '04_merge_reads_contigs_stats: Threads.'),
        'contigs_taxonomy_classify_method': (
            ['Params_c_sample', '01_taxonomy', 'contigs_taxonomy_classify_method'], 'kraken2', str,
            '01_taxonomy: Contigs taxonomy classifier using one of [CAT_pack/kraken2].'
        ),
        'contigs_taxonomy_threads': (
            ['Params_c_sample', '01_taxonomy', 'threads'], 60, int,
            '01_taxonomy: Threads for [CAT_pack/kraken2].'
        ),
        'contigs_taxonomy_cat_range': (
            ['Params_c_sample', '01_taxonomy_CAT_pack', 'range'], 10, float,
            '\b\n01_taxonomy_CAT_pack: Percentage range [0-100] of DIAMOND hits relative to top bit-score to include in LCA calculation.\n'
            'Higher values yield more conservative classifications by incorporating divergent hits, while lower values enable specific assignments.\n'
            'Applicable when --contigs-taxonomy-classify-method is set to CAT_pack.'
        ),
        'contigs_taxonomy_cat_fraction': (
            ['Params_c_sample', '01_taxonomy_CAT_pack', 'fraction'], 0.50, float,
            '\b\n01_taxonomy_CAT_pack: Minimum bit-score support fraction [0-0.99].\n'
            'Higher values enforce stricter consensus (one classification per sequence), while lower values allow multiple lineage detections.\n'
            'Applicable when --contigs-taxonomy-classify-method is set to CAT_pack.'
        ),
        'contigs_taxonomy_cat_pack_extra_args': (
            ['Params_c_sample', '01_taxonomy_CAT_pack', 'extra_args'], '', str,
            '01_taxonomy_CAT_pack: Extra arguments forwarded to CAT_pack contigs (e.g., "--block_size 12.0").'),
        'contigs_taxonomy_kraken2_confidence': (
            ['Params_c_sample', '01_taxonomy_kraken2', 'confidence'], 0.7, float,
            '\b\n01_taxonomy_kraken2: Confidence threshold [0-1] for taxonomic classification, applicable when --contigs-taxonomy-classify-method is set to kraken2.\n'
            'Default value was tuned on the CAMI2 marine dataset. For samples with higher complexity, a lower confidence threshold is usually recommended.'
        ),
        'contigs_taxonomy_kraken2_extra_args': (
            ['Params_c_sample', '01_taxonomy_kraken2', 'extra_args'], '', str,
            '01_taxonomy_kraken2: Extra arguments forwarded to kraken2 (e.g., "--memory-mapping").'),
        'contigs_function_classify_method': (
            ['Params_c_sample', '02_function', 'contigs_function_classify_method'], 'kraken2', str,
            '02_function: Contigs function classifier using one of [kraken2/emapper].'
        ),
        'contigs_function_threads': (
            ['Params_c_sample', '02_function', 'threads'], 60, int,
            '02_function: Threads for [kraken2/emapper].'),
        'contigs_function_emapper_evalue': (
            ['Params_c_sample', '02_function_emapper', 'seed_ortholog_evalue'], 1e-10, float,
            '02_function_emapper: Report only alignments with e-value below or equal to the threshold, applicable when --contigs-function-classify-method is set to emapper.'
        ),
        'contigs_function_emapper_method': (
            ['Params_c_sample', '02_function_emapper', 'method'], 'diamond', str,
            '\b\n02_function_emapper: search seed orthologs using one of [diamond/mmseqs], applicable when --contigs-function-classify-method is set to emapper.\n'
            'diamond: using default parameters of diamond mode of emapper.\n'
            'mmseqs: using default parameters of mmseqs mode of emapper'
        ),
        'contigs_function_emapper_query_cover': (
            ['Params_c_sample', '02_function_emapper', 'query_cover'], 30, float,
            '\b\n02_function_emapper: Report only alignments with query cover above or equal to the given percentage [0-100],.\n'
            'applicable when --contigs-function-classify-method is set to emapper'
        ),
        'contigs_function_emapper_extra_args': (
            ['Params_c_sample', '02_function_emapper', 'extra_args'], '', str,
            '02_function_emapper: Extra arguments forwarded to eggNOG-mapper (e.g., "--pident 30").'),
        'contigs_function_kraken2_confidence': (
            ['Params_c_sample', '02_function_kraken2', 'confidence'], 0.05, float,
            '02_function_kraken2: Confidence threshold [0-1] for functional classification, applicable when --contigs-function-classify-method is set to kraken2.'
        ),
        'contigs_function_kraken2_extra_args': (
            ['Params_c_sample', '02_function_kraken2', 'extra_args'], '', str,
            '02_function_kraken2: Extra arguments forwarded to kraken2 (e.g., "--memory-mapping").'),
        'cr_unmapped_reads_merge_threads': (
            ['Params_cr_sample', '01_unmapped_reads_trim_merge', 'threads'], 16, int,
            '01_unmapped_reads_trim_merge_fastp: Threads for fastp, samtools, and seqkit (paired-end), or samtools and seqkit only (single-end).'),
        'cr_unmapped_reads_merge_length_required': (
            ['Params_cr_sample', '01_unmapped_reads_trim_merge', 'length_required'], 15, int,
            '01_unmapped_reads_trim_merge_fastp: Reads shorter than length_required will be discarded. Only applicable when read_type is paired.'),
        'cr_unmapped_reads_merge_n_base_limit': (
            ['Params_cr_sample', '01_unmapped_reads_trim_merge', 'n_base_limit'], 5, int,
            "01_unmapped_reads_trim_merge_fastp: if one read's number of N base is >n_base_limit, then this read/pair is discarded. Only applicable when read_type is paired."),
        'cr_unmapped_reads_merge_compression': (
            ['Params_cr_sample', '01_unmapped_reads_trim_merge', 'compression'], 4, int,
            '01_unmapped_reads_trim_merge_fastp: Compression level for gzip output [1-9]. 1 is fastest, 9 is smallest. Only applicable when read_type is paired.'),
        'cr_unmapped_reads_fastp_extra_args': (
            ['Params_cr_sample', '01_unmapped_reads_trim_merge', 'extra_args'], '', str,
            '01_unmapped_reads_trim_merge: Extra arguments forwarded to fastp for unmapped reads trimming (e.g., "--reads_to_process 0 -Q"). Only applicable when read_type is paired.'),
        'cr_fraggenescanrs_model': (
            ['Params_cr_sample', '02_FragGeneScanRs', 'model'], 'illumina_1', str,
            '\b\n02_FragGeneScanRs: Model parameter used for short-read open reading frame prediction.\n'
            '[complete] for complete genomic sequences or short sequence reads without sequencing error\n'
            '[sanger_5] for Sanger sequencing reads with about 0.5% error rate\n'
            '[sanger_10] for Sanger sequencing reads with about 1% error rate\n'
            '[454_5] for 454 pyrosequencing reads with about 0.5% error rate\n'
            '[454_10] for 454 pyrosequencing reads with about 1% error rate\n'
            '[454_30] for 454 pyrosequencing reads with about 3% error rate\n'
            '[illumina_1] for Illumina sequencing reads with about 0.1% error rate\n'
            '[illumina_5] for Illumina sequencing reads with about 0.5% error rate\n'
            '[illumina_10] for Illumina sequencing reads with about 1% error rate'
        ),
        'cr_fraggenescanrs_threads': (
            ['Params_cr_sample', '02_FragGeneScanRs', 'threads'], 20, int,
            '02_FragGeneScanRs: Threads.'),
        'cr_fraggenescanrs_stats_threads': (
            ['Params_cr_sample', '02_FragGeneScanRs_stats', 'threads'], 10, int,
            '02_FragGeneScanRs_stats: Threads.'),
        'cr_taxonomy_threads': (
            ['Params_cr_sample', '03_taxonomy', 'threads'], 60, int,
            '03_taxonomy: Threads.'),
        'cr_taxonomy_kraken2_confidence': (
            ['Params_cr_sample', '03_taxonomy_kraken2', 'confidence'], 0.7, float,
            '\b\n03_taxonomy_kraken2: Confidence threshold [0-1] for taxonomic classification of unmapped reads with kraken2.\n'
            'Default value was tuned on the CAMI2 marine dataset. For samples with higher complexity, a lower confidence threshold is usually recommended.'
        ),
        'cr_taxonomy_kraken2_extra_args': (
            ['Params_cr_sample', '03_taxonomy_kraken2', 'extra_args'], '', str,
            '03_taxonomy_kraken2: Extra arguments forwarded to kraken2 (e.g., "--memory-mapping").'),
        'cr_function_classify_method': (
            ['Params_cr_sample', '04_function', 'reads_function_classify_method'], 'kraken2', str,
            '04_function: reads functional classifier using one of [kraken2/emapper].'
        ),
        'cr_function_threads': (
            ['Params_cr_sample', '04_function', 'threads'], 60, int,
            '04_function: Threads for [kraken2/emapper].'),
        'cr_function_kraken2_confidence': (
            ['Params_cr_sample', '04_function_kraken2', 'confidence'], 0.05, float,
            '04_function_kraken2: Confidence threshold [0-1] for functional classification for reads when --cr-function-classify-method is set to kraken2.'
        ),
        'cr_function_kraken2_extra_args': (
            ['Params_cr_sample', '04_function_kraken2', 'extra_args'], '', str,
            '04_function_kraken2: Extra arguments forwarded to kraken2 (e.g., "--memory-mapping").'),
        'cr_function_emapper_evalue': (
            ['Params_cr_sample', '04_function_emapper', 'seed_ortholog_evalue'], 1e-10, float,
            '04_function_emapper: Report only alignments with e-value below or equal to the threshold, applicable when --cr-function-classify-method is set to emapper.'
        ),
        'cr_function_emapper_method': (
            ['Params_cr_sample', '04_function_emapper', 'method'], 'diamond', str,
            '\b\n04_function_emapper: search seed orthologs using one of [diamond/mmseqs], applicable when --cr-function-classify-method is set to emapper.\n'
            'diamond: using default parameters of diamond mode of emapper.\n'
            'mmseqs: using default parameters of mmseqs mode of emapper'
        ),
        'cr_function_emapper_query_cover': (
            ['Params_cr_sample', '04_function_emapper', 'query_cover'], 30, float,
            '\b\n04_function_emapper: Report only alignments with query cover above or equal to the given percentage [0-100],.\n'
            'applicable when --cr-function-classify-method is set to emapper'
        ),
        'cr_function_emapper_extra_args': (
            ['Params_cr_sample', '04_function_emapper', 'extra_args'], '', str,
            '04_function_emapper: Extra arguments forwarded to eggNOG-mapper (e.g., "--pident 30").'),
        'merge_min_orf_length': (
            ['Params_cr_sample', '05_merge_tax_fun_depth', 'minORFLen'], 50, int,
            '05_merge_tax_fun_depth: Minimum ORF length for finally used ORFs of contigs and reads.'
        ),
        'merge_multiple_ko_handling_mode': (
            ['Params_cr_sample', '05_merge_tax_fun_depth', 'KO_mode'], 'ignore_multiple_KOs', str,
            '\b\n05_merge_tax_fun_depth: Method using one of [split_multiple_KOs/keep_multiple_KOs/ignore_multiple_KOs] when an ORF has multiple KO annotations:.\n'
            'split_multiple_KOs: evenly distributing the total abundance across all assigned KOs;\n'
            'keep_multiple_KOs: retain original multiple KO annotations;\n'
            'ignore_multiple_KOs: ignore the ORFs with multiple KO annotations and only use the ORFs with single KO annotation.'
        ),
        'merge_min_ko_num_per_taxon': (
            ['Params_cr_sample', '05_merge_tax_fun_depth', 'min_KO_num_in_one_taxon'], 100, int,
            '\b\n05_merge_tax_fun_depth: Minimum number of KOs per taxon, used to filter outputs in KO_Taxa_RelativeAbundance and KO_Taxa_SumCoverage.\n'
            'This threshold filters out taxa with insufficient KOs,\n'
            'retaining only KO-rich taxa to ensure robust metabolic reconstruction\n'
            'while preventing excessive output file sizes caused by sparse matrices.'
        ),
        'merge_tax_fun_depth_threads': (
            ['Params_cr_sample', '05_merge_tax_fun_depth', 'threads'], 20, int,
            '05_merge_tax_fun_depth: Threads.'),
        'merge_cami_format_output': (
            ['Params_cr_sample', '05_merge_tax_fun_depth', 'CAMI_format_output'], False, bool,
            '05_merge_tax_fun_depth: Report a taxonomic profile using the CAMI output format (flag only).'
        ),
        'merge_output_detail': (
            ['Params_cr_sample', '05_merge_tax_fun_depth', 'output_detail'], False, bool,
            '05_merge_tax_fun_depth: Whether to output detail tables (flag only).'),
        'mclink_sample_min_ko_num_per_taxon': (
            ['Params_cr_sample', '06_mclink_taxon', 'min_KO_num_in_one_taxon'], 100, int,
            '06_mclink_taxon: Minimum number of KOs per taxon for calculating taxon-level metabolic pathway/module completeness or abundance per sample. Taxa with fewer than this number of detected KOs in a given sample will be filtered out. Calculated independently for each sample.'
            ),
        'mclink_sample_min_ko_abundance_per_taxon': (
            ['Params_cr_sample', '06_mclink_taxon', 'min_mean_KO_RA_in_one_taxon'], 0, float,
            '06_mclink_taxon: Minimum mean KO abundance [0-1e6] per taxon for calculating taxon-level metabolic pathway/module completeness or abundance per sample. Taxa with mean KO abundance below this threshold in a given sample will be filtered out. Note: GPM/TPM sum of all taxa in each sample is 1e6.'
            ),
        'mclink_sample_topn_taxon': (
            ['Params_cr_sample', '06_mclink_taxon', 'TopN'], 100, int,
            '06_mclink_taxon: Number of top taxa (ranked by mean KO abundance) retained for calculating taxon-level metabolic pathway/module completeness or abundance per sample. Set to 0 to disable filtering.'
        ),
        'mclink_sample_completeness_per_taxon': (
            ['Params_cr_sample', '06_mclink_taxon', 'completeness_threshold'], 0.5, float,
            '06_mclink_taxon: Completeness threshold [0-1] of modules for calculating taxon-level metabolic pathway/module abundance. Higher values retain only more complete profiles.'),
        'mclink_sample_plus_scale_method_per_taxon': (
            ['Params_cr_sample', '06_mclink_taxon', 'plus_scale_method'], 'mean', str,
             "06_mclink_taxon: Aggregation method for '+'-separated KOs in a module [mean/max/min]. 'mean': average abundance across all required KOs; 'max': highest KO abundance; 'min': lowest KO abundance."),
        'mclink_sample_comma_scale_method_per_taxon': (
            ['Params_cr_sample', '06_mclink_taxon', 'comma_scale_method'], 'sum', str,
            "06_mclink_taxon: Aggregation method for ','-separated alternative KOs in a module [sum/max]. 'sum': cumulative abundance of all alternative KOs; 'max': abundance of the highest single KO among alternatives. Note: this applies to abundance calculation only; module completeness is always computed with 'max'."),
        'barplot_topn_taxon': (
            ['Params_cr_project', '01_Taxa_visualize', 'TopN'], 10, int,
            '01_Taxa_visualize: Number of top taxa to display in stacked barplots. Remaining categories are merged into Others. Maximum 29.'
        ),
        'barplot_topn_fun': (
            ['Params_cr_project', '03_Fun_visualize', 'TopN'], 10, int,
            '03_Fun_visualize: Number of top functional categories to display in stacked barplots. Remaining categories are merged into Others. Maximum 29.'
        ),
        'mclink_min_ko_num_per_sample': (
            ['Params_cr_project', '03_mclink_sample', 'min_sample_num_per_KO'], 1, int,
            '03_mclink_sample: Minimum number of samples in which a KO must be detected for calculating community-level metabolic pathway/module completeness or abundance.'
        ),
        'mclink_min_ko_prevalence_per_sample': (
            ['Params_cr_project', '03_mclink_sample', 'preva'], 0, float,
            '03_mclink_sample: Minimum detection rate [0-1] of a KO across all samples for calculating community-level metabolic pathway/module completeness or abundance.'
        ),
        'mclink_min_ko_abundance_per_sample': (
            ['Params_cr_project', '03_mclink_sample', 'abund'], 0, float,
            '03_mclink_sample: Minimum mean relative abundance [0-1e6] of a KO across all samples for calculating community-level metabolic pathway/module completeness or abundance. KOs with mean abundance below this threshold are excluded. Note: GPM/TPM sum of a sample is 1e6.'
        ),
        'mclink_completeness_per_sample': (
            ['Params_cr_project', '03_mclink_sample', 'completeness_threshold'], 0.5, float,
            '03_mclink_sample: Completeness threshold [0-1] of modules for calculating community-level metabolic pathway/module abundance. Higher values retain only more complete profiles.'
        ),
        'mclink_plus_scale_method_per_sample': (
            ['Params_cr_project', '03_mclink_sample', 'plus_scale_method'], 'mean', str,
            "03_mclink_sample: Aggregation method for '+'-separated KOs in a module [mean/max/min]. 'mean': average abundance across all required KOs; 'max': highest KO abundance; 'min': lowest KO abundance."
        ),
        'mclink_comma_scale_method_per_sample': (
            ['Params_cr_project', '03_mclink_sample', 'comma_scale_method'], 'sum', str,
            "03_mclink_sample: Aggregation method for ','-separated alternative KOs in a module [sum/max]. 'sum': cumulative abundance of all alternative KOs; 'max': abundance of the highest single KO among alternatives. Note: this applies to abundance calculation only; module completeness is always computed with 'max'."
        ),
        'merge_min_taxa_ko_prevalence_across_samples': (
            ['Params_cr_project', '05_merge_tax_KO_depth', 'preva'], 0.10, float,
            '05_merge_tax_KO_depth: Minimum prevalence [0-1] of a Taxon-KO across samples. Taxa-KOs detected in fewer samples than this proportion are excluded.'
        ),
        'merge_min_taxa_ko_abundance_across_samples': (
            ['Params_cr_project', '05_merge_tax_KO_depth', 'abund'], 0.0, float,
            '05_merge_tax_KO_depth: Minimum mean relative abundance [0-1e6] of a Taxon-KO across samples. Taxa-KOs with mean abundance below this threshold are excluded. Note: GPM/TPM sum of a sample is 1e6.'
        ),
        'merge_min_taxa_module_prevalence_across_samples': (
            ['Params_cr_project', '05_merge_tax_Module_RelativeAbundance', 'preva'], 0.10, float,
            '05_merge_tax_Module_RelativeAbundance: Minimum prevalence [0-1] of a Taxon-Module across samples. Taxa-Modules detected in fewer samples than this proportion are excluded.'
        ),
        'merge_min_taxa_module_abundance_across_samples': (
            ['Params_cr_project', '05_merge_tax_Module_RelativeAbundance', 'abund'], 0.0, float,
            '05_merge_tax_Module_RelativeAbundance: Minimum mean relative abundance [0-1e6] of a Taxon-Module across samples. Taxa-Modules with mean abundance below this threshold are excluded. Note: GPM/TPM sum of a sample is 1e6.'
        ),
        'merge_min_taxa_module_sc_prevalence_across_samples': (
            ['Params_cr_project', '05_merge_tax_Module_SumCoverage', 'preva'], 0.10, float,
            '05_merge_tax_Module_SumCoverage: Minimum prevalence [0-1] of a Taxon-Module across samples. Taxa-Modules detected in fewer samples than this proportion are excluded.'
        ),
        'merge_min_taxa_module_sc_abundance_across_samples': (
            ['Params_cr_project', '05_merge_tax_Module_SumCoverage', 'abund'], 0.0, float,
            '05_merge_tax_Module_SumCoverage: Minimum mean sum coverage [0-1e6] of a Taxon-Module across samples. Taxa-Modules with mean sum coverage below this threshold are excluded. Note: GPM/TPM sum of a sample is 1e6.'
        ),
        'mclink_project_min_ko_num_per_taxon': (
            ['Params_cr_project', '06_mclink_taxon', 'min_KO_num_in_one_taxon'], 100, int,
            '06_mclink_taxon: Minimum number of KOs per taxon for calculating taxon-level metabolic pathway/module completeness on the project-averaged taxon-function matrix. Taxa with fewer than this number of detected KOs in the project-averaged matrix will be filtered out. KO abundances are first averaged across all samples to generate a single project-level matrix before calculation.'
        ),
        'mclink_project_min_ko_abundance_per_taxon': (
            ['Params_cr_project', '06_mclink_taxon', 'min_mean_KO_RA_in_one_taxon'], 0, float,
            '06_mclink_taxon: Minimum mean KO abundance [0-1e6] per taxon for calculating taxon-level metabolic pathway/module completeness on the project-averaged taxon-function matrix. Taxa with mean KO abundance below this threshold in the project-averaged matrix will be filtered out. Note: GPM/TPM sum of all taxa is 1e6 per sample; values are averaged across samples prior to filtering and calculation.'
        ),
        'mclink_project_topn_taxon': (
            ['Params_cr_project', '06_mclink_taxon', 'TopN'], 100, int,
            '06_mclink_taxon: Number of top taxa (ranked by mean KO abundance in the project-averaged taxon-function matrix) retained for calculating taxon-level metabolic pathway/module completeness. Set to 0 to disable filtering.'
        ),
        'mclink_project_completeness_per_taxon': (
            ['Params_cr_project', '06_mclink_taxon', 'completeness_threshold'], 0.5, float,
            '06_mclink_taxon: Completeness threshold [0-1] of modules for calculating taxon-level metabolic pathway/module abundance on the project-averaged taxon-function matrix. Higher values retain only more complete profiles.'),
        'mclink_project_plus_scale_method_per_taxon': (
            ['Params_cr_project', '06_mclink_taxon', 'plus_scale_method'], 'mean', str,
            "06_mclink_taxon: Aggregation method for '+'-separated KOs in a module [mean/max/min] on the project-averaged taxon-function matrix. 'mean': average abundance across all required KOs; 'max': highest KO abundance; 'min': lowest KO abundance."),
        'mclink_project_comma_scale_method_per_taxon': (
            ['Params_cr_project', '06_mclink_taxon', 'comma_scale_method'], 'sum', str,
            "06_mclink_taxon: Aggregation method for ','-separated alternative KOs in a module [sum/max] on the project-averaged taxon-function matrix. 'sum': cumulative abundance of all alternative KOs; 'max': abundance of the highest single KO among alternatives. Note: this applies to abundance calculation only; module completeness is always computed with 'max'."),
    },
    'r': {
        'raw_stats_sample_threads': (
            ['Params_b_sample', '00_raw_stats', 'threads'], 4, int,
            '00_raw_stats: Threads for seqkit.'),
        'trim_sample_threads': (
            ['Params_b_sample', '01_trim', 'threads'], 16, int,
            '01_trim: Threads for fastp (next-generation short reads) or fastplong (third-generation long reads).'),
        'trim_length_required': (
            ['Params_b_sample', '01_trim', 'length_required'], 15, int,
            '01_trim: Reads shorter than length_required will be discarded'),
        'trim_n_base_limit': (
            ['Params_b_sample', '01_trim', 'n_base_limit'], 5, int,
            "01_trim: if one read's number of N base is >n_base_limit, then this read/pair is discarded."),
        'trim_compression': (
            ['Params_b_sample', '01_trim', 'compression'], 4, int,
            '01_trim: Compression level for gzip output [1-9]. 1 is fastest, 9 is smallest.'),
        'trim_extra_args': (
            ['Params_b_sample', '01_trim', 'extra_args'], '', str,
            '01_trim: Extra arguments forwarded to fastp or fastplong for raw reads trimming (e.g., "--reads_to_process 0 -Q").'),
        'trim_stats_sample_threads': (
            ['Params_b_sample', '01_trim_stats', 'threads'], 16, int,
            '01_trim_stats: Threads.'),
        'raw_stats_project_threads': (
            ['Params_b_project', '01_cat_raw_stats', 'threads'], 4, int,
            '01_cat_raw_stats: Threads.'),
        'trimmed_stats_project_threads': (
            ['Params_b_project', '02_cat_trimmed_stats', 'threads'], 4, int,
            '02_cat_trimmed_stats: Threads.'),
        'target_stage': (
            ['Params_r_sample', '00_target_stage'], 'all', str,
            '\b\nTarget execution stage for reads-level analysis, using one of [all/taxonomy/function/taxonomy_and_function].\n'
            'taxonomy: only run taxonomy classification (kraken2 and bracken).\n'
            'function: only run function classification (FragGeneScanRs + kraken2/emapper).\n'
            'taxonomy_and_function: run taxonomy and function classification (kraken2/bracken + FragGeneScanRs + kraken2/emapper) (skip merge step).\n'
            'all: run taxonomy and function classification, then merge results.'
        ),
        'reads_sample_num': (
            ['Params_r_sample', '00_sample_reads', 'target_reads'], 0, int,
            '\b\n00_sample_reads: Target number of reads for downsampling. Set to 0 to disable sampling.\n'
            'For example, 10000000 reads typically yield ~1 GB compressed FASTQ files.\n'
            'If the actual read count exceeds this value, seqkit sample will randomly extract the target number;\n'
            'otherwise, the original files are used directly via symbolic links.'
        ),
        'reads_sample_seed': (
            ['Params_r_sample', '00_sample_reads', 'sample_seed'], 6, int,
            '00_sample_reads: Random seed for seqkit sample, ensuring reproducible read downsampling.'
        ),
        'reads_trim_merge_threads': (
            ['Params_r_sample', '01_sampled_reads_trim_merge', 'threads'], 16, int,
            '01_sampled_reads_trim_merge_fastp: Threads. Only applicable when read_type is paired.'),
        'reads_trim_merge_length_required': (
            ['Params_r_sample', '01_sampled_reads_trim_merge', 'length_required'], 15, int,
            '01_sampled_reads_trim_merge_fastp: Reads shorter than length_required will be discarded. Only applicable when read_type is paired.'),
        'reads_trim_merge_n_base_limit': (
            ['Params_r_sample', '01_sampled_reads_trim_merge', 'n_base_limit'], 5, int,
            "01_sampled_reads_trim_merge_fastp: if one read's number of N base is >n_base_limit, then this read/pair is discarded. Only applicable when read_type is paired."),
        'reads_trim_merge_compression': (
            ['Params_r_sample', '01_sampled_reads_trim_merge', 'compression'], 4, int,
            '01_sampled_reads_trim_merge_fastp: Compression level for gzip output [1-9]. 1 is fastest, 9 is smallest. Only applicable when read_type is paired.'),
        'reads_unmapped_reads_fastp_extra_args': (
            ['Params_r_sample', '01_sampled_reads_trim_merge', 'extra_args'], '', str,
            '01_sampled_reads_trim_merge: Extra arguments forwarded to fastp for sampled reads trimming (e.g., "--reads_to_process 0 -Q"). Only applicable when read_type is paired.'),
        'reads_fraggenescanrs_model': (
            ['Params_r_sample', '02_FragGeneScanRs', 'model'], 'illumina_1', str,
            '\b\n02_FragGeneScanRs: Model parameter used for short-read open reading frame prediction.\n'
            '[complete] for complete genomic sequences or short sequence reads without sequencing error\n'
            '[sanger_5] for Sanger sequencing reads with about 0.5% error rate\n'
            '[sanger_10] for Sanger sequencing reads with about 1% error rate\n'
            '[454_5] for 454 pyrosequencing reads with about 0.5% error rate\n'
            '[454_10] for 454 pyrosequencing reads with about 1% error rate\n'
            '[454_30] for 454 pyrosequencing reads with about 3% error rate\n'
            '[illumina_1] for Illumina sequencing reads with about 0.1% error rate\n'
            '[illumina_5] for Illumina sequencing reads with about 0.5% error rate\n'
            '[illumina_10] for Illumina sequencing reads with about 1% error rate'
        ),
        'reads_fraggenescanrs_threads': (
            ['Params_r_sample', '02_FragGeneScanRs', 'threads'], 20, int,
            '02_FragGeneScanRs: Threads.'),
        'reads_fraggenescanrs_stats_threads': (
            ['Params_r_sample', '02_FragGeneScanRs_stats', 'threads'], 10, int,
            '02_FragGeneScanRs_stats: Threads.'),
        'reads_taxonomy_threads': (
            ['Params_r_sample', '03_taxonomy', 'threads'], 60, int,
            '03_taxonomy: Threads for [kraken2].'),
        'reads_taxonomy_kraken2_confidence': (
            ['Params_r_sample', '03_taxonomy_kraken2', 'confidence'], 0.6, float,
            '\b\n03_taxonomy_kraken2: Confidence threshold [0-1] for taxonomic classification for reads using kraken2 as taxonomy classifier.\n'
            'Default value was tuned on the CAMI2 marine dataset. For samples with higher complexity, a lower confidence threshold is usually recommended.'
        ),
        'reads_taxonomy_kraken2_extra_args': (
            ['Params_r_sample', '03_taxonomy_kraken2', 'extra_args'], '', str,
            '03_taxonomy_kraken2: Extra arguments forwarded to kraken2 (e.g., "--memory-mapping").'),
        'reads_taxonomy_bracken_level': (
            ['Params_r_sample', '03_taxonomy', 'bracken_level'], 'S', str,
            '\b\n03_taxonomy: Taxonomic level at which bracken re-estimates abundances, using one of [D/P/C/O/F/G/S].'
        ),
        'reads_taxonomy_bracken_threshold': (
            ['Params_r_sample', '03_taxonomy', 'bracken_threshold'], 0, int,
            '\b\n03_taxonomy: Minimum number of clade reads a taxon must have at the bracken level\n'
            'to be included in abundance re-estimation (bracken -t; 0 = keep everything).'
        ),
        'reads_function_classify_method': (
            ['Params_r_sample', '04_function', 'reads_function_classify_method'], 'kraken2', str,
            '04_function: reads functional classifier using one of [kraken2/emapper].'
        ),
        'reads_function_threads': (
            ['Params_r_sample', '04_function', 'threads'], 60, int,
            '04_function: Threads for [kraken2/emapper].'),
        'reads_function_kraken2_confidence': (
            ['Params_r_sample', '04_function_kraken2', 'confidence'], 0.00, float,
            '04_function_kraken2: Confidence threshold [0-1] for functional classification for reads when --reads-function-classify-method is set to kraken2.'
        ),
        'reads_function_kraken2_extra_args': (
            ['Params_r_sample', '04_function_kraken2', 'extra_args'], '', str,
            '04_function_kraken2: Extra arguments forwarded to kraken2 (e.g., "--memory-mapping").'),
        'reads_function_emapper_evalue': (
            ['Params_r_sample', '04_function_emapper', 'seed_ortholog_evalue'], 1e-10, float,
            '04_function_emapper: Report only alignments with e-value below or equal to the threshold, applicable when --reads-function-classify-method is set to emapper.'
        ),
        'reads_function_emapper_method': (
            ['Params_r_sample', '04_function_emapper', 'method'], 'diamond', str,
            '\b\n04_function_emapper: search seed orthologs using one of [diamond/mmseqs], applicable when --reads-function-classify-method is set to emapper.\n'
            'diamond: using default parameters of diamond mode of emapper.\n'
            'mmseqs: using default parameters of mmseqs mode of emapper'
        ),
        'reads_function_emapper_query_cover': (
            ['Params_r_sample', '04_function_emapper', 'query_cover'], 30, float,
            '\b\n04_function_emapper: Report only alignments with query cover above or equal to the given percentage [0-100],.\n'
            'applicable when --read-function-classify-method is set to emapper'
        ),
        'reads_function_emapper_extra_args': (
            ['Params_r_sample', '04_function_emapper', 'extra_args'], '', str,
            '04_function_emapper: Extra arguments forwarded to eggNOG-mapper (e.g., "--pident 30").'),
        'merge_min_orf_length': (
            ['Params_r_sample', '05_merge_tax_fun_depth', 'minORFLen'], 50, int,
            '05_merge_tax_fun_depth: Minimum ORF length for finally used ORFs of reads.'
        ),
        'merge_multiple_ko_handling_mode': (
            ['Params_r_sample', '05_merge_tax_fun_depth', 'KO_mode'], 'ignore_multiple_KOs', str,
            '\b\n05_merge_tax_fun_depth: Method using one of [split_multiple_KOs/keep_multiple_KOs/ignore_multiple_KOs] when an ORF has multiple KO annotations:.\n'
            'split_multiple_KOs: evenly distributing the total abundance across all assigned KOs;\n'
            'keep_multiple_KOs: retain original multiple KO annotations;\n'
            'ignore_multiple_KOs: ignore the ORFs with multiple KO annotations and only use the ORFs with single KO annotation.'
        ),
        'merge_min_ko_num_per_taxon': (
            ['Params_r_sample', '05_merge_tax_fun_depth', 'min_KO_num_in_one_taxon'], 100, int,
            '\b\n05_merge_tax_fun_depth: Minimum number of KOs per taxon, used to filter outputs in KO_Taxa_RelativeAbundance and KO_Taxa_SumCoverage.\n'
            'This threshold filters out taxa with insufficient KOs,\n'
            'retaining only KO-rich taxa to ensure robust metabolic reconstruction\n'
            'while preventing excessive output file sizes caused by sparse matrices.'
        ),
        'merge_tax_fun_depth_threads': (
            ['Params_r_sample', '05_merge_tax_fun_depth', 'threads'], 60, int,
            '05_merge_tax_fun_depth: Threads.'),
        'merge_cami_format_output': (
            ['Params_r_sample', '05_merge_tax_fun_depth', 'CAMI_format_output'], False, bool,
            '05_merge_tax_fun_depth: Report a taxonomic profile using the CAMI output format (flag only).'
        ),
        'merge_output_detail': (
            ['Params_r_sample', '05_merge_tax_fun_depth', 'output_detail'], False, bool,
            '05_merge_tax_fun_depth: Whether to output detail tables (flag only).'),
        'mclink_sample_min_ko_num_per_taxon': (
            ['Params_r_sample', '06_mclink_taxon', 'min_KO_num_in_one_taxon'], 100, int,
            '06_mclink_taxon: Minimum number of KOs per taxon for calculating taxon-level metabolic pathway/module completeness or abundance per sample. Taxa with fewer than this number of detected KOs in a given sample will be filtered out. Calculated independently for each sample.'
            ),
        'mclink_sample_min_ko_abundance_per_taxon': (
            ['Params_r_sample', '06_mclink_taxon', 'min_mean_KO_RA_in_one_taxon'], 0, float,
            '06_mclink_taxon: Minimum mean KO abundance [0-1e6] per taxon for calculating taxon-level metabolic pathway/module completeness or abundance per sample. Taxa with mean KO abundance below this threshold in a given sample will be filtered out. Note: GPM/TPM sum of all taxa in each sample is 1e6.'
            ),
        'mclink_sample_topn_taxon': (
            ['Params_r_sample', '06_mclink_taxon', 'TopN'], 100, int,
            '06_mclink_taxon: Number of top taxa (ranked by mean KO abundance) retained for calculating taxon-level metabolic pathway/module completeness or abundance per sample. Set to 0 to disable filtering.'
        ),
        'mclink_sample_completeness_per_taxon': (
            ['Params_r_sample', '06_mclink_taxon', 'completeness_threshold'], 0.5, float,
            '06_mclink_taxon: Completeness threshold [0-1] of modules for calculating taxon-level metabolic pathway/module abundance. Higher values retain only more complete profiles.'),
        'mclink_sample_plus_scale_method_per_taxon': (
            ['Params_r_sample', '06_mclink_taxon', 'plus_scale_method'], 'mean', str,
             "06_mclink_taxon: Aggregation method for '+'-separated KOs in a module [mean/max/min]. 'mean': average abundance across all required KOs; 'max': highest KO abundance; 'min': lowest KO abundance."),
        'mclink_sample_comma_scale_method_per_taxon': (
            ['Params_r_sample', '06_mclink_taxon', 'comma_scale_method'], 'sum', str,
            "06_mclink_taxon: Aggregation method for ','-separated alternative KOs in a module [sum/max]. 'sum': cumulative abundance of all alternative KOs; 'max': abundance of the highest single KO among alternatives. Note: this applies to abundance calculation only; module completeness is always computed with 'max'."),
        'barplot_topn_taxon': (
            ['Params_r_project', '01_Taxa_visualize', 'TopN'], 10, int,
            '01_Taxa_visualize: Number of top taxa to display in stacked barplots. Remaining categories are merged into Others. Maximum 29.'
        ),
        'barplot_topn_fun': (
            ['Params_r_project', '03_Fun_visualize', 'TopN'], 10, int,
            '03_Fun_visualize: Number of top functional categories to display in stacked barplots. Remaining categories are merged into Others. Maximum 29.'
        ),
        'mclink_min_ko_num_per_sample': (
            ['Params_r_project', '03_mclink_sample', 'min_sample_num_per_KO'], 1, int,
            '03_mclink_sample: Minimum number of samples in which a KO must be detected for calculating community-level metabolic pathway/module completeness or abundance.'
        ),
        'mclink_min_ko_prevalence_per_sample': (
            ['Params_r_project', '03_mclink_sample', 'preva'], 0, float,
            '03_mclink_sample: Minimum detection rate [0-1] of a KO across all samples for calculating community-level metabolic pathway/module completeness or abundance.'
        ),
        'mclink_min_ko_abundance_per_sample': (
            ['Params_r_project', '03_mclink_sample', 'abund'], 0, float,
            '03_mclink_sample: Minimum mean relative abundance [0-1e6] of a KO across all samples for calculating community-level metabolic pathway/module completeness or abundance. KOs with mean abundance below this threshold are excluded. Note: GPM/TPM sum of a sample is 1e6.'
        ),
        'mclink_completeness_per_sample': (
            ['Params_r_project', '03_mclink_sample', 'completeness_threshold'], 0.5, float,
            '03_mclink_sample: Completeness threshold [0-1] of modules for calculating community-level metabolic pathway/module abundance. Higher values retain only more complete profiles.'
        ),
        'mclink_plus_scale_method_per_sample': (
            ['Params_r_project', '03_mclink_sample', 'plus_scale_method'], 'mean', str,
            "03_mclink_sample: Aggregation method for '+'-separated KOs in a module [mean/max/min]. 'mean': average abundance across all required KOs; 'max': highest KO abundance; 'min': lowest KO abundance."
        ),
        'mclink_comma_scale_method_per_sample': (
            ['Params_r_project', '03_mclink_sample', 'comma_scale_method'], 'sum', str,
            "03_mclink_sample: Aggregation method for ','-separated alternative KOs in a module [sum/max]. 'sum': cumulative abundance of all alternative KOs; 'max': abundance of the highest single KO among alternatives. Note: this applies to abundance calculation only; module completeness is always computed with 'max'."
        ),
        'merge_min_taxa_ko_prevalence_across_samples': (
            ['Params_r_project', '05_merge_tax_KO_depth', 'preva'], 0.10, float,
            '05_merge_tax_KO_depth: Minimum prevalence [0-1] of a Taxon-KO across samples. Taxa-KOs detected in fewer samples than this proportion are excluded.'
        ),
        'merge_min_taxa_ko_abundance_across_samples': (
            ['Params_r_project', '05_merge_tax_KO_depth', 'abund'], 0.0, float,
            '05_merge_tax_KO_depth: Minimum mean relative abundance [0-1e6] of a Taxon-KO across samples. Taxa-KOs with mean abundance below this threshold are excluded. Note: GPM/TPM sum of a sample is 1e6.'
        ),
        'merge_min_taxa_module_prevalence_across_samples': (
            ['Params_r_project', '05_merge_tax_Module_RelativeAbundance', 'preva'], 0.10, float,
            '05_merge_tax_Module_RelativeAbundance: Minimum prevalence [0-1] of a Taxon-Module across samples. Taxa-Modules detected in fewer samples than this proportion are excluded.'
        ),
        'merge_min_taxa_module_abundance_across_samples': (
            ['Params_r_project', '05_merge_tax_Module_RelativeAbundance', 'abund'], 0.0, float,
            '05_merge_tax_Module_RelativeAbundance: Minimum mean relative abundance [0-1e6] of a Taxon-Module across samples. Taxa-Modules with mean abundance below this threshold are excluded. Note: GPM/TPM sum of a sample is 1e6.'
        ),
        'merge_min_taxa_module_sc_prevalence_across_samples': (
            ['Params_r_project', '05_merge_tax_Module_SumCoverage', 'preva'], 0.10, float,
            '05_merge_tax_Module_SumCoverage: Minimum prevalence [0-1] of a Taxon-Module across samples. Taxa-Modules detected in fewer samples than this proportion are excluded.'
        ),
        'merge_min_taxa_module_sc_abundance_across_samples': (
            ['Params_r_project', '05_merge_tax_Module_SumCoverage', 'abund'], 0.0, float,
            '05_merge_tax_Module_SumCoverage: Minimum mean sum coverage [0-1e6] of a Taxon-Module across samples. Taxa-Modules with mean sum coverage below this threshold are excluded. Note: GPM/TPM sum of a sample is 1e6.'
        ),
        'mclink_project_min_ko_num_per_taxon': (
            ['Params_r_project', '06_mclink_taxon', 'min_KO_num_in_one_taxon'], 100, int,
            '06_mclink_taxon: Minimum number of KOs per taxon for calculating taxon-level metabolic pathway/module completeness on the project-averaged taxon-function matrix. Taxa with fewer than this number of detected KOs in the project-averaged matrix will be filtered out. KO abundances are first averaged across all samples to generate a single project-level matrix before calculation.'
        ),
        'mclink_project_min_ko_abundance_per_taxon': (
            ['Params_r_project', '06_mclink_taxon', 'min_mean_KO_RA_in_one_taxon'], 0, float,
            '06_mclink_taxon: Minimum mean KO abundance [0-1e6] per taxon for calculating taxon-level metabolic pathway/module completeness on the project-averaged taxon-function matrix. Taxa with mean KO abundance below this threshold in the project-averaged matrix will be filtered out. Note: GPM/TPM sum of all taxa is 1e6 per sample; values are averaged across samples prior to filtering and calculation.'
        ),
        'mclink_project_topn_taxon': (
            ['Params_r_project', '06_mclink_taxon', 'TopN'], 100, int,
            '06_mclink_taxon: Number of top taxa (ranked by mean KO abundance in the project-averaged taxon-function matrix) retained for calculating taxon-level metabolic pathway/module completeness. Set to 0 to disable filtering.'
        ),
        'mclink_project_completeness_per_taxon': (
            ['Params_r_project', '06_mclink_taxon', 'completeness_threshold'], 0.5, float,
            '06_mclink_taxon: Completeness threshold [0-1] of modules for calculating taxon-level metabolic pathway/module abundance on the project-averaged taxon-function matrix. Higher values retain only more complete profiles.'),
        'mclink_project_plus_scale_method_per_taxon': (
            ['Params_r_project', '06_mclink_taxon', 'plus_scale_method'], 'mean', str,
            "06_mclink_taxon: Aggregation method for '+'-separated KOs in a module [mean/max/min] on the project-averaged taxon-function matrix. 'mean': average abundance across all required KOs; 'max': highest KO abundance; 'min': lowest KO abundance."),
        'mclink_project_comma_scale_method_per_taxon': (
            ['Params_r_project', '06_mclink_taxon', 'comma_scale_method'], 'sum', str,
            "06_mclink_taxon: Aggregation method for ','-separated alternative KOs in a module [sum/max] on the project-averaged taxon-function matrix. 'sum': cumulative abundance of all alternative KOs; 'max': abundance of the highest single KO among alternatives. Note: this applies to abundance calculation only; module completeness is always computed with 'max'."),
    },
    'm': {
        'raw_stats_sample_threads': (
            ['Params_b_sample', '00_raw_stats', 'threads'], 4, int,
            '00_raw_stats: Threads for seqkit.'),
        'trim_sample_threads': (
            ['Params_b_sample', '01_trim', 'threads'], 16, int,
            '01_trim: Threads for fastp (next-generation short reads) or fastplong (third-generation long reads).'),
        'trim_length_required': (
            ['Params_b_sample', '01_trim', 'length_required'], 15, int,
            '01_trim: Reads shorter than length_required will be discarded'),
        'trim_n_base_limit': (
            ['Params_b_sample', '01_trim', 'n_base_limit'], 5, int,
            "01_trim: if one read's number of N base is >n_base_limit, then this read/pair is discarded."),
        'trim_compression': (
            ['Params_b_sample', '01_trim', 'compression'], 4, int,
            '01_trim: Compression level for gzip output [1-9]. 1 is fastest, 9 is smallest.'),
        'trim_extra_args': (
            ['Params_b_sample', '01_trim', 'extra_args'], '', str,
            '01_trim: Extra arguments forwarded to fastp or fastplong for raw reads trimming (e.g., "--reads_to_process 0 -Q").'),
        'trim_stats_sample_threads': (
            ['Params_b_sample', '01_trim_stats', 'threads'], 16, int,
            '01_trim_stats: Threads.'),
        'assembly_method': (
            ['Params_b_sample', '02_assembly', 'method'], 'megahit', str,
            '\b\n02_assembly: Underlying assembly software using one of [megahit/spades/metaflye].\n'
            'Must be compatible with --read-type (config key read_type):\n'
            'megahit: MEGAHIT (second-generation sequencing; for paired or single DNA reads).\n'
            'spades: SPAdes (second-generation sequencing; metaSPAdes for paired DNA reads, rnaSPAdes for paired-rna or single-rna RNA reads).\n'
            'metaflye: metaFlye (third-generation sequencing; required for pacbio-raw/corr/hifi or nano-raw/corr/hq long reads).'),
        'assembly_threads': (
            ['Params_b_sample', '02_assembly', 'threads'], 40, int,
            '02_assembly: Threads for [megahit/spades/metaflye].'),
        'megahit_presets': (
            ['Params_b_sample', '02_megahit', 'presets'], 'meta-sensitive', str,
            '\b\n02_megahit: MEGAHIT presets using one of [meta-sensitive/meta-large].\n'
            "meta-sensitive: '--min-count 1 --k-list 21,29,39,49,...,129,141'.\n"
            "meta-large: '--k-min 27 --k-max 127 --k-step 10'.\n"
            '(large & complex metagenomes, like soil).'),
        'megahit_k_list': (
            ['Params_b_sample', '02_megahit', 'k_list'], '', str,
            '02_megahit: MEGAHIT --k-list (optional). Explicit comma-separated k-mer sizes (e.g. "21,29,39,59,79,99,119,141"). Overrides --megahit-presets. If set, --megahit-k-min, --megahit-k-max, --megahit-k-step must NOT be set. Optional --megahit-min-count can be used together.'),
        'megahit_min_count': (
            ['Params_b_sample', '02_megahit', 'min_count'], 2, int,
            '02_megahit: MEGAHIT --min-count (optional). Minimum multiplicity for filtering (k_min+1)-mers. Can be used with --megahit-k-list or --megahit-k-min/--k-max/--k-step.'),
        'megahit_k_min': (
            ['Params_b_sample', '02_megahit', 'k_min'], '', str,
            '02_megahit: MEGAHIT --k-min (optional). Minimum kmer size (<= 255), must be odd number. If set with --k-max and --k-step, overrides --megahit-presets and uses custom k-mer mode.'),
        'megahit_k_max': (
            ['Params_b_sample', '02_megahit', 'k_max'], '', str,
            '02_megahit: MEGAHIT --k-max (optional). Maximum kmer size (<= 255), must be odd number.'),
        'megahit_k_step': (
            ['Params_b_sample', '02_megahit', 'k_step'], '', str,
            '02_megahit: MEGAHIT --k-step (optional). Increment of kmer size of each iteration (<= 28), must be even number.'),
        'megahit_extra_args': (
            ['Params_b_sample', '02_megahit', 'extra_args'], '', str,
            '02_megahit: Extra arguments forwarded to MEGAHIT (e.g. "--no-mercy --tmp-dir /path/to/your/tmp/dir"). Do NOT use this for k-mer settings (--k-list, --k-min, --presets, etc.).'),
        'spades_extra_args': (
            ['Params_b_sample', '02_spades', 'extra_args'], '', str,
            '02_spades: Extra arguments forwarded to SPAdes. Do NOT use this for input/output or thread settings.'),
        'metaflye_extra_args': (
            ['Params_b_sample', '02_metaflye', 'extra_args'], '', str,
            '02_metaflye: Extra arguments forwarded to Flye for third-generation long reads assembly (e.g., "--iterations 2"). Do NOT use this for read type or input/output settings.'),
        'mapping_threads': (
            ['Params_b_sample', '03_mapping', 'threads'], 20, int,
            '03_mapping: Threads for [minibwa/bwamem/bwamem2/strobealign].'),
        'mapping_method': (
            ['Params_b_sample', '03_mapping', 'method'], 'minibwa', str,
            '\b\n03_mapping: Underlying mapping software used [minibwa/bwamem/bwamem2/strobealign].\n'
            'minibwa: minibwa using default parameters.\n'
            'bwamem: bwa mem using default parameters.\n'
            'bwamem2: bwa-mem2 using default parameters.\n'
            'strobealign: strobealign using default parameters.'),
        'mapping_samtools_sort_mem': (
            ['Params_b_sample', '03_mapping', 'samtools_sort_mem'], '768M', str,
            '03_mapping: Set maximum memory per thread; suffix K/M/G recognized, e.g. 2G, 768M.'),
        'mapping_extra_args': (
            ['Params_b_sample', '03_mapping', 'extra_args'], '', str,
            '03_mapping: Extra arguments forwarded to the active mapping aligner (minibwa map, bwa mem, bwa-mem2 mem, or strobealign --use-index), (e.g. "-B 8" for minibwa map).'),
        'binning_metawrap_depth_threads': (
            ['Params_b_sample', '05_metawrap_binning_depth', 'threads'], 8, int,
            '05_metawrap_binning_depth: Threads.'),
        'binning_min_read_percent_identity': (
            ['Params_b_sample', '05_metawrap_binning_depth', 'min_read_percent_identity'], 97, float,
            '05_metawrap_binning_depth: Exclude reads by overall percent identity [0-100] during the binning process, e.g. 97 for 97%.'),
        'binning_contig_end_exclusion': (
            ['Params_b_sample', '05_metawrap_binning_depth', 'contig_end_exclusion'], 75, int,
            '05_metawrap_binning_depth: Exclude bases at the ends of contig sequences used for depth calculation and bining.'),
        'binning_min_contig_length': (
            ['Params_b_sample', '05_metawrap_binning_depth', 'min_contig_length'], 1500, int,
            '05_metawrap_binning_depth: Minimum contig length for binning.'),
        'binning_tools': (
            ['Params_b_sample', '05_metawrap_binning_depth', 'binning_tools'], 'metabat2,maxbin2,concoct', str,
            '\b\n05_metawrap_binning_depth: Comma-separated list of binning tools, choose 2 or 3 from [metabat2,maxbin2,concoct].\n'
            'At least 2 tools recommended for bin_refinement.'),
        'binning_max_memory_gb': (
            ['Params_b_sample', '05_metawrap_binning_depth', 'max_memory_gb'], 512, int,
            '05_metawrap_binning_depth: Maximum memory (GB) for MetaWRAP binning.'),
        'binning_metawrap_extra_args': (
            ['Params_b_sample', '05_metawrap_binning_depth', 'extra_args'], '', str,
            '05_metawrap_binning_depth: Extra arguments forwarded to metaWRAP binning (e.g., "-m 5").'),
        'metawrap_refinement_max_memory_gb': (
            ['Params_b_sample', '06_metawrap_bin_refinement', 'max_memory_gb'], 512, int,
            '06_metawrap_bin_refinement: Maximum memory in GB for MetaWRAP bin_refinement.'),
        'metawrap_refinement_extra_args': (
            ['Params_b_sample', '06_metawrap_bin_refinement', 'extra_args'], '', str,
            '06_metawrap_bin_refinement: Extra arguments forwarded to metaWRAP bin_refinement (e.g., "--remove-ambiguous").'),
        'bin_checkm_threads': (
            ['Params_b_sample', '06_bin_checkm', 'threads'], 8, int,
            '06_bin_checkm or 06_metawrap_bin_refinement: Threads.'),
        'bin_min_completeness': (
            ['Params_b_sample', '06_bin_checkm', 'min_completeness'], 50, float,
            '06_bin_checkm or 06_metawrap_bin_refinement: Minimum completeness (%) bin filtering, e.g. 50 for 50%.'),
        'bin_max_contamination': (
            ['Params_b_sample', '06_bin_checkm', 'max_contamination'], 10, float,
            '06_bin_checkm or 06_metawrap_bin_refinement: Maximum contamination (%) bin filtering, e.g. 10 for 10%.'),
        'bin_min_quality_score': (
            ['Params_b_sample', '06_bin_checkm', 'min_quality_score'], 50, float,
            '06_bin_checkm: Minimum quality score for bin filtering (GTDB-like), QS = completeness - 5*contamination, e.g. 50 for QS>=50.'),
        'contigs_depth_sample_threads': (
            ['Params_b_sample', '07_contigs_depth', 'threads'], 8, int,
            '07_contigs_depth: Threads.'),
        'depth_min_read_percent_identity': (
            ['Params_b_sample', '07_contigs_depth', 'min_read_percent_identity'], 95, float,
            '\b\n07_contigs_depth: Exclude reads by overall percent identity [0-100] for contig depth estimation\n'
            'using jgi_summarize_bam_contig_depths --percentIdentity, e.g. 95 for 95%.\n'
            'Reads with percent identity below this threshold are excluded from depth calculation.\n'
            'Higher values increase stringency, retaining only high-quality alignments.'),
        'depth_contig_end_exclusion': (
            ['Params_b_sample', '07_contigs_depth', 'contig_end_exclusion'], 75, int,
            '\b\n07_contigs_depth: Exclude bases at the ends of contig sequences for depth estimation using jgi_summarize_bam_contig_depths --maxEdgeBases.\n'
            'Bases within this distance from either end of a contig are ignored when calculating depth,\n'
            'preventing edge effects from partial alignments.'),
        'depth_min_contig_length': (
            ['Params_b_sample', '07_contigs_depth', 'min_contig_length'], 300, int,
            '\b\n07_contigs_depth: Minimum contig length for depth estimation [bp], e.g. 300 for 300 bp.\n'
            'Contigs shorter than this threshold are excluded from depth calculation.'),
        'contigs_prodigal_extra_args': (
            ['Params_b_sample', '07_contigs_depth', 'extra_args'], '', str,
            '07_contigs_depth: Extra arguments forwarded to prodigal for contig gene prediction (e.g., "-m -n").'),
        'bin_prodigal_threads': (
            ['Params_b_sample', '08_bin_prodigal', 'threads'], 8, int,
            '08_bin_prodigal: Threads.'),
        'bin_prodigal_extra_args': (
            ['Params_b_sample', '08_bin_prodigal', 'extra_args'], '', str,
            '08_bin_prodigal: Extra arguments forwarded to prodigal for MAG gene prediction (e.g., "-m -n").'),
        'raw_stats_project_threads': (
            ['Params_b_project', '01_cat_raw_stats', 'threads'], 4, int,
            '01_cat_raw_stats: Threads.'),
        'trimmed_stats_project_threads': (
            ['Params_b_project', '02_cat_trimmed_stats', 'threads'], 4, int,
            '02_cat_trimmed_stats: Threads.'),
        'contigs_stats_project_threads': (
            ['Params_b_project', '03_contigs_stats', 'threads'], 4, int,
            '03_contigs_stats: Threads.'),
        'merge_stats_project_threads': (
            ['Params_b_project', '04_merge_reads_contigs_stats', 'threads'], 4, int,
            '04_merge_reads_contigs_stats: Threads.'),
        'gather_checkm_project_threads': (
            ['Params_b_project', '05_gather_project_bin_checkm_stats', 'threads'], 4, int,
            '05_gather_project_bin_checkm_stats: Threads.'),
        'mmags_checkm2_threads': (
            ['Params_m_sample', '01_checkm2', 'threads'], 60, int,
            '01_checkm2: Threads.'),
        'mmags_checkm2_extra_args': (
            ['Params_m_sample', '01_checkm2', 'extra_args'], '', str,
            '01_checkm2: Extra arguments forwarded to checkm2 predict (e.g., "--remove_intermediates").'),
        'mmags_taxonomy_classify_method': (
            ['Params_m_sample', '02_taxonomy', 'MAGs_taxonomy_classify_method'], 'gtdbtk', str,
            '02_taxonomy: mMAGs taxonomy classifier, should be one of the [gtdbtk/BAT].'
        ),
        'mmags_taxonomy_threads': (
            ['Params_m_sample', '02_taxonomy', 'threads'], 60, int,
            '02_taxonomy_gtdbtk: Threads for [gtdbtk/BAT].'),
        'mmags_taxonomy_gtdbtk_extra_args': (
            ['Params_m_sample', '02_taxonomy_gtdbtk', 'extra_args'], '', str,
            '02_taxonomy_gtdbtk: Extra arguments forwarded to gtdbtk classify_wf (e.g., "--min_af 0.5").'),
        'mmags_taxonomy_bat_range': (
            ['Params_m_sample', '02_taxonomy_BAT', 'range'], 10, float,
            '\b\n02_taxonomy_BAT: Percentage range [0-100] of DIAMOND hits relative to top bit-score to include in LCA calculation.\n'
            'Higher values yield more conservative classifications by incorporating divergent hits, while lower values enable specific assignments.\n'
            'Applicable when --mMAGs-taxonomy-classify-method is set to BAT.'
        ),
        'mmags_taxonomy_bat_fraction': (
            ['Params_m_sample', '02_taxonomy_BAT', 'fraction'], 0.50, float,
            '\b\n02_taxonomy_BAT: Minimum bit-score support fraction [0-0.99].\n'
            'Higher values enforce stricter consensus (one classification per sequence), while lower values allow multiple lineage detections.\n'
            'Applicable when --mMAGs-taxonomy-classify-method is set to BAT.'
        ),
        'mmags_taxonomy_bat_extra_args': (
            ['Params_m_sample', '02_taxonomy_BAT', 'extra_args'], '', str,
            '02_taxonomy_BAT: Extra arguments forwarded to CAT_pack bins (e.g., "--block_size 12.0").'),
        'mmags_drep99_completeness_per_sample': (
            ['Params_m_sample', '05_drep_eachsample_99', 'comp'], 50, float,
            '05_drep_eachsample_99: The first level dereplication of mMAGs in a sample using completeness threshold [0-100], e.g., 50 for 50%.'
        ),
        'mmags_drep99_contamination_per_sample': (
            ['Params_m_sample', '05_drep_eachsample_99', 'con'], 10, float,
            '05_drep_eachsample_99: The first level dereplication of mMAGs in a sample using contamination threshold [0-100], e.g., 10 for 10%.'
        ),
        'mmags_drep99_ani_per_sample': (
            ['Params_m_sample', '05_drep_eachsample_99', 'sa'], 0.99, float,
            '05_drep_eachsample_99: The first level dereplication of mMAGs in a sample using ANI threshold [0-1], e.g., 0.99 for 99%.'
        ),
        'mmags_drep99_sample_threads': (
            ['Params_m_sample', '05_drep_eachsample_99', 'threads'], 20, int,
            '05_drep_eachsample_99: Threads.'),
        'mmags_drep99_sample_extra_args': (
            ['Params_m_sample', '05_drep_eachsample_99', 'extra_args'], '', str,
            '05_drep_eachsample: Extra arguments forwarded to dRep dereplicate (e.g., "--completeness_weight 1 --contamination_weight 5").'),
        'mmags_drep95_completeness_per_sample': (
            ['Params_m_sample', '06_drep_eachsample_95', 'comp'], 50, float,
            '06_drep_eachsample_95: The second level dereplication of mMAGs in a sample using completeness threshold [0-100], e.g., 50 for 50%.'
        ),
        'mmags_drep95_contamination_per_sample': (
            ['Params_m_sample', '06_drep_eachsample_95', 'con'], 10, float,
            '06_drep_eachsample_95: The second level dereplication of mMAGs in a sample using contamination threshold [0-100], e.g., 10 for 10%.'
        ),
        'mmags_drep95_ani_per_sample': (
            ['Params_m_sample', '06_drep_eachsample_95', 'sa'], 0.95, float,
            '\b\n06_drep_eachsample_95: The second level dereplication of mMAGs in a sample using ANI threshold [0-1], e.g., 0.95 for 95%.\n'
            'This value should be smaller than --drep99-ANI-per-sample'
        ),
        'mmags_drep95_sample_threads': (
            ['Params_m_sample', '06_drep_eachsample_95', 'threads'], 20, int,
            '06_drep_eachsample_95: Threads.'),
        'mmags_drep95_sample_extra_args': (
            ['Params_m_sample', '06_drep_eachsample_95', 'extra_args'], '', str,
            '06_drep_eachsample: Extra arguments forwarded to dRep dereplicate (e.g., "--completeness_weight 1 --contamination_weight 5").'),
        'mmags_coverm_method': (
            ['Params_m_sample', '07_coverm_eachsample', 'method'], 'trimmed_mean', str,
            '\b\n07_coverm_eachsample: Method(s) for coverm-based coverage calculation, must be one of [trimmed_mean/mean].\n'
            'A more detailed description of the different methods is available at https://github.com/wwood/CoverM#calculation-methods'
        ),
        'mmags_coverm_min_read_percent_identity': (
            ['Params_m_sample', '07_coverm_eachsample', 'min_read_percent_identity'], 95, float,
            '07_coverm_eachsample: Exclude reads by overall percent identity [0-100], e.g. 95 for 95%.'
        ),
        'mmags_coverm_min_read_aligned_length': (
            ['Params_m_sample', '07_coverm_eachsample', 'min_read_aligned_length'], 45, int,
            '07_coverm_eachsample: Exclude reads with fewer aligned bases.'
        ),
        'mmags_coverm_contig_end_exclusion': (
            ['Params_m_sample', '07_coverm_eachsample', 'contig_end_exclusion'], 75, int,
            '07_coverm_eachsample: Exclude bases at the ends of contig sequences for coverm calculation.'
        ),
        'mmags_coverm_threads': (
            ['Params_m_sample', '07_coverm_eachsample', 'threads'], 20, int,
            '07_coverm_eachsample: Threads.'),
        'mmags_coverm_sample_extra_args': (
            ['Params_m_sample', '07_coverm_eachsample', 'extra_args'], '', str,
            '07_coverm_eachsample: Extra arguments forwarded to coverm genome (e.g., "--quiet").'),
        'mmags_function_classify_method': (
            ['Params_m_sample', '08_function', 'MAGs_function_classify_method'], 'emapper', str,
            '08_function: mMAGs function classifier using one of [kraken2/emapper].'
        ),
        'mmags_function_threads': (
            ['Params_m_sample', '08_function', 'threads'], 60, int,
            '08_function: Threads for [kraken2/emapper].'),
        'mmags_function_kraken2_confidence': (
            ['Params_m_sample', '08_function_kraken2', 'confidence'], 0.05, float,
            '08_function_kraken2: Confidence threshold [0-1] for functional classification for mMAGs when --mMAGs-function-classify-method is set to kraken2.'
        ),
        'mmags_function_kraken2_extra_args': (
            ['Params_m_sample', '08_function_kraken2', 'extra_args'], '', str,
            '08_function_kraken2: Extra arguments forwarded to kraken2 (e.g., "--memory-mapping").'),
        'mmags_function_emapper_evalue': (
            ['Params_m_sample', '08_function_emapper', 'seed_ortholog_evalue'], 1e-10, float,
            '08_function_emapper: Report only alignments with e-value below or equal to the threshold, applicable when --mMAGs-function-classify-method is set to emapper.'
        ),
        'mmags_function_emapper_method': (
            ['Params_m_sample', '08_function_emapper', 'method'], 'diamond', str,
            '\b\n08_function_emapper: search seed orthologs using one of [diamond/mmseqs], applicable when --mMAGs-function-classify-method is set to emapper.\n'
            'diamond: using default parameters of diamond mode of emapper.\n'
            'mmseqs: using default parameters of mmseqs mode of emapper'
        ),
        'mmags_function_emapper_query_cover': (
            ['Params_m_sample', '08_function_emapper', 'query_cover'], 30, float,
            '\b\n08_function_emapper: Report only alignments with query cover above or equal to the given percentage [0-100],.\n'
            'applicable when --mMAGs-function-classify-method is set to emapper'
        ),
        'mmags_function_emapper_extra_args': (
            ['Params_m_sample', '08_function_emapper', 'extra_args'], '', str,
            '08_function_emapper: Extra arguments forwarded to eggNOG-mapper (e.g., "--pident 30").'),
        'merge_mmags_min_completeness1': (
            ['Params_m_sample', '05_merge_tax_fun_depth', 'MAGs_min_completeness1'], 50, float,
            '02_quality_gate: Minimum completeness threshold (%) for CheckM1-based mMAG quality filtering [0-100], e.g., 50 for 50%.'
        ),
        'merge_mmags_max_contamination1': (
            ['Params_m_sample', '05_merge_tax_fun_depth', 'MAGs_max_contamination1'], 10, float,
            '02_quality_gate: Maximum contamination threshold (%) for CheckM1-based mMAG quality filtering [0-100], e.g., 10 for 10%.'
        ),
        'merge_mmags_min_completeness2': (
            ['Params_m_sample', '05_merge_tax_fun_depth', 'MAGs_min_completeness2'], 50, float,
            '02_quality_gate: Minimum completeness threshold (%) for CheckM2-based mMAG quality filtering [0-100], e.g., 50 for 50%.'
        ),
        'merge_mmags_max_contamination2': (
            ['Params_m_sample', '05_merge_tax_fun_depth', 'MAGs_max_contamination2'], 10, float,
            '02_quality_gate: Maximum contamination threshold (%) for CheckM2-based mMAG quality filtering [0-100], e.g., 10 for 10%.'
        ),
        'merge_min_orf_length': (
            ['Params_m_sample', '05_merge_tax_fun_depth', 'minORFLen'], 50, int,
            '05_merge_tax_fun_depth: Minimum ORF length for finally used ORFs of mMAGs.'
        ),
        'merge_multiple_ko_handling_mode': (
            ['Params_m_sample', '05_merge_tax_fun_depth', 'KO_mode'], 'ignore_multiple_KOs', str,
            '\b\n05_merge_tax_fun_depth: Method using one of [split_multiple_KOs/keep_multiple_KOs/ignore_multiple_KOs] when an ORF has multiple KO annotations:.\n'
            'split_multiple_KOs: evenly distributing the total abundance across all assigned KOs;\n'
            'keep_multiple_KOs: retain original multiple KO annotations;\n'
            'ignore_multiple_KOs: ignore the ORFs with multiple KO annotations and only use the ORFs with single KO annotation.'
        ),
        'merge_min_ko_num_per_taxon': (
            ['Params_m_sample', '05_merge_tax_fun_depth', 'min_KO_num_in_one_taxon'], 100, int,
            '\b\n05_merge_tax_fun_depth: Minimum number of KOs per taxon, used to filter outputs in KO_Taxa_RelativeAbundance and KO_Taxa_SumCoverage.\n'
            'This threshold filters out taxa with insufficient KOs,\n'
            'retaining only KO-rich taxa to ensure robust metabolic reconstruction\n'
            'while preventing excessive output file sizes caused by sparse matrices.'
        ),
        'merge_tax_fun_depth_threads': (
            ['Params_m_sample', '05_merge_tax_fun_depth', 'threads'], 20, int,
            '05_merge_tax_fun_depth: Threads.'),
        'merge_cami_format_output': (
            ['Params_m_sample', '05_merge_tax_fun_depth', 'CAMI_format_output'], False, bool,
            '05_merge_tax_fun_depth: Report a taxonomic profile using the CAMI output format (flag only). Disabled when mMAG mode uses gtdbtk.'
        ),
        'merge_output_detail': (
            ['Params_m_sample', '05_merge_tax_fun_depth', 'output_detail'], False, bool,
            '05_merge_tax_fun_depth: Whether to output detail tables (flag only).'),
        'mclink_sample_min_ko_num_per_taxon': (
            ['Params_m_sample', '06_mclink_taxon', 'min_KO_num_in_one_taxon'], 100, int,
            '06_mclink_taxon: Minimum number of KOs per taxon for calculating taxon-level metabolic pathway/module completeness or abundance per sample. Taxa with fewer than this number of detected KOs in a given sample will be filtered out. Calculated independently for each sample.'
        ),
        'mclink_sample_min_ko_abundance_per_taxon': (
            ['Params_m_sample', '06_mclink_taxon', 'min_mean_KO_RA_in_one_taxon'], 0, float,
            '06_mclink_taxon: Minimum mean KO abundance [0-1e6] per taxon for calculating taxon-level metabolic pathway/module completeness or abundance per sample. Taxa with mean KO abundance below this threshold in a given sample will be filtered out. Note: GPM/TPM sum of all taxa in each sample is 1e6.'
        ),
        'mclink_sample_topn_taxon': (
            ['Params_m_sample', '06_mclink_taxon', 'TopN'], 100, int,
            '06_mclink_taxon: Number of top taxa (ranked by mean KO abundance) retained for calculating taxon-level metabolic pathway/module completeness or abundance per sample. Set to 0 to disable filtering.'
        ),
        'mclink_sample_completeness_per_taxon': (
            ['Params_m_sample', '06_mclink_taxon', 'completeness_threshold'], 0.5, float,
            '06_mclink_taxon: Completeness threshold [0-1] of modules for calculating taxon-level metabolic pathway/module abundance. Higher values retain only more complete profiles.'),
        'mclink_sample_plus_scale_method_per_taxon': (
            ['Params_m_sample', '06_mclink_taxon', 'plus_scale_method'], 'mean', str,
            "06_mclink_taxon: Aggregation method for '+'-separated KOs in a module [mean/max/min]. 'mean': average abundance across all required KOs; 'max': highest KO abundance; 'min': lowest KO abundance."),
        'mclink_sample_comma_scale_method_per_taxon': (
            ['Params_m_sample', '06_mclink_taxon', 'comma_scale_method'], 'sum', str,
            "06_mclink_taxon: Aggregation method for ','-separated alternative KOs in a module [sum/max]. 'sum': cumulative abundance of all alternative KOs; 'max': abundance of the highest single KO among alternatives. Note: this applies to abundance calculation only; module completeness is always computed with 'max'."),
        'mmags_drep99_completeness_per_project': (
            ['Params_m_project', '01_drep_project_99', 'comp'], 50, float,
            '01_drep_project_99: Dereplicate mMAGs from all samples in a project with completeness threshold [0-100], e.g., 50 for 50%.'
        ),
        'mmags_drep99_contamination_per_project': (
            ['Params_m_project', '01_drep_project_99', 'con'], 10, float,
            '01_drep_project_99: Dereplicate mMAGs from all samples in a project with contamination threshold [0-100], e.g., 5 for 5%.'
        ),
        'mmags_drep99_project_threads': (
            ['Params_m_project', '01_drep_project_99', 'threads'], 60, int,
            '01_drep_project_99: Threads.'),
        'mmags_drep99_project_extra_args': (
            ['Params_m_project', '01_drep_project_99', 'extra_args'], '', str,
            '01_drep_project_99: Extra arguments forwarded to dRep dereplicate ("--completeness_weight 1 --contamination_weight 5").'),
        'mmags_drep95_completeness_per_project': (
            ['Params_m_project', '01_drep_project_95', 'comp'], 50, float,
            '01_drep_project_95: Dereplicate mMAGs from all samples in a project with completeness threshold [0-100], e.g., 50 for 50%.'
        ),
        'mmags_drep95_contamination_per_project': (
            ['Params_m_project', '01_drep_project_95', 'con'], 10, float,
            '01_drep_project_95: Dereplicate mMAGs from all samples in a project with contamination threshold [0-100], e.g., 5 for 5%.'
        ),
        'mmags_drep95_project_threads': (
            ['Params_m_project', '01_drep_project_95', 'threads'], 60, int,
            '01_drep_project_95: Threads.'),
        'mmags_drep95_project_extra_args': (
            ['Params_m_project', '01_drep_project_95', 'extra_args'], '', str,
            '01_drep_project_95: Extra arguments forwarded to dRep dereplicate (e.g., "--completeness_weight 1 --contamination_weight 5").'),
        'barplot_topn_taxon': (
            ['Params_m_project', '01_Taxa_visualize', 'TopN'], 10, int,
            '01_Taxa_visualize: Number of top taxa to display in stacked barplots. Remaining categories are merged into Others. Maximum 29.'
        ),
        'barplot_topn_fun': (
            ['Params_m_project', '03_Fun_visualize', 'TopN'], 10, int,
            '03_Fun_visualize: Number of top functional categories to display in stacked barplots. Remaining categories are merged into Others. Maximum 29.'
        ),
        'mclink_min_ko_num_per_sample': (
            ['Params_m_project', '03_mclink_sample', 'min_sample_num_per_KO'], 1, int,
            '03_mclink_sample: Minimum number of samples in which a KO must be detected for calculating community-level metabolic pathway/module completeness or abundance.'
        ),
        'mclink_min_ko_prevalence_per_sample': (
            ['Params_m_project', '03_mclink_sample', 'preva'], 0, float,
            '03_mclink_sample: Minimum detection rate [0-1] of a KO across all samples for calculating community-level metabolic pathway/module completeness or abundance.'
        ),
        'mclink_min_ko_abundance_per_sample': (
            ['Params_m_project', '03_mclink_sample', 'abund'], 0, float,
            '03_mclink_sample: Minimum mean relative abundance [0-1e6] of a KO across all samples for calculating community-level metabolic pathway/module completeness or abundance. KOs with mean abundance below this threshold are excluded. Note: GPM/TPM sum of a sample is 1e6.'
        ),
        'mclink_completeness_per_sample': (
            ['Params_m_project', '03_mclink_sample', 'completeness_threshold'], 0.5, float,
            '03_mclink_sample: Completeness threshold [0-1] of modules for calculating community-level metabolic pathway/module abundance. Higher values retain only more complete profiles.'
        ),
        'mclink_plus_scale_method_per_sample': (
            ['Params_m_project', '03_mclink_sample', 'plus_scale_method'], 'mean', str,
            "03_mclink_sample: Aggregation method for '+'-separated KOs in a module [mean/max/min]. 'mean': average abundance across all required KOs; 'max': highest KO abundance; 'min': lowest KO abundance."
        ),
        'mclink_comma_scale_method_per_sample': (
            ['Params_m_project', '03_mclink_sample', 'comma_scale_method'], 'sum', str,
            "03_mclink_sample: Aggregation method for ','-separated alternative KOs in a module [sum/max]. 'sum': cumulative abundance of all alternative KOs; 'max': abundance of the highest single KO among alternatives. Note: this applies to abundance calculation only; module completeness is always computed with 'max'."
        ),
        'merge_min_taxa_ko_prevalence_across_samples': (
            ['Params_m_project', '05_merge_tax_KO_depth', 'preva'], 0.10, float,
            '05_merge_tax_KO_depth: Minimum prevalence [0-1] of a Taxon-KO across samples. Taxa-KOs detected in fewer samples than this proportion are excluded.'
        ),
        'merge_min_taxa_ko_abundance_across_samples': (
            ['Params_m_project', '05_merge_tax_KO_depth', 'abund'], 0.0, float,
            '05_merge_tax_KO_depth: Minimum mean relative abundance [0-1e6] of a Taxon-KO across samples. Taxa-KOs with mean abundance below this threshold are excluded. Note: GPM/TPM sum of a sample is 1e6.'
        ),
        'merge_min_taxa_module_prevalence_across_samples': (
            ['Params_m_project', '05_merge_tax_Module_RelativeAbundance', 'preva'], 0.10, float,
            '05_merge_tax_Module_RelativeAbundance: Minimum prevalence [0-1] of a Taxon-Module across samples. Taxa-Modules detected in fewer samples than this proportion are excluded.'
        ),
        'merge_min_taxa_module_abundance_across_samples': (
            ['Params_m_project', '05_merge_tax_Module_RelativeAbundance', 'abund'], 0.0, float,
            '05_merge_tax_Module_RelativeAbundance: Minimum mean relative abundance [0-1e6] of a Taxon-Module across samples. Taxa-Modules with mean abundance below this threshold are excluded. Note: GPM/TPM sum of a sample is 1e6.'
        ),
        'merge_min_taxa_module_sc_prevalence_across_samples': (
            ['Params_m_project', '05_merge_tax_Module_SumCoverage', 'preva'], 0.10, float,
            '05_merge_tax_Module_SumCoverage: Minimum prevalence [0-1] of a Taxon-Module across samples. Taxa-Modules detected in fewer samples than this proportion are excluded.'
        ),
        'merge_min_taxa_module_sc_abundance_across_samples': (
            ['Params_m_project', '05_merge_tax_Module_SumCoverage', 'abund'], 0.0, float,
            '05_merge_tax_Module_SumCoverage: Minimum mean sum coverage [0-1e6] of a Taxon-Module across samples. Taxa-Modules with mean sum coverage below this threshold are excluded. Note: GPM/TPM sum of a sample is 1e6.'
        ),
        'mclink_project_min_ko_num_per_taxon': (
            ['Params_m_project', '06_mclink_taxon', 'min_KO_num_in_one_taxon'], 100, int,
            '06_mclink_taxon: Minimum number of KOs per taxon for calculating taxon-level metabolic pathway/module completeness on the project-averaged taxon-function matrix. Taxa with fewer than this number of detected KOs in the project-averaged matrix will be filtered out. KO abundances are first averaged across all samples to generate a single project-level matrix before calculation.'
        ),
        'mclink_project_min_ko_abundance_per_taxon': (
            ['Params_m_project', '06_mclink_taxon', 'min_mean_KO_RA_in_one_taxon'], 0, float,
            '06_mclink_taxon: Minimum mean KO abundance [0-1e6] per taxon for calculating taxon-level metabolic pathway/module completeness on the project-averaged taxon-function matrix. Taxa with mean KO abundance below this threshold in the project-averaged matrix will be filtered out. Note: GPM/TPM sum of all taxa is 1e6 per sample; values are averaged across samples prior to filtering and calculation.'
        ),
        'mclink_project_topn_taxon': (
            ['Params_m_project', '06_mclink_taxon', 'TopN'], 100, int,
            '06_mclink_taxon: Number of top taxa (ranked by mean KO abundance in the project-averaged taxon-function matrix) retained for calculating taxon-level metabolic pathway/module completeness. Set to 0 to disable filtering.'
        ),
        'mclink_project_completeness_per_taxon': (
            ['Params_m_project', '06_mclink_taxon', 'completeness_threshold'], 0.5, float,
            '06_mclink_taxon: Completeness threshold [0-1] of modules for calculating taxon-level metabolic pathway/module abundance on the project-averaged taxon-function matrix. Higher values retain only more complete profiles.'),
        'mclink_project_plus_scale_method_per_taxon': (
            ['Params_m_project', '06_mclink_taxon', 'plus_scale_method'], 'mean', str,
            "06_mclink_taxon: Aggregation method for '+'-separated KOs in a module [mean/max/min] on the project-averaged taxon-function matrix. 'mean': average abundance across all required KOs; 'max': highest KO abundance; 'min': lowest KO abundance."),
        'mclink_project_comma_scale_method_per_taxon': (
            ['Params_m_project', '06_mclink_taxon', 'comma_scale_method'], 'sum', str,
            "06_mclink_taxon: Aggregation method for ','-separated alternative KOs in a module [sum/max] on the project-averaged taxon-function matrix. 'sum': cumulative abundance of all alternative KOs; 'max': abundance of the highest single KO among alternatives. Note: this applies to abundance calculation only; module completeness is always computed with 'max'."),
    },
    'mc': {
        'raw_stats_sample_threads': (
            ['Params_b_sample', '00_raw_stats', 'threads'], 4, int,
            '00_raw_stats: Threads for seqkit.'),
        'trim_sample_threads': (
            ['Params_b_sample', '01_trim', 'threads'], 16, int,
            '01_trim: Threads for fastp (next-generation short reads) or fastplong (third-generation long reads).'),
        'trim_length_required': (
            ['Params_b_sample', '01_trim', 'length_required'], 15, int,
            '01_trim: Reads shorter than length_required will be discarded'),
        'trim_n_base_limit': (
            ['Params_b_sample', '01_trim', 'n_base_limit'], 5, int,
            "01_trim: if one read's number of N base is >n_base_limit, then this read/pair is discarded."),
        'trim_compression': (
            ['Params_b_sample', '01_trim', 'compression'], 4, int,
            '01_trim: Compression level for gzip output [1-9]. 1 is fastest, 9 is smallest.'),
        'trim_extra_args': (
            ['Params_b_sample', '01_trim', 'extra_args'], '', str,
            '01_trim: Extra arguments forwarded to fastp or fastplong for raw reads trimming (e.g., "--reads_to_process 0 -Q").'),
        'trim_stats_sample_threads': (
            ['Params_b_sample', '01_trim_stats', 'threads'], 16, int,
            '01_trim_stats: Threads.'),
        'assembly_method': (
            ['Params_b_sample', '02_assembly', 'method'], 'megahit', str,
            '\b\n02_assembly: Underlying assembly software using one of [megahit/spades/metaflye].\n'
            'Must be compatible with --read-type (config key read_type):\n'
            'megahit: MEGAHIT (second-generation sequencing; for paired or single DNA reads).\n'
            'spades: SPAdes (second-generation sequencing; metaSPAdes for paired DNA reads, rnaSPAdes for paired-rna or single-rna RNA reads).\n'
            'metaflye: metaFlye (third-generation sequencing; required for pacbio-raw/corr/hifi or nano-raw/corr/hq long reads).'),
        'assembly_threads': (
            ['Params_b_sample', '02_assembly', 'threads'], 40, int,
            '02_assembly: Threads for [megahit/spades/metaflye].'),
        'megahit_presets': (
            ['Params_b_sample', '02_megahit', 'presets'], 'meta-sensitive', str,
            '\b\n02_megahit: MEGAHIT presets using one of [meta-sensitive/meta-large].\n'
            "meta-sensitive: '--min-count 1 --k-list 21,29,39,49,...,129,141'.\n"
            "meta-large: '--k-min 27 --k-max 127 --k-step 10'.\n"
            '(large & complex metagenomes, like soil).'),
        'megahit_k_list': (
            ['Params_b_sample', '02_megahit', 'k_list'], '', str,
            '02_megahit: MEGAHIT --k-list (optional). Explicit comma-separated k-mer sizes (e.g. "21,29,39,59,79,99,119,141"). Overrides --megahit-presets. If set, --megahit-k-min, --megahit-k-max, --megahit-k-step must NOT be set. Optional --megahit-min-count can be used together.'),
        'megahit_min_count': (
            ['Params_b_sample', '02_megahit', 'min_count'], 2, int,
            '02_megahit: MEGAHIT --min-count (optional). Minimum multiplicity for filtering (k_min+1)-mers. Can be used with --megahit-k-list or --megahit-k-min/--k-max/--k-step.'),
        'megahit_k_min': (
            ['Params_b_sample', '02_megahit', 'k_min'], '', str,
            '02_megahit: MEGAHIT --k-min (optional). Minimum kmer size (<= 255), must be odd number. If set with --k-max and --k-step, overrides --megahit-presets and uses custom k-mer mode.'),
        'megahit_k_max': (
            ['Params_b_sample', '02_megahit', 'k_max'], '', str,
            '02_megahit: MEGAHIT --k-max (optional). Maximum kmer size (<= 255), must be odd number.'),
        'megahit_k_step': (
            ['Params_b_sample', '02_megahit', 'k_step'], '', str,
            '02_megahit: MEGAHIT --k-step (optional). Increment of kmer size of each iteration (<= 28), must be even number.'),
        'megahit_extra_args': (
            ['Params_b_sample', '02_megahit', 'extra_args'], '', str,
            '02_megahit: Extra arguments forwarded to MEGAHIT (e.g. "--no-mercy --tmp-dir /path/to/your/tmp/dir"). Do NOT use this for k-mer settings (--k-list, --k-min, --presets, etc.).'),
        'spades_extra_args': (
            ['Params_b_sample', '02_spades', 'extra_args'], '', str,
            '02_spades: Extra arguments forwarded to SPAdes. Do NOT use this for input/output or thread settings.'),
        'metaflye_extra_args': (
            ['Params_b_sample', '02_metaflye', 'extra_args'], '', str,
            '02_metaflye: Extra arguments forwarded to Flye for third-generation long reads assembly (e.g., "--iterations 2"). Do NOT use this for read type or input/output settings.'),
        'mapping_threads': (
            ['Params_b_sample', '03_mapping', 'threads'], 20, int,
            '03_mapping: Threads for [minibwa/bwamem/bwamem2/strobealign].'),
        'mapping_method': (
            ['Params_b_sample', '03_mapping', 'method'], 'minibwa', str,
            '\b\n03_mapping: Underlying mapping software used [minibwa/bwamem/bwamem2/strobealign].\n'
            'minibwa: minibwa using default parameters.\n'
            'bwamem: bwa mem using default parameters.\n'
            'bwamem2: bwa-mem2 using default parameters.\n'
            'strobealign: strobealign using default parameters.'),
        'mapping_samtools_sort_mem': (
            ['Params_b_sample', '03_mapping', 'samtools_sort_mem'], '768M', str,
            '03_mapping: Set maximum memory per thread; suffix K/M/G recognized, e.g. 2G, 768M.'),
        'mapping_extra_args': (
            ['Params_b_sample', '03_mapping', 'extra_args'], '', str,
            '03_mapping: Extra arguments forwarded to the active mapping aligner (minibwa map, bwa mem, bwa-mem2 mem, or strobealign --use-index), (e.g. "-B 8" for minibwa map).'),
        'binning_metawrap_depth_threads': (
            ['Params_b_sample', '05_metawrap_binning_depth', 'threads'], 8, int,
            '05_metawrap_binning_depth: Threads.'),
        'binning_min_read_percent_identity': (
            ['Params_b_sample', '05_metawrap_binning_depth', 'min_read_percent_identity'], 97, float,
            '05_metawrap_binning_depth: Exclude reads by overall percent identity [0-100] during the binning process, e.g. 97 for 97%.'),
        'binning_contig_end_exclusion': (
            ['Params_b_sample', '05_metawrap_binning_depth', 'contig_end_exclusion'], 75, int,
            '05_metawrap_binning_depth: Exclude bases at the ends of contig sequences used for depth calculation and bining.'),
        'binning_min_contig_length': (
            ['Params_b_sample', '05_metawrap_binning_depth', 'min_contig_length'], 1500, int,
            '05_metawrap_binning_depth: Minimum contig length for binning.'),
        'binning_tools': (
            ['Params_b_sample', '05_metawrap_binning_depth', 'binning_tools'], 'metabat2,maxbin2,concoct', str,
            '\b\n05_metawrap_binning_depth: Comma-separated list of binning tools, choose 2 or 3 from [metabat2,maxbin2,concoct].\n'
            'At least 2 tools recommended for bin_refinement.'),
        'binning_max_memory_gb': (
            ['Params_b_sample', '05_metawrap_binning_depth', 'max_memory_gb'], 512, int,
            '05_metawrap_binning_depth: Maximum memory (GB) for MetaWRAP binning.'),
        'binning_metawrap_extra_args': (
            ['Params_b_sample', '05_metawrap_binning_depth', 'extra_args'], '', str,
            '05_metawrap_binning_depth: Extra arguments forwarded to metaWRAP binning (e.g., "-m 5").'),
        'metawrap_refinement_max_memory_gb': (
            ['Params_b_sample', '06_metawrap_bin_refinement', 'max_memory_gb'], 512, int,
            '06_metawrap_bin_refinement: Maximum memory in GB for MetaWRAP bin_refinement.'),
        'metawrap_refinement_extra_args': (
            ['Params_b_sample', '06_metawrap_bin_refinement', 'extra_args'], '', str,
            '06_metawrap_bin_refinement: Extra arguments forwarded to metaWRAP bin_refinement (e.g., "--remove-ambiguous").'),
        'bin_checkm_threads': (
            ['Params_b_sample', '06_bin_checkm', 'threads'], 8, int,
            '06_bin_checkm or 06_metawrap_bin_refinement: Threads.'),
        'bin_min_completeness': (
            ['Params_b_sample', '06_bin_checkm', 'min_completeness'], 50, float,
            '06_bin_checkm or 06_metawrap_bin_refinement: Minimum completeness (%) bin filtering, e.g. 50 for 50%.'),
        'bin_max_contamination': (
            ['Params_b_sample', '06_bin_checkm', 'max_contamination'], 10, float,
            '06_bin_checkm or 06_metawrap_bin_refinement: Maximum contamination (%) bin filtering, e.g. 10 for 10%.'),
        'bin_min_quality_score': (
            ['Params_b_sample', '06_bin_checkm', 'min_quality_score'], 50, float,
            '06_bin_checkm: Minimum quality score for bin filtering (GTDB-like), QS = completeness - 5*contamination, e.g. 50 for QS>=50.'),
        'contigs_depth_sample_threads': (
            ['Params_b_sample', '07_contigs_depth', 'threads'], 8, int,
            '07_contigs_depth: Threads.'),
        'depth_min_read_percent_identity': (
            ['Params_b_sample', '07_contigs_depth', 'min_read_percent_identity'], 95, float,
            '\b\n07_contigs_depth: Exclude reads by overall percent identity [0-100] for contig depth estimation\n'
            'using jgi_summarize_bam_contig_depths --percentIdentity, e.g. 95 for 95%.\n'
            'Reads with percent identity below this threshold are excluded from depth calculation.\n'
            'Higher values increase stringency, retaining only high-quality alignments.'),
        'depth_contig_end_exclusion': (
            ['Params_b_sample', '07_contigs_depth', 'contig_end_exclusion'], 75, int,
            '\b\n07_contigs_depth: Exclude bases at the ends of contig sequences for depth estimation using jgi_summarize_bam_contig_depths --maxEdgeBases.\n'
            'Bases within this distance from either end of a contig are ignored when calculating depth,\n'
            'preventing edge effects from partial alignments.'),
        'depth_min_contig_length': (
            ['Params_b_sample', '07_contigs_depth', 'min_contig_length'], 300, int,
            '\b\n07_contigs_depth: Minimum contig length for depth estimation [bp], e.g. 300 for 300 bp.\n'
            'Contigs shorter than this threshold are excluded from depth calculation.'),
        'contigs_prodigal_extra_args': (
            ['Params_b_sample', '07_contigs_depth', 'extra_args'], '', str,
            '07_contigs_depth: Extra arguments forwarded to prodigal for contig gene prediction (e.g., "-m -n").'),
        'bin_prodigal_threads': (
            ['Params_b_sample', '08_bin_prodigal', 'threads'], 8, int,
            '08_bin_prodigal: Threads.'),
        'bin_prodigal_extra_args': (
            ['Params_b_sample', '08_bin_prodigal', 'extra_args'], '', str,
            '08_bin_prodigal: Extra arguments forwarded to prodigal for MAG gene prediction (e.g., "-m -n").'),
        'raw_stats_project_threads': (
            ['Params_b_project', '01_cat_raw_stats', 'threads'], 4, int,
            '01_cat_raw_stats: Threads.'),
        'trimmed_stats_project_threads': (
            ['Params_b_project', '02_cat_trimmed_stats', 'threads'], 4, int,
            '02_cat_trimmed_stats: Threads.'),
        'contigs_stats_project_threads': (
            ['Params_b_project', '03_contigs_stats', 'threads'], 4, int,
            '03_contigs_stats: Threads.'),
        'merge_stats_project_threads': (
            ['Params_b_project', '04_merge_reads_contigs_stats', 'threads'], 4, int,
            '04_merge_reads_contigs_stats: Threads.'),
        'gather_checkm_project_threads': (
            ['Params_b_project', '05_gather_project_bin_checkm_stats', 'threads'], 4, int,
            '05_gather_project_bin_checkm_stats: Threads.'),
        'mmags_checkm2_threads': (
            ['Params_m_sample', '01_checkm2', 'threads'], 60, int,
            '01_checkm2: Threads.'),
        'mmags_checkm2_extra_args': (
            ['Params_m_sample', '01_checkm2', 'extra_args'], '', str,
            '01_checkm2: Extra arguments forwarded to checkm2 predict (e.g., "--remove_intermediates").'),
        'mmags_taxonomy_classify_method': (
            ['Params_m_sample', '02_taxonomy', 'MAGs_taxonomy_classify_method'], 'gtdbtk', str,
            '02_taxonomy: mMAGs taxonomy classifier, should be one of the [gtdbtk/BAT].'
        ),
        'mmags_taxonomy_threads': (
            ['Params_m_sample', '02_taxonomy', 'threads'], 60, int,
            '02_taxonomy_gtdbtk: Threads for [gtdbtk/BAT].'),
        'mmags_taxonomy_gtdbtk_extra_args': (
            ['Params_m_sample', '02_taxonomy_gtdbtk', 'extra_args'], '', str,
            '02_taxonomy_gtdbtk: Extra arguments forwarded to gtdbtk classify_wf (e.g., "--min_af 0.5").'),
        'mmags_taxonomy_bat_range': (
            ['Params_m_sample', '02_taxonomy_BAT', 'range'], 10, float,
            '\b\n02_taxonomy_BAT: Percentage range [0-100] of DIAMOND hits relative to top bit-score to include in LCA calculation.\n'
            'Higher values yield more conservative classifications by incorporating divergent hits, while lower values enable specific assignments.\n'
            'Applicable when --mMAGs-taxonomy-classify-method is set to BAT.'
        ),
        'mmags_taxonomy_bat_fraction': (
            ['Params_m_sample', '02_taxonomy_BAT', 'fraction'], 0.50, float,
            '\b\n02_taxonomy_BAT: Minimum bit-score support fraction [0-0.99].\n'
            'Higher values enforce stricter consensus (one classification per sequence), while lower values allow multiple lineage detections.\n'
            'Applicable when --mMAGs-taxonomy-classify-method is set to BAT.'
        ),
        'mmags_taxonomy_bat_extra_args': (
            ['Params_m_sample', '02_taxonomy_BAT', 'extra_args'], '', str,
            '02_taxonomy_BAT: Extra arguments forwarded to CAT_pack bins (e.g., "--block_size 12.0").'),
        'mmags_drep99_completeness_per_sample': (
            ['Params_m_sample', '05_drep_eachsample_99', 'comp'], 50, float,
            '05_drep_eachsample_99: The first level dereplication of mMAGs in a sample using completeness threshold [0-100], e.g., 50 for 50%.'
        ),
        'mmags_drep99_contamination_per_sample': (
            ['Params_m_sample', '05_drep_eachsample_99', 'con'], 10, float,
            '05_drep_eachsample_99: The first level dereplication of mMAGs in a sample using contamination threshold [0-100], e.g., 10 for 10%.'
        ),
        'mmags_drep99_ani_per_sample': (
            ['Params_m_sample', '05_drep_eachsample_99', 'sa'], 0.99, float,
            '05_drep_eachsample_99: The first level dereplication of mMAGs in a sample using ANI threshold [0-1], e.g., 0.99 for 99%.'
        ),
        'mmags_drep99_sample_threads': (
            ['Params_m_sample', '05_drep_eachsample_99', 'threads'], 20, int,
            '05_drep_eachsample_99: Threads.'),
        'mmags_drep99_sample_extra_args': (
            ['Params_m_sample', '05_drep_eachsample_99', 'extra_args'], '', str,
            '05_drep_eachsample: Extra arguments forwarded to dRep dereplicate (e.g., "--completeness_weight 1 --contamination_weight 5").'),
        'mmags_drep95_completeness_per_sample': (
            ['Params_m_sample', '06_drep_eachsample_95', 'comp'], 50, float,
            '06_drep_eachsample_95: The second level dereplication of mMAGs in a sample using completeness threshold [0-100], e.g., 50 for 50%.'
        ),
        'mmags_drep95_contamination_per_sample': (
            ['Params_m_sample', '06_drep_eachsample_95', 'con'], 10, float,
            '06_drep_eachsample_95: The second level dereplication of mMAGs in a sample using contamination threshold [0-100], e.g., 10 for 10%.'
        ),
        'mmags_drep95_ani_per_sample': (
            ['Params_m_sample', '06_drep_eachsample_95', 'sa'], 0.95, float,
            '\b\n06_drep_eachsample_95: The second level dereplication of mMAGs in a sample using ANI threshold [0-1], e.g., 0.95 for 95%.\n'
            'This value should be smaller than --drep99-ANI-per-sample'
        ),
        'mmags_drep95_sample_threads': (
            ['Params_m_sample', '06_drep_eachsample_95', 'threads'], 20, int,
            '06_drep_eachsample_95: Threads.'),
        'mmags_drep95_sample_extra_args': (
            ['Params_m_sample', '06_drep_eachsample_95', 'extra_args'], '', str,
            '06_drep_eachsample: Extra arguments forwarded to dRep dereplicate (e.g., "--completeness_weight 1 --contamination_weight 5").'),
        'mmags_coverm_method': (
            ['Params_m_sample', '07_coverm_eachsample', 'method'], 'trimmed_mean', str,
            '\b\n07_coverm_eachsample: Method(s) for coverm-based coverage calculation, must be one of [trimmed_mean/mean].\n'
            'A more detailed description of the different methods is available at https://github.com/wwood/CoverM#calculation-methods'
        ),
        'mmags_coverm_min_read_percent_identity': (
            ['Params_m_sample', '07_coverm_eachsample', 'min_read_percent_identity'], 95, float,
            '07_coverm_eachsample: Exclude reads by overall percent identity [0-100], e.g. 95 for 95%.'
        ),
        'mmags_coverm_min_read_aligned_length': (
            ['Params_m_sample', '07_coverm_eachsample', 'min_read_aligned_length'], 45, int,
            '07_coverm_eachsample: Exclude reads with fewer aligned bases.'
        ),
        'mmags_coverm_contig_end_exclusion': (
            ['Params_m_sample', '07_coverm_eachsample', 'contig_end_exclusion'], 75, int,
            '07_coverm_eachsample: Exclude bases at the ends of contig sequences for coverm calculation.'
        ),
        'mmags_coverm_threads': (
            ['Params_m_sample', '07_coverm_eachsample', 'threads'], 20, int,
            '07_coverm_eachsample: Threads.'),
        'mmags_coverm_sample_extra_args': (
            ['Params_m_sample', '07_coverm_eachsample', 'extra_args'], '', str,
            '07_coverm_eachsample: Extra arguments forwarded to coverm genome (e.g., "--quiet").'),
        'mmags_function_classify_method': (
            ['Params_m_sample', '08_function', 'MAGs_function_classify_method'], 'emapper', str,
            '08_function: mMAGs function classifier using one of [kraken2/emapper].'
        ),
        'mmags_function_threads': (
            ['Params_m_sample', '08_function', 'threads'], 60, int,
            '08_function: Threads for [kraken2/emapper].'),
        'mmags_function_kraken2_confidence': (
            ['Params_m_sample', '08_function_kraken2', 'confidence'], 0.05, float,
            '08_function_kraken2: Confidence threshold [0-1] for functional classification for mMAGs when --mMAGs-function-classify-method is set to kraken2.'
        ),
        'mmags_function_kraken2_extra_args': (
            ['Params_m_sample', '08_function_kraken2', 'extra_args'], '', str,
            '08_function_kraken2: Extra arguments forwarded to kraken2 (e.g., "--memory-mapping").'),
        'mmags_function_emapper_evalue': (
            ['Params_m_sample', '08_function_emapper', 'seed_ortholog_evalue'], 1e-10, float,
            '08_function_emapper: Report only alignments with e-value below or equal to the threshold, applicable when --mMAGs-function-classify-method is set to emapper.'
        ),
        'mmags_function_emapper_method': (
            ['Params_m_sample', '08_function_emapper', 'method'], 'diamond', str,
            '\b\n08_function_emapper: search seed orthologs using one of [diamond/mmseqs], applicable when --mMAGs-function-classify-method is set to emapper.\n'
            'diamond: using default parameters of diamond mode of emapper.\n'
            'mmseqs: using default parameters of mmseqs mode of emapper'
        ),
        'mmags_function_emapper_query_cover': (
            ['Params_m_sample', '08_function_emapper', 'query_cover'], 30, float,
            '\b\n08_function_emapper: Report only alignments with query cover above or equal to the given percentage [0-100],.\n'
            'applicable when --mMAGs-function-classify-method is set to emapper'
        ),
        'mmags_function_emapper_extra_args': (
            ['Params_m_sample', '08_function_emapper', 'extra_args'], '', str,
            '08_function_emapper: Extra arguments forwarded to eggNOG-mapper (e.g., "--pident 30").'),
        'contigs_taxonomy_classify_method': (
            ['Params_c_sample', '01_taxonomy', 'contigs_taxonomy_classify_method'], 'kraken2', str,
            '01_taxonomy: Contigs taxonomy classifier using one of [CAT_pack/kraken2].'
        ),
        'contigs_taxonomy_threads': (
            ['Params_c_sample', '01_taxonomy', 'threads'], 60, int,
            '01_taxonomy: Threads for [CAT_pack/kraken2].'
        ),
        'contigs_taxonomy_cat_range': (
            ['Params_c_sample', '01_taxonomy_CAT_pack', 'range'], 10, float,
            '\b\n01_taxonomy_CAT_pack: Percentage range [0-100] of DIAMOND hits relative to top bit-score to include in LCA calculation.\n'
            'Higher values yield more conservative classifications by incorporating divergent hits, while lower values enable specific assignments.\n'
            'Applicable when --contigs-taxonomy-classify-method is set to CAT_pack.'
        ),
        'contigs_taxonomy_cat_fraction': (
            ['Params_c_sample', '01_taxonomy_CAT_pack', 'fraction'], 0.50, float,
            '\b\n01_taxonomy_CAT_pack: Minimum bit-score support fraction [0-0.99].\n'
            'Higher values enforce stricter consensus (one classification per sequence), while lower values allow multiple lineage detections.\n'
            'Applicable when --contigs-taxonomy-classify-method is set to CAT_pack.'
        ),
        'contigs_taxonomy_cat_pack_extra_args': (
            ['Params_c_sample', '01_taxonomy_CAT_pack', 'extra_args'], '', str,
            '01_taxonomy_CAT_pack: Extra arguments forwarded to CAT_pack contigs (e.g., "--block_size 12.0").'),
        'contigs_taxonomy_kraken2_confidence': (
            ['Params_c_sample', '01_taxonomy_kraken2', 'confidence'], 0.7, float,
            '\b\n01_taxonomy_kraken2: Confidence threshold [0-1] for taxonomic classification, applicable when --contigs-taxonomy-classify-method is set to kraken2.\n'
            'Default value was tuned on the CAMI2 marine dataset. For samples with higher complexity, a lower confidence threshold is usually recommended.'
        ),
        'contigs_taxonomy_kraken2_extra_args': (
            ['Params_c_sample', '01_taxonomy_kraken2', 'extra_args'], '', str,
            '01_taxonomy_kraken2: Extra arguments forwarded to kraken2 (e.g., "--memory-mapping").'),
        'contigs_function_classify_method': (
            ['Params_c_sample', '02_function', 'contigs_function_classify_method'], 'kraken2', str,
            '02_function: Contigs function classifier using one of [kraken2/emapper].'
        ),
        'contigs_function_threads': (
            ['Params_c_sample', '02_function', 'threads'], 60, int,
            '02_function: Threads for [kraken2/emapper].'),
        'contigs_function_emapper_evalue': (
            ['Params_c_sample', '02_function_emapper', 'seed_ortholog_evalue'], 1e-10, float,
            '02_function_emapper: Report only alignments with e-value below or equal to the threshold, applicable when --contigs-function-classify-method is set to emapper.'
        ),
        'contigs_function_emapper_method': (
            ['Params_c_sample', '02_function_emapper', 'method'], 'diamond', str,
            '\b\n02_function_emapper: search seed orthologs using one of [diamond/mmseqs], applicable when --contigs-function-classify-method is set to emapper.\n'
            'diamond: using default parameters of diamond mode of emapper.\n'
            'mmseqs: using default parameters of mmseqs mode of emapper'
        ),
        'contigs_function_emapper_query_cover': (
            ['Params_c_sample', '02_function_emapper', 'query_cover'], 30, float,
            '\b\n02_function_emapper: Report only alignments with query cover above or equal to the given percentage [0-100],.\n'
            'applicable when --contigs-function-classify-method is set to emapper'
        ),
        'contigs_function_emapper_extra_args': (
            ['Params_c_sample', '02_function_emapper', 'extra_args'], '', str,
            '02_function_emapper: Extra arguments forwarded to eggNOG-mapper (e.g., "--pident 30").'),
        'contigs_function_kraken2_confidence': (
            ['Params_c_sample', '02_function_kraken2', 'confidence'], 0.05, float,
            '02_function_kraken2: Confidence threshold [0-1] for functional classification, applicable when --contigs-function-classify-method is set to kraken2.'
        ),
        'contigs_function_kraken2_extra_args': (
            ['Params_c_sample', '02_function_kraken2', 'extra_args'], '', str,
            '02_function_kraken2: Extra arguments forwarded to kraken2 (e.g., "--memory-mapping").'),
        'merge_mmags_min_completeness1': (
            ['Params_mc_sample', '05_merge_tax_fun_depth', 'MAGs_min_completeness1'], 50, float,
            '02_quality_gate: Minimum completeness threshold (%) for CheckM1-based mMAG quality filtering [0-100], e.g., 50 for 50%.'
        ),
        'merge_mmags_max_contamination1': (
            ['Params_mc_sample', '05_merge_tax_fun_depth', 'MAGs_max_contamination1'], 1, float,
            '02_quality_gate: Maximum contamination threshold (%) for CheckM1-based mMAG quality filtering [0-100], e.g., 5 for 5%.'
        ),
        'merge_mmags_min_completeness2': (
            ['Params_mc_sample', '05_merge_tax_fun_depth', 'MAGs_min_completeness2'], 50, float,
            '02_quality_gate: Minimum completeness threshold (%) for CheckM2-based mMAG quality filtering [0-100], e.g., 50 for 50%.'
        ),
        'merge_mmags_max_contamination2': (
            ['Params_mc_sample', '05_merge_tax_fun_depth', 'MAGs_max_contamination2'], 1, float,
            '02_quality_gate: Maximum contamination threshold (%) for CheckM2-based mMAG quality filtering [0-100], e.g., 5 for 5%.'
        ),
        'merge_min_orf_length': (
            ['Params_mc_sample', '05_merge_tax_fun_depth', 'minORFLen'], 50, int,
            '05_merge_tax_fun_depth: Minimum ORF length for finally used ORFs of mMAGs and contigs.'
        ),
        'merge_multiple_ko_handling_mode': (
            ['Params_mc_sample', '05_merge_tax_fun_depth', 'KO_mode'], 'ignore_multiple_KOs', str,
            '\b\n05_merge_tax_fun_depth: Method using one of [split_multiple_KOs/keep_multiple_KOs/ignore_multiple_KOs] when an ORF has multiple KO annotations:.\n'
            'split_multiple_KOs: evenly distributing the total abundance across all assigned KOs;\n'
            'keep_multiple_KOs: retain original multiple KO annotations;\n'
            'ignore_multiple_KOs: ignore the ORFs with multiple KO annotations and only use the ORFs with single KO annotation.'
        ),
        'merge_min_ko_num_per_taxon': (
            ['Params_mc_sample', '05_merge_tax_fun_depth', 'min_KO_num_in_one_taxon'], 100, int,
            '\b\n05_merge_tax_fun_depth: Minimum number of KOs per taxon, used to filter outputs in KO_Taxa_RelativeAbundance and KO_Taxa_SumCoverage.\n'
            'This threshold filters out taxa with insufficient KOs,\n'
            'retaining only KO-rich taxa to ensure robust metabolic reconstruction\n'
            'while preventing excessive output file sizes caused by sparse matrices.'
        ),
        'merge_tax_fun_depth_threads': (
            ['Params_mc_sample', '05_merge_tax_fun_depth', 'threads'], 20, int,
            '05_merge_tax_fun_depth: Threads.'),
        'merge_cami_format_output': (
            ['Params_mc_sample', '05_merge_tax_fun_depth', 'CAMI_format_output'], False, bool,
            '05_merge_tax_fun_depth: Report a taxonomic profile using the CAMI output format (flag only). Disabled when mMAG mode uses gtdbtk.'
        ),
        'merge_output_detail': (
            ['Params_mc_sample', '05_merge_tax_fun_depth', 'output_detail'], False, bool,
            '05_merge_tax_fun_depth: Whether to output detail tables (flag only).'),
        'mclink_sample_min_ko_num_per_taxon': (
            ['Params_mc_sample', '06_mclink_taxon', 'min_KO_num_in_one_taxon'], 100, int,
            '06_mclink_taxon: Minimum number of KOs per taxon for calculating taxon-level metabolic pathway/module completeness or abundance per sample. Taxa with fewer than this number of detected KOs in a given sample will be filtered out. Calculated independently for each sample.'
        ),
        'mclink_sample_min_ko_abundance_per_taxon': (
            ['Params_mc_sample', '06_mclink_taxon', 'min_mean_KO_RA_in_one_taxon'], 0, float,
            '06_mclink_taxon: Minimum mean KO abundance [0-1e6] per taxon for calculating taxon-level metabolic pathway/module completeness or abundance per sample. Taxa with mean KO abundance below this threshold in a given sample will be filtered out. Note: GPM/TPM sum of all taxa in each sample is 1e6.'
        ),
        'mclink_sample_topn_taxon': (
            ['Params_mc_sample', '06_mclink_taxon', 'TopN'], 100, int,
            '06_mclink_taxon: Number of top taxa (ranked by mean KO abundance) retained for calculating taxon-level metabolic pathway/module completeness or abundance per sample. Set to 0 to disable filtering.'
        ),
        'mclink_sample_completeness_per_taxon': (
            ['Params_mc_sample', '06_mclink_taxon', 'completeness_threshold'], 0.5, float,
            '06_mclink_taxon: Completeness threshold [0-1] of modules for calculating taxon-level metabolic pathway/module abundance. Higher values retain only more complete profiles.'),
        'mclink_sample_plus_scale_method_per_taxon': (
            ['Params_mc_sample', '06_mclink_taxon', 'plus_scale_method'], 'mean', str,
            "06_mclink_taxon: Aggregation method for '+'-separated KOs in a module [mean/max/min]. 'mean': average abundance across all required KOs; 'max': highest KO abundance; 'min': lowest KO abundance."),
        'mclink_sample_comma_scale_method_per_taxon': (
            ['Params_mc_sample', '06_mclink_taxon', 'comma_scale_method'], 'sum', str,
            "06_mclink_taxon: Aggregation method for ','-separated alternative KOs in a module [sum/max]. 'sum': cumulative abundance of all alternative KOs; 'max': abundance of the highest single KO among alternatives. Note: this applies to abundance calculation only; module completeness is always computed with 'max'."),
        'barplot_topn_taxon': (
            ['Params_mc_project', '01_Taxa_visualize', 'TopN'], 10, int,
            '01_Taxa_visualize: Number of top taxa to display in stacked barplots. Remaining categories are merged into Others. Maximum 29.'
        ),
        'barplot_topn_fun': (
            ['Params_mc_project', '03_Fun_visualize', 'TopN'], 10, int,
            '03_Fun_visualize: Number of top functional categories to display in stacked barplots. Remaining categories are merged into Others. Maximum 29.'
        ),
        'mclink_min_ko_num_per_sample': (
            ['Params_mc_project', '03_mclink_sample', 'min_sample_num_per_KO'], 1, int,
            '03_mclink_sample: Minimum number of samples in which a KO must be detected for calculating community-level metabolic pathway/module completeness or abundance.'
        ),
        'mclink_min_ko_prevalence_per_sample': (
            ['Params_mc_project', '03_mclink_sample', 'preva'], 0, float,
            '03_mclink_sample: Minimum detection rate [0-1] of a KO across all samples for calculating community-level metabolic pathway/module completeness or abundance.'
        ),
        'mclink_min_ko_abundance_per_sample': (
            ['Params_mc_project', '03_mclink_sample', 'abund'], 0, float,
            '03_mclink_sample: Minimum mean relative abundance [0-1e6] of a KO across all samples for calculating community-level metabolic pathway/module completeness or abundance. KOs with mean abundance below this threshold are excluded. Note: GPM/TPM sum of a sample is 1e6.'
        ),
        'mclink_completeness_per_sample': (
            ['Params_mc_project', '03_mclink_sample', 'completeness_threshold'], 0.5, float,
            '03_mclink_sample: Completeness threshold [0-1] of modules for calculating community-level metabolic pathway/module abundance. Higher values retain only more complete profiles.'
        ),
        'mclink_plus_scale_method_per_sample': (
            ['Params_mc_project', '03_mclink_sample', 'plus_scale_method'], 'mean', str,
            "03_mclink_sample: Aggregation method for '+'-separated KOs in a module [mean/max/min]. 'mean': average abundance across all required KOs; 'max': highest KO abundance; 'min': lowest KO abundance."
        ),
        'mclink_comma_scale_method_per_sample': (
            ['Params_mc_project', '03_mclink_sample', 'comma_scale_method'], 'sum', str,
            "03_mclink_sample: Aggregation method for ','-separated alternative KOs in a module [sum/max]. 'sum': cumulative abundance of all alternative KOs; 'max': abundance of the highest single KO among alternatives. Note: this applies to abundance calculation only; module completeness is always computed with 'max'."
        ),
        'merge_min_taxa_ko_prevalence_across_samples': (
            ['Params_mc_project', '05_merge_tax_KO_depth', 'preva'], 0.10, float,
            '05_merge_tax_KO_depth: Minimum prevalence [0-1] of a Taxon-KO across samples. Taxa-KOs detected in fewer samples than this proportion are excluded.'
        ),
        'merge_min_taxa_ko_abundance_across_samples': (
            ['Params_mc_project', '05_merge_tax_KO_depth', 'abund'], 0.0, float,
            '05_merge_tax_KO_depth: Minimum mean relative abundance [0-1e6] of a Taxon-KO across samples. Taxa-KOs with mean abundance below this threshold are excluded. Note: GPM/TPM sum of a sample is 1e6.'
        ),
        'merge_min_taxa_module_prevalence_across_samples': (
            ['Params_mc_project', '05_merge_tax_Module_RelativeAbundance', 'preva'], 0.10, float,
            '05_merge_tax_Module_RelativeAbundance: Minimum prevalence [0-1] of a Taxon-Module across samples. Taxa-Modules detected in fewer samples than this proportion are excluded.'
        ),
        'merge_min_taxa_module_abundance_across_samples': (
            ['Params_mc_project', '05_merge_tax_Module_RelativeAbundance', 'abund'], 0.0, float,
            '05_merge_tax_Module_RelativeAbundance: Minimum mean relative abundance [0-1e6] of a Taxon-Module across samples. Taxa-Modules with mean abundance below this threshold are excluded. Note: GPM/TPM sum of a sample is 1e6.'
        ),
        'merge_min_taxa_module_sc_prevalence_across_samples': (
            ['Params_mc_project', '05_merge_tax_Module_SumCoverage', 'preva'], 0.10, float,
            '05_merge_tax_Module_SumCoverage: Minimum prevalence [0-1] of a Taxon-Module across samples. Taxa-Modules detected in fewer samples than this proportion are excluded.'
        ),
        'merge_min_taxa_module_sc_abundance_across_samples': (
            ['Params_mc_project', '05_merge_tax_Module_SumCoverage', 'abund'], 0.0, float,
            '05_merge_tax_Module_SumCoverage: Minimum mean sum coverage [0-1e6] of a Taxon-Module across samples. Taxa-Modules with mean sum coverage below this threshold are excluded. Note: GPM/TPM sum of a sample is 1e6.'
        ),
        'mclink_project_min_ko_num_per_taxon': (
            ['Params_mc_project', '06_mclink_taxon', 'min_KO_num_in_one_taxon'], 100, int,
            '06_mclink_taxon: Minimum number of KOs per taxon for calculating taxon-level metabolic pathway/module completeness on the project-averaged taxon-function matrix. Taxa with fewer than this number of detected KOs in the project-averaged matrix will be filtered out. KO abundances are first averaged across all samples to generate a single project-level matrix before calculation.'
        ),
        'mclink_project_min_ko_abundance_per_taxon': (
            ['Params_mc_project', '06_mclink_taxon', 'min_mean_KO_RA_in_one_taxon'], 0, float,
            '06_mclink_taxon: Minimum mean KO abundance [0-1e6] per taxon for calculating taxon-level metabolic pathway/module completeness on the project-averaged taxon-function matrix. Taxa with mean KO abundance below this threshold in the project-averaged matrix will be filtered out. Note: GPM/TPM sum of all taxa is 1e6 per sample; values are averaged across samples prior to filtering and calculation.'
        ),
        'mclink_project_topn_taxon': (
            ['Params_mc_project', '06_mclink_taxon', 'TopN'], 100, int,
            '06_mclink_taxon: Number of top taxa (ranked by mean KO abundance in the project-averaged taxon-function matrix) retained for calculating taxon-level metabolic pathway/module completeness. Set to 0 to disable filtering.'
        ),
        'mclink_project_completeness_per_taxon': (
            ['Params_mc_project', '06_mclink_taxon', 'completeness_threshold'], 0.5, float,
            '06_mclink_taxon: Completeness threshold [0-1] of modules for calculating taxon-level metabolic pathway/module abundance on the project-averaged taxon-function matrix. Higher values retain only more complete profiles.'),
        'mclink_project_plus_scale_method_per_taxon': (
            ['Params_mc_project', '06_mclink_taxon', 'plus_scale_method'], 'mean', str,
            "06_mclink_taxon: Aggregation method for '+'-separated KOs in a module [mean/max/min] on the project-averaged taxon-function matrix. 'mean': average abundance across all required KOs; 'max': highest KO abundance; 'min': lowest KO abundance."),
        'mclink_project_comma_scale_method_per_taxon': (
            ['Params_mc_project', '06_mclink_taxon', 'comma_scale_method'], 'sum', str,
            "06_mclink_taxon: Aggregation method for ','-separated alternative KOs in a module [sum/max] on the project-averaged taxon-function matrix. 'sum': cumulative abundance of all alternative KOs; 'max': abundance of the highest single KO among alternatives. Note: this applies to abundance calculation only; module completeness is always computed with 'max'."),
    },
    'mcr': {
        'raw_stats_sample_threads': (
            ['Params_b_sample', '00_raw_stats', 'threads'], 4, int,
            '00_raw_stats: Threads for seqkit.'),
        'trim_sample_threads': (
            ['Params_b_sample', '01_trim', 'threads'], 16, int,
            '01_trim: Threads for fastp (next-generation short reads) or fastplong (third-generation long reads).'),
        'trim_length_required': (
            ['Params_b_sample', '01_trim', 'length_required'], 15, int,
            '01_trim: Reads shorter than length_required will be discarded'),
        'trim_n_base_limit': (
            ['Params_b_sample', '01_trim', 'n_base_limit'], 5, int,
            "01_trim: if one read's number of N base is >n_base_limit, then this read/pair is discarded."),
        'trim_compression': (
            ['Params_b_sample', '01_trim', 'compression'], 4, int,
            '01_trim: Compression level for gzip output [1-9]. 1 is fastest, 9 is smallest.'),
        'trim_extra_args': (
            ['Params_b_sample', '01_trim', 'extra_args'], '', str,
            '01_trim: Extra arguments forwarded to fastp or fastplong for raw reads trimming (e.g., "--reads_to_process 0 -Q").'),
        'trim_stats_sample_threads': (
            ['Params_b_sample', '01_trim_stats', 'threads'], 16, int,
            '01_trim_stats: Threads.'),
        'assembly_method': (
            ['Params_b_sample', '02_assembly', 'method'], 'megahit', str,
            '\b\n02_assembly: Underlying assembly software using one of [megahit/spades/metaflye].\n'
            'Must be compatible with --read-type (config key read_type):\n'
            'megahit: MEGAHIT (second-generation sequencing; for paired or single DNA reads).\n'
            'spades: SPAdes (second-generation sequencing; metaSPAdes for paired DNA reads, rnaSPAdes for paired-rna or single-rna RNA reads).\n'
            'metaflye: metaFlye (third-generation sequencing; required for pacbio-raw/corr/hifi or nano-raw/corr/hq long reads).'),
        'assembly_threads': (
            ['Params_b_sample', '02_assembly', 'threads'], 40, int,
            '02_assembly: Threads for [megahit/spades/metaflye].'),
        'megahit_presets': (
            ['Params_b_sample', '02_megahit', 'presets'], 'meta-sensitive', str,
            '\b\n02_megahit: MEGAHIT presets using one of [meta-sensitive/meta-large].\n'
            "meta-sensitive: '--min-count 1 --k-list 21,29,39,49,...,129,141'.\n"
            "meta-large: '--k-min 27 --k-max 127 --k-step 10'.\n"
            '(large & complex metagenomes, like soil).'),
        'megahit_k_list': (
            ['Params_b_sample', '02_megahit', 'k_list'], '', str,
            '02_megahit: MEGAHIT --k-list (optional). Explicit comma-separated k-mer sizes (e.g. "21,29,39,59,79,99,119,141"). Overrides --megahit-presets. If set, --megahit-k-min, --megahit-k-max, --megahit-k-step must NOT be set. Optional --megahit-min-count can be used together.'),
        'megahit_min_count': (
            ['Params_b_sample', '02_megahit', 'min_count'], 2, int,
            '02_megahit: MEGAHIT --min-count (optional). Minimum multiplicity for filtering (k_min+1)-mers. Can be used with --megahit-k-list or --megahit-k-min/--k-max/--k-step.'),
        'megahit_k_min': (
            ['Params_b_sample', '02_megahit', 'k_min'], '', str,
            '02_megahit: MEGAHIT --k-min (optional). Minimum kmer size (<= 255), must be odd number. If set with --k-max and --k-step, overrides --megahit-presets and uses custom k-mer mode.'),
        'megahit_k_max': (
            ['Params_b_sample', '02_megahit', 'k_max'], '', str,
            '02_megahit: MEGAHIT --k-max (optional). Maximum kmer size (<= 255), must be odd number.'),
        'megahit_k_step': (
            ['Params_b_sample', '02_megahit', 'k_step'], '', str,
            '02_megahit: MEGAHIT --k-step (optional). Increment of kmer size of each iteration (<= 28), must be even number.'),
        'megahit_extra_args': (
            ['Params_b_sample', '02_megahit', 'extra_args'], '', str,
            '02_megahit: Extra arguments forwarded to MEGAHIT (e.g. "--no-mercy --tmp-dir /path/to/your/tmp/dir"). Do NOT use this for k-mer settings (--k-list, --k-min, --presets, etc.).'),
        'spades_extra_args': (
            ['Params_b_sample', '02_spades', 'extra_args'], '', str,
            '02_spades: Extra arguments forwarded to SPAdes. Do NOT use this for input/output or thread settings.'),
        'metaflye_extra_args': (
            ['Params_b_sample', '02_metaflye', 'extra_args'], '', str,
            '02_metaflye: Extra arguments forwarded to Flye for third-generation long reads assembly (e.g., "--iterations 2"). Do NOT use this for read type or input/output settings.'),
        'mapping_threads': (
            ['Params_b_sample', '03_mapping', 'threads'], 20, int,
            '03_mapping: Threads for [minibwa/bwamem/bwamem2/strobealign].'),
        'mapping_method': (
            ['Params_b_sample', '03_mapping', 'method'], 'minibwa', str,
            '\b\n03_mapping: Underlying mapping software used [minibwa/bwamem/bwamem2/strobealign].\n'
            'minibwa: minibwa using default parameters.\n'
            'bwamem: bwa mem using default parameters.\n'
            'bwamem2: bwa-mem2 using default parameters.\n'
            'strobealign: strobealign using default parameters.'),
        'mapping_samtools_sort_mem': (
            ['Params_b_sample', '03_mapping', 'samtools_sort_mem'], '768M', str,
            '03_mapping: Set maximum memory per thread; suffix K/M/G recognized, e.g. 2G, 768M.'),
        'mapping_extra_args': (
            ['Params_b_sample', '03_mapping', 'extra_args'], '', str,
            '03_mapping: Extra arguments forwarded to the active mapping aligner (minibwa map, bwa mem, bwa-mem2 mem, or strobealign --use-index), (e.g. "-B 8" for minibwa map).'),
        'binning_metawrap_depth_threads': (
            ['Params_b_sample', '05_metawrap_binning_depth', 'threads'], 8, int,
            '05_metawrap_binning_depth: Threads.'),
        'binning_min_read_percent_identity': (
            ['Params_b_sample', '05_metawrap_binning_depth', 'min_read_percent_identity'], 97, float,
            '05_metawrap_binning_depth: Exclude reads by overall percent identity [0-100] during the binning process, e.g. 97 for 97%.'),
        'binning_contig_end_exclusion': (
            ['Params_b_sample', '05_metawrap_binning_depth', 'contig_end_exclusion'], 75, int,
            '05_metawrap_binning_depth: Exclude bases at the ends of contig sequences used for depth calculation and bining.'),
        'binning_min_contig_length': (
            ['Params_b_sample', '05_metawrap_binning_depth', 'min_contig_length'], 1500, int,
            '05_metawrap_binning_depth: Minimum contig length for binning.'),
        'binning_tools': (
            ['Params_b_sample', '05_metawrap_binning_depth', 'binning_tools'], 'metabat2,maxbin2,concoct', str,
            '\b\n05_metawrap_binning_depth: Comma-separated list of binning tools, choose 2 or 3 from [metabat2,maxbin2,concoct].\n'
            'At least 2 tools recommended for bin_refinement.'),
        'binning_max_memory_gb': (
            ['Params_b_sample', '05_metawrap_binning_depth', 'max_memory_gb'], 512, int,
            '05_metawrap_binning_depth: Maximum memory (GB) for MetaWRAP binning.'),
        'binning_metawrap_extra_args': (
            ['Params_b_sample', '05_metawrap_binning_depth', 'extra_args'], '', str,
            '05_metawrap_binning_depth: Extra arguments forwarded to metaWRAP binning (e.g., "-m 5").'),
        'metawrap_refinement_max_memory_gb': (
            ['Params_b_sample', '06_metawrap_bin_refinement', 'max_memory_gb'], 512, int,
            '06_metawrap_bin_refinement: Maximum memory in GB for MetaWRAP bin_refinement.'),
        'metawrap_refinement_extra_args': (
            ['Params_b_sample', '06_metawrap_bin_refinement', 'extra_args'], '', str,
            '06_metawrap_bin_refinement: Extra arguments forwarded to metaWRAP bin_refinement (e.g., "--remove-ambiguous").'),
        'bin_checkm_threads': (
            ['Params_b_sample', '06_bin_checkm', 'threads'], 8, int,
            '06_bin_checkm or 06_metawrap_bin_refinement: Threads.'),
        'bin_min_completeness': (
            ['Params_b_sample', '06_bin_checkm', 'min_completeness'], 50, float,
            '06_bin_checkm or 06_metawrap_bin_refinement: Minimum completeness (%) bin filtering, e.g. 50 for 50%.'),
        'bin_max_contamination': (
            ['Params_b_sample', '06_bin_checkm', 'max_contamination'], 10, float,
            '06_bin_checkm or 06_metawrap_bin_refinement: Maximum contamination (%) bin filtering, e.g. 10 for 10%.'),
        'bin_min_quality_score': (
            ['Params_b_sample', '06_bin_checkm', 'min_quality_score'], 50, float,
            '06_bin_checkm: Minimum quality score for bin filtering (GTDB-like), QS = completeness - 5*contamination, e.g. 50 for QS>=50.'),
        'contigs_depth_sample_threads': (
            ['Params_b_sample', '07_contigs_depth', 'threads'], 8, int,
            '07_contigs_depth: Threads.'),
        'depth_min_read_percent_identity': (
            ['Params_b_sample', '07_contigs_depth', 'min_read_percent_identity'], 95, float,
            '\b\n07_contigs_depth: Exclude reads by overall percent identity [0-100] for contig depth estimation\n'
            'using jgi_summarize_bam_contig_depths --percentIdentity, e.g. 95 for 95%.\n'
            'Reads with percent identity below this threshold are excluded from depth calculation.\n'
            'Higher values increase stringency, retaining only high-quality alignments.'),
        'depth_contig_end_exclusion': (
            ['Params_b_sample', '07_contigs_depth', 'contig_end_exclusion'], 75, int,
            '\b\n07_contigs_depth: Exclude bases at the ends of contig sequences for depth estimation using jgi_summarize_bam_contig_depths --maxEdgeBases.\n'
            'Bases within this distance from either end of a contig are ignored when calculating depth,\n'
            'preventing edge effects from partial alignments.'),
        'depth_min_contig_length': (
            ['Params_b_sample', '07_contigs_depth', 'min_contig_length'], 300, int,
            '\b\n07_contigs_depth: Minimum contig length for depth estimation [bp], e.g. 300 for 300 bp.\n'
            'Contigs shorter than this threshold are excluded from depth calculation.'),
        'contigs_prodigal_extra_args': (
            ['Params_b_sample', '07_contigs_depth', 'extra_args'], '', str,
            '07_contigs_depth: Extra arguments forwarded to prodigal for contig gene prediction (e.g., "-m -n").'),
        'bin_prodigal_threads': (
            ['Params_b_sample', '08_bin_prodigal', 'threads'], 8, int,
            '08_bin_prodigal: Threads.'),
        'bin_prodigal_extra_args': (
            ['Params_b_sample', '08_bin_prodigal', 'extra_args'], '', str,
            '08_bin_prodigal: Extra arguments forwarded to prodigal for MAG gene prediction (e.g., "-m -n").'),
        'raw_stats_project_threads': (
            ['Params_b_project', '01_cat_raw_stats', 'threads'], 4, int,
            '01_cat_raw_stats: Threads.'),
        'trimmed_stats_project_threads': (
            ['Params_b_project', '02_cat_trimmed_stats', 'threads'], 4, int,
            '02_cat_trimmed_stats: Threads.'),
        'contigs_stats_project_threads': (
            ['Params_b_project', '03_contigs_stats', 'threads'], 4, int,
            '03_contigs_stats: Threads.'),
        'merge_stats_project_threads': (
            ['Params_b_project', '04_merge_reads_contigs_stats', 'threads'], 4, int,
            '04_merge_reads_contigs_stats: Threads.'),
        'gather_checkm_project_threads': (
            ['Params_b_project', '05_gather_project_bin_checkm_stats', 'threads'], 4, int,
            '05_gather_project_bin_checkm_stats: Threads.'),
        'mmags_checkm2_threads': (
            ['Params_m_sample', '01_checkm2', 'threads'], 60, int,
            '01_checkm2: Threads.'),
        'mmags_checkm2_extra_args': (
            ['Params_m_sample', '01_checkm2', 'extra_args'], '', str,
            '01_checkm2: Extra arguments forwarded to checkm2 predict (e.g., "--remove_intermediates").'),
        'mmags_taxonomy_classify_method': (
            ['Params_m_sample', '02_taxonomy', 'MAGs_taxonomy_classify_method'], 'gtdbtk', str,
            '02_taxonomy: mMAGs taxonomy classifier, should be one of the [gtdbtk/BAT].'
        ),
        'mmags_taxonomy_threads': (
            ['Params_m_sample', '02_taxonomy', 'threads'], 60, int,
            '02_taxonomy_gtdbtk: Threads for [gtdbtk/BAT].'),
        'mmags_taxonomy_gtdbtk_extra_args': (
            ['Params_m_sample', '02_taxonomy_gtdbtk', 'extra_args'], '', str,
            '02_taxonomy_gtdbtk: Extra arguments forwarded to gtdbtk classify_wf (e.g., "--min_af 0.5").'),
        'mmags_taxonomy_bat_range': (
            ['Params_m_sample', '02_taxonomy_BAT', 'range'], 10, float,
            '\b\n02_taxonomy_BAT: Percentage range [0-100] of DIAMOND hits relative to top bit-score to include in LCA calculation.\n'
            'Higher values yield more conservative classifications by incorporating divergent hits, while lower values enable specific assignments.\n'
            'Applicable when --mMAGs-taxonomy-classify-method is set to BAT.'
        ),
        'mmags_taxonomy_bat_fraction': (
            ['Params_m_sample', '02_taxonomy_BAT', 'fraction'], 0.50, float,
            '\b\n02_taxonomy_BAT: Minimum bit-score support fraction [0-0.99].\n'
            'Higher values enforce stricter consensus (one classification per sequence), while lower values allow multiple lineage detections.\n'
            'Applicable when --mMAGs-taxonomy-classify-method is set to BAT.'
        ),
        'mmags_taxonomy_bat_extra_args': (
            ['Params_m_sample', '02_taxonomy_BAT', 'extra_args'], '', str,
            '02_taxonomy_BAT: Extra arguments forwarded to CAT_pack bins (e.g., "--block_size 12.0").'),
        'mmags_drep99_completeness_per_sample': (
            ['Params_m_sample', '05_drep_eachsample_99', 'comp'], 50, float,
            '05_drep_eachsample_99: The first level dereplication of mMAGs in a sample using completeness threshold [0-100], e.g., 50 for 50%.'
        ),
        'mmags_drep99_contamination_per_sample': (
            ['Params_m_sample', '05_drep_eachsample_99', 'con'], 10, float,
            '05_drep_eachsample_99: The first level dereplication of mMAGs in a sample using contamination threshold [0-100], e.g., 10 for 10%.'
        ),
        'mmags_drep99_ani_per_sample': (
            ['Params_m_sample', '05_drep_eachsample_99', 'sa'], 0.99, float,
            '05_drep_eachsample_99: The first level dereplication of mMAGs in a sample using ANI threshold [0-1], e.g., 0.99 for 99%.'
        ),
        'mmags_drep99_sample_threads': (
            ['Params_m_sample', '05_drep_eachsample_99', 'threads'], 20, int,
            '05_drep_eachsample_99: Threads.'),
        'mmags_drep99_sample_extra_args': (
            ['Params_m_sample', '05_drep_eachsample_99', 'extra_args'], '', str,
            '05_drep_eachsample: Extra arguments forwarded to dRep dereplicate (e.g., "--completeness_weight 1 --contamination_weight 5").'),
        'mmags_drep95_completeness_per_sample': (
            ['Params_m_sample', '06_drep_eachsample_95', 'comp'], 50, float,
            '06_drep_eachsample_95: The second level dereplication of mMAGs in a sample using completeness threshold [0-100], e.g., 50 for 50%.'
        ),
        'mmags_drep95_contamination_per_sample': (
            ['Params_m_sample', '06_drep_eachsample_95', 'con'], 10, float,
            '06_drep_eachsample_95: The second level dereplication of mMAGs in a sample using contamination threshold [0-100], e.g., 10 for 10%.'
        ),
        'mmags_drep95_ani_per_sample': (
            ['Params_m_sample', '06_drep_eachsample_95', 'sa'], 0.95, float,
            '\b\n06_drep_eachsample_95: The second level dereplication of mMAGs in a sample using ANI threshold [0-1], e.g., 0.95 for 95%.\n'
            'This value should be smaller than --drep99-ANI-per-sample'
        ),
        'mmags_drep95_sample_threads': (
            ['Params_m_sample', '06_drep_eachsample_95', 'threads'], 20, int,
            '06_drep_eachsample_95: Threads.'),
        'mmags_drep95_sample_extra_args': (
            ['Params_m_sample', '06_drep_eachsample_95', 'extra_args'], '', str,
            '06_drep_eachsample: Extra arguments forwarded to dRep dereplicate (e.g., "--completeness_weight 1 --contamination_weight 5").'),
        'mmags_coverm_method': (
            ['Params_m_sample', '07_coverm_eachsample', 'method'], 'trimmed_mean', str,
            '\b\n07_coverm_eachsample: Method(s) for coverm-based coverage calculation, must be one of [trimmed_mean/mean].\n'
            'A more detailed description of the different methods is available at https://github.com/wwood/CoverM#calculation-methods'
        ),
        'mmags_coverm_min_read_percent_identity': (
            ['Params_m_sample', '07_coverm_eachsample', 'min_read_percent_identity'], 95, float,
            '07_coverm_eachsample: Exclude reads by overall percent identity [0-100], e.g. 95 for 95%.'
        ),
        'mmags_coverm_min_read_aligned_length': (
            ['Params_m_sample', '07_coverm_eachsample', 'min_read_aligned_length'], 45, int,
            '07_coverm_eachsample: Exclude reads with fewer aligned bases.'
        ),
        'mmags_coverm_contig_end_exclusion': (
            ['Params_m_sample', '07_coverm_eachsample', 'contig_end_exclusion'], 75, int,
            '07_coverm_eachsample: Exclude bases at the ends of contig sequences for coverm calculation.'
        ),
        'mmags_coverm_threads': (
            ['Params_m_sample', '07_coverm_eachsample', 'threads'], 20, int,
            '07_coverm_eachsample: Threads.'),
        'mmags_coverm_sample_extra_args': (
            ['Params_m_sample', '07_coverm_eachsample', 'extra_args'], '', str,
            '07_coverm_eachsample: Extra arguments forwarded to coverm genome (e.g., "--quiet").'),
        'mmags_function_classify_method': (
            ['Params_m_sample', '08_function', 'MAGs_function_classify_method'], 'emapper', str,
            '08_function: mMAGs function classifier using one of [kraken2/emapper].'
        ),
        'mmags_function_threads': (
            ['Params_m_sample', '08_function', 'threads'], 60, int,
            '08_function: Threads for [kraken2/emapper].'),
        'mmags_function_kraken2_confidence': (
            ['Params_m_sample', '08_function_kraken2', 'confidence'], 0.05, float,
            '08_function_kraken2: Confidence threshold [0-1] for functional classification for mMAGs when --mMAGs-function-classify-method is set to kraken2.'
        ),
        'mmags_function_kraken2_extra_args': (
            ['Params_m_sample', '08_function_kraken2', 'extra_args'], '', str,
            '08_function_kraken2: Extra arguments forwarded to kraken2 (e.g., "--memory-mapping").'),
        'mmags_function_emapper_evalue': (
            ['Params_m_sample', '08_function_emapper', 'seed_ortholog_evalue'], 1e-10, float,
            '08_function_emapper: Report only alignments with e-value below or equal to the threshold, applicable when --mMAGs-function-classify-method is set to emapper.'
        ),
        'mmags_function_emapper_method': (
            ['Params_m_sample', '08_function_emapper', 'method'], 'diamond', str,
            '\b\n08_function_emapper: search seed orthologs using one of [diamond/mmseqs], applicable when --mMAGs-function-classify-method is set to emapper.\n'
            'diamond: using default parameters of diamond mode of emapper.\n'
            'mmseqs: using default parameters of mmseqs mode of emapper'
        ),
        'mmags_function_emapper_query_cover': (
            ['Params_m_sample', '08_function_emapper', 'query_cover'], 30, float,
            '\b\n08_function_emapper: Report only alignments with query cover above or equal to the given percentage [0-100],.\n'
            'applicable when --mMAGs-function-classify-method is set to emapper'
        ),
        'mmags_function_emapper_extra_args': (
            ['Params_m_sample', '08_function_emapper', 'extra_args'], '', str,
            '08_function_emapper: Extra arguments forwarded to eggNOG-mapper (e.g., "--pident 30").'),
        'contigs_taxonomy_classify_method': (
            ['Params_c_sample', '01_taxonomy', 'contigs_taxonomy_classify_method'], 'kraken2', str,
            '01_taxonomy: Contigs taxonomy classifier using one of [CAT_pack/kraken2].'
        ),
        'contigs_taxonomy_threads': (
            ['Params_c_sample', '01_taxonomy', 'threads'], 60, int,
            '01_taxonomy: Threads for [CAT_pack/kraken2].'
        ),
        'contigs_taxonomy_cat_range': (
            ['Params_c_sample', '01_taxonomy_CAT_pack', 'range'], 10, float,
            '\b\n01_taxonomy_CAT_pack: Percentage range [0-100] of DIAMOND hits relative to top bit-score to include in LCA calculation.\n'
            'Higher values yield more conservative classifications by incorporating divergent hits, while lower values enable specific assignments.\n'
            'Applicable when --contigs-taxonomy-classify-method is set to CAT_pack.'
        ),
        'contigs_taxonomy_cat_fraction': (
            ['Params_c_sample', '01_taxonomy_CAT_pack', 'fraction'], 0.50, float,
            '\b\n01_taxonomy_CAT_pack: Minimum bit-score support fraction [0-0.99].\n'
            'Higher values enforce stricter consensus (one classification per sequence), while lower values allow multiple lineage detections.\n'
            'Applicable when --contigs-taxonomy-classify-method is set to CAT_pack.'
        ),
        'contigs_taxonomy_cat_pack_extra_args': (
            ['Params_c_sample', '01_taxonomy_CAT_pack', 'extra_args'], '', str,
            '01_taxonomy_CAT_pack: Extra arguments forwarded to CAT_pack contigs (e.g., "--block_size 12.0").'),
        'contigs_taxonomy_kraken2_confidence': (
            ['Params_c_sample', '01_taxonomy_kraken2', 'confidence'], 0.7, float,
            '\b\n01_taxonomy_kraken2: Confidence threshold [0-1] for taxonomic classification, applicable when --contigs-taxonomy-classify-method is set to kraken2.\n'
            'Default value was tuned on the CAMI2 marine dataset. For samples with higher complexity, a lower confidence threshold is usually recommended.'
        ),
        'contigs_taxonomy_kraken2_extra_args': (
            ['Params_c_sample', '01_taxonomy_kraken2', 'extra_args'], '', str,
            '01_taxonomy_kraken2: Extra arguments forwarded to kraken2 (e.g., "--memory-mapping").'),
        'contigs_function_classify_method': (
            ['Params_c_sample', '02_function', 'contigs_function_classify_method'], 'kraken2', str,
            '02_function: Contigs function classifier using one of [kraken2/emapper].'
        ),
        'contigs_function_threads': (
            ['Params_c_sample', '02_function', 'threads'], 60, int,
            '02_function: Threads for [kraken2/emapper].'),
        'contigs_function_emapper_evalue': (
            ['Params_c_sample', '02_function_emapper', 'seed_ortholog_evalue'], 1e-10, float,
            '02_function_emapper: Report only alignments with e-value below or equal to the threshold, applicable when --contigs-function-classify-method is set to emapper.'
        ),
        'contigs_function_emapper_method': (
            ['Params_c_sample', '02_function_emapper', 'method'], 'diamond', str,
            '\b\n02_function_emapper: search seed orthologs using one of [diamond/mmseqs], applicable when --contigs-function-classify-method is set to emapper.\n'
            'diamond: using default parameters of diamond mode of emapper.\n'
            'mmseqs: using default parameters of mmseqs mode of emapper'
        ),
        'contigs_function_emapper_query_cover': (
            ['Params_c_sample', '02_function_emapper', 'query_cover'], 30, float,
            '\b\n02_function_emapper: Report only alignments with query cover above or equal to the given percentage [0-100],.\n'
            'applicable when --contigs-function-classify-method is set to emapper'
        ),
        'contigs_function_emapper_extra_args': (
            ['Params_c_sample', '02_function_emapper', 'extra_args'], '', str,
            '02_function_emapper: Extra arguments forwarded to eggNOG-mapper (e.g., "--pident 30").'),
        'contigs_function_kraken2_confidence': (
            ['Params_c_sample', '02_function_kraken2', 'confidence'], 0.05, float,
            '02_function_kraken2: Confidence threshold [0-1] for functional classification, applicable when --contigs-function-classify-method is set to kraken2.'
        ),
        'contigs_function_kraken2_extra_args': (
            ['Params_c_sample', '02_function_kraken2', 'extra_args'], '', str,
            '02_function_kraken2: Extra arguments forwarded to kraken2 (e.g., "--memory-mapping").'),
        'mcr_unmapped_reads_merge_threads': (
            ['Params_mcr_sample', '01_unmapped_reads_trim_merge', 'threads'], 16, int,
            '01_unmapped_reads_trim_merge_fastp: Threads for fastp, samtools, and seqkit (paired-end), or samtools and seqkit only (single-end).'),
        'mcr_unmapped_reads_merge_length_required': (
            ['Params_mcr_sample', '01_unmapped_reads_trim_merge', 'length_required'], 15, int,
            '01_unmapped_reads_trim_merge_fastp: Reads shorter than length_required will be discarded. Only applicable when read_type is paired.'),
        'mcr_unmapped_reads_merge_n_base_limit': (
            ['Params_mcr_sample', '01_unmapped_reads_trim_merge', 'n_base_limit'], 5, int,
            "01_unmapped_reads_trim_merge_fastp: if one read's number of N base is >n_base_limit, then this read/pair is discarded. Only applicable when read_type is paired."),
        'mcr_unmapped_reads_merge_compression': (
            ['Params_mcr_sample', '01_unmapped_reads_trim_merge', 'compression'], 4, int,
            '01_unmapped_reads_trim_merge_fastp: Compression level for gzip output [1-9]. 1 is fastest, 9 is smallest. Only applicable when read_type is paired.'),
        'mcr_unmapped_reads_fastp_extra_args': (
            ['Params_mcr_sample', '01_unmapped_reads_trim_merge', 'extra_args'], '', str,
            '01_unmapped_reads_trim_merge: Extra arguments forwarded to fastp for unmapped reads trimming (e.g., "--reads_to_process 0 -Q"). Only applicable when read_type is paired.'),
        'mcr_fraggenescanrs_model': (
            ['Params_mcr_sample', '02_FragGeneScanRs', 'model'], 'illumina_1', str,
            '\b\n02_FragGeneScanRs: Model parameter used for short-read open reading frame prediction.\n'
            '[complete] for complete genomic sequences or short sequence reads without sequencing error\n'
            '[sanger_5] for Sanger sequencing reads with about 0.5% error rate\n'
            '[sanger_10] for Sanger sequencing reads with about 1% error rate\n'
            '[454_5] for 454 pyrosequencing reads with about 0.5% error rate\n'
            '[454_10] for 454 pyrosequencing reads with about 1% error rate\n'
            '[454_30] for 454 pyrosequencing reads with about 3% error rate\n'
            '[illumina_1] for Illumina sequencing reads with about 0.1% error rate\n'
            '[illumina_5] for Illumina sequencing reads with about 0.5% error rate\n'
            '[illumina_10] for Illumina sequencing reads with about 1% error rate'
        ),
        'mcr_fraggenescanrs_threads': (
            ['Params_mcr_sample', '02_FragGeneScanRs', 'threads'], 20, int,
            '02_FragGeneScanRs: Threads.'),
        'mcr_fraggenescanrs_stats_threads': (
            ['Params_mcr_sample', '02_FragGeneScanRs_stats', 'threads'], 10, int,
            '02_FragGeneScanRs_stats: Threads.'),
        'mcr_taxonomy_threads': (
            ['Params_mcr_sample', '03_taxonomy', 'threads'], 60, int,
            '03_taxonomy: Threads.'),
        'mcr_taxonomy_kraken2_confidence': (
            ['Params_mcr_sample', '03_taxonomy_kraken2', 'confidence'], 0.7, float,
            '\b\n03_taxonomy_kraken2: Confidence threshold [0-1] for taxonomic classification of unmapped reads with kraken2.\n'
            'Default value was tuned on the CAMI2 marine dataset. For samples with higher complexity, a lower confidence threshold is usually recommended.'
        ),
        'mcr_taxonomy_kraken2_extra_args': (
            ['Params_mcr_sample', '03_taxonomy_kraken2', 'extra_args'], '', str,
            '03_taxonomy_kraken2: Extra arguments forwarded to kraken2 (e.g., "--memory-mapping").'),
        'mcr_function_classify_method': (
            ['Params_mcr_sample', '04_function', 'reads_function_classify_method'], 'kraken2', str,
            '04_function: reads functional classifier using one of [kraken2/emapper].'
        ),
        'mcr_function_threads': (
            ['Params_mcr_sample', '04_function', 'threads'], 60, int,
            '04_function: Threads for [kraken2/emapper].'),
        'mcr_function_kraken2_confidence': (
            ['Params_mcr_sample', '04_function_kraken2', 'confidence'], 0.05, float,
            '04_function_kraken2: Confidence threshold [0-1] for functional classification for reads when --mcr-function-classify-method is set to kraken2.'
        ),
        'mcr_function_kraken2_extra_args': (
            ['Params_mcr_sample', '04_function_kraken2', 'extra_args'], '', str,
            '04_function_kraken2: Extra arguments forwarded to kraken2 (e.g., "--memory-mapping").'),
        'mcr_function_emapper_evalue': (
            ['Params_mcr_sample', '04_function_emapper', 'seed_ortholog_evalue'], 1e-10, float,
            '04_function_emapper: Report only alignments with e-value below or equal to the threshold, applicable when --mcr-function-classify-method is set to emapper.'
        ),
        'mcr_function_emapper_method': (
            ['Params_mcr_sample', '04_function_emapper', 'method'], 'diamond', str,
            '\b\n04_function_emapper: search seed orthologs using one of [diamond/mmseqs], applicable when --mcr-function-classify-method is set to emapper.\n'
            'diamond: using default parameters of diamond mode of emapper.\n'
            'mmseqs: using default parameters of mmseqs mode of emapper'
        ),
        'mcr_function_emapper_query_cover': (
            ['Params_mcr_sample', '04_function_emapper', 'query_cover'], 30, float,
            '\b\n04_function_emapper: Report only alignments with query cover above or equal to the given percentage [0-100],.\n'
            'applicable when --mcr-function-classify-method is set to emapper'
        ),
        'mcr_function_emapper_extra_args': (
            ['Params_mcr_sample', '04_function_emapper', 'extra_args'], '', str,
            '04_function_emapper: Extra arguments forwarded to eggNOG-mapper (e.g., "--pident 30").'),
        'merge_mmags_min_completeness1': (
            ['Params_mcr_sample', '05_merge_tax_fun_depth', 'MAGs_min_completeness1'], 50, float,
            '02_quality_gate: Minimum completeness threshold (%) for CheckM1-based mMAG quality filtering [0-100], e.g., 50 for 50%.'
        ),
        'merge_mmags_max_contamination1': (
            ['Params_mcr_sample', '05_merge_tax_fun_depth', 'MAGs_max_contamination1'], 1, float,
            '02_quality_gate: Maximum contamination threshold (%) for CheckM1-based mMAG quality filtering [0-100], e.g., 5 for 5%.'
        ),
        'merge_mmags_min_completeness2': (
            ['Params_mcr_sample', '05_merge_tax_fun_depth', 'MAGs_min_completeness2'], 50, float,
            '02_quality_gate: Minimum completeness threshold (%) for CheckM2-based mMAG quality filtering [0-100], e.g., 50 for 50%.'
        ),
        'merge_mmags_max_contamination2': (
            ['Params_mcr_sample', '05_merge_tax_fun_depth', 'MAGs_max_contamination2'], 1, float,
            '02_quality_gate: Maximum contamination threshold (%) for CheckM2-based mMAG quality filtering [0-100], e.g., 5 for 5%.'
        ),
        'merge_min_orf_length': (
            ['Params_mcr_sample', '05_merge_tax_fun_depth', 'minORFLen'], 50, int,
            '05_merge_tax_fun_depth: Minimum ORF length for finally used ORFs of mMAGs, contigs, and reads.'
        ),
        'merge_multiple_ko_handling_mode': (
            ['Params_mcr_sample', '05_merge_tax_fun_depth', 'KO_mode'], 'ignore_multiple_KOs', str,
            '\b\n05_merge_tax_fun_depth: Method using one of [split_multiple_KOs/keep_multiple_KOs/ignore_multiple_KOs] when an ORF has multiple KO annotations:.\n'
            'split_multiple_KOs: evenly distributing the total abundance across all assigned KOs;\n'
            'keep_multiple_KOs: retain original multiple KO annotations;\n'
            'ignore_multiple_KOs: ignore the ORFs with multiple KO annotations and only use the ORFs with single KO annotation.'
        ),
        'merge_min_ko_num_per_taxon': (
            ['Params_mcr_sample', '05_merge_tax_fun_depth', 'min_KO_num_in_one_taxon'], 100, int,
            '\b\n05_merge_tax_fun_depth: Minimum number of KOs per taxon, used to filter outputs in KO_Taxa_RelativeAbundance and KO_Taxa_SumCoverage.\n'
            'This threshold filters out taxa with insufficient KOs,\n'
            'retaining only KO-rich taxa to ensure robust metabolic reconstruction\n'
            'while preventing excessive output file sizes caused by sparse matrices.'
        ),
        'merge_tax_fun_depth_threads': (
            ['Params_mcr_sample', '05_merge_tax_fun_depth', 'threads'], 20, int,
            '05_merge_tax_fun_depth: Threads.'),
        'merge_cami_format_output': (
            ['Params_mcr_sample', '05_merge_tax_fun_depth', 'CAMI_format_output'], False, bool,
            '05_merge_tax_fun_depth: Report a taxonomic profile using the CAMI output format (flag only). Disabled when mMAG mode uses gtdbtk.'
        ),
        'merge_output_detail': (
            ['Params_mcr_sample', '05_merge_tax_fun_depth', 'output_detail'], False, bool,
            '05_merge_tax_fun_depth: Whether to output detail tables (flag only).'),
        'mclink_sample_min_ko_num_per_taxon': (
            ['Params_mcr_sample', '06_mclink_taxon', 'min_KO_num_in_one_taxon'], 100, int,
            '06_mclink_taxon: Minimum number of KOs per taxon for calculating taxon-level metabolic pathway/module completeness or abundance per sample. Taxa with fewer than this number of detected KOs in a given sample will be filtered out. Calculated independently for each sample.'
            ),
        'mclink_sample_min_ko_abundance_per_taxon': (
            ['Params_mcr_sample', '06_mclink_taxon', 'min_mean_KO_RA_in_one_taxon'], 0, float,
            '06_mclink_taxon: Minimum mean KO abundance [0-1e6] per taxon for calculating taxon-level metabolic pathway/module completeness or abundance per sample. Taxa with mean KO abundance below this threshold in a given sample will be filtered out. Note: GPM/TPM sum of all taxa in each sample is 1e6.'
            ),
        'mclink_sample_topn_taxon': (
            ['Params_mcr_sample', '06_mclink_taxon', 'TopN'], 100, int,
            '06_mclink_taxon: Number of top taxa (ranked by mean KO abundance) retained for calculating taxon-level metabolic pathway/module completeness or abundance per sample. Set to 0 to disable filtering.'
        ),
        'mclink_sample_completeness_per_taxon': (
            ['Params_mcr_sample', '06_mclink_taxon', 'completeness_threshold'], 0.5, float,
            '06_mclink_taxon: Completeness threshold [0-1] of modules for calculating taxon-level metabolic pathway/module abundance. Higher values retain only more complete profiles.'),
        'mclink_sample_plus_scale_method_per_taxon': (
            ['Params_mcr_sample', '06_mclink_taxon', 'plus_scale_method'], 'mean', str,
            "06_mclink_taxon: Aggregation method for '+'-separated KOs in a module [mean/max/min]. 'mean': average abundance across all required KOs; 'max': highest KO abundance; 'min': lowest KO abundance."),
        'mclink_sample_comma_scale_method_per_taxon': (
            ['Params_mcr_sample', '06_mclink_taxon', 'comma_scale_method'], 'sum', str,
            "06_mclink_taxon: Aggregation method for ','-separated alternative KOs in a module [sum/max]. 'sum': cumulative abundance of all alternative KOs; 'max': abundance of the highest single KO among alternatives. Note: this applies to abundance calculation only; module completeness is always computed with 'max'."),
        'barplot_topn_taxon': (
            ['Params_mcr_project', '01_Taxa_visualize', 'TopN'], 10, int,
            '01_Taxa_visualize: Number of top taxa to display in stacked barplots. Remaining categories are merged into Others. Maximum 29.'
        ),
        'barplot_topn_fun': (
            ['Params_mcr_project', '03_Fun_visualize', 'TopN'], 10, int,
            '03_Fun_visualize: Number of top functional categories to display in stacked barplots. Remaining categories are merged into Others. Maximum 29.'
        ),
        'mclink_min_ko_num_per_sample': (
            ['Params_mcr_project', '03_mclink_sample', 'min_sample_num_per_KO'], 1, int,
            '03_mclink_sample: Minimum number of samples in which a KO must be detected for calculating community-level metabolic pathway/module completeness or abundance.'
        ),
        'mclink_min_ko_prevalence_per_sample': (
            ['Params_mcr_project', '03_mclink_sample', 'preva'], 0, float,
            '03_mclink_sample: Minimum detection rate [0-1] of a KO across all samples for calculating community-level metabolic pathway/module completeness or abundance.'
        ),
        'mclink_min_ko_abundance_per_sample': (
            ['Params_mcr_project', '03_mclink_sample', 'abund'], 0, float,
            '03_mclink_sample: Minimum mean relative abundance [0-1e6] of a KO across all samples for calculating community-level metabolic pathway/module completeness or abundance. KOs with mean abundance below this threshold are excluded. Note: GPM/TPM sum of a sample is 1e6.'
        ),
        'mclink_completeness_per_sample': (
            ['Params_mcr_project', '03_mclink_sample', 'completeness_threshold'], 0.5, float,
            '03_mclink_sample: Completeness threshold [0-1] of modules for calculating community-level metabolic pathway/module abundance. Higher values retain only more complete profiles.'
        ),
        'mclink_plus_scale_method_per_sample': (
            ['Params_mcr_project', '03_mclink_sample', 'plus_scale_method'], 'mean', str,
            "03_mclink_sample: Aggregation method for '+'-separated KOs in a module [mean/max/min]. 'mean': average abundance across all required KOs; 'max': highest KO abundance; 'min': lowest KO abundance."
        ),
        'mclink_comma_scale_method_per_sample': (
            ['Params_mcr_project', '03_mclink_sample', 'comma_scale_method'], 'sum', str,
            "03_mclink_sample: Aggregation method for ','-separated alternative KOs in a module [sum/max]. 'sum': cumulative abundance of all alternative KOs; 'max': abundance of the highest single KO among alternatives. Note: this applies to abundance calculation only; module completeness is always computed with 'max'."
        ),
        'merge_min_taxa_ko_prevalence_across_samples': (
            ['Params_mcr_project', '05_merge_tax_KO_depth', 'preva'], 0.10, float,
            '05_merge_tax_KO_depth: Minimum prevalence [0-1] of a Taxon-KO across samples. Taxa-KOs detected in fewer samples than this proportion are excluded.'
        ),
        'merge_min_taxa_ko_abundance_across_samples': (
            ['Params_mcr_project', '05_merge_tax_KO_depth', 'abund'], 0.0, float,
            '05_merge_tax_KO_depth: Minimum mean relative abundance [0-1e6] of a Taxon-KO across samples. Taxa-KOs with mean abundance below this threshold are excluded. Note: GPM/TPM sum of a sample is 1e6.'
        ),
        'merge_min_taxa_module_prevalence_across_samples': (
            ['Params_mcr_project', '05_merge_tax_Module_RelativeAbundance', 'preva'], 0.10, float,
            '05_merge_tax_Module_RelativeAbundance: Minimum prevalence [0-1] of a Taxon-Module across samples. Taxa-Modules detected in fewer samples than this proportion are excluded.'
        ),
        'merge_min_taxa_module_abundance_across_samples': (
            ['Params_mcr_project', '05_merge_tax_Module_RelativeAbundance', 'abund'], 0.0, float,
            '05_merge_tax_Module_RelativeAbundance: Minimum mean relative abundance [0-1e6] of a Taxon-Module across samples. Taxa-Modules with mean abundance below this threshold are excluded. Note: GPM/TPM sum of a sample is 1e6.'
        ),
        'merge_min_taxa_module_sc_prevalence_across_samples': (
            ['Params_mcr_project', '05_merge_tax_Module_SumCoverage', 'preva'], 0.10, float,
            '05_merge_tax_Module_SumCoverage: Minimum prevalence [0-1] of a Taxon-Module across samples. Taxa-Modules detected in fewer samples than this proportion are excluded.'
        ),
        'merge_min_taxa_module_sc_abundance_across_samples': (
            ['Params_mcr_project', '05_merge_tax_Module_SumCoverage', 'abund'], 0.0, float,
            '05_merge_tax_Module_SumCoverage: Minimum mean sum coverage [0-1e6] of a Taxon-Module across samples. Taxa-Modules with mean sum coverage below this threshold are excluded. Note: GPM/TPM sum of a sample is 1e6.'
        ),
        'mclink_project_min_ko_num_per_taxon': (
            ['Params_mcr_project', '06_mclink_taxon', 'min_KO_num_in_one_taxon'], 100, int,
            '06_mclink_taxon: Minimum number of KOs per taxon for calculating taxon-level metabolic pathway/module completeness on the project-averaged taxon-function matrix. Taxa with fewer than this number of detected KOs in the project-averaged matrix will be filtered out. KO abundances are first averaged across all samples to generate a single project-level matrix before calculation.'
        ),
        'mclink_project_min_ko_abundance_per_taxon': (
            ['Params_mcr_project', '06_mclink_taxon', 'min_mean_KO_RA_in_one_taxon'], 0, float,
            '06_mclink_taxon: Minimum mean KO abundance [0-1e6] per taxon for calculating taxon-level metabolic pathway/module completeness on the project-averaged taxon-function matrix. Taxa with mean KO abundance below this threshold in the project-averaged matrix will be filtered out. Note: GPM/TPM sum of all taxa is 1e6 per sample; values are averaged across samples prior to filtering and calculation.'
        ),
        'mclink_project_topn_taxon': (
            ['Params_mcr_project', '06_mclink_taxon', 'TopN'], 100, int,
            '06_mclink_taxon: Number of top taxa (ranked by mean KO abundance in the project-averaged taxon-function matrix) retained for calculating taxon-level metabolic pathway/module completeness. Set to 0 to disable filtering.'
        ),
        'mclink_project_completeness_per_taxon': (
            ['Params_mcr_project', '06_mclink_taxon', 'completeness_threshold'], 0.5, float,
            '06_mclink_taxon: Completeness threshold [0-1] of modules for calculating taxon-level metabolic pathway/module abundance on the project-averaged taxon-function matrix. Higher values retain only more complete profiles.'),
        'mclink_project_plus_scale_method_per_taxon': (
            ['Params_mcr_project', '06_mclink_taxon', 'plus_scale_method'], 'mean', str,
            "06_mclink_taxon: Aggregation method for '+'-separated KOs in a module [mean/max/min] on the project-averaged taxon-function matrix. 'mean': average abundance across all required KOs; 'max': highest KO abundance; 'min': lowest KO abundance."),
        'mclink_project_comma_scale_method_per_taxon': (
            ['Params_mcr_project', '06_mclink_taxon', 'comma_scale_method'], 'sum', str,
            "06_mclink_taxon: Aggregation method for ','-separated alternative KOs in a module [sum/max] on the project-averaged taxon-function matrix. 'sum': cumulative abundance of all alternative KOs; 'max': abundance of the highest single KO among alternatives. Note: this applies to abundance calculation only; module completeness is always computed with 'max'."),
    },
    'pmr': {
        'mags': (
            ['mags'], 'm', str,
            '\b\nMAG source for the pmr chains, using one of [m/v/vm].\n'
            'm: project + mMAGs + reads (b + m-s + m-p + pmr_m).\n'
            'v: project + vMAGs + reads (b + v-s + v-p + pmr_v).\n'
            'vm: project + vMAGs + mMAGs + reads; vMAGs and mMAGs are combined into one reference for a single mapping pass (b + v-s + v-p + m-s + m-p + pmr_vm).'),
        'raw_stats_sample_threads': (
            ['Params_b_sample', '00_raw_stats', 'threads'], 4, int,
            '00_raw_stats: Threads for seqkit.'),
        'trim_sample_threads': (
            ['Params_b_sample', '01_trim', 'threads'], 16, int,
            '01_trim: Threads for fastp (next-generation short reads) or fastplong (third-generation long reads).'),
        'trim_length_required': (
            ['Params_b_sample', '01_trim', 'length_required'], 15, int,
            '01_trim: Reads shorter than length_required will be discarded'),
        'trim_n_base_limit': (
            ['Params_b_sample', '01_trim', 'n_base_limit'], 5, int,
            "01_trim: if one read's number of N base is >n_base_limit, then this read/pair is discarded."),
        'trim_compression': (
            ['Params_b_sample', '01_trim', 'compression'], 4, int,
            '01_trim: Compression level for gzip output [1-9]. 1 is fastest, 9 is smallest.'),
        'trim_extra_args': (
            ['Params_b_sample', '01_trim', 'extra_args'], '', str,
            '01_trim: Extra arguments forwarded to fastp or fastplong for raw reads trimming (e.g., "--reads_to_process 0 -Q").'),
        'trim_stats_sample_threads': (
            ['Params_b_sample', '01_trim_stats', 'threads'], 16, int,
            '01_trim_stats: Threads.'),
        'assembly_method': (
            ['Params_b_sample', '02_assembly', 'method'], 'megahit', str,
            '\b\n02_assembly: Underlying assembly software using one of [megahit/spades/metaflye].\n'
            'Must be compatible with --read-type (config key read_type):\n'
            'megahit: MEGAHIT (second-generation sequencing; for paired or single DNA reads).\n'
            'spades: SPAdes (second-generation sequencing; metaSPAdes for paired DNA reads, rnaSPAdes for paired-rna or single-rna RNA reads).\n'
            'metaflye: metaFlye (third-generation sequencing; required for pacbio-raw/corr/hifi or nano-raw/corr/hq long reads).'),
        'assembly_threads': (
            ['Params_b_sample', '02_assembly', 'threads'], 40, int,
            '02_assembly: Threads for [megahit/spades/metaflye].'),
        'megahit_presets': (
            ['Params_b_sample', '02_megahit', 'presets'], 'meta-sensitive', str,
            '\b\n02_megahit: MEGAHIT presets using one of [meta-sensitive/meta-large].\n'
            "meta-sensitive: '--min-count 1 --k-list 21,29,39,49,...,129,141'.\n"
            "meta-large: '--k-min 27 --k-max 127 --k-step 10'.\n"
            '(large & complex metagenomes, like soil).'),
        'megahit_k_list': (
            ['Params_b_sample', '02_megahit', 'k_list'], '', str,
            '02_megahit: MEGAHIT --k-list (optional). Explicit comma-separated k-mer sizes (e.g. "21,29,39,59,79,99,119,141"). Overrides --megahit-presets. If set, --megahit-k-min, --megahit-k-max, --megahit-k-step must NOT be set. Optional --megahit-min-count can be used together.'),
        'megahit_min_count': (
            ['Params_b_sample', '02_megahit', 'min_count'], 2, int,
            '02_megahit: MEGAHIT --min-count (optional). Minimum multiplicity for filtering (k_min+1)-mers. Can be used with --megahit-k-list or --megahit-k-min/--k-max/--k-step.'),
        'megahit_k_min': (
            ['Params_b_sample', '02_megahit', 'k_min'], '', str,
            '02_megahit: MEGAHIT --k-min (optional). Minimum kmer size (<= 255), must be odd number. If set with --k-max and --k-step, overrides --megahit-presets and uses custom k-mer mode.'),
        'megahit_k_max': (
            ['Params_b_sample', '02_megahit', 'k_max'], '', str,
            '02_megahit: MEGAHIT --k-max (optional). Maximum kmer size (<= 255), must be odd number.'),
        'megahit_k_step': (
            ['Params_b_sample', '02_megahit', 'k_step'], '', str,
            '02_megahit: MEGAHIT --k-step (optional). Increment of kmer size of each iteration (<= 28), must be even number.'),
        'megahit_extra_args': (
            ['Params_b_sample', '02_megahit', 'extra_args'], '', str,
            '02_megahit: Extra arguments forwarded to MEGAHIT (e.g. "--no-mercy --tmp-dir /path/to/your/tmp/dir"). Do NOT use this for k-mer settings (--k-list, --k-min, --presets, etc.).'),
        'spades_extra_args': (
            ['Params_b_sample', '02_spades', 'extra_args'], '', str,
            '02_spades: Extra arguments forwarded to SPAdes. Do NOT use this for input/output or thread settings.'),
        'metaflye_extra_args': (
            ['Params_b_sample', '02_metaflye', 'extra_args'], '', str,
            '02_metaflye: Extra arguments forwarded to Flye for third-generation long reads assembly (e.g., "--iterations 2"). Do NOT use this for read type or input/output settings.'),
        'mapping_threads': (
            ['Params_b_sample', '03_mapping', 'threads'], 20, int,
            '03_mapping: Threads for [minibwa/bwamem/bwamem2/strobealign].'),
        'mapping_method': (
            ['Params_b_sample', '03_mapping', 'method'], 'minibwa', str,
            '\b\n03_mapping: Underlying mapping software used [minibwa/bwamem/bwamem2/strobealign].\n'
            'minibwa: minibwa using default parameters.\n'
            'bwamem: bwa mem using default parameters.\n'
            'bwamem2: bwa-mem2 using default parameters.\n'
            'strobealign: strobealign using default parameters.'),
        'mapping_samtools_sort_mem': (
            ['Params_b_sample', '03_mapping', 'samtools_sort_mem'], '768M', str,
            '03_mapping: Set maximum memory per thread; suffix K/M/G recognized, e.g. 2G, 768M.'),
        'mapping_extra_args': (
            ['Params_b_sample', '03_mapping', 'extra_args'], '', str,
            '03_mapping: Extra arguments forwarded to the active mapping aligner (minibwa map, bwa mem, bwa-mem2 mem, or strobealign --use-index), (e.g. "-B 8" for minibwa map).'),
        'binning_metawrap_depth_threads': (
            ['Params_b_sample', '05_metawrap_binning_depth', 'threads'], 8, int,
            '05_metawrap_binning_depth: Threads.'),
        'binning_min_read_percent_identity': (
            ['Params_b_sample', '05_metawrap_binning_depth', 'min_read_percent_identity'], 97, float,
            '05_metawrap_binning_depth: Exclude reads by overall percent identity [0-100] during the binning process, e.g. 97 for 97%.'),
        'binning_contig_end_exclusion': (
            ['Params_b_sample', '05_metawrap_binning_depth', 'contig_end_exclusion'], 75, int,
            '05_metawrap_binning_depth: Exclude bases at the ends of contig sequences used for depth calculation and bining.'),
        'binning_min_contig_length': (
            ['Params_b_sample', '05_metawrap_binning_depth', 'min_contig_length'], 1500, int,
            '05_metawrap_binning_depth: Minimum contig length for binning.'),
        'binning_tools': (
            ['Params_b_sample', '05_metawrap_binning_depth', 'binning_tools'], 'metabat2,maxbin2,concoct', str,
            '\b\n05_metawrap_binning_depth: Comma-separated list of binning tools, choose 2 or 3 from [metabat2,maxbin2,concoct].\n'
            'At least 2 tools recommended for bin_refinement.'),
        'binning_max_memory_gb': (
            ['Params_b_sample', '05_metawrap_binning_depth', 'max_memory_gb'], 512, int,
            '05_metawrap_binning_depth: Maximum memory (GB) for MetaWRAP binning.'),
        'binning_metawrap_extra_args': (
            ['Params_b_sample', '05_metawrap_binning_depth', 'extra_args'], '', str,
            '05_metawrap_binning_depth: Extra arguments forwarded to metaWRAP binning (e.g., "-m 5").'),
        'metawrap_refinement_max_memory_gb': (
            ['Params_b_sample', '06_metawrap_bin_refinement', 'max_memory_gb'], 512, int,
            '06_metawrap_bin_refinement: Maximum memory in GB for MetaWRAP bin_refinement.'),
        'metawrap_refinement_extra_args': (
            ['Params_b_sample', '06_metawrap_bin_refinement', 'extra_args'], '', str,
            '06_metawrap_bin_refinement: Extra arguments forwarded to metaWRAP bin_refinement (e.g., "--remove-ambiguous").'),
        'bin_checkm_threads': (
            ['Params_b_sample', '06_bin_checkm', 'threads'], 8, int,
            '06_bin_checkm or 06_metawrap_bin_refinement: Threads.'),
        'bin_min_completeness': (
            ['Params_b_sample', '06_bin_checkm', 'min_completeness'], 50, float,
            '06_bin_checkm or 06_metawrap_bin_refinement: Minimum completeness (%) bin filtering, e.g. 50 for 50%.'),
        'bin_max_contamination': (
            ['Params_b_sample', '06_bin_checkm', 'max_contamination'], 10, float,
            '06_bin_checkm or 06_metawrap_bin_refinement: Maximum contamination (%) bin filtering, e.g. 10 for 10%.'),
        'bin_min_quality_score': (
            ['Params_b_sample', '06_bin_checkm', 'min_quality_score'], 50, float,
            '06_bin_checkm: Minimum quality score for bin filtering (GTDB-like), QS = completeness - 5*contamination, e.g. 50 for QS>=50.'),
        'contigs_depth_sample_threads': (
            ['Params_b_sample', '07_contigs_depth', 'threads'], 8, int,
            '07_contigs_depth: Threads.'),
        'depth_min_read_percent_identity': (
            ['Params_b_sample', '07_contigs_depth', 'min_read_percent_identity'], 95, float,
            '\b\n07_contigs_depth: Exclude reads by overall percent identity [0-100] for contig depth estimation\n'
            'using jgi_summarize_bam_contig_depths --percentIdentity, e.g. 95 for 95%.\n'
            'Reads with percent identity below this threshold are excluded from depth calculation.\n'
            'Higher values increase stringency, retaining only high-quality alignments.'),
        'depth_contig_end_exclusion': (
            ['Params_b_sample', '07_contigs_depth', 'contig_end_exclusion'], 75, int,
            '\b\n07_contigs_depth: Exclude bases at the ends of contig sequences for depth estimation using jgi_summarize_bam_contig_depths --maxEdgeBases.\n'
            'Bases within this distance from either end of a contig are ignored when calculating depth,\n'
            'preventing edge effects from partial alignments.'),
        'depth_min_contig_length': (
            ['Params_b_sample', '07_contigs_depth', 'min_contig_length'], 300, int,
            '\b\n07_contigs_depth: Minimum contig length for depth estimation [bp], e.g. 300 for 300 bp.\n'
            'Contigs shorter than this threshold are excluded from depth calculation.'),
        'contigs_prodigal_extra_args': (
            ['Params_b_sample', '07_contigs_depth', 'extra_args'], '', str,
            '07_contigs_depth: Extra arguments forwarded to prodigal for contig gene prediction (e.g., "-m -n").'),
        'bin_prodigal_threads': (
            ['Params_b_sample', '08_bin_prodigal', 'threads'], 8, int,
            '08_bin_prodigal: Threads.'),
        'bin_prodigal_extra_args': (
            ['Params_b_sample', '08_bin_prodigal', 'extra_args'], '', str,
            '08_bin_prodigal: Extra arguments forwarded to prodigal for MAG gene prediction (e.g., "-m -n").'),
        'raw_stats_project_threads': (
            ['Params_b_project', '01_cat_raw_stats', 'threads'], 4, int,
            '01_cat_raw_stats: Threads.'),
        'trimmed_stats_project_threads': (
            ['Params_b_project', '02_cat_trimmed_stats', 'threads'], 4, int,
            '02_cat_trimmed_stats: Threads.'),
        'contigs_stats_project_threads': (
            ['Params_b_project', '03_contigs_stats', 'threads'], 4, int,
            '03_contigs_stats: Threads.'),
        'merge_stats_project_threads': (
            ['Params_b_project', '04_merge_reads_contigs_stats', 'threads'], 4, int,
            '04_merge_reads_contigs_stats: Threads.'),
        'gather_checkm_project_threads': (
            ['Params_b_project', '05_gather_project_bin_checkm_stats', 'threads'], 4, int,
            '05_gather_project_bin_checkm_stats: Threads.'),
        'mmags_checkm2_threads': (
            ['Params_m_sample', '01_checkm2', 'threads'], 60, int,
            '01_checkm2: Threads.'),
        'mmags_checkm2_extra_args': (
            ['Params_m_sample', '01_checkm2', 'extra_args'], '', str,
            '01_checkm2: Extra arguments forwarded to checkm2 predict (e.g., "--remove_intermediates").'),
        'mmags_taxonomy_classify_method': (
            ['Params_m_sample', '02_taxonomy', 'MAGs_taxonomy_classify_method'], 'gtdbtk', str,
            '02_taxonomy: mMAGs taxonomy classifier, should be one of the [gtdbtk/BAT].'
        ),
        'mmags_taxonomy_threads': (
            ['Params_m_sample', '02_taxonomy', 'threads'], 60, int,
            '02_taxonomy_gtdbtk: Threads for [gtdbtk/BAT].'),
        'mmags_taxonomy_gtdbtk_extra_args': (
            ['Params_m_sample', '02_taxonomy_gtdbtk', 'extra_args'], '', str,
            '02_taxonomy_gtdbtk: Extra arguments forwarded to gtdbtk classify_wf (e.g., "--min_af 0.5").'),
        'mmags_taxonomy_bat_range': (
            ['Params_m_sample', '02_taxonomy_BAT', 'range'], 10, float,
            '\b\n02_taxonomy_BAT: Percentage range [0-100] of DIAMOND hits relative to top bit-score to include in LCA calculation.\n'
            'Higher values yield more conservative classifications by incorporating divergent hits, while lower values enable specific assignments.\n'
            'Applicable when --mMAGs-taxonomy-classify-method is set to BAT.'
        ),
        'mmags_taxonomy_bat_fraction': (
            ['Params_m_sample', '02_taxonomy_BAT', 'fraction'], 0.50, float,
            '\b\n02_taxonomy_BAT: Minimum bit-score support fraction [0-0.99].\n'
            'Higher values enforce stricter consensus (one classification per sequence), while lower values allow multiple lineage detections.\n'
            'Applicable when --mMAGs-taxonomy-classify-method is set to BAT.'
        ),
        'mmags_taxonomy_bat_extra_args': (
            ['Params_m_sample', '02_taxonomy_BAT', 'extra_args'], '', str,
            '02_taxonomy_BAT: Extra arguments forwarded to CAT_pack bins (e.g., "--block_size 12.0").'),
        'mmags_drep99_completeness_per_sample': (
            ['Params_m_sample', '05_drep_eachsample_99', 'comp'], 50, float,
            '05_drep_eachsample_99: The first level dereplication of mMAGs in a sample using completeness threshold [0-100], e.g., 50 for 50%.'
        ),
        'mmags_drep99_contamination_per_sample': (
            ['Params_m_sample', '05_drep_eachsample_99', 'con'], 10, float,
            '05_drep_eachsample_99: The first level dereplication of mMAGs in a sample using contamination threshold [0-100], e.g., 10 for 10%.'
        ),
        'mmags_drep99_ani_per_sample': (
            ['Params_m_sample', '05_drep_eachsample_99', 'sa'], 0.99, float,
            '05_drep_eachsample_99: The first level dereplication of mMAGs in a sample using ANI threshold [0-1], e.g., 0.99 for 99%.'
        ),
        'mmags_drep99_sample_threads': (
            ['Params_m_sample', '05_drep_eachsample_99', 'threads'], 20, int,
            '05_drep_eachsample_99: Threads.'),
        'mmags_drep99_sample_extra_args': (
            ['Params_m_sample', '05_drep_eachsample_99', 'extra_args'], '', str,
            '05_drep_eachsample: Extra arguments forwarded to dRep dereplicate (e.g., "--completeness_weight 1 --contamination_weight 5").'),
        'mmags_drep95_completeness_per_sample': (
            ['Params_m_sample', '06_drep_eachsample_95', 'comp'], 50, float,
            '06_drep_eachsample_95: The second level dereplication of mMAGs in a sample using completeness threshold [0-100], e.g., 50 for 50%.'
        ),
        'mmags_drep95_contamination_per_sample': (
            ['Params_m_sample', '06_drep_eachsample_95', 'con'], 10, float,
            '06_drep_eachsample_95: The second level dereplication of mMAGs in a sample using contamination threshold [0-100], e.g., 10 for 10%.'
        ),
        'mmags_drep95_ani_per_sample': (
            ['Params_m_sample', '06_drep_eachsample_95', 'sa'], 0.95, float,
            '\b\n06_drep_eachsample_95: The second level dereplication of mMAGs in a sample using ANI threshold [0-1], e.g., 0.95 for 95%.\n'
            'This value should be smaller than --drep99-ANI-per-sample'
        ),
        'mmags_drep95_sample_threads': (
            ['Params_m_sample', '06_drep_eachsample_95', 'threads'], 20, int,
            '06_drep_eachsample_95: Threads.'),
        'mmags_drep95_sample_extra_args': (
            ['Params_m_sample', '06_drep_eachsample_95', 'extra_args'], '', str,
            '06_drep_eachsample: Extra arguments forwarded to dRep dereplicate (e.g., "--completeness_weight 1 --contamination_weight 5").'),
        'mmags_coverm_method': (
            ['Params_m_sample', '07_coverm_eachsample', 'method'], 'trimmed_mean', str,
            '\b\n07_coverm_eachsample: Method(s) for coverm-based coverage calculation, must be one of [trimmed_mean/mean].\n'
            'A more detailed description of the different methods is available at https://github.com/wwood/CoverM#calculation-methods'
        ),
        'mmags_coverm_min_read_percent_identity': (
            ['Params_m_sample', '07_coverm_eachsample', 'min_read_percent_identity'], 95, float,
            '07_coverm_eachsample: Exclude reads by overall percent identity [0-100], e.g. 95 for 95%.'
        ),
        'mmags_coverm_min_read_aligned_length': (
            ['Params_m_sample', '07_coverm_eachsample', 'min_read_aligned_length'], 45, int,
            '07_coverm_eachsample: Exclude reads with fewer aligned bases.'
        ),
        'mmags_coverm_contig_end_exclusion': (
            ['Params_m_sample', '07_coverm_eachsample', 'contig_end_exclusion'], 75, int,
            '07_coverm_eachsample: Exclude bases at the ends of contig sequences for coverm calculation.'
        ),
        'mmags_coverm_threads': (
            ['Params_m_sample', '07_coverm_eachsample', 'threads'], 20, int,
            '07_coverm_eachsample: Threads.'),
        'mmags_coverm_sample_extra_args': (
            ['Params_m_sample', '07_coverm_eachsample', 'extra_args'], '', str,
            '07_coverm_eachsample: Extra arguments forwarded to coverm genome (e.g., "--quiet").'),
        'mmags_function_classify_method': (
            ['Params_m_sample', '08_function', 'MAGs_function_classify_method'], 'emapper', str,
            '08_function: mMAGs function classifier using one of [kraken2/emapper].'
        ),
        'mmags_function_threads': (
            ['Params_m_sample', '08_function', 'threads'], 60, int,
            '08_function: Threads for [kraken2/emapper].'),
        'mmags_function_kraken2_confidence': (
            ['Params_m_sample', '08_function_kraken2', 'confidence'], 0.05, float,
            '08_function_kraken2: Confidence threshold [0-1] for functional classification for mMAGs when --mMAGs-function-classify-method is set to kraken2.'
        ),
        'mmags_function_kraken2_extra_args': (
            ['Params_m_sample', '08_function_kraken2', 'extra_args'], '', str,
            '08_function_kraken2: Extra arguments forwarded to kraken2 (e.g., "--memory-mapping").'),
        'mmags_function_emapper_evalue': (
            ['Params_m_sample', '08_function_emapper', 'seed_ortholog_evalue'], 1e-10, float,
            '08_function_emapper: Report only alignments with e-value below or equal to the threshold, applicable when --mMAGs-function-classify-method is set to emapper.'
        ),
        'mmags_function_emapper_method': (
            ['Params_m_sample', '08_function_emapper', 'method'], 'diamond', str,
            '\b\n08_function_emapper: search seed orthologs using one of [diamond/mmseqs], applicable when --mMAGs-function-classify-method is set to emapper.\n'
            'diamond: using default parameters of diamond mode of emapper.\n'
            'mmseqs: using default parameters of mmseqs mode of emapper'
        ),
        'mmags_function_emapper_query_cover': (
            ['Params_m_sample', '08_function_emapper', 'query_cover'], 30, float,
            '\b\n08_function_emapper: Report only alignments with query cover above or equal to the given percentage [0-100],.\n'
            'applicable when --mMAGs-function-classify-method is set to emapper'
        ),
        'mmags_function_emapper_extra_args': (
            ['Params_m_sample', '08_function_emapper', 'extra_args'], '', str,
            '08_function_emapper: Extra arguments forwarded to eggNOG-mapper (e.g., "--pident 30").'),
        'mmags_drep99_completeness_per_project': (
            ['Params_m_project', '01_drep_project_99', 'comp'], 50, float,
            '01_drep_project_99: Dereplicate mMAGs from all samples in a project with completeness threshold [0-100], e.g., 50 for 50%.'
        ),
        'mmags_drep99_contamination_per_project': (
            ['Params_m_project', '01_drep_project_99', 'con'], 10, float,
            '01_drep_project_99: Dereplicate mMAGs from all samples in a project with contamination threshold [0-100], e.g., 5 for 5%.'
        ),
        'mmags_drep99_project_threads': (
            ['Params_m_project', '01_drep_project_99', 'threads'], 60, int,
            '01_drep_project_99: Threads.'),
        'mmags_drep99_project_extra_args': (
            ['Params_m_project', '01_drep_project_99', 'extra_args'], '', str,
            '01_drep_project_99: Extra arguments forwarded to dRep dereplicate ("--completeness_weight 1 --contamination_weight 5").'),
        'mmags_drep95_completeness_per_project': (
            ['Params_m_project', '01_drep_project_95', 'comp'], 50, float,
            '01_drep_project_95: Dereplicate mMAGs from all samples in a project with completeness threshold [0-100], e.g., 50 for 50%.'
        ),
        'mmags_drep95_contamination_per_project': (
            ['Params_m_project', '01_drep_project_95', 'con'], 10, float,
            '01_drep_project_95: Dereplicate mMAGs from all samples in a project with contamination threshold [0-100], e.g., 5 for 5%.'
        ),
        'mmags_drep95_project_threads': (
            ['Params_m_project', '01_drep_project_95', 'threads'], 60, int,
            '01_drep_project_95: Threads.'),
        'mmags_drep95_project_extra_args': (
            ['Params_m_project', '01_drep_project_95', 'extra_args'], '', str,
            '01_drep_project_95: Extra arguments forwarded to dRep dereplicate (e.g., "--completeness_weight 1 --contamination_weight 5").'),
        'virus_identify_method': (
            ['Params_v_sample', '01_virus_identify', 'identify_method'], 'genomad_vs_vb', str,
            '\b\n01_virus_identify: Virus identification method, should be one of the [genomad/vs_vb/genomad_vs_vb].\n'
            'genomad: Use geNomad;\n'
            'vs_vb: Use VirSorter2 and CheckV, then use VIBRANT to check the low quality viruses from VirSorter2;\n'
            'genomad_vs_vb: Combine genomad and vs_vb virus to get the union results.'
        ),
        'virus_length_limit': (
            ['Params_v_sample', '01_virus_identify', 'virus_length_limit'], 2000, int,
            '01_virus_identify: Minimum length (bp) for viral contigs.'
        ),
        'virus_identify_threads': (
            ['Params_v_sample', '01_virus_identify', 'threads'], 8, int,
            '01_virus_identify: Threads.'),
        'virus_contig_taxonomy_cat_range': (
            ['Params_v_sample', '03_taxonomy_CAT_pack', 'range'], 10, float,
            '\b\n03_taxonomy_CAT_pack: Percentage range [0-100] of DIAMOND hits relative to top bit-score to include in LCA calculation.\n'
            'Higher values yield more conservative classifications by incorporating divergent hits, while lower values enable specific assignments.'
        ),
        'virus_contig_taxonomy_cat_fraction': (
            ['Params_v_sample', '03_taxonomy_CAT_pack', 'fraction'], 0.50, float,
            '\b\n03_taxonomy_CAT_pack: Minimum bit-score support fraction [0-0.99].\n'
            'Higher values enforce stricter consensus (one classification per sequence), while lower values allow multiple lineage detections.'
        ),
        'virus_contig_taxonomy_cat_pack_threads': (
            ['Params_v_sample', '03_taxonomy_CAT_pack', 'threads'], 8, int,
            '03_taxonomy_CAT_pack: Threads.'),
        'virus_contig_taxonomy_cat_pack_extra_args': (
            ['Params_v_sample', '03_taxonomy_CAT_pack', 'extra_args'], '', str,
            '03_taxonomy_CAT_pack: Extra arguments forwarded to CAT_pack contigs (e.g., "--block_size 12.0").'),
        'virus_mapping_threads': (
            ['Params_v_sample', '03_mapping', 'threads'], 20, int,
            '03_mapping: Threads for [minibwa/bwamem/bwamem2/strobealign].'),
        'virus_mapping_method': (
            ['Params_v_sample', '03_mapping', 'method'], 'minibwa', str,
            '\b\n03_mapping: Underlying mapping software used [minibwa/bwamem/bwamem2/strobealign].\n'
            'minibwa: minibwa using default parameters.\n'
            'bwamem: bwa mem using default parameters.\n'
            'bwamem2: bwa-mem2 using default parameters.\n'
            'strobealign: strobealign using default parameters.'),
        'virus_mapping_samtools_sort_mem': (
            ['Params_v_sample', '03_mapping', 'samtools_sort_mem'], '768M', str,
            '03_mapping: Set maximum memory per thread; suffix K/M/G recognized, e.g. 2G, 768M.'),
        'virus_mapping_extra_args': (
            ['Params_v_sample', '03_mapping', 'extra_args'], '', str,
            '03_mapping: Extra arguments forwarded to the active mapping aligner (minibwa map, bwa mem, bwa-mem2 mem, or strobealign --use-index), (e.g. "-B 8" for minibwa map).'),
        'virus_binning_threads': (
            ['Params_v_sample', '04_binning', 'threads'], 8, int,
            '04_binning_vRhyme: Threads.'),
        'virus_checkv_threads': (
            ['Params_v_sample', '05_checkv', 'threads'], 8, int,
            '05_checkv: Threads.'),
        'virus_vclust_ani_sample': (
            ['Params_v_sample', '06_vclust', 'ani'], 0.95, float,
            '06_vclust: ANI threshold [0-1] of vclust for sample-level virus clustering.'
        ),
        'virus_vclust_threads_sample': (
            ['Params_v_sample', '06_vclust', 'threads'], 8, int,
            '06_vclust: Threads.'),
        'vclust_sample_prefilter_extra_args': (
            ['Params_v_sample', '06_vclust', 'prefilter_extra_args'], '', str,
            '06_vclust: Extra arguments forwarded to vclust prefilter (e.g., "-k 25").'),
        'vclust_sample_align_extra_args': (
            ['Params_v_sample', '06_vclust', 'align_extra_args'], '', str,
            '06_vclust: Extra arguments forwarded to vclust align (e.g., "--mal 11").'),
        'vclust_sample_cluster_extra_args': (
            ['Params_v_sample', '06_vclust', 'cluster_extra_args'], '', str,
            '06_vclust: Extra arguments forwarded to vclust cluster (e.g., "--len_ratio 0").'),
        'virus_iphop_threads_sample': (
            ['Params_v_sample', '07_iPHoP', 'threads'], 20, int,
            '07_iPHoP: Threads.'),
        'virus_vmags_taxonomy_bat_range': (
            ['Params_v_sample', '08_taxonomy_BAT', 'range'], 10, float,
            '\b\n08_taxonomy_BAT: Percentage range [0-100] of DIAMOND hits relative to top bit-score to include in LCA calculation.\n'
            'Higher values yield more conservative classifications by incorporating divergent hits, while lower values enable specific assignments.'
        ),
        'virus_vmags_taxonomy_bat_fraction': (
            ['Params_v_sample', '08_taxonomy_BAT', 'fraction'], 0.50, float,
            '\b\n08_taxonomy_BAT: Minimum bit-score support fraction [0-0.99].\n'
            'Higher values enforce stricter consensus (one classification per sequence), while lower values allow multiple lineage detections.'
        ),
        'virus_vmags_taxonomy_bat_threads': (
            ['Params_v_sample', '08_taxonomy_BAT', 'threads'], 8, int,
            '08_taxonomy_BAT: Threads.'),
        'virus_vmags_taxonomy_bat_extra_args': (
            ['Params_v_sample', '08_taxonomy_BAT', 'extra_args'], '', str,
            '08_taxonomy_BAT: Extra arguments forwarded to CAT_pack bins (e.g., "--block_size 12.0").'),
        'virus_coverm_method_sample': (
            ['Params_v_sample', '09_coverm_eachsample_vMAGs', 'method'], 'trimmed_mean', str,
            '\b\n09_coverm_eachsample_vMAGs: Method(s) for coverm-based coverage calculation, must be one of [trimmed_mean/mean].\n'
            'A more detailed description of the different methods is available at https://github.com/wwood/CoverM#calculation-methods'
        ),
        'virus_coverm_min_read_percent_identity_sample': (
            ['Params_v_sample', '09_coverm_eachsample_vMAGs', 'min_read_percent_identity'], 95, float,
            '09_coverm_eachsample_vMAGs: Exclude reads by overall percent identity [0-100], e.g. 95 for 95%.'
        ),
        'virus_coverm_min_read_aligned_length_sample': (
            ['Params_v_sample', '09_coverm_eachsample_vMAGs', 'min_read_aligned_length'], 45, int,
            '09_coverm_eachsample_vMAGs: Exclude reads with fewer aligned bases.'
        ),
        'virus_coverm_contig_end_exclusion_sample': (
            ['Params_v_sample', '09_coverm_eachsample_vMAGs', 'contig_end_exclusion'], 75, int,
            '09_coverm_eachsample_vMAGs: Exclude bases at the ends of contig sequences for coverm calculation.'
        ),
        'virus_coverm_threads_sample': (
            ['Params_v_sample', '09_coverm_eachsample_vMAGs', 'threads'], 8, int,
            '09_coverm_eachsample_vMAGs: Threads.'),
        'virus_coverm_sample_extra_args': (
            ['Params_v_sample', '09_coverm_eachsample_vMAGs', 'extra_args'], '', str,
            '09_coverm_eachsample_vMAGs: Extra arguments forwarded to coverm genome (e.g., "--quiet").'),
        'virus_vclust_ani_project': (
            ['Params_v_project', '06_vclust', 'ani'], 0.95, float,
            '06_vclust: ANI threshold [0-1] of vclust for project-level virus clustering.'
        ),
        'virus_vclust_threads_project': (
            ['Params_v_project', '06_vclust', 'threads'], 8, int,
            '06_vclust: Threads.'),
        'vclust_project_prefilter_extra_args': (
            ['Params_v_project', '06_vclust', 'prefilter_extra_args'], '', str,
            '06_vclust: Extra arguments forwarded to vclust prefilter (e.g., "-k 25".'),
        'vclust_project_align_extra_args': (
            ['Params_v_project', '06_vclust', 'align_extra_args'], '', str,
            '06_vclust: Extra arguments forwarded to vclust align (e.g., "--mal 11").'),
        'vclust_project_cluster_extra_args': (
            ['Params_v_project', '06_vclust', 'cluster_extra_args'], '', str,
            '06_vclust: Extra arguments forwarded to vclust cluster (e.g., "--len_ratio 0").'),
        'virus_iphop_custom_mags_dir_project': (
            ['Params_v_project', '07_iPHoP', 'custom_MAGs_dir'], 'none', str,
            '\b\n07_iPHoP: Custom mMAGs directory for iPHoP host prediction, using one of ["none"/user-defined-path].\n'
            'Use "none" to disable custom mMAGs.\n'
            'Use an absolute path to a custom MAGs directory (including *.fna files).\n'
            'NOTE: "auto_drep99" and "auto_drep95" (automatically using dereplicated mMAGs with metalink default path) are only effective when the m-p chain runs (vm mode or pmr --mags vm).\n'
            'This parameter is not applicable in vc/vcr/vmc/vmcr modes (iPHoP host prediction is skipped).'
        ),
        'virus_iphop_threads_project': (
            ['Params_v_project', '07_iPHoP', 'threads'], 60, int,
            '07_iPHoP: Threads.'),
        'pmr_mapping_threads': (
            ['Params_pmr_{mags}_sample', '03_mapping', 'threads'], 20, int,
            '03_mapping: Threads for [minibwa/bwamem/bwamem2/strobealign].'),
        'pmr_mapping_method': (
            ['Params_pmr_{mags}_sample', '03_mapping', 'method'], 'minibwa', str,
            '\b\n03_mapping: Underlying mapping software used [minibwa/bwamem/bwamem2/strobealign].\n'
            'minibwa: minibwa using default parameters.\n'
            'bwamem: bwa mem using default parameters.\n'
            'bwamem2: bwa-mem2 using default parameters.\n'
            'strobealign: strobealign using default parameters.'),
        'pmr_mapping_samtools_sort_mem': (
            ['Params_pmr_{mags}_sample', '03_mapping', 'samtools_sort_mem'], '768M', str,
            '03_mapping: Set maximum memory per thread; suffix K/M/G recognized, e.g. 2G, 768M.'),
        'pmr_mapping_extra_args': (
            ['Params_pmr_{mags}_sample', '03_mapping', 'extra_args'], '', str,
            '03_mapping: Extra arguments forwarded to the active mapping aligner (minibwa map, bwa mem, bwa-mem2 mem, or strobealign --use-index), (e.g. "-B 8" for minibwa map).'),
        'coverm_mags_method': (
            ['Params_pmr_{mags}_sample', '03_coverm_MAGs', 'method'], 'trimmed_mean', str,
            '\b\n03_coverm_MAGs: Method(s) for coverm-based coverage calculation, must be one of [trimmed_mean/mean].\n'
            'A more detailed description of the different methods is available at https://github.com/wwood/CoverM#calculation-methods'
        ),
        'coverm_mags_read_percent_identity': (
            ['Params_pmr_{mags}_sample', '03_coverm_MAGs', 'min_read_percent_identity'], 95, float,
            '03_coverm_MAGs: Exclude reads by overall percent identity [0-100], e.g. 95 for 95%.'
        ),
        'coverm_mags_min_read_aligned_length': (
            ['Params_pmr_{mags}_sample', '03_coverm_MAGs', 'min_read_aligned_length'], 45, int,
            '03_coverm_MAGs: Exclude reads with fewer aligned bases.'
        ),
        'coverm_mags_contig_end_exclusion': (
            ['Params_pmr_{mags}_sample', '03_coverm_MAGs', 'contig_end_exclusion'], 75, int,
            '03_coverm_MAGs: Exclude bases at the ends of contig sequences for coverm calculation.'
        ),
        'coverm_mags_threads': (
            ['Params_pmr_{mags}_sample', '03_coverm_MAGs', 'threads'], 8, int,
            '03_coverm_MAGs: Threads.'),
        'coverm_mags_extra_args': (
            ['Params_pmr_{mags}_sample', '03_coverm_MAGs', 'extra_args'], '', str,
            '03_coverm_MAGs: Extra arguments forwarded to coverm genome (e.g., "--quiet").'),
        'coverm_contigs_method': (
            ['Params_pmr_{mags}_sample', '03_coverm_contigs', 'method'], 'trimmed_mean', str,
            '\b\n03_coverm_contigs: Method(s) for coverm-based coverage calculation, must be one of [trimmed_mean/mean].\n'
            'A more detailed description of the different methods is available at https://github.com/wwood/CoverM#calculation-methods'
        ),
        'coverm_contigs_read_percent_identity': (
            ['Params_pmr_{mags}_sample', '03_coverm_contigs', 'min_read_percent_identity'], 95, float,
            '03_coverm_contigs: Exclude reads by overall percent identity [0-100], e.g. 95 for 95%.'
        ),
        'coverm_contigs_min_read_aligned_length': (
            ['Params_pmr_{mags}_sample', '03_coverm_contigs', 'min_read_aligned_length'], 45, int,
            '03_coverm_contigs: Exclude reads with fewer aligned bases.'
        ),
        'coverm_contigs_contig_end_exclusion': (
            ['Params_pmr_{mags}_sample', '03_coverm_contigs', 'contig_end_exclusion'], 75, int,
            '03_coverm_contigs: Exclude bases at the ends of contig sequences for coverm calculation.'
        ),
        'coverm_contigs_threads': (
            ['Params_pmr_{mags}_sample', '03_coverm_contigs', 'threads'], 20, int,
            '03_coverm_contigs: Threads.'),
        'coverm_contigs_extra_args': (
            ['Params_pmr_{mags}_sample', '03_coverm_contigs', 'extra_args'], '', str,
            '03_coverm_contigs: Extra arguments forwarded to coverm contig (e.g., "--quiet").'),
        'mmags_min_completeness1': (
            ['Params_pmr_{mags}_project', '02_quality_gate', 'MAGs_min_completeness1'], 0, float,
            '02_quality_gate: Minimum completeness threshold (%) for CheckM1-based mMAG quality filtering [0-100], e.g., 50 for 50%.'
        ),
        'mmags_max_contamination1': (
            ['Params_pmr_{mags}_project', '02_quality_gate', 'MAGs_max_contamination1'], 100, float,
            '02_quality_gate: Maximum contamination threshold (%) for CheckM1-based mMAG quality filtering [0-100], e.g., 10 for 10%.'
        ),
        'mmags_min_completeness2': (
            ['Params_pmr_{mags}_project', '02_quality_gate', 'MAGs_min_completeness2'], 0, float,
            '02_quality_gate: Minimum completeness threshold (%) for CheckM2-based mMAG quality filtering [0-100], e.g., 50 for 50%.'
        ),
        'mmags_max_contamination2': (
            ['Params_pmr_{mags}_project', '02_quality_gate', 'MAGs_max_contamination2'], 100, float,
            '02_quality_gate: Maximum contamination threshold (%) for CheckM2-based mMAG quality filtering [0-100], e.g., 10 for 10%.'
        ),
        'vmags_checkv_quality': (
            ['Params_pmr_{mags}_project', '02_quality_gate', 'checkv_quality'], 'Not-determined', str,
            '\b\n02_quality_gate: Minimum CheckV quality tier for vMAG quality filtering, using one of\n'
            '[Complete/High-quality/Medium-quality/Low-quality/Not-determined].\n'
            'vMAGs below the given tier are removed (MAGs_v_quality_filter.cmd.R). Only applies to --mags v/vm.'
        ),
        'merge_multiple_ko_handling_mode': (
            ['Params_pmr_{mags}_sample', '05_merge_tax_fun_depth', 'KO_mode'], 'ignore_multiple_KOs', str,
            '\b\n05_merge_tax_fun_depth: Method using one of [split_multiple_KOs/keep_multiple_KOs/ignore_multiple_KOs] when an ORF has multiple KO annotations:.\n'
            'split_multiple_KOs: evenly distributing the total abundance across all assigned KOs;\n'
            'keep_multiple_KOs: retain original multiple KO annotations;\n'
            'ignore_multiple_KOs: ignore the ORFs with multiple KO annotations and only use the ORFs with single KO annotation.'
        ),
        'merge_min_orf_length': (
            ['Params_pmr_{mags}_sample', '05_merge_tax_fun_depth', 'minORFLen'], 50, int,
            '05_merge_tax_fun_depth: Minimum ORF length for finally used ORFs of mMAGs and reads.'
        ),
        'merge_cami_format_output': (
            ['Params_pmr_{mags}_sample', '05_merge_tax_fun_depth', 'CAMI_format_output'], False, bool,
            '05_merge_tax_fun_depth: Report a taxonomic profile using the CAMI output format (flag only). Disabled when mMAG mode uses gtdbtk.'
        ),
        'merge_output_detail': (
            ['Params_pmr_{mags}_sample', '05_merge_tax_fun_depth', 'output_detail'], False, bool,
            '05_merge_tax_fun_depth: Whether to output detail tables (flag only).'),
        'merge_min_ko_num_per_taxon': (
            ['Params_pmr_{mags}_sample', '05_merge_tax_fun_depth', 'min_KO_num_in_one_taxon'], 100, int,
            '\b\n05_merge_tax_fun_depth: Minimum number of KOs per taxon, used to filter outputs in KO_Taxa_RelativeAbundance and KO_Taxa_SumCoverage.\n'
            'This threshold filters out taxa with insufficient KOs,\n'
            'retaining only KO-rich taxa to ensure robust metabolic reconstruction\n'
            'while preventing excessive output file sizes caused by sparse matrices.'
        ),
        'merge_tax_fun_depth_threads': (
            ['Params_pmr_{mags}_sample', '05_merge_tax_fun_depth', 'threads'], 20, int,
            '05_merge_tax_fun_depth: Threads.'),
        'mclink_sample_min_ko_num_per_taxon': (
            ['Params_pmr_{mags}_sample', '06_mclink_taxon', 'min_KO_num_in_one_taxon'], 100, int,
            '06_mclink_taxon: Minimum number of KOs per taxon for calculating taxon-level metabolic pathway/module completeness or abundance per sample. Taxa with fewer than this number of detected KOs in a given sample will be filtered out. Calculated independently for each sample.'
            ),
        'mclink_sample_min_ko_abundance_per_taxon': (
            ['Params_pmr_{mags}_sample', '06_mclink_taxon', 'min_mean_KO_RA_in_one_taxon'], 0, float,
            '06_mclink_taxon: Minimum mean KO abundance [0-1e6] per taxon for calculating taxon-level metabolic pathway/module completeness or abundance per sample. Taxa with mean KO abundance below this threshold in a given sample will be filtered out. Note: GPM/TPM sum of all taxa in each sample is 1e6.'
            ),
        'mclink_sample_topn_taxon': (
            ['Params_pmr_{mags}_sample', '06_mclink_taxon', 'TopN'], 100, int,
            '06_mclink_taxon: Number of top taxa (ranked by mean KO abundance) retained for calculating taxon-level metabolic pathway/module completeness or abundance per sample. Set to 0 to disable filtering.'
        ),
        'mclink_sample_completeness_per_taxon': (
            ['Params_pmr_{mags}_sample', '06_mclink_taxon', 'completeness_threshold'], 0.5, float,
            '06_mclink_taxon: Completeness threshold [0-1] of modules for calculating taxon-level metabolic pathway/module abundance. Higher values retain only more complete profiles.'),
        'mclink_sample_plus_scale_method_per_taxon': (
            ['Params_pmr_{mags}_sample', '06_mclink_taxon', 'plus_scale_method'], 'mean', str,
            "06_mclink_taxon: Aggregation method for '+'-separated KOs in a module [mean/max/min]. 'mean': average abundance across all required KOs; 'max': highest KO abundance; 'min': lowest KO abundance."),
        'mclink_sample_comma_scale_method_per_taxon': (
            ['Params_pmr_{mags}_sample', '06_mclink_taxon', 'comma_scale_method'], 'sum', str,
            "06_mclink_taxon: Aggregation method for ','-separated alternative KOs in a module [sum/max]. 'sum': cumulative abundance of all alternative KOs; 'max': abundance of the highest single KO among alternatives. Note: this applies to abundance calculation only; module completeness is always computed with 'max'."),
        'barplot_topn_taxon': (
            ['Params_pmr_{mags}_project', '01_Taxa_visualize', 'TopN'], 10, int,
            '01_Taxa_visualize: Number of top taxa to display in stacked barplots. Remaining categories are merged into Others. Maximum 29.'
        ),
        'barplot_topn_fun': (
            ['Params_pmr_{mags}_project', '03_Fun_visualize', 'TopN'], 10, int,
            '03_Fun_visualize: Number of top functional categories to display in stacked barplots. Remaining categories are merged into Others. Maximum 29.'
        ),
        'mclink_min_ko_num_per_sample': (
            ['Params_pmr_{mags}_project', '03_mclink_sample', 'min_sample_num_per_KO'], 1, int,
            '03_mclink_sample: Minimum number of samples in which a KO must be detected for calculating community-level metabolic pathway/module completeness or abundance.'
        ),
        'mclink_min_ko_prevalence_per_sample': (
            ['Params_pmr_{mags}_project', '03_mclink_sample', 'preva'], 0, float,
            '03_mclink_sample: Minimum detection rate [0-1] of a KO across all samples for calculating community-level metabolic pathway/module completeness or abundance.'
        ),
        'mclink_min_ko_abundance_per_sample': (
            ['Params_pmr_{mags}_project', '03_mclink_sample', 'abund'], 0, float,
            '03_mclink_sample: Minimum mean relative abundance [0-1e6] of a KO across all samples for calculating community-level metabolic pathway/module completeness or abundance. KOs with mean abundance below this threshold are excluded. Note: GPM/TPM sum of a sample is 1e6.'
        ),
        'mclink_completeness_per_sample': (
            ['Params_pmr_{mags}_project', '03_mclink_sample', 'completeness_threshold'], 0.5, float,
            '03_mclink_sample: Completeness threshold [0-1] of modules for calculating community-level metabolic pathway/module abundance. Higher values retain only more complete profiles.'
        ),
        'mclink_plus_scale_method_per_sample': (
            ['Params_pmr_{mags}_project', '03_mclink_sample', 'plus_scale_method'], 'mean', str,
            "03_mclink_sample: Aggregation method for '+'-separated KOs in a module [mean/max/min]. 'mean': average abundance across all required KOs; 'max': highest KO abundance; 'min': lowest KO abundance."
        ),
        'mclink_comma_scale_method_per_sample': (
            ['Params_pmr_{mags}_project', '03_mclink_sample', 'comma_scale_method'], 'sum', str,
            "03_mclink_sample: Aggregation method for ','-separated alternative KOs in a module [sum/max]. 'sum': cumulative abundance of all alternative KOs; 'max': abundance of the highest single KO among alternatives. Note: this applies to abundance calculation only; module completeness is always computed with 'max'."
        ),
        'merge_min_taxa_ko_prevalence_across_samples': (
            ['Params_pmr_{mags}_project', '05_merge_tax_KO_depth', 'preva'], 0.10, float,
            '05_merge_tax_KO_depth: Minimum prevalence [0-1] of a Taxon-KO across samples. Taxa-KOs detected in fewer samples than this proportion are excluded.'
        ),
        'merge_min_taxa_ko_abundance_across_samples': (
            ['Params_pmr_{mags}_project', '05_merge_tax_KO_depth', 'abund'], 0.0, float,
            '05_merge_tax_KO_depth: Minimum mean relative abundance [0-1e6] of a Taxon-KO across samples. Taxa-KOs with mean abundance below this threshold are excluded. Note: GPM/TPM sum of a sample is 1e6.'
        ),
        'merge_min_taxa_module_prevalence_across_samples': (
            ['Params_pmr_{mags}_project', '05_merge_tax_Module_RelativeAbundance', 'preva'], 0.10, float,
            '05_merge_tax_Module_RelativeAbundance: Minimum prevalence [0-1] of a Taxon-Module across samples. Taxa-Modules detected in fewer samples than this proportion are excluded.'
        ),
        'merge_min_taxa_module_abundance_across_samples': (
            ['Params_pmr_{mags}_project', '05_merge_tax_Module_RelativeAbundance', 'abund'], 0.0, float,
            '05_merge_tax_Module_RelativeAbundance: Minimum mean relative abundance [0-1e6] of a Taxon-Module across samples. Taxa-Modules with mean abundance below this threshold are excluded. Note: GPM/TPM sum of a sample is 1e6.'
        ),
        'merge_min_taxa_module_sc_prevalence_across_samples': (
            ['Params_pmr_{mags}_project', '05_merge_tax_Module_SumCoverage', 'preva'], 0.10, float,
            '05_merge_tax_Module_SumCoverage: Minimum prevalence [0-1] of a Taxon-Module across samples. Taxa-Modules detected in fewer samples than this proportion are excluded.'
        ),
        'merge_min_taxa_module_sc_abundance_across_samples': (
            ['Params_pmr_{mags}_project', '05_merge_tax_Module_SumCoverage', 'abund'], 0.0, float,
            '05_merge_tax_Module_SumCoverage: Minimum mean sum coverage [0-1e6] of a Taxon-Module across samples. Taxa-Modules with mean sum coverage below this threshold are excluded. Note: GPM/TPM sum of a sample is 1e6.'
        ),
        'mclink_project_min_ko_num_per_taxon': (
            ['Params_pmr_{mags}_project', '06_mclink_taxon', 'min_KO_num_in_one_taxon'], 100, int,
            '06_mclink_taxon: Minimum number of KOs per taxon for calculating taxon-level metabolic pathway/module completeness on the project-averaged taxon-function matrix. Taxa with fewer than this number of detected KOs in the project-averaged matrix will be filtered out. KO abundances are first averaged across all samples to generate a single project-level matrix before calculation.'
        ),
        'mclink_project_min_ko_abundance_per_taxon': (
            ['Params_pmr_{mags}_project', '06_mclink_taxon', 'min_mean_KO_RA_in_one_taxon'], 0, float,
            '06_mclink_taxon: Minimum mean KO abundance [0-1e6] per taxon for calculating taxon-level metabolic pathway/module completeness on the project-averaged taxon-function matrix. Taxa with mean KO abundance below this threshold in the project-averaged matrix will be filtered out. Note: GPM/TPM sum of all taxa is 1e6 per sample; values are averaged across samples prior to filtering and calculation.'
        ),
        'mclink_project_topn_taxon': (
            ['Params_pmr_{mags}_project', '06_mclink_taxon', 'TopN'], 100, int,
            '06_mclink_taxon: Number of top taxa (ranked by mean KO abundance in the project-averaged taxon-function matrix) retained for calculating taxon-level metabolic pathway/module completeness. Set to 0 to disable filtering.'
        ),
        'mclink_project_completeness_per_taxon': (
            ['Params_pmr_{mags}_project', '06_mclink_taxon', 'completeness_threshold'], 0.5, float,
            '06_mclink_taxon: Completeness threshold [0-1] of modules for calculating taxon-level metabolic pathway/module abundance on the project-averaged taxon-function matrix. Higher values retain only more complete profiles.'),
        'mclink_project_plus_scale_method_per_taxon': (
            ['Params_pmr_{mags}_project', '06_mclink_taxon', 'plus_scale_method'], 'mean', str,
            "06_mclink_taxon: Aggregation method for '+'-separated KOs in a module [mean/max/min] on the project-averaged taxon-function matrix. 'mean': average abundance across all required KOs; 'max': highest KO abundance; 'min': lowest KO abundance."),
        'mclink_project_comma_scale_method_per_taxon': (
            ['Params_pmr_{mags}_project', '06_mclink_taxon', 'comma_scale_method'], 'sum', str,
            "06_mclink_taxon: Aggregation method for ','-separated alternative KOs in a module [sum/max] on the project-averaged taxon-function matrix. 'sum': cumulative abundance of all alternative KOs; 'max': abundance of the highest single KO among alternatives. Note: this applies to abundance calculation only; module completeness is always computed with 'max'."),
    },
    'v': {
        'raw_stats_sample_threads': (
            ['Params_b_sample', '00_raw_stats', 'threads'], 4, int,
            '00_raw_stats: Threads for seqkit.'),
        'trim_sample_threads': (
            ['Params_b_sample', '01_trim', 'threads'], 16, int,
            '01_trim: Threads for fastp (next-generation short reads) or fastplong (third-generation long reads).'),
        'trim_length_required': (
            ['Params_b_sample', '01_trim', 'length_required'], 15, int,
            '01_trim: Reads shorter than length_required will be discarded'),
        'trim_n_base_limit': (
            ['Params_b_sample', '01_trim', 'n_base_limit'], 5, int,
            "01_trim: if one read's number of N base is >n_base_limit, then this read/pair is discarded."),
        'trim_compression': (
            ['Params_b_sample', '01_trim', 'compression'], 4, int,
            '01_trim: Compression level for gzip output [1-9]. 1 is fastest, 9 is smallest.'),
        'trim_extra_args': (
            ['Params_b_sample', '01_trim', 'extra_args'], '', str,
            '01_trim: Extra arguments forwarded to fastp or fastplong for raw reads trimming (e.g., "--reads_to_process 0 -Q").'),
        'trim_stats_sample_threads': (
            ['Params_b_sample', '01_trim_stats', 'threads'], 16, int,
            '01_trim_stats: Threads.'),
        'assembly_method': (
            ['Params_b_sample', '02_assembly', 'method'], 'megahit', str,
            '\b\n02_assembly: Underlying assembly software using one of [megahit/spades/metaflye].\n'
            'Must be compatible with --read-type (config key read_type):\n'
            'megahit: MEGAHIT (second-generation sequencing; for paired or single DNA reads).\n'
            'spades: SPAdes (second-generation sequencing; metaSPAdes for paired DNA reads, rnaSPAdes for paired-rna or single-rna RNA reads).\n'
            'metaflye: metaFlye (third-generation sequencing; required for pacbio-raw/corr/hifi or nano-raw/corr/hq long reads).'),
        'assembly_threads': (
            ['Params_b_sample', '02_assembly', 'threads'], 40, int,
            '02_assembly: Threads for [megahit/spades/metaflye].'),
        'megahit_presets': (
            ['Params_b_sample', '02_megahit', 'presets'], 'meta-sensitive', str,
            '\b\n02_megahit: MEGAHIT presets using one of [meta-sensitive/meta-large].\n'
            "meta-sensitive: '--min-count 1 --k-list 21,29,39,49,...,129,141'.\n"
            "meta-large: '--k-min 27 --k-max 127 --k-step 10'.\n"
            '(large & complex metagenomes, like soil).'),
        'megahit_k_list': (
            ['Params_b_sample', '02_megahit', 'k_list'], '', str,
            '02_megahit: MEGAHIT --k-list (optional). Explicit comma-separated k-mer sizes (e.g. "21,29,39,59,79,99,119,141"). Overrides --megahit-presets. If set, --megahit-k-min, --megahit-k-max, --megahit-k-step must NOT be set. Optional --megahit-min-count can be used together.'),
        'megahit_min_count': (
            ['Params_b_sample', '02_megahit', 'min_count'], 2, int,
            '02_megahit: MEGAHIT --min-count (optional). Minimum multiplicity for filtering (k_min+1)-mers. Can be used with --megahit-k-list or --megahit-k-min/--k-max/--k-step.'),
        'megahit_k_min': (
            ['Params_b_sample', '02_megahit', 'k_min'], '', str,
            '02_megahit: MEGAHIT --k-min (optional). Minimum kmer size (<= 255), must be odd number. If set with --k-max and --k-step, overrides --megahit-presets and uses custom k-mer mode.'),
        'megahit_k_max': (
            ['Params_b_sample', '02_megahit', 'k_max'], '', str,
            '02_megahit: MEGAHIT --k-max (optional). Maximum kmer size (<= 255), must be odd number.'),
        'megahit_k_step': (
            ['Params_b_sample', '02_megahit', 'k_step'], '', str,
            '02_megahit: MEGAHIT --k-step (optional). Increment of kmer size of each iteration (<= 28), must be even number.'),
        'megahit_extra_args': (
            ['Params_b_sample', '02_megahit', 'extra_args'], '', str,
            '02_megahit: Extra arguments forwarded to MEGAHIT (e.g. "--no-mercy --tmp-dir /path/to/your/tmp/dir"). Do NOT use this for k-mer settings (--k-list, --k-min, --presets, etc.).'),
        'spades_extra_args': (
            ['Params_b_sample', '02_spades', 'extra_args'], '', str,
            '02_spades: Extra arguments forwarded to SPAdes. Do NOT use this for input/output or thread settings.'),
        'metaflye_extra_args': (
            ['Params_b_sample', '02_metaflye', 'extra_args'], '', str,
            '02_metaflye: Extra arguments forwarded to Flye for third-generation long reads assembly (e.g., "--iterations 2"). Do NOT use this for read type or input/output settings.'),
        'mapping_threads': (
            ['Params_b_sample', '03_mapping', 'threads'], 20, int,
            '03_mapping: Threads for [minibwa/bwamem/bwamem2/strobealign].'),
        'mapping_method': (
            ['Params_b_sample', '03_mapping', 'method'], 'minibwa', str,
            '\b\n03_mapping: Underlying mapping software used [minibwa/bwamem/bwamem2/strobealign].\n'
            'minibwa: minibwa using default parameters.\n'
            'bwamem: bwa mem using default parameters.\n'
            'bwamem2: bwa-mem2 using default parameters.\n'
            'strobealign: strobealign using default parameters.'),
        'mapping_samtools_sort_mem': (
            ['Params_b_sample', '03_mapping', 'samtools_sort_mem'], '768M', str,
            '03_mapping: Set maximum memory per thread; suffix K/M/G recognized, e.g. 2G, 768M.'),
        'mapping_extra_args': (
            ['Params_b_sample', '03_mapping', 'extra_args'], '', str,
            '03_mapping: Extra arguments forwarded to the active mapping aligner (minibwa map, bwa mem, bwa-mem2 mem, or strobealign --use-index), (e.g. "-B 8" for minibwa map).'),
        'contigs_depth_sample_threads': (
            ['Params_b_sample', '07_contigs_depth', 'threads'], 8, int,
            '07_contigs_depth: Threads.'),
        'depth_min_read_percent_identity': (
            ['Params_b_sample', '07_contigs_depth', 'min_read_percent_identity'], 95, float,
            '\b\n07_contigs_depth: Exclude reads by overall percent identity [0-100] for contig depth estimation\n'
            'using jgi_summarize_bam_contig_depths --percentIdentity, e.g. 95 for 95%.\n'
            'Reads with percent identity below this threshold are excluded from depth calculation.\n'
            'Higher values increase stringency, retaining only high-quality alignments.'),
        'depth_contig_end_exclusion': (
            ['Params_b_sample', '07_contigs_depth', 'contig_end_exclusion'], 75, int,
            '\b\n07_contigs_depth: Exclude bases at the ends of contig sequences for depth estimation using jgi_summarize_bam_contig_depths --maxEdgeBases.\n'
            'Bases within this distance from either end of a contig are ignored when calculating depth,\n'
            'preventing edge effects from partial alignments.'),
        'depth_min_contig_length': (
            ['Params_b_sample', '07_contigs_depth', 'min_contig_length'], 300, int,
            '\b\n07_contigs_depth: Minimum contig length for depth estimation [bp], e.g. 300 for 300 bp.\n'
            'Contigs shorter than this threshold are excluded from depth calculation.'),
        'contigs_prodigal_extra_args': (
            ['Params_b_sample', '07_contigs_depth', 'extra_args'], '', str,
            '07_contigs_depth: Extra arguments forwarded to prodigal for contig gene prediction (e.g., "-m -n").'),
        'raw_stats_project_threads': (
            ['Params_b_project', '01_cat_raw_stats', 'threads'], 4, int,
            '01_cat_raw_stats: Threads.'),
        'trimmed_stats_project_threads': (
            ['Params_b_project', '02_cat_trimmed_stats', 'threads'], 4, int,
            '02_cat_trimmed_stats: Threads.'),
        'contigs_stats_project_threads': (
            ['Params_b_project', '03_contigs_stats', 'threads'], 4, int,
            '03_contigs_stats: Threads.'),
        'merge_stats_project_threads': (
            ['Params_b_project', '04_merge_reads_contigs_stats', 'threads'], 4, int,
            '04_merge_reads_contigs_stats: Threads.'),
        'virus_identify_method': (
            ['Params_v_sample', '01_virus_identify', 'identify_method'], 'genomad_vs_vb', str,
            '\b\n01_virus_identify: Virus identification method, should be one of the [genomad/vs_vb/genomad_vs_vb].\n'
            'genomad: Use geNomad;\n'
            'vs_vb: Use VirSorter2 and CheckV, then use VIBRANT to check the low quality viruses from VirSorter2;\n'
            'genomad_vs_vb: Combine genomad and vs_vb virus to get the union results.'
        ),
        'virus_length_limit': (
            ['Params_v_sample', '01_virus_identify', 'virus_length_limit'], 2000, int,
            '01_virus_identify: Minimum length (bp) for viral contigs.'
        ),
        'virus_identify_threads': (
            ['Params_v_sample', '01_virus_identify', 'threads'], 8, int,
            '01_virus_identify: Threads.'),
        'virus_contig_taxonomy_cat_range': (
            ['Params_v_sample', '03_taxonomy_CAT_pack', 'range'], 10, float,
            '\b\n03_taxonomy_CAT_pack: Percentage range [0-100] of DIAMOND hits relative to top bit-score to include in LCA calculation.\n'
            'Higher values yield more conservative classifications by incorporating divergent hits, while lower values enable specific assignments.'
        ),
        'virus_contig_taxonomy_cat_fraction': (
            ['Params_v_sample', '03_taxonomy_CAT_pack', 'fraction'], 0.50, float,
            '\b\n03_taxonomy_CAT_pack: Minimum bit-score support fraction [0-0.99].\n'
            'Higher values enforce stricter consensus (one classification per sequence), while lower values allow multiple lineage detections.'
        ),
        'virus_contig_taxonomy_cat_pack_threads': (
            ['Params_v_sample', '03_taxonomy_CAT_pack', 'threads'], 8, int,
            '03_taxonomy_CAT_pack: Threads.'),
        'virus_contig_taxonomy_cat_pack_extra_args': (
            ['Params_v_sample', '03_taxonomy_CAT_pack', 'extra_args'], '', str,
            '03_taxonomy_CAT_pack: Extra arguments forwarded to CAT_pack contigs (e.g., "--block_size 12.0").'),
        'virus_mapping_threads': (
            ['Params_v_sample', '03_mapping', 'threads'], 20, int,
            '03_mapping: Threads for [minibwa/bwamem/bwamem2/strobealign].'),
        'virus_mapping_method': (
            ['Params_v_sample', '03_mapping', 'method'], 'minibwa', str,
            '\b\n03_mapping: Underlying mapping software used [minibwa/bwamem/bwamem2/strobealign].\n'
            'minibwa: minibwa using default parameters.\n'
            'bwamem: bwa mem using default parameters.\n'
            'bwamem2: bwa-mem2 using default parameters.\n'
            'strobealign: strobealign using default parameters.'),
        'virus_mapping_samtools_sort_mem': (
            ['Params_v_sample', '03_mapping', 'samtools_sort_mem'], '768M', str,
            '03_mapping: Set maximum memory per thread; suffix K/M/G recognized, e.g. 2G, 768M.'),
        'virus_mapping_extra_args': (
            ['Params_v_sample', '03_mapping', 'extra_args'], '', str,
            '03_mapping: Extra arguments forwarded to the active mapping aligner (minibwa map, bwa mem, bwa-mem2 mem, or strobealign --use-index), (e.g. "-B 8" for minibwa map).'),
        'virus_binning_threads': (
            ['Params_v_sample', '04_binning', 'threads'], 8, int,
            '04_binning_vRhyme: Threads.'),
        'virus_checkv_threads': (
            ['Params_v_sample', '05_checkv', 'threads'], 8, int,
            '05_checkv: Threads.'),
        'virus_vclust_ani_sample': (
            ['Params_v_sample', '06_vclust', 'ani'], 0.95, float,
            '06_vclust: ANI threshold [0-1] of vclust for sample-level virus clustering.'
        ),
        'virus_vclust_threads_sample': (
            ['Params_v_sample', '06_vclust', 'threads'], 8, int,
            '06_vclust: Threads.'),
        'vclust_sample_prefilter_extra_args': (
            ['Params_v_sample', '06_vclust', 'prefilter_extra_args'], '', str,
            '06_vclust: Extra arguments forwarded to vclust prefilter (e.g., "-k 25").'),
        'vclust_sample_align_extra_args': (
            ['Params_v_sample', '06_vclust', 'align_extra_args'], '', str,
            '06_vclust: Extra arguments forwarded to vclust align (e.g., "--mal 11").'),
        'vclust_sample_cluster_extra_args': (
            ['Params_v_sample', '06_vclust', 'cluster_extra_args'], '', str,
            '06_vclust: Extra arguments forwarded to vclust cluster (e.g., "--len_ratio 0").'),
        'virus_iphop_threads_sample': (
            ['Params_v_sample', '07_iPHoP', 'threads'], 20, int,
            '07_iPHoP: Threads.'),
        'virus_vmags_taxonomy_bat_range': (
            ['Params_v_sample', '08_taxonomy_BAT', 'range'], 10, float,
            '\b\n08_taxonomy_BAT: Percentage range [0-100] of DIAMOND hits relative to top bit-score to include in LCA calculation.\n'
            'Higher values yield more conservative classifications by incorporating divergent hits, while lower values enable specific assignments.'
        ),
        'virus_vmags_taxonomy_bat_fraction': (
            ['Params_v_sample', '08_taxonomy_BAT', 'fraction'], 0.50, float,
            '\b\n08_taxonomy_BAT: Minimum bit-score support fraction [0-0.99].\n'
            'Higher values enforce stricter consensus (one classification per sequence), while lower values allow multiple lineage detections.'
        ),
        'virus_vmags_taxonomy_bat_threads': (
            ['Params_v_sample', '08_taxonomy_BAT', 'threads'], 8, int,
            '08_taxonomy_BAT: Threads.'),
        'virus_vmags_taxonomy_bat_extra_args': (
            ['Params_v_sample', '08_taxonomy_BAT', 'extra_args'], '', str,
            '08_taxonomy_BAT: Extra arguments forwarded to CAT_pack bins (e.g., "--block_size 12.0").'),
        'virus_coverm_method_sample': (
            ['Params_v_sample', '09_coverm_eachsample_vMAGs', 'method'], 'trimmed_mean', str,
            '\b\n09_coverm_eachsample_vMAGs: Method(s) for coverm-based coverage calculation, must be one of [trimmed_mean/mean].\n'
            'A more detailed description of the different methods is available at https://github.com/wwood/CoverM#calculation-methods'
        ),
        'virus_coverm_min_read_percent_identity_sample': (
            ['Params_v_sample', '09_coverm_eachsample_vMAGs', 'min_read_percent_identity'], 95, float,
            '09_coverm_eachsample_vMAGs: Exclude reads by overall percent identity [0-100], e.g. 95 for 95%.'
        ),
        'virus_coverm_min_read_aligned_length_sample': (
            ['Params_v_sample', '09_coverm_eachsample_vMAGs', 'min_read_aligned_length'], 45, int,
            '09_coverm_eachsample_vMAGs: Exclude reads with fewer aligned bases.'
        ),
        'virus_coverm_contig_end_exclusion_sample': (
            ['Params_v_sample', '09_coverm_eachsample_vMAGs', 'contig_end_exclusion'], 75, int,
            '09_coverm_eachsample_vMAGs: Exclude bases at the ends of contig sequences for coverm calculation.'
        ),
        'virus_coverm_threads_sample': (
            ['Params_v_sample', '09_coverm_eachsample_vMAGs', 'threads'], 8, int,
            '09_coverm_eachsample_vMAGs: Threads.'),
        'virus_coverm_sample_extra_args': (
            ['Params_v_sample', '09_coverm_eachsample_vMAGs', 'extra_args'], '', str,
            '09_coverm_eachsample_vMAGs: Extra arguments forwarded to coverm genome (e.g., "--quiet").'),
        'virus_vclust_ani_project': (
            ['Params_v_project', '06_vclust', 'ani'], 0.95, float,
            '06_vclust: ANI threshold [0-1] of vclust for project-level virus clustering.'
        ),
        'virus_vclust_threads_project': (
            ['Params_v_project', '06_vclust', 'threads'], 8, int,
            '06_vclust: Threads.'),
        'vclust_project_prefilter_extra_args': (
            ['Params_v_project', '06_vclust', 'prefilter_extra_args'], '', str,
            '06_vclust: Extra arguments forwarded to vclust prefilter (e.g., "-k 25".'),
        'vclust_project_align_extra_args': (
            ['Params_v_project', '06_vclust', 'align_extra_args'], '', str,
            '06_vclust: Extra arguments forwarded to vclust align (e.g., "--mal 11").'),
        'vclust_project_cluster_extra_args': (
            ['Params_v_project', '06_vclust', 'cluster_extra_args'], '', str,
            '06_vclust: Extra arguments forwarded to vclust cluster (e.g., "--len_ratio 0").'),
        'virus_iphop_custom_mags_dir_project': (
            ['Params_v_project', '07_iPHoP', 'custom_MAGs_dir'], 'none', str,
            '\b\n07_iPHoP: Custom mMAGs directory for iPHoP host prediction, using one of ["none"/user-defined-path].\n'
            'Use "none" to disable custom mMAGs.\n'
            'Use an absolute path to a custom MAGs directory (including *.fna files).\n'
            'NOTE: "auto_drep99" and "auto_drep95" (automatically using dereplicated mMAGs with metalink default path) are only effective in vm mode.\n'
            'This parameter is not applicable in vc/vcr/vmc/vmcr modes (iPHoP host prediction is skipped).'
        ),
        'virus_iphop_threads_project': (
            ['Params_v_project', '07_iPHoP', 'threads'], 60, int,
            '07_iPHoP: Threads.'),
    },
    'vmcr': {
        'raw_stats_sample_threads': (
            ['Params_b_sample', '00_raw_stats', 'threads'], 4, int,
            '00_raw_stats: Threads for seqkit.'),
        'trim_sample_threads': (
            ['Params_b_sample', '01_trim', 'threads'], 16, int,
            '01_trim: Threads for fastp (next-generation short reads) or fastplong (third-generation long reads).'),
        'trim_length_required': (
            ['Params_b_sample', '01_trim', 'length_required'], 15, int,
            '01_trim: Reads shorter than length_required will be discarded'),
        'trim_n_base_limit': (
            ['Params_b_sample', '01_trim', 'n_base_limit'], 5, int,
            "01_trim: if one read's number of N base is >n_base_limit, then this read/pair is discarded."),
        'trim_compression': (
            ['Params_b_sample', '01_trim', 'compression'], 4, int,
            '01_trim: Compression level for gzip output [1-9]. 1 is fastest, 9 is smallest.'),
        'trim_extra_args': (
            ['Params_b_sample', '01_trim', 'extra_args'], '', str,
            '01_trim: Extra arguments forwarded to fastp or fastplong for raw reads trimming (e.g., "--reads_to_process 0 -Q").'),
        'trim_stats_sample_threads': (
            ['Params_b_sample', '01_trim_stats', 'threads'], 16, int,
            '01_trim_stats: Threads.'),
        'assembly_method': (
            ['Params_b_sample', '02_assembly', 'method'], 'megahit', str,
            '\b\n02_assembly: Underlying assembly software using one of [megahit/spades/metaflye].\n'
            'Must be compatible with --read-type (config key read_type):\n'
            'megahit: MEGAHIT (second-generation sequencing; for paired or single DNA reads).\n'
            'spades: SPAdes (second-generation sequencing; metaSPAdes for paired DNA reads, rnaSPAdes for paired-rna or single-rna RNA reads).\n'
            'metaflye: metaFlye (third-generation sequencing; required for pacbio-raw/corr/hifi or nano-raw/corr/hq long reads).'),
        'assembly_threads': (
            ['Params_b_sample', '02_assembly', 'threads'], 40, int,
            '02_assembly: Threads for [megahit/spades/metaflye].'),
        'megahit_presets': (
            ['Params_b_sample', '02_megahit', 'presets'], 'meta-sensitive', str,
            '\b\n02_megahit: MEGAHIT presets using one of [meta-sensitive/meta-large].\n'
            "meta-sensitive: '--min-count 1 --k-list 21,29,39,49,...,129,141'.\n"
            "meta-large: '--k-min 27 --k-max 127 --k-step 10'.\n"
            '(large & complex metagenomes, like soil).'),
        'megahit_k_list': (
            ['Params_b_sample', '02_megahit', 'k_list'], '', str,
            '02_megahit: MEGAHIT --k-list (optional). Explicit comma-separated k-mer sizes (e.g. "21,29,39,59,79,99,119,141"). Overrides --megahit-presets. If set, --megahit-k-min, --megahit-k-max, --megahit-k-step must NOT be set. Optional --megahit-min-count can be used together.'),
        'megahit_min_count': (
            ['Params_b_sample', '02_megahit', 'min_count'], 2, int,
            '02_megahit: MEGAHIT --min-count (optional). Minimum multiplicity for filtering (k_min+1)-mers. Can be used with --megahit-k-list or --megahit-k-min/--k-max/--k-step.'),
        'megahit_k_min': (
            ['Params_b_sample', '02_megahit', 'k_min'], '', str,
            '02_megahit: MEGAHIT --k-min (optional). Minimum kmer size (<= 255), must be odd number. If set with --k-max and --k-step, overrides --megahit-presets and uses custom k-mer mode.'),
        'megahit_k_max': (
            ['Params_b_sample', '02_megahit', 'k_max'], '', str,
            '02_megahit: MEGAHIT --k-max (optional). Maximum kmer size (<= 255), must be odd number.'),
        'megahit_k_step': (
            ['Params_b_sample', '02_megahit', 'k_step'], '', str,
            '02_megahit: MEGAHIT --k-step (optional). Increment of kmer size of each iteration (<= 28), must be even number.'),
        'megahit_extra_args': (
            ['Params_b_sample', '02_megahit', 'extra_args'], '', str,
            '02_megahit: Extra arguments forwarded to MEGAHIT (e.g. "--no-mercy --tmp-dir /path/to/your/tmp/dir"). Do NOT use this for k-mer settings (--k-list, --k-min, --presets, etc.).'),
        'spades_extra_args': (
            ['Params_b_sample', '02_spades', 'extra_args'], '', str,
            '02_spades: Extra arguments forwarded to SPAdes. Do NOT use this for input/output or thread settings.'),
        'metaflye_extra_args': (
            ['Params_b_sample', '02_metaflye', 'extra_args'], '', str,
            '02_metaflye: Extra arguments forwarded to Flye for third-generation long reads assembly (e.g., "--iterations 2"). Do NOT use this for read type or input/output settings.'),
        'mapping_threads': (
            ['Params_b_sample', '03_mapping', 'threads'], 20, int,
            '03_mapping: Threads for [minibwa/bwamem/bwamem2/strobealign].'),
        'mapping_method': (
            ['Params_b_sample', '03_mapping', 'method'], 'minibwa', str,
            '\b\n03_mapping: Underlying mapping software used [minibwa/bwamem/bwamem2/strobealign].\n'
            'minibwa: minibwa using default parameters.\n'
            'bwamem: bwa mem using default parameters.\n'
            'bwamem2: bwa-mem2 using default parameters.\n'
            'strobealign: strobealign using default parameters.'),
        'mapping_samtools_sort_mem': (
            ['Params_b_sample', '03_mapping', 'samtools_sort_mem'], '768M', str,
            '03_mapping: Set maximum memory per thread; suffix K/M/G recognized, e.g. 2G, 768M.'),
        'mapping_extra_args': (
            ['Params_b_sample', '03_mapping', 'extra_args'], '', str,
            '03_mapping: Extra arguments forwarded to the active mapping aligner (minibwa map, bwa mem, bwa-mem2 mem, or strobealign --use-index), (e.g. "-B 8" for minibwa map).'),
        'binning_metawrap_depth_threads': (
            ['Params_b_sample', '05_metawrap_binning_depth', 'threads'], 8, int,
            '05_metawrap_binning_depth: Threads.'),
        'binning_min_read_percent_identity': (
            ['Params_b_sample', '05_metawrap_binning_depth', 'min_read_percent_identity'], 97, float,
            '05_metawrap_binning_depth: Exclude reads by overall percent identity [0-100] during the binning process, e.g. 97 for 97%.'),
        'binning_contig_end_exclusion': (
            ['Params_b_sample', '05_metawrap_binning_depth', 'contig_end_exclusion'], 75, int,
            '05_metawrap_binning_depth: Exclude bases at the ends of contig sequences used for depth calculation and bining.'),
        'binning_min_contig_length': (
            ['Params_b_sample', '05_metawrap_binning_depth', 'min_contig_length'], 1500, int,
            '05_metawrap_binning_depth: Minimum contig length for binning.'),
        'binning_tools': (
            ['Params_b_sample', '05_metawrap_binning_depth', 'binning_tools'], 'metabat2,maxbin2,concoct', str,
            '\b\n05_metawrap_binning_depth: Comma-separated list of binning tools, choose 2 or 3 from [metabat2,maxbin2,concoct].\n'
            'At least 2 tools recommended for bin_refinement.'),
        'binning_max_memory_gb': (
            ['Params_b_sample', '05_metawrap_binning_depth', 'max_memory_gb'], 512, int,
            '05_metawrap_binning_depth: Maximum memory (GB) for MetaWRAP binning.'),
        'binning_metawrap_extra_args': (
            ['Params_b_sample', '05_metawrap_binning_depth', 'extra_args'], '', str,
            '05_metawrap_binning_depth: Extra arguments forwarded to metaWRAP binning (e.g., "-m 5").'),
        'metawrap_refinement_max_memory_gb': (
            ['Params_b_sample', '06_metawrap_bin_refinement', 'max_memory_gb'], 512, int,
            '06_metawrap_bin_refinement: Maximum memory in GB for MetaWRAP bin_refinement.'),
        'metawrap_refinement_extra_args': (
            ['Params_b_sample', '06_metawrap_bin_refinement', 'extra_args'], '', str,
            '06_metawrap_bin_refinement: Extra arguments forwarded to metaWRAP bin_refinement (e.g., "--remove-ambiguous").'),
        'bin_checkm_threads': (
            ['Params_b_sample', '06_bin_checkm', 'threads'], 8, int,
            '06_bin_checkm or 06_metawrap_bin_refinement: Threads.'),
        'bin_min_completeness': (
            ['Params_b_sample', '06_bin_checkm', 'min_completeness'], 50, float,
            '06_bin_checkm or 06_metawrap_bin_refinement: Minimum completeness (%) bin filtering, e.g. 50 for 50%.'),
        'bin_max_contamination': (
            ['Params_b_sample', '06_bin_checkm', 'max_contamination'], 10, float,
            '06_bin_checkm or 06_metawrap_bin_refinement: Maximum contamination (%) bin filtering, e.g. 10 for 10%.'),
        'bin_min_quality_score': (
            ['Params_b_sample', '06_bin_checkm', 'min_quality_score'], 50, float,
            '06_bin_checkm: Minimum quality score for bin filtering (GTDB-like), QS = completeness - 5*contamination, e.g. 50 for QS>=50.'),
        'contigs_depth_sample_threads': (
            ['Params_b_sample', '07_contigs_depth', 'threads'], 8, int,
            '07_contigs_depth: Threads.'),
        'depth_min_read_percent_identity': (
            ['Params_b_sample', '07_contigs_depth', 'min_read_percent_identity'], 95, float,
            '\b\n07_contigs_depth: Exclude reads by overall percent identity [0-100] for contig depth estimation\n'
            'using jgi_summarize_bam_contig_depths --percentIdentity, e.g. 95 for 95%.\n'
            'Reads with percent identity below this threshold are excluded from depth calculation.\n'
            'Higher values increase stringency, retaining only high-quality alignments.'),
        'depth_contig_end_exclusion': (
            ['Params_b_sample', '07_contigs_depth', 'contig_end_exclusion'], 75, int,
            '\b\n07_contigs_depth: Exclude bases at the ends of contig sequences for depth estimation using jgi_summarize_bam_contig_depths --maxEdgeBases.\n'
            'Bases within this distance from either end of a contig are ignored when calculating depth,\n'
            'preventing edge effects from partial alignments.'),
        'depth_min_contig_length': (
            ['Params_b_sample', '07_contigs_depth', 'min_contig_length'], 300, int,
            '\b\n07_contigs_depth: Minimum contig length for depth estimation [bp], e.g. 300 for 300 bp.\n'
            'Contigs shorter than this threshold are excluded from depth calculation.'),
        'contigs_prodigal_extra_args': (
            ['Params_b_sample', '07_contigs_depth', 'extra_args'], '', str,
            '07_contigs_depth: Extra arguments forwarded to prodigal for contig gene prediction (e.g., "-m -n").'),
        'bin_prodigal_threads': (
            ['Params_b_sample', '08_bin_prodigal', 'threads'], 8, int,
            '08_bin_prodigal: Threads.'),
        'bin_prodigal_extra_args': (
            ['Params_b_sample', '08_bin_prodigal', 'extra_args'], '', str,
            '08_bin_prodigal: Extra arguments forwarded to prodigal for MAG gene prediction (e.g., "-m -n").'),
        'raw_stats_project_threads': (
            ['Params_b_project', '01_cat_raw_stats', 'threads'], 4, int,
            '01_cat_raw_stats: Threads.'),
        'trimmed_stats_project_threads': (
            ['Params_b_project', '02_cat_trimmed_stats', 'threads'], 4, int,
            '02_cat_trimmed_stats: Threads.'),
        'contigs_stats_project_threads': (
            ['Params_b_project', '03_contigs_stats', 'threads'], 4, int,
            '03_contigs_stats: Threads.'),
        'merge_stats_project_threads': (
            ['Params_b_project', '04_merge_reads_contigs_stats', 'threads'], 4, int,
            '04_merge_reads_contigs_stats: Threads.'),
        'gather_checkm_project_threads': (
            ['Params_b_project', '05_gather_project_bin_checkm_stats', 'threads'], 4, int,
            '05_gather_project_bin_checkm_stats: Threads.'),
        'virus_identify_method': (
            ['Params_v_sample', '01_virus_identify', 'identify_method'], 'genomad_vs_vb', str,
            '\b\n01_virus_identify: Virus identification method, should be one of the [genomad/vs_vb/genomad_vs_vb].\n'
            'genomad: Use geNomad;\n'
            'vs_vb: Use VirSorter2 and CheckV, then use VIBRANT to check the low quality viruses from VirSorter2;\n'
            'genomad_vs_vb: Combine genomad and vs_vb virus to get the union results.'
        ),
        'virus_length_limit': (
            ['Params_v_sample', '01_virus_identify', 'virus_length_limit'], 2000, int,
            '01_virus_identify: Minimum length (bp) for viral contigs.'
        ),
        'virus_identify_threads': (
            ['Params_v_sample', '01_virus_identify', 'threads'], 8, int,
            '01_virus_identify: Threads.'),
        'virus_contig_taxonomy_cat_range': (
            ['Params_v_sample', '03_taxonomy_CAT_pack', 'range'], 10, float,
            '\b\n03_taxonomy_CAT_pack: Percentage range [0-100] of DIAMOND hits relative to top bit-score to include in LCA calculation.\n'
            'Higher values yield more conservative classifications by incorporating divergent hits, while lower values enable specific assignments.'
        ),
        'virus_contig_taxonomy_cat_fraction': (
            ['Params_v_sample', '03_taxonomy_CAT_pack', 'fraction'], 0.50, float,
            '\b\n03_taxonomy_CAT_pack: Minimum bit-score support fraction [0-0.99].\n'
            'Higher values enforce stricter consensus (one classification per sequence), while lower values allow multiple lineage detections.'
        ),
        'virus_contig_taxonomy_cat_pack_threads': (
            ['Params_v_sample', '03_taxonomy_CAT_pack', 'threads'], 8, int,
            '03_taxonomy_CAT_pack: Threads.'),
        'virus_contig_taxonomy_cat_pack_extra_args': (
            ['Params_v_sample', '03_taxonomy_CAT_pack', 'extra_args'], '', str,
            '03_taxonomy_CAT_pack: Extra arguments forwarded to CAT_pack contigs (e.g., "--block_size 12.0").'),
        'virus_mapping_threads': (
            ['Params_v_sample', '03_mapping', 'threads'], 20, int,
            '03_mapping: Threads for [minibwa/bwamem/bwamem2/strobealign].'),
        'virus_mapping_method': (
            ['Params_v_sample', '03_mapping', 'method'], 'minibwa', str,
            '\b\n03_mapping: Underlying mapping software used [minibwa/bwamem/bwamem2/strobealign].\n'
            'minibwa: minibwa using default parameters.\n'
            'bwamem: bwa mem using default parameters.\n'
            'bwamem2: bwa-mem2 using default parameters.\n'
            'strobealign: strobealign using default parameters.'),
        'virus_mapping_samtools_sort_mem': (
            ['Params_v_sample', '03_mapping', 'samtools_sort_mem'], '768M', str,
            '03_mapping: Set maximum memory per thread; suffix K/M/G recognized, e.g. 2G, 768M.'),
        'virus_mapping_extra_args': (
            ['Params_v_sample', '03_mapping', 'extra_args'], '', str,
            '03_mapping: Extra arguments forwarded to the active mapping aligner (minibwa map, bwa mem, bwa-mem2 mem, or strobealign --use-index), (e.g. "-B 8" for minibwa map).'),
        'virus_binning_threads': (
            ['Params_v_sample', '04_binning', 'threads'], 8, int,
            '04_binning_vRhyme: Threads.'),
        'virus_checkv_threads': (
            ['Params_v_sample', '05_checkv', 'threads'], 8, int,
            '05_checkv: Threads.'),
        'virus_vclust_ani_sample': (
            ['Params_v_sample', '06_vclust', 'ani'], 0.95, float,
            '06_vclust: ANI threshold [0-1] of vclust for sample-level virus clustering.'
        ),
        'virus_vclust_threads_sample': (
            ['Params_v_sample', '06_vclust', 'threads'], 8, int,
            '06_vclust: Threads.'),
        'vclust_sample_prefilter_extra_args': (
            ['Params_v_sample', '06_vclust', 'prefilter_extra_args'], '', str,
            '06_vclust: Extra arguments forwarded to vclust prefilter (e.g., "-k 25").'),
        'vclust_sample_align_extra_args': (
            ['Params_v_sample', '06_vclust', 'align_extra_args'], '', str,
            '06_vclust: Extra arguments forwarded to vclust align (e.g., "--mal 11").'),
        'vclust_sample_cluster_extra_args': (
            ['Params_v_sample', '06_vclust', 'cluster_extra_args'], '', str,
            '06_vclust: Extra arguments forwarded to vclust cluster (e.g., "--len_ratio 0").'),
        'virus_iphop_threads_sample': (
            ['Params_v_sample', '07_iPHoP', 'threads'], 20, int,
            '07_iPHoP: Threads.'),
        'virus_vmags_taxonomy_bat_range': (
            ['Params_v_sample', '08_taxonomy_BAT', 'range'], 10, float,
            '\b\n08_taxonomy_BAT: Percentage range [0-100] of DIAMOND hits relative to top bit-score to include in LCA calculation.\n'
            'Higher values yield more conservative classifications by incorporating divergent hits, while lower values enable specific assignments.'
        ),
        'virus_vmags_taxonomy_bat_fraction': (
            ['Params_v_sample', '08_taxonomy_BAT', 'fraction'], 0.50, float,
            '\b\n08_taxonomy_BAT: Minimum bit-score support fraction [0-0.99].\n'
            'Higher values enforce stricter consensus (one classification per sequence), while lower values allow multiple lineage detections.'
        ),
        'virus_vmags_taxonomy_bat_threads': (
            ['Params_v_sample', '08_taxonomy_BAT', 'threads'], 8, int,
            '08_taxonomy_BAT: Threads.'),
        'virus_vmags_taxonomy_bat_extra_args': (
            ['Params_v_sample', '08_taxonomy_BAT', 'extra_args'], '', str,
            '08_taxonomy_BAT: Extra arguments forwarded to CAT_pack bins (e.g., "--block_size 12.0").'),
        'virus_coverm_method_sample': (
            ['Params_v_sample', '09_coverm_eachsample_vMAGs', 'method'], 'trimmed_mean', str,
            '\b\n09_coverm_eachsample_vMAGs: Method(s) for coverm-based coverage calculation, must be one of [trimmed_mean/mean].\n'
            'A more detailed description of the different methods is available at https://github.com/wwood/CoverM#calculation-methods'
        ),
        'virus_coverm_min_read_percent_identity_sample': (
            ['Params_v_sample', '09_coverm_eachsample_vMAGs', 'min_read_percent_identity'], 95, float,
            '09_coverm_eachsample_vMAGs: Exclude reads by overall percent identity [0-100], e.g. 95 for 95%.'
        ),
        'virus_coverm_min_read_aligned_length_sample': (
            ['Params_v_sample', '09_coverm_eachsample_vMAGs', 'min_read_aligned_length'], 45, int,
            '09_coverm_eachsample_vMAGs: Exclude reads with fewer aligned bases.'
        ),
        'virus_coverm_contig_end_exclusion_sample': (
            ['Params_v_sample', '09_coverm_eachsample_vMAGs', 'contig_end_exclusion'], 75, int,
            '09_coverm_eachsample_vMAGs: Exclude bases at the ends of contig sequences for coverm calculation.'
        ),
        'virus_coverm_threads_sample': (
            ['Params_v_sample', '09_coverm_eachsample_vMAGs', 'threads'], 8, int,
            '09_coverm_eachsample_vMAGs: Threads.'),
        'virus_coverm_sample_extra_args': (
            ['Params_v_sample', '09_coverm_eachsample_vMAGs', 'extra_args'], '', str,
            '09_coverm_eachsample_vMAGs: Extra arguments forwarded to coverm genome (e.g., "--quiet").'),
        'mmags_checkm2_threads': (
            ['Params_m_sample', '01_checkm2', 'threads'], 60, int,
            '01_checkm2: Threads.'),
        'mmags_checkm2_extra_args': (
            ['Params_m_sample', '01_checkm2', 'extra_args'], '', str,
            '01_checkm2: Extra arguments forwarded to checkm2 predict (e.g., "--remove_intermediates").'),
        'mmags_taxonomy_classify_method': (
            ['Params_m_sample', '02_taxonomy', 'MAGs_taxonomy_classify_method'], 'gtdbtk', str,
            '02_taxonomy: mMAGs taxonomy classifier, should be one of the [gtdbtk/BAT].'
        ),
        'mmags_taxonomy_threads': (
            ['Params_m_sample', '02_taxonomy', 'threads'], 60, int,
            '02_taxonomy_gtdbtk: Threads for [gtdbtk/BAT].'),
        'mmags_taxonomy_gtdbtk_extra_args': (
            ['Params_m_sample', '02_taxonomy_gtdbtk', 'extra_args'], '', str,
            '02_taxonomy_gtdbtk: Extra arguments forwarded to gtdbtk classify_wf (e.g., "--min_af 0.5").'),
        'mmags_taxonomy_bat_range': (
            ['Params_m_sample', '02_taxonomy_BAT', 'range'], 10, float,
            '\b\n02_taxonomy_BAT: Percentage range [0-100] of DIAMOND hits relative to top bit-score to include in LCA calculation.\n'
            'Higher values yield more conservative classifications by incorporating divergent hits, while lower values enable specific assignments.\n'
            'Applicable when --mMAGs-taxonomy-classify-method is set to BAT.'
        ),
        'mmags_taxonomy_bat_fraction': (
            ['Params_m_sample', '02_taxonomy_BAT', 'fraction'], 0.50, float,
            '\b\n02_taxonomy_BAT: Minimum bit-score support fraction [0-0.99].\n'
            'Higher values enforce stricter consensus (one classification per sequence), while lower values allow multiple lineage detections.\n'
            'Applicable when --mMAGs-taxonomy-classify-method is set to BAT.'
        ),
        'mmags_taxonomy_bat_extra_args': (
            ['Params_m_sample', '02_taxonomy_BAT', 'extra_args'], '', str,
            '02_taxonomy_BAT: Extra arguments forwarded to CAT_pack bins (e.g., "--block_size 12.0").'),
        'mmags_drep99_completeness_per_sample': (
            ['Params_m_sample', '05_drep_eachsample_99', 'comp'], 50, float,
            '05_drep_eachsample_99: The first level dereplication of mMAGs in a sample using completeness threshold [0-100], e.g., 50 for 50%.'
        ),
        'mmags_drep99_contamination_per_sample': (
            ['Params_m_sample', '05_drep_eachsample_99', 'con'], 10, float,
            '05_drep_eachsample_99: The first level dereplication of mMAGs in a sample using contamination threshold [0-100], e.g., 10 for 10%.'
        ),
        'mmags_drep99_ani_per_sample': (
            ['Params_m_sample', '05_drep_eachsample_99', 'sa'], 0.99, float,
            '05_drep_eachsample_99: The first level dereplication of mMAGs in a sample using ANI threshold [0-1], e.g., 0.99 for 99%.'
        ),
        'mmags_drep99_sample_threads': (
            ['Params_m_sample', '05_drep_eachsample_99', 'threads'], 20, int,
            '05_drep_eachsample_99: Threads.'),
        'mmags_drep99_sample_extra_args': (
            ['Params_m_sample', '05_drep_eachsample_99', 'extra_args'], '', str,
            '05_drep_eachsample: Extra arguments forwarded to dRep dereplicate (e.g., "--completeness_weight 1 --contamination_weight 5").'),
        'mmags_drep95_completeness_per_sample': (
            ['Params_m_sample', '06_drep_eachsample_95', 'comp'], 50, float,
            '06_drep_eachsample_95: The second level dereplication of mMAGs in a sample using completeness threshold [0-100], e.g., 50 for 50%.'
        ),
        'mmags_drep95_contamination_per_sample': (
            ['Params_m_sample', '06_drep_eachsample_95', 'con'], 10, float,
            '06_drep_eachsample_95: The second level dereplication of mMAGs in a sample using contamination threshold [0-100], e.g., 10 for 10%.'
        ),
        'mmags_drep95_ani_per_sample': (
            ['Params_m_sample', '06_drep_eachsample_95', 'sa'], 0.95, float,
            '\b\n06_drep_eachsample_95: The second level dereplication of mMAGs in a sample using ANI threshold [0-1], e.g., 0.95 for 95%.\n'
            'This value should be smaller than --drep99-ANI-per-sample'
        ),
        'mmags_drep95_sample_threads': (
            ['Params_m_sample', '06_drep_eachsample_95', 'threads'], 20, int,
            '06_drep_eachsample_95: Threads.'),
        'mmags_drep95_sample_extra_args': (
            ['Params_m_sample', '06_drep_eachsample_95', 'extra_args'], '', str,
            '06_drep_eachsample: Extra arguments forwarded to dRep dereplicate (e.g., "--completeness_weight 1 --contamination_weight 5").'),
        'mmags_coverm_method': (
            ['Params_m_sample', '07_coverm_eachsample', 'method'], 'trimmed_mean', str,
            '\b\n07_coverm_eachsample: Method(s) for coverm-based coverage calculation, must be one of [trimmed_mean/mean].\n'
            'A more detailed description of the different methods is available at https://github.com/wwood/CoverM#calculation-methods'
        ),
        'mmags_coverm_min_read_percent_identity': (
            ['Params_m_sample', '07_coverm_eachsample', 'min_read_percent_identity'], 95, float,
            '07_coverm_eachsample: Exclude reads by overall percent identity [0-100], e.g. 95 for 95%.'
        ),
        'mmags_coverm_min_read_aligned_length': (
            ['Params_m_sample', '07_coverm_eachsample', 'min_read_aligned_length'], 45, int,
            '07_coverm_eachsample: Exclude reads with fewer aligned bases.'
        ),
        'mmags_coverm_contig_end_exclusion': (
            ['Params_m_sample', '07_coverm_eachsample', 'contig_end_exclusion'], 75, int,
            '07_coverm_eachsample: Exclude bases at the ends of contig sequences for coverm calculation.'
        ),
        'mmags_coverm_threads': (
            ['Params_m_sample', '07_coverm_eachsample', 'threads'], 20, int,
            '07_coverm_eachsample: Threads.'),
        'mmags_coverm_sample_extra_args': (
            ['Params_m_sample', '07_coverm_eachsample', 'extra_args'], '', str,
            '07_coverm_eachsample: Extra arguments forwarded to coverm genome (e.g., "--quiet").'),
        'mmags_function_classify_method': (
            ['Params_m_sample', '08_function', 'MAGs_function_classify_method'], 'emapper', str,
            '08_function: mMAGs function classifier using one of [kraken2/emapper].'
        ),
        'mmags_function_threads': (
            ['Params_m_sample', '08_function', 'threads'], 60, int,
            '08_function: Threads for [kraken2/emapper].'),
        'mmags_function_kraken2_confidence': (
            ['Params_m_sample', '08_function_kraken2', 'confidence'], 0.05, float,
            '08_function_kraken2: Confidence threshold [0-1] for functional classification for mMAGs when --mMAGs-function-classify-method is set to kraken2.'
        ),
        'mmags_function_kraken2_extra_args': (
            ['Params_m_sample', '08_function_kraken2', 'extra_args'], '', str,
            '08_function_kraken2: Extra arguments forwarded to kraken2 (e.g., "--memory-mapping").'),
        'mmags_function_emapper_evalue': (
            ['Params_m_sample', '08_function_emapper', 'seed_ortholog_evalue'], 1e-10, float,
            '08_function_emapper: Report only alignments with e-value below or equal to the threshold, applicable when --mMAGs-function-classify-method is set to emapper.'
        ),
        'mmags_function_emapper_method': (
            ['Params_m_sample', '08_function_emapper', 'method'], 'diamond', str,
            '\b\n08_function_emapper: search seed orthologs using one of [diamond/mmseqs], applicable when --mMAGs-function-classify-method is set to emapper.\n'
            'diamond: using default parameters of diamond mode of emapper.\n'
            'mmseqs: using default parameters of mmseqs mode of emapper'
        ),
        'mmags_function_emapper_query_cover': (
            ['Params_m_sample', '08_function_emapper', 'query_cover'], 30, float,
            '\b\n08_function_emapper: Report only alignments with query cover above or equal to the given percentage [0-100],.\n'
            'applicable when --mMAGs-function-classify-method is set to emapper'
        ),
        'mmags_function_emapper_extra_args': (
            ['Params_m_sample', '08_function_emapper', 'extra_args'], '', str,
            '08_function_emapper: Extra arguments forwarded to eggNOG-mapper (e.g., "--pident 30").'),
        'contigs_taxonomy_classify_method': (
            ['Params_c_sample', '01_taxonomy', 'contigs_taxonomy_classify_method'], 'kraken2', str,
            '01_taxonomy: Contigs taxonomy classifier using one of [CAT_pack/kraken2].'
        ),
        'contigs_taxonomy_threads': (
            ['Params_c_sample', '01_taxonomy', 'threads'], 60, int,
            '01_taxonomy: Threads for [CAT_pack/kraken2].'
        ),
        'contigs_taxonomy_cat_range': (
            ['Params_c_sample', '01_taxonomy_CAT_pack', 'range'], 10, float,
            '\b\n01_taxonomy_CAT_pack: Percentage range [0-100] of DIAMOND hits relative to top bit-score to include in LCA calculation.\n'
            'Higher values yield more conservative classifications by incorporating divergent hits, while lower values enable specific assignments.\n'
            'Applicable when --contigs-taxonomy-classify-method is set to CAT_pack.'
        ),
        'contigs_taxonomy_cat_fraction': (
            ['Params_c_sample', '01_taxonomy_CAT_pack', 'fraction'], 0.50, float,
            '\b\n01_taxonomy_CAT_pack: Minimum bit-score support fraction [0-0.99].\n'
            'Higher values enforce stricter consensus (one classification per sequence), while lower values allow multiple lineage detections.\n'
            'Applicable when --contigs-taxonomy-classify-method is set to CAT_pack.'
        ),
        'contigs_taxonomy_cat_pack_extra_args': (
            ['Params_c_sample', '01_taxonomy_CAT_pack', 'extra_args'], '', str,
            '01_taxonomy_CAT_pack: Extra arguments forwarded to CAT_pack contigs (e.g., "--block_size 12.0").'),
        'contigs_taxonomy_kraken2_confidence': (
            ['Params_c_sample', '01_taxonomy_kraken2', 'confidence'], 0.7, float,
            '\b\n01_taxonomy_kraken2: Confidence threshold [0-1] for taxonomic classification, applicable when --contigs-taxonomy-classify-method is set to kraken2.\n'
            'Default value was tuned on the CAMI2 marine dataset. For samples with higher complexity, a lower confidence threshold is usually recommended.'
        ),
        'contigs_taxonomy_kraken2_extra_args': (
            ['Params_c_sample', '01_taxonomy_kraken2', 'extra_args'], '', str,
            '01_taxonomy_kraken2: Extra arguments forwarded to kraken2 (e.g., "--memory-mapping").'),
        'contigs_function_classify_method': (
            ['Params_c_sample', '02_function', 'contigs_function_classify_method'], 'kraken2', str,
            '02_function: Contigs function classifier using one of [kraken2/emapper].'
        ),
        'contigs_function_threads': (
            ['Params_c_sample', '02_function', 'threads'], 60, int,
            '02_function: Threads for [kraken2/emapper].'),
        'contigs_function_emapper_evalue': (
            ['Params_c_sample', '02_function_emapper', 'seed_ortholog_evalue'], 1e-10, float,
            '02_function_emapper: Report only alignments with e-value below or equal to the threshold, applicable when --contigs-function-classify-method is set to emapper.'
        ),
        'contigs_function_emapper_method': (
            ['Params_c_sample', '02_function_emapper', 'method'], 'diamond', str,
            '\b\n02_function_emapper: search seed orthologs using one of [diamond/mmseqs], applicable when --contigs-function-classify-method is set to emapper.\n'
            'diamond: using default parameters of diamond mode of emapper.\n'
            'mmseqs: using default parameters of mmseqs mode of emapper'
        ),
        'contigs_function_emapper_query_cover': (
            ['Params_c_sample', '02_function_emapper', 'query_cover'], 30, float,
            '\b\n02_function_emapper: Report only alignments with query cover above or equal to the given percentage [0-100],.\n'
            'applicable when --contigs-function-classify-method is set to emapper'
        ),
        'contigs_function_emapper_extra_args': (
            ['Params_c_sample', '02_function_emapper', 'extra_args'], '', str,
            '02_function_emapper: Extra arguments forwarded to eggNOG-mapper (e.g., "--pident 30").'),
        'contigs_function_kraken2_confidence': (
            ['Params_c_sample', '02_function_kraken2', 'confidence'], 0.05, float,
            '02_function_kraken2: Confidence threshold [0-1] for functional classification, applicable when --contigs-function-classify-method is set to kraken2.'
        ),
        'contigs_function_kraken2_extra_args': (
            ['Params_c_sample', '02_function_kraken2', 'extra_args'], '', str,
            '02_function_kraken2: Extra arguments forwarded to kraken2 (e.g., "--memory-mapping").'),
        'vmcr_unmapped_reads_merge_threads': (
            ['Params_vmcr_sample', '01_unmapped_reads_trim_merge', 'threads'], 16, int,
            '01_unmapped_reads_trim_merge_fastp: Threads for fastp, samtools, and seqkit (paired-end), or samtools and seqkit only (single-end).'),
        'vmcr_unmapped_reads_merge_length_required': (
            ['Params_vmcr_sample', '01_unmapped_reads_trim_merge', 'length_required'], 15, int,
            '01_unmapped_reads_trim_merge_fastp: Reads shorter than length_required will be discarded. Only applicable when read_type is paired.'),
        'vmcr_unmapped_reads_merge_n_base_limit': (
            ['Params_vmcr_sample', '01_unmapped_reads_trim_merge', 'n_base_limit'], 5, int,
            "01_unmapped_reads_trim_merge_fastp: if one read's number of N base is >n_base_limit, then this read/pair is discarded. Only applicable when read_type is paired."),
        'vmcr_unmapped_reads_merge_compression': (
            ['Params_vmcr_sample', '01_unmapped_reads_trim_merge', 'compression'], 4, int,
            '01_unmapped_reads_trim_merge_fastp: Compression level for gzip output [1-9]. 1 is fastest, 9 is smallest. Only applicable when read_type is paired.'),
        'vmcr_unmapped_reads_fastp_extra_args': (
            ['Params_vmcr_sample', '01_unmapped_reads_trim_merge', 'extra_args'], '', str,
            '01_unmapped_reads_trim_merge: Extra arguments forwarded to fastp for unmapped reads trimming (e.g., "--reads_to_process 0 -Q"). Only applicable when read_type is paired.'),
        'vmcr_fraggenescanrs_model': (
            ['Params_vmcr_sample', '02_FragGeneScanRs', 'model'], 'illumina_1', str,
            '\b\n02_FragGeneScanRs: Model parameter used for short-read open reading frame prediction.\n'
            '[complete] for complete genomic sequences or short sequence reads without sequencing error\n'
            '[sanger_5] for Sanger sequencing reads with about 0.5% error rate\n'
            '[sanger_10] for Sanger sequencing reads with about 1% error rate\n'
            '[454_5] for 454 pyrosequencing reads with about 0.5% error rate\n'
            '[454_10] for 454 pyrosequencing reads with about 1% error rate\n'
            '[454_30] for 454 pyrosequencing reads with about 3% error rate\n'
            '[illumina_1] for Illumina sequencing reads with about 0.1% error rate\n'
            '[illumina_5] for Illumina sequencing reads with about 0.5% error rate\n'
            '[illumina_10] for Illumina sequencing reads with about 1% error rate'
        ),
        'vmcr_fraggenescanrs_threads': (
            ['Params_vmcr_sample', '02_FragGeneScanRs', 'threads'], 20, int,
            '02_FragGeneScanRs: Threads.'),
        'vmcr_fraggenescanrs_stats_threads': (
            ['Params_vmcr_sample', '02_FragGeneScanRs_stats', 'threads'], 10, int,
            '02_FragGeneScanRs_stats: Threads.'),
        'vmcr_taxonomy_threads': (
            ['Params_vmcr_sample', '03_taxonomy', 'threads'], 60, int,
            '03_taxonomy: Threads.'),
        'vmcr_taxonomy_kraken2_confidence': (
            ['Params_vmcr_sample', '03_taxonomy_kraken2', 'confidence'], 0.7, float,
            '\b\n03_taxonomy_kraken2: Confidence threshold [0-1] for taxonomic classification of unmapped reads with kraken2.\n'
            'Default value was tuned on the CAMI2 marine dataset. For samples with higher complexity, a lower confidence threshold is usually recommended.'
        ),
        'vmcr_taxonomy_kraken2_extra_args': (
            ['Params_vmcr_sample', '03_taxonomy_kraken2', 'extra_args'], '', str,
            '03_taxonomy_kraken2: Extra arguments forwarded to kraken2 (e.g., "--memory-mapping").'),
        'vmcr_function_classify_method': (
            ['Params_vmcr_sample', '04_function', 'reads_function_classify_method'], 'kraken2', str,
            '04_function: reads functional classifier using one of [kraken2/emapper].'
        ),
        'vmcr_function_threads': (
            ['Params_vmcr_sample', '04_function', 'threads'], 60, int,
            '04_function: Threads for [kraken2/emapper].'),
        'vmcr_function_kraken2_confidence': (
            ['Params_vmcr_sample', '04_function_kraken2', 'confidence'], 0.05, float,
            '04_function_kraken2: Confidence threshold [0-1] for functional classification for reads when --vmcr-function-classify-method is set to kraken2.'
        ),
        'vmcr_function_kraken2_extra_args': (
            ['Params_vmcr_sample', '04_function_kraken2', 'extra_args'], '', str,
            '04_function_kraken2: Extra arguments forwarded to kraken2 (e.g., "--memory-mapping").'),
        'vmcr_function_emapper_evalue': (
            ['Params_vmcr_sample', '04_function_emapper', 'seed_ortholog_evalue'], 1e-10, float,
            '04_function_emapper: Report only alignments with e-value below or equal to the threshold, applicable when --vmcr-function-classify-method is set to emapper.'
        ),
        'vmcr_function_emapper_method': (
            ['Params_vmcr_sample', '04_function_emapper', 'method'], 'diamond', str,
            '\b\n04_function_emapper: search seed orthologs using one of [diamond/mmseqs], applicable when --vmcr-function-classify-method is set to emapper.\n'
            'diamond: using default parameters of diamond mode of emapper.\n'
            'mmseqs: using default parameters of mmseqs mode of emapper'
        ),
        'vmcr_function_emapper_query_cover': (
            ['Params_vmcr_sample', '04_function_emapper', 'query_cover'], 30, float,
            '\b\n04_function_emapper: Report only alignments with query cover above or equal to the given percentage [0-100],.\n'
            'applicable when --vmcr-function-classify-method is set to emapper'
        ),
        'vmcr_function_emapper_extra_args': (
            ['Params_vmcr_sample', '04_function_emapper', 'extra_args'], '', str,
            '04_function_emapper: Extra arguments forwarded to eggNOG-mapper (e.g., "--pident 30").'),
        'merge_mmags_min_completeness1': (
            ['Params_vmcr_sample', '05_merge_tax_fun_depth', 'MAGs_min_completeness1'], 50, float,
            '02_quality_gate: Minimum completeness threshold (%) for CheckM1-based mMAG quality filtering [0-100], e.g., 50 for 50%.'
        ),
        'merge_mmags_max_contamination1': (
            ['Params_vmcr_sample', '05_merge_tax_fun_depth', 'MAGs_max_contamination1'], 1, float,
            '02_quality_gate: Maximum contamination threshold (%) for CheckM1-based mMAG quality filtering [0-100], e.g., 5 for 5%.'
        ),
        'merge_mmags_min_completeness2': (
            ['Params_vmcr_sample', '05_merge_tax_fun_depth', 'MAGs_min_completeness2'], 50, float,
            '02_quality_gate: Minimum completeness threshold (%) for CheckM2-based mMAG quality filtering [0-100], e.g., 50 for 50%.'
        ),
        'merge_mmags_max_contamination2': (
            ['Params_vmcr_sample', '05_merge_tax_fun_depth', 'MAGs_max_contamination2'], 1, float,
            '02_quality_gate: Maximum contamination threshold (%) for CheckM2-based mMAG quality filtering [0-100], e.g., 5 for 5%.'
        ),
        'merge_min_orf_length': (
            ['Params_vmcr_sample', '05_merge_tax_fun_depth', 'minORFLen'], 50, int,
            '05_merge_tax_fun_depth: Minimum ORF length for finally used ORFs of vMAGs, mMAGs, contigs, and reads.'
        ),
        'merge_vmags_checkv_quality': (
            ['Params_vmcr_sample', '05_merge_tax_fun_depth', 'checkv_quality'], 'High-quality', str,
            '\b\n05_merge_tax_fun_depth: CheckV quality filter for viral MAGs (lower bound filter).\n'
            'viral MAGs with quality equal to or higher than the specified threshold will be retained.\n'
            'Options: Complete, High-quality, Medium-quality, Low-quality, Not-determined.\n'
            'Complete: complete genome (highest quality)\n'
            'High-quality: genome completeness ≥90%\n'
            'Medium-quality: genome completeness 50-90%\n'
            'Low-quality: genome completeness <50%\n'
            'Not-determined: completeness could not be determined (lowest quality)'
        ),
        'merge_multiple_ko_handling_mode': (
            ['Params_vmcr_sample', '05_merge_tax_fun_depth', 'KO_mode'], 'ignore_multiple_KOs', str,
            '\b\n05_merge_tax_fun_depth: Method using one of [split_multiple_KOs/keep_multiple_KOs/ignore_multiple_KOs] when an ORF has multiple KO annotations:.\n'
            'split_multiple_KOs: evenly distributing the total abundance across all assigned KOs;\n'
            'keep_multiple_KOs: retain original multiple KO annotations;\n'
            'ignore_multiple_KOs: ignore the ORFs with multiple KO annotations and only use the ORFs with single KO annotation.'
        ),
        'merge_min_ko_num_per_taxon': (
            ['Params_vmcr_sample', '05_merge_tax_fun_depth', 'min_KO_num_in_one_taxon'], 100, int,
            '\b\n05_merge_tax_fun_depth: Minimum number of KOs per taxon, used to filter outputs in KO_Taxa_RelativeAbundance and KO_Taxa_SumCoverage.\n'
            'This threshold filters out taxa with insufficient KOs,\n'
            'retaining only KO-rich taxa to ensure robust metabolic reconstruction\n'
            'while preventing excessive output file sizes caused by sparse matrices.'
        ),
        'merge_tax_fun_depth_threads': (
            ['Params_vmcr_sample', '05_merge_tax_fun_depth', 'threads'], 20, int,
            '05_merge_tax_fun_depth: Threads.'),
        'merge_cami_format_output': (
            ['Params_vmcr_sample', '05_merge_tax_fun_depth', 'CAMI_format_output'], False, bool,
            '05_merge_tax_fun_depth: Report a taxonomic profile using the CAMI output format (flag only). Disabled when mMAG mode uses gtdbtk.'
        ),
        'merge_output_detail': (
            ['Params_vmcr_sample', '05_merge_tax_fun_depth', 'output_detail'], False, bool,
            '05_merge_tax_fun_depth: Whether to output detail tables (flag only).'),
        'mclink_sample_min_ko_num_per_taxon': (
            ['Params_vmcr_sample', '06_mclink_taxon', 'min_KO_num_in_one_taxon'], 100, int,
            '06_mclink_taxon: Minimum number of KOs per taxon for calculating taxon-level metabolic pathway/module completeness or abundance per sample. Taxa with fewer than this number of detected KOs in a given sample will be filtered out. Calculated independently for each sample.'
            ),
        'mclink_sample_min_ko_abundance_per_taxon': (
            ['Params_vmcr_sample', '06_mclink_taxon', 'min_mean_KO_RA_in_one_taxon'], 0, float,
            '06_mclink_taxon: Minimum mean KO abundance [0-1e6] per taxon for calculating taxon-level metabolic pathway/module completeness or abundance per sample. Taxa with mean KO abundance below this threshold in a given sample will be filtered out. Note: GPM/TPM sum of all taxa in each sample is 1e6.'
            ),
        'mclink_sample_topn_taxon': (
            ['Params_vmcr_sample', '06_mclink_taxon', 'TopN'], 100, int,
            '06_mclink_taxon: Number of top taxa (ranked by mean KO abundance) retained for calculating taxon-level metabolic pathway/module completeness or abundance per sample. Set to 0 to disable filtering.'
        ),
        'mclink_sample_completeness_per_taxon': (
            ['Params_vmcr_sample', '06_mclink_taxon', 'completeness_threshold'], 0.5, float,
            '06_mclink_taxon: Completeness threshold [0-1] of modules for calculating taxon-level metabolic pathway/module abundance. Higher values retain only more complete profiles.'),
        'mclink_sample_plus_scale_method_per_taxon': (
            ['Params_vmcr_sample', '06_mclink_taxon', 'plus_scale_method'], 'mean', str,
            "06_mclink_taxon: Aggregation method for '+'-separated KOs in a module [mean/max/min]. 'mean': average abundance across all required KOs; 'max': highest KO abundance; 'min': lowest KO abundance."),
        'mclink_sample_comma_scale_method_per_taxon': (
            ['Params_vmcr_sample', '06_mclink_taxon', 'comma_scale_method'], 'sum', str,
            "06_mclink_taxon: Aggregation method for ','-separated alternative KOs in a module [sum/max]. 'sum': cumulative abundance of all alternative KOs; 'max': abundance of the highest single KO among alternatives. Note: this applies to abundance calculation only; module completeness is always computed with 'max'."),
        'barplot_topn_taxon': (
            ['Params_vmcr_project', '01_Taxa_visualize', 'TopN'], 10, int,
            '01_Taxa_visualize: Number of top taxa to display in stacked barplots. Remaining categories are merged into Others. Maximum 29.'
        ),
        'barplot_topn_fun': (
            ['Params_vmcr_project', '03_Fun_visualize', 'TopN'], 10, int,
            '03_Fun_visualize: Number of top functional categories to display in stacked barplots. Remaining categories are merged into Others. Maximum 29.'
        ),
        'mclink_min_ko_num_per_sample': (
            ['Params_vmcr_project', '03_mclink_sample', 'min_sample_num_per_KO'], 1, int,
            '03_mclink_sample: Minimum number of samples in which a KO must be detected for calculating community-level metabolic pathway/module completeness or abundance.'
        ),
        'mclink_min_ko_prevalence_per_sample': (
            ['Params_vmcr_project', '03_mclink_sample', 'preva'], 0, float,
            '03_mclink_sample: Minimum detection rate [0-1] of a KO across all samples for calculating community-level metabolic pathway/module completeness or abundance.'
        ),
        'mclink_min_ko_abundance_per_sample': (
            ['Params_vmcr_project', '03_mclink_sample', 'abund'], 0, float,
            '03_mclink_sample: Minimum mean relative abundance [0-1e6] of a KO across all samples for calculating community-level metabolic pathway/module completeness or abundance. KOs with mean abundance below this threshold are excluded. Note: GPM/TPM sum of a sample is 1e6.'
        ),
        'mclink_completeness_per_sample': (
            ['Params_vmcr_project', '03_mclink_sample', 'completeness_threshold'], 0.5, float,
            '03_mclink_sample: Completeness threshold [0-1] of modules for calculating community-level metabolic pathway/module abundance. Higher values retain only more complete profiles.'
        ),
        'mclink_plus_scale_method_per_sample': (
            ['Params_vmcr_project', '03_mclink_sample', 'plus_scale_method'], 'mean', str,
            "03_mclink_sample: Aggregation method for '+'-separated KOs in a module [mean/max/min]. 'mean': average abundance across all required KOs; 'max': highest KO abundance; 'min': lowest KO abundance."
        ),
        'mclink_comma_scale_method_per_sample': (
            ['Params_vmcr_project', '03_mclink_sample', 'comma_scale_method'], 'sum', str,
            "03_mclink_sample: Aggregation method for ','-separated alternative KOs in a module [sum/max]. 'sum': cumulative abundance of all alternative KOs; 'max': abundance of the highest single KO among alternatives. Note: this applies to abundance calculation only; module completeness is always computed with 'max'."
        ),
        'merge_min_taxa_ko_prevalence_across_samples': (
            ['Params_vmcr_project', '05_merge_tax_KO_depth', 'preva'], 0.10, float,
            '05_merge_tax_KO_depth: Minimum prevalence [0-1] of a Taxon-KO across samples. Taxa-KOs detected in fewer samples than this proportion are excluded.'
        ),
        'merge_min_taxa_ko_abundance_across_samples': (
            ['Params_vmcr_project', '05_merge_tax_KO_depth', 'abund'], 0.0, float,
            '05_merge_tax_KO_depth: Minimum mean relative abundance [0-1e6] of a Taxon-KO across samples. Taxa-KOs with mean abundance below this threshold are excluded. Note: GPM/TPM sum of a sample is 1e6.'
        ),
        'merge_min_taxa_module_prevalence_across_samples': (
            ['Params_vmcr_project', '05_merge_tax_Module_RelativeAbundance', 'preva'], 0.10, float,
            '05_merge_tax_Module_RelativeAbundance: Minimum prevalence [0-1] of a Taxon-Module across samples. Taxa-Modules detected in fewer samples than this proportion are excluded.'
        ),
        'merge_min_taxa_module_abundance_across_samples': (
            ['Params_vmcr_project', '05_merge_tax_Module_RelativeAbundance', 'abund'], 0.0, float,
            '05_merge_tax_Module_RelativeAbundance: Minimum mean relative abundance [0-1e6] of a Taxon-Module across samples. Taxa-Modules with mean abundance below this threshold are excluded. Note: GPM/TPM sum of a sample is 1e6.'
        ),
        'merge_min_taxa_module_sc_prevalence_across_samples': (
            ['Params_vmcr_project', '05_merge_tax_Module_SumCoverage', 'preva'], 0.10, float,
            '05_merge_tax_Module_SumCoverage: Minimum prevalence [0-1] of a Taxon-Module across samples. Taxa-Modules detected in fewer samples than this proportion are excluded.'
        ),
        'merge_min_taxa_module_sc_abundance_across_samples': (
            ['Params_vmcr_project', '05_merge_tax_Module_SumCoverage', 'abund'], 0.0, float,
            '05_merge_tax_Module_SumCoverage: Minimum mean sum coverage [0-1e6] of a Taxon-Module across samples. Taxa-Modules with mean sum coverage below this threshold are excluded. Note: GPM/TPM sum of a sample is 1e6.'
        ),
        'mclink_project_min_ko_num_per_taxon': (
            ['Params_vmcr_project', '06_mclink_taxon', 'min_KO_num_in_one_taxon'], 100, int,
            '06_mclink_taxon: Minimum number of KOs per taxon for calculating taxon-level metabolic pathway/module completeness on the project-averaged taxon-function matrix. Taxa with fewer than this number of detected KOs in the project-averaged matrix will be filtered out. KO abundances are first averaged across all samples to generate a single project-level matrix before calculation.'
        ),
        'mclink_project_min_ko_abundance_per_taxon': (
            ['Params_vmcr_project', '06_mclink_taxon', 'min_mean_KO_RA_in_one_taxon'], 0, float,
            '06_mclink_taxon: Minimum mean KO abundance [0-1e6] per taxon for calculating taxon-level metabolic pathway/module completeness on the project-averaged taxon-function matrix. Taxa with mean KO abundance below this threshold in the project-averaged matrix will be filtered out. Note: GPM/TPM sum of all taxa is 1e6 per sample; values are averaged across samples prior to filtering and calculation.'
        ),
        'mclink_project_topn_taxon': (
            ['Params_vmcr_project', '06_mclink_taxon', 'TopN'], 100, int,
            '06_mclink_taxon: Number of top taxa (ranked by mean KO abundance in the project-averaged taxon-function matrix) retained for calculating taxon-level metabolic pathway/module completeness. Set to 0 to disable filtering.'
        ),
        'mclink_project_completeness_per_taxon': (
            ['Params_vmcr_project', '06_mclink_taxon', 'completeness_threshold'], 0.5, float,
            '06_mclink_taxon: Completeness threshold [0-1] of modules for calculating taxon-level metabolic pathway/module abundance on the project-averaged taxon-function matrix. Higher values retain only more complete profiles.'),
        'mclink_project_plus_scale_method_per_taxon': (
            ['Params_vmcr_project', '06_mclink_taxon', 'plus_scale_method'], 'mean', str,
            "06_mclink_taxon: Aggregation method for '+'-separated KOs in a module [mean/max/min] on the project-averaged taxon-function matrix. 'mean': average abundance across all required KOs; 'max': highest KO abundance; 'min': lowest KO abundance."),
        'mclink_project_comma_scale_method_per_taxon': (
            ['Params_vmcr_project', '06_mclink_taxon', 'comma_scale_method'], 'sum', str,
            "06_mclink_taxon: Aggregation method for ','-separated alternative KOs in a module [sum/max] on the project-averaged taxon-function matrix. 'sum': cumulative abundance of all alternative KOs; 'max': abundance of the highest single KO among alternatives. Note: this applies to abundance calculation only; module completeness is always computed with 'max'."),
    },
    'vc': {
        'raw_stats_sample_threads': (
            ['Params_b_sample', '00_raw_stats', 'threads'], 4, int,
            '00_raw_stats: Threads for seqkit.'),
        'trim_sample_threads': (
            ['Params_b_sample', '01_trim', 'threads'], 16, int,
            '01_trim: Threads for fastp (next-generation short reads) or fastplong (third-generation long reads).'),
        'trim_length_required': (
            ['Params_b_sample', '01_trim', 'length_required'], 15, int,
            '01_trim: Reads shorter than length_required will be discarded'),
        'trim_n_base_limit': (
            ['Params_b_sample', '01_trim', 'n_base_limit'], 5, int,
            "01_trim: if one read's number of N base is >n_base_limit, then this read/pair is discarded."),
        'trim_compression': (
            ['Params_b_sample', '01_trim', 'compression'], 4, int,
            '01_trim: Compression level for gzip output [1-9]. 1 is fastest, 9 is smallest.'),
        'trim_extra_args': (
            ['Params_b_sample', '01_trim', 'extra_args'], '', str,
            '01_trim: Extra arguments forwarded to fastp or fastplong for raw reads trimming (e.g., "--reads_to_process 0 -Q").'),
        'trim_stats_sample_threads': (
            ['Params_b_sample', '01_trim_stats', 'threads'], 16, int,
            '01_trim_stats: Threads.'),
        'assembly_method': (
            ['Params_b_sample', '02_assembly', 'method'], 'megahit', str,
            '\b\n02_assembly: Underlying assembly software using one of [megahit/spades/metaflye].\n'
            'Must be compatible with --read-type (config key read_type):\n'
            'megahit: MEGAHIT (second-generation sequencing; for paired or single DNA reads).\n'
            'spades: SPAdes (second-generation sequencing; metaSPAdes for paired DNA reads, rnaSPAdes for paired-rna or single-rna RNA reads).\n'
            'metaflye: metaFlye (third-generation sequencing; required for pacbio-raw/corr/hifi or nano-raw/corr/hq long reads).'),
        'assembly_threads': (
            ['Params_b_sample', '02_assembly', 'threads'], 40, int,
            '02_assembly: Threads for [megahit/spades/metaflye].'),
        'megahit_presets': (
            ['Params_b_sample', '02_megahit', 'presets'], 'meta-sensitive', str,
            '\b\n02_megahit: MEGAHIT presets using one of [meta-sensitive/meta-large].\n'
            "meta-sensitive: '--min-count 1 --k-list 21,29,39,49,...,129,141'.\n"
            "meta-large: '--k-min 27 --k-max 127 --k-step 10'.\n"
            '(large & complex metagenomes, like soil).'),
        'megahit_k_list': (
            ['Params_b_sample', '02_megahit', 'k_list'], '', str,
            '02_megahit: MEGAHIT --k-list (optional). Explicit comma-separated k-mer sizes (e.g. "21,29,39,59,79,99,119,141"). Overrides --megahit-presets. If set, --megahit-k-min, --megahit-k-max, --megahit-k-step must NOT be set. Optional --megahit-min-count can be used together.'),
        'megahit_min_count': (
            ['Params_b_sample', '02_megahit', 'min_count'], 2, int,
            '02_megahit: MEGAHIT --min-count (optional). Minimum multiplicity for filtering (k_min+1)-mers. Can be used with --megahit-k-list or --megahit-k-min/--k-max/--k-step.'),
        'megahit_k_min': (
            ['Params_b_sample', '02_megahit', 'k_min'], '', str,
            '02_megahit: MEGAHIT --k-min (optional). Minimum kmer size (<= 255), must be odd number. If set with --k-max and --k-step, overrides --megahit-presets and uses custom k-mer mode.'),
        'megahit_k_max': (
            ['Params_b_sample', '02_megahit', 'k_max'], '', str,
            '02_megahit: MEGAHIT --k-max (optional). Maximum kmer size (<= 255), must be odd number.'),
        'megahit_k_step': (
            ['Params_b_sample', '02_megahit', 'k_step'], '', str,
            '02_megahit: MEGAHIT --k-step (optional). Increment of kmer size of each iteration (<= 28), must be even number.'),
        'megahit_extra_args': (
            ['Params_b_sample', '02_megahit', 'extra_args'], '', str,
            '02_megahit: Extra arguments forwarded to MEGAHIT (e.g. "--no-mercy --tmp-dir /path/to/your/tmp/dir"). Do NOT use this for k-mer settings (--k-list, --k-min, --presets, etc.).'),
        'spades_extra_args': (
            ['Params_b_sample', '02_spades', 'extra_args'], '', str,
            '02_spades: Extra arguments forwarded to SPAdes. Do NOT use this for input/output or thread settings.'),
        'metaflye_extra_args': (
            ['Params_b_sample', '02_metaflye', 'extra_args'], '', str,
            '02_metaflye: Extra arguments forwarded to Flye for third-generation long reads assembly (e.g., "--iterations 2"). Do NOT use this for read type or input/output settings.'),
        'mapping_threads': (
            ['Params_b_sample', '03_mapping', 'threads'], 20, int,
            '03_mapping: Threads for [minibwa/bwamem/bwamem2/strobealign].'),
        'mapping_method': (
            ['Params_b_sample', '03_mapping', 'method'], 'minibwa', str,
            '\b\n03_mapping: Underlying mapping software used [minibwa/bwamem/bwamem2/strobealign].\n'
            'minibwa: minibwa using default parameters.\n'
            'bwamem: bwa mem using default parameters.\n'
            'bwamem2: bwa-mem2 using default parameters.\n'
            'strobealign: strobealign using default parameters.'),
        'mapping_samtools_sort_mem': (
            ['Params_b_sample', '03_mapping', 'samtools_sort_mem'], '768M', str,
            '03_mapping: Set maximum memory per thread; suffix K/M/G recognized, e.g. 2G, 768M.'),
        'mapping_extra_args': (
            ['Params_b_sample', '03_mapping', 'extra_args'], '', str,
            '03_mapping: Extra arguments forwarded to the active mapping aligner (minibwa map, bwa mem, bwa-mem2 mem, or strobealign --use-index), (e.g. "-B 8" for minibwa map).'),
        'contigs_depth_sample_threads': (
            ['Params_b_sample', '07_contigs_depth', 'threads'], 8, int,
            '07_contigs_depth: Threads.'),
        'depth_min_read_percent_identity': (
            ['Params_b_sample', '07_contigs_depth', 'min_read_percent_identity'], 95, float,
            '\b\n07_contigs_depth: Exclude reads by overall percent identity [0-100] for contig depth estimation\n'
            'using jgi_summarize_bam_contig_depths --percentIdentity, e.g. 95 for 95%.\n'
            'Reads with percent identity below this threshold are excluded from depth calculation.\n'
            'Higher values increase stringency, retaining only high-quality alignments.'),
        'depth_contig_end_exclusion': (
            ['Params_b_sample', '07_contigs_depth', 'contig_end_exclusion'], 75, int,
            '\b\n07_contigs_depth: Exclude bases at the ends of contig sequences for depth estimation using jgi_summarize_bam_contig_depths --maxEdgeBases.\n'
            'Bases within this distance from either end of a contig are ignored when calculating depth,\n'
            'preventing edge effects from partial alignments.'),
        'depth_min_contig_length': (
            ['Params_b_sample', '07_contigs_depth', 'min_contig_length'], 300, int,
            '\b\n07_contigs_depth: Minimum contig length for depth estimation [bp], e.g. 300 for 300 bp.\n'
            'Contigs shorter than this threshold are excluded from depth calculation.'),
        'contigs_prodigal_extra_args': (
            ['Params_b_sample', '07_contigs_depth', 'extra_args'], '', str,
            '07_contigs_depth: Extra arguments forwarded to prodigal for contig gene prediction (e.g., "-m -n").'),
        'raw_stats_project_threads': (
            ['Params_b_project', '01_cat_raw_stats', 'threads'], 4, int,
            '01_cat_raw_stats: Threads.'),
        'trimmed_stats_project_threads': (
            ['Params_b_project', '02_cat_trimmed_stats', 'threads'], 4, int,
            '02_cat_trimmed_stats: Threads.'),
        'contigs_stats_project_threads': (
            ['Params_b_project', '03_contigs_stats', 'threads'], 4, int,
            '03_contigs_stats: Threads.'),
        'merge_stats_project_threads': (
            ['Params_b_project', '04_merge_reads_contigs_stats', 'threads'], 4, int,
            '04_merge_reads_contigs_stats: Threads.'),
        'virus_identify_method': (
            ['Params_v_sample', '01_virus_identify', 'identify_method'], 'genomad_vs_vb', str,
            '\b\n01_virus_identify: Virus identification method, should be one of the [genomad/vs_vb/genomad_vs_vb].\n'
            'genomad: Use geNomad;\n'
            'vs_vb: Use VirSorter2 and CheckV, then use VIBRANT to check the low quality viruses from VirSorter2;\n'
            'genomad_vs_vb: Combine genomad and vs_vb virus to get the union results.'
        ),
        'virus_length_limit': (
            ['Params_v_sample', '01_virus_identify', 'virus_length_limit'], 2000, int,
            '01_virus_identify: Minimum length (bp) for viral contigs.'
        ),
        'virus_identify_threads': (
            ['Params_v_sample', '01_virus_identify', 'threads'], 8, int,
            '01_virus_identify: Threads.'),
        'virus_contig_taxonomy_cat_range': (
            ['Params_v_sample', '03_taxonomy_CAT_pack', 'range'], 10, float,
            '\b\n03_taxonomy_CAT_pack: Percentage range [0-100] of DIAMOND hits relative to top bit-score to include in LCA calculation.\n'
            'Higher values yield more conservative classifications by incorporating divergent hits, while lower values enable specific assignments.'
        ),
        'virus_contig_taxonomy_cat_fraction': (
            ['Params_v_sample', '03_taxonomy_CAT_pack', 'fraction'], 0.50, float,
            '\b\n03_taxonomy_CAT_pack: Minimum bit-score support fraction [0-0.99].\n'
            'Higher values enforce stricter consensus (one classification per sequence), while lower values allow multiple lineage detections.'
        ),
        'virus_contig_taxonomy_cat_pack_threads': (
            ['Params_v_sample', '03_taxonomy_CAT_pack', 'threads'], 8, int,
            '03_taxonomy_CAT_pack: Threads.'),
        'virus_contig_taxonomy_cat_pack_extra_args': (
            ['Params_v_sample', '03_taxonomy_CAT_pack', 'extra_args'], '', str,
            '03_taxonomy_CAT_pack: Extra arguments forwarded to CAT_pack contigs (e.g., "--block_size 12.0").'),
        'virus_mapping_threads': (
            ['Params_v_sample', '03_mapping', 'threads'], 20, int,
            '03_mapping: Threads for [minibwa/bwamem/bwamem2/strobealign].'),
        'virus_mapping_method': (
            ['Params_v_sample', '03_mapping', 'method'], 'minibwa', str,
            '\b\n03_mapping: Underlying mapping software used [minibwa/bwamem/bwamem2/strobealign].\n'
            'minibwa: minibwa using default parameters.\n'
            'bwamem: bwa mem using default parameters.\n'
            'bwamem2: bwa-mem2 using default parameters.\n'
            'strobealign: strobealign using default parameters.'),
        'virus_mapping_samtools_sort_mem': (
            ['Params_v_sample', '03_mapping', 'samtools_sort_mem'], '768M', str,
            '03_mapping: Set maximum memory per thread; suffix K/M/G recognized, e.g. 2G, 768M.'),
        'virus_mapping_extra_args': (
            ['Params_v_sample', '03_mapping', 'extra_args'], '', str,
            '03_mapping: Extra arguments forwarded to the active mapping aligner (minibwa map, bwa mem, bwa-mem2 mem, or strobealign --use-index), (e.g. "-B 8" for minibwa map).'),
        'virus_binning_threads': (
            ['Params_v_sample', '04_binning', 'threads'], 8, int,
            '04_binning_vRhyme: Threads.'),
        'virus_checkv_threads': (
            ['Params_v_sample', '05_checkv', 'threads'], 8, int,
            '05_checkv: Threads.'),
        'virus_vclust_ani_sample': (
            ['Params_v_sample', '06_vclust', 'ani'], 0.95, float,
            '06_vclust: ANI threshold [0-1] of vclust for sample-level virus clustering.'
        ),
        'virus_vclust_threads_sample': (
            ['Params_v_sample', '06_vclust', 'threads'], 8, int,
            '06_vclust: Threads.'),
        'vclust_sample_prefilter_extra_args': (
            ['Params_v_sample', '06_vclust', 'prefilter_extra_args'], '', str,
            '06_vclust: Extra arguments forwarded to vclust prefilter (e.g., "-k 25").'),
        'vclust_sample_align_extra_args': (
            ['Params_v_sample', '06_vclust', 'align_extra_args'], '', str,
            '06_vclust: Extra arguments forwarded to vclust align (e.g., "--mal 11").'),
        'vclust_sample_cluster_extra_args': (
            ['Params_v_sample', '06_vclust', 'cluster_extra_args'], '', str,
            '06_vclust: Extra arguments forwarded to vclust cluster (e.g., "--len_ratio 0").'),
        'virus_iphop_threads_sample': (
            ['Params_v_sample', '07_iPHoP', 'threads'], 20, int,
            '07_iPHoP: Threads.'),
        'virus_vmags_taxonomy_bat_range': (
            ['Params_v_sample', '08_taxonomy_BAT', 'range'], 10, float,
            '\b\n08_taxonomy_BAT: Percentage range [0-100] of DIAMOND hits relative to top bit-score to include in LCA calculation.\n'
            'Higher values yield more conservative classifications by incorporating divergent hits, while lower values enable specific assignments.'
        ),
        'virus_vmags_taxonomy_bat_fraction': (
            ['Params_v_sample', '08_taxonomy_BAT', 'fraction'], 0.50, float,
            '\b\n08_taxonomy_BAT: Minimum bit-score support fraction [0-0.99].\n'
            'Higher values enforce stricter consensus (one classification per sequence), while lower values allow multiple lineage detections.'
        ),
        'virus_vmags_taxonomy_bat_threads': (
            ['Params_v_sample', '08_taxonomy_BAT', 'threads'], 8, int,
            '08_taxonomy_BAT: Threads.'),
        'virus_vmags_taxonomy_bat_extra_args': (
            ['Params_v_sample', '08_taxonomy_BAT', 'extra_args'], '', str,
            '08_taxonomy_BAT: Extra arguments forwarded to CAT_pack bins (e.g., "--block_size 12.0").'),
        'virus_coverm_method_sample': (
            ['Params_v_sample', '09_coverm_eachsample_vMAGs', 'method'], 'trimmed_mean', str,
            '\b\n09_coverm_eachsample_vMAGs: Method(s) for coverm-based coverage calculation, must be one of [trimmed_mean/mean].\n'
            'A more detailed description of the different methods is available at https://github.com/wwood/CoverM#calculation-methods'
        ),
        'virus_coverm_min_read_percent_identity_sample': (
            ['Params_v_sample', '09_coverm_eachsample_vMAGs', 'min_read_percent_identity'], 95, float,
            '09_coverm_eachsample_vMAGs: Exclude reads by overall percent identity [0-100], e.g. 95 for 95%.'
        ),
        'virus_coverm_min_read_aligned_length_sample': (
            ['Params_v_sample', '09_coverm_eachsample_vMAGs', 'min_read_aligned_length'], 45, int,
            '09_coverm_eachsample_vMAGs: Exclude reads with fewer aligned bases.'
        ),
        'virus_coverm_contig_end_exclusion_sample': (
            ['Params_v_sample', '09_coverm_eachsample_vMAGs', 'contig_end_exclusion'], 75, int,
            '09_coverm_eachsample_vMAGs: Exclude bases at the ends of contig sequences for coverm calculation.'
        ),
        'virus_coverm_threads_sample': (
            ['Params_v_sample', '09_coverm_eachsample_vMAGs', 'threads'], 8, int,
            '09_coverm_eachsample_vMAGs: Threads.'),
        'virus_coverm_sample_extra_args': (
            ['Params_v_sample', '09_coverm_eachsample_vMAGs', 'extra_args'], '', str,
            '09_coverm_eachsample_vMAGs: Extra arguments forwarded to coverm genome (e.g., "--quiet").'),
        'contigs_taxonomy_classify_method': (
            ['Params_c_sample', '01_taxonomy', 'contigs_taxonomy_classify_method'], 'kraken2', str,
            '01_taxonomy: Contigs taxonomy classifier using one of [CAT_pack/kraken2].'
        ),
        'contigs_taxonomy_threads': (
            ['Params_c_sample', '01_taxonomy', 'threads'], 60, int,
            '01_taxonomy: Threads for [CAT_pack/kraken2].'
        ),
        'contigs_taxonomy_cat_range': (
            ['Params_c_sample', '01_taxonomy_CAT_pack', 'range'], 10, float,
            '\b\n01_taxonomy_CAT_pack: Percentage range [0-100] of DIAMOND hits relative to top bit-score to include in LCA calculation.\n'
            'Higher values yield more conservative classifications by incorporating divergent hits, while lower values enable specific assignments.\n'
            'Applicable when --contigs-taxonomy-classify-method is set to CAT_pack.'
        ),
        'contigs_taxonomy_cat_fraction': (
            ['Params_c_sample', '01_taxonomy_CAT_pack', 'fraction'], 0.50, float,
            '\b\n01_taxonomy_CAT_pack: Minimum bit-score support fraction [0-0.99].\n'
            'Higher values enforce stricter consensus (one classification per sequence), while lower values allow multiple lineage detections.\n'
            'Applicable when --contigs-taxonomy-classify-method is set to CAT_pack.'
        ),
        'contigs_taxonomy_cat_pack_extra_args': (
            ['Params_c_sample', '01_taxonomy_CAT_pack', 'extra_args'], '', str,
            '01_taxonomy_CAT_pack: Extra arguments forwarded to CAT_pack contigs (e.g., "--block_size 12.0").'),
        'contigs_taxonomy_kraken2_confidence': (
            ['Params_c_sample', '01_taxonomy_kraken2', 'confidence'], 0.7, float,
            '\b\n01_taxonomy_kraken2: Confidence threshold [0-1] for taxonomic classification, applicable when --contigs-taxonomy-classify-method is set to kraken2.\n'
            'Default value was tuned on the CAMI2 marine dataset. For samples with higher complexity, a lower confidence threshold is usually recommended.'
        ),
        'contigs_taxonomy_kraken2_extra_args': (
            ['Params_c_sample', '01_taxonomy_kraken2', 'extra_args'], '', str,
            '01_taxonomy_kraken2: Extra arguments forwarded to kraken2 (e.g., "--memory-mapping").'),
        'contigs_function_classify_method': (
            ['Params_c_sample', '02_function', 'contigs_function_classify_method'], 'kraken2', str,
            '02_function: Contigs function classifier using one of [kraken2/emapper].'
        ),
        'contigs_function_threads': (
            ['Params_c_sample', '02_function', 'threads'], 60, int,
            '02_function: Threads for [kraken2/emapper].'),
        'contigs_function_emapper_evalue': (
            ['Params_c_sample', '02_function_emapper', 'seed_ortholog_evalue'], 1e-10, float,
            '02_function_emapper: Report only alignments with e-value below or equal to the threshold, applicable when --contigs-function-classify-method is set to emapper.'
        ),
        'contigs_function_emapper_method': (
            ['Params_c_sample', '02_function_emapper', 'method'], 'diamond', str,
            '\b\n02_function_emapper: search seed orthologs using one of [diamond/mmseqs], applicable when --contigs-function-classify-method is set to emapper.\n'
            'diamond: using default parameters of diamond mode of emapper.\n'
            'mmseqs: using default parameters of mmseqs mode of emapper'
        ),
        'contigs_function_emapper_query_cover': (
            ['Params_c_sample', '02_function_emapper', 'query_cover'], 30, float,
            '\b\n02_function_emapper: Report only alignments with query cover above or equal to the given percentage [0-100],.\n'
            'applicable when --contigs-function-classify-method is set to emapper'
        ),
        'contigs_function_emapper_extra_args': (
            ['Params_c_sample', '02_function_emapper', 'extra_args'], '', str,
            '02_function_emapper: Extra arguments forwarded to eggNOG-mapper (e.g., "--pident 30").'),
        'contigs_function_kraken2_confidence': (
            ['Params_c_sample', '02_function_kraken2', 'confidence'], 0.05, float,
            '02_function_kraken2: Confidence threshold [0-1] for functional classification, applicable when --contigs-function-classify-method is set to kraken2.'
        ),
        'contigs_function_kraken2_extra_args': (
            ['Params_c_sample', '02_function_kraken2', 'extra_args'], '', str,
            '02_function_kraken2: Extra arguments forwarded to kraken2 (e.g., "--memory-mapping").'),
        'merge_min_orf_length': (
            ['Params_vc_sample', '05_merge_tax_fun_depth', 'minORFLen'], 50, int,
            '05_merge_tax_fun_depth: Minimum ORF length for finally used ORFs of vMAGs and contigs.'
        ),
        'merge_vmags_checkv_quality': (
            ['Params_vc_sample', '05_merge_tax_fun_depth', 'checkv_quality'], 'High-quality', str,
            '\b\n05_merge_tax_fun_depth: CheckV quality filter for viral MAGs (lower bound filter).\n'
            'viral MAGs with quality equal to or higher than the specified threshold will be retained.\n'
            'Options: Complete, High-quality, Medium-quality, Low-quality, Not-determined.\n'
            'Complete: complete genome (highest quality)\n'
            'High-quality: genome completeness ≥90%\n'
            'Medium-quality: genome completeness 50-90%\n'
            'Low-quality: genome completeness <50%\n'
            'Not-determined: completeness could not be determined (lowest quality)'
        ),
        'merge_multiple_ko_handling_mode': (
            ['Params_vc_sample', '05_merge_tax_fun_depth', 'KO_mode'], 'ignore_multiple_KOs', str,
            '\b\n05_merge_tax_fun_depth: Method using one of [split_multiple_KOs/keep_multiple_KOs/ignore_multiple_KOs] when an ORF has multiple KO annotations:.\n'
            'split_multiple_KOs: evenly distributing the total abundance across all assigned KOs;\n'
            'keep_multiple_KOs: retain original multiple KO annotations;\n'
            'ignore_multiple_KOs: ignore the ORFs with multiple KO annotations and only use the ORFs with single KO annotation.'
        ),
        'merge_min_ko_num_per_taxon': (
            ['Params_vc_sample', '05_merge_tax_fun_depth', 'min_KO_num_in_one_taxon'], 100, int,
            '\b\n05_merge_tax_fun_depth: Minimum number of KOs per taxon, used to filter outputs in KO_Taxa_RelativeAbundance and KO_Taxa_SumCoverage.\n'
            'This threshold filters out taxa with insufficient KOs,\n'
            'retaining only KO-rich taxa to ensure robust metabolic reconstruction\n'
            'while preventing excessive output file sizes caused by sparse matrices.'
        ),
        'merge_tax_fun_depth_threads': (
            ['Params_vc_sample', '05_merge_tax_fun_depth', 'threads'], 20, int,
            '05_merge_tax_fun_depth: Threads.'),
        'merge_cami_format_output': (
            ['Params_vc_sample', '05_merge_tax_fun_depth', 'CAMI_format_output'], False, bool,
            '05_merge_tax_fun_depth: Report a taxonomic profile using the CAMI output format (flag only).'
        ),
        'merge_output_detail': (
            ['Params_vc_sample', '05_merge_tax_fun_depth', 'output_detail'], False, bool,
            '05_merge_tax_fun_depth: Whether to output detail tables (flag only).'),
        'mclink_sample_min_ko_num_per_taxon': (
            ['Params_vc_sample', '06_mclink_taxon', 'min_KO_num_in_one_taxon'], 100, int,
            '06_mclink_taxon: Minimum number of KOs per taxon for calculating taxon-level metabolic pathway/module completeness or abundance per sample. Taxa with fewer than this number of detected KOs in a given sample will be filtered out. Calculated independently for each sample.'
        ),
        'mclink_sample_min_ko_abundance_per_taxon': (
            ['Params_vc_sample', '06_mclink_taxon', 'min_mean_KO_RA_in_one_taxon'], 0, float,
            '06_mclink_taxon: Minimum mean KO abundance [0-1e6] per taxon for calculating taxon-level metabolic pathway/module completeness or abundance per sample. Taxa with mean KO abundance below this threshold in a given sample will be filtered out. Note: GPM/TPM sum of all taxa in each sample is 1e6.'
        ),
        'mclink_sample_topn_taxon': (
            ['Params_vc_sample', '06_mclink_taxon', 'TopN'], 100, int,
            '06_mclink_taxon: Number of top taxa (ranked by mean KO abundance) retained for calculating taxon-level metabolic pathway/module completeness or abundance per sample. Set to 0 to disable filtering.'
        ),
        'mclink_sample_completeness_per_taxon': (
            ['Params_vc_sample', '06_mclink_taxon', 'completeness_threshold'], 0.5, float,
            '06_mclink_taxon: Completeness threshold [0-1] of modules for calculating taxon-level metabolic pathway/module abundance. Higher values retain only more complete profiles.'),
        'mclink_sample_plus_scale_method_per_taxon': (
            ['Params_vc_sample', '06_mclink_taxon', 'plus_scale_method'], 'mean', str,
            "06_mclink_taxon: Aggregation method for '+'-separated KOs in a module [mean/max/min]. 'mean': average abundance across all required KOs; 'max': highest KO abundance; 'min': lowest KO abundance."),
        'mclink_sample_comma_scale_method_per_taxon': (
            ['Params_vc_sample', '06_mclink_taxon', 'comma_scale_method'], 'sum', str,
            "06_mclink_taxon: Aggregation method for ','-separated alternative KOs in a module [sum/max]. 'sum': cumulative abundance of all alternative KOs; 'max': abundance of the highest single KO among alternatives. Note: this applies to abundance calculation only; module completeness is always computed with 'max'."),
        'barplot_topn_taxon': (
            ['Params_vc_project', '01_Taxa_visualize', 'TopN'], 10, int,
            '01_Taxa_visualize: Number of top taxa to display in stacked barplots. Remaining categories are merged into Others. Maximum 29.'
        ),
        'barplot_topn_fun': (
            ['Params_vc_project', '03_Fun_visualize', 'TopN'], 10, int,
            '03_Fun_visualize: Number of top functional categories to display in stacked barplots. Remaining categories are merged into Others. Maximum 29.'
        ),
        'mclink_min_ko_num_per_sample': (
            ['Params_vc_project', '03_mclink_sample', 'min_sample_num_per_KO'], 1, int,
            '03_mclink_sample: Minimum number of samples in which a KO must be detected for calculating community-level metabolic pathway/module completeness or abundance.'
        ),
        'mclink_min_ko_prevalence_per_sample': (
            ['Params_vc_project', '03_mclink_sample', 'preva'], 0, float,
            '03_mclink_sample: Minimum detection rate [0-1] of a KO across all samples for calculating community-level metabolic pathway/module completeness or abundance.'
        ),
        'mclink_min_ko_abundance_per_sample': (
            ['Params_vc_project', '03_mclink_sample', 'abund'], 0, float,
            '03_mclink_sample: Minimum mean relative abundance [0-1e6] of a KO across all samples for calculating community-level metabolic pathway/module completeness or abundance. KOs with mean abundance below this threshold are excluded. Note: GPM/TPM sum of a sample is 1e6.'
        ),
        'mclink_completeness_per_sample': (
            ['Params_vc_project', '03_mclink_sample', 'completeness_threshold'], 0.5, float,
            '03_mclink_sample: Completeness threshold [0-1] of modules for calculating community-level metabolic pathway/module abundance. Higher values retain only more complete profiles.'
        ),
        'mclink_plus_scale_method_per_sample': (
            ['Params_vc_project', '03_mclink_sample', 'plus_scale_method'], 'mean', str,
            "03_mclink_sample: Aggregation method for '+'-separated KOs in a module [mean/max/min]. 'mean': average abundance across all required KOs; 'max': highest KO abundance; 'min': lowest KO abundance."
        ),
        'mclink_comma_scale_method_per_sample': (
            ['Params_vc_project', '03_mclink_sample', 'comma_scale_method'], 'sum', str,
            "03_mclink_sample: Aggregation method for ','-separated alternative KOs in a module [sum/max]. 'sum': cumulative abundance of all alternative KOs; 'max': abundance of the highest single KO among alternatives. Note: this applies to abundance calculation only; module completeness is always computed with 'max'."
        ),
        'merge_min_taxa_ko_prevalence_across_samples': (
            ['Params_vc_project', '05_merge_tax_KO_depth', 'preva'], 0.10, float,
            '05_merge_tax_KO_depth: Minimum prevalence [0-1] of a Taxon-KO across samples. Taxa-KOs detected in fewer samples than this proportion are excluded.'
        ),
        'merge_min_taxa_ko_abundance_across_samples': (
            ['Params_vc_project', '05_merge_tax_KO_depth', 'abund'], 0.0, float,
            '05_merge_tax_KO_depth: Minimum mean relative abundance [0-1e6] of a Taxon-KO across samples. Taxa-KOs with mean abundance below this threshold are excluded. Note: GPM/TPM sum of a sample is 1e6.'
        ),
        'merge_min_taxa_module_prevalence_across_samples': (
            ['Params_vc_project', '05_merge_tax_Module_RelativeAbundance', 'preva'], 0.10, float,
            '05_merge_tax_Module_RelativeAbundance: Minimum prevalence [0-1] of a Taxon-Module across samples. Taxa-Modules detected in fewer samples than this proportion are excluded.'
        ),
        'merge_min_taxa_module_abundance_across_samples': (
            ['Params_vc_project', '05_merge_tax_Module_RelativeAbundance', 'abund'], 0.0, float,
            '05_merge_tax_Module_RelativeAbundance: Minimum mean relative abundance [0-1e6] of a Taxon-Module across samples. Taxa-Modules with mean abundance below this threshold are excluded. Note: GPM/TPM sum of a sample is 1e6.'
        ),
        'merge_min_taxa_module_sc_prevalence_across_samples': (
            ['Params_vc_project', '05_merge_tax_Module_SumCoverage', 'preva'], 0.10, float,
            '05_merge_tax_Module_SumCoverage: Minimum prevalence [0-1] of a Taxon-Module across samples. Taxa-Modules detected in fewer samples than this proportion are excluded.'
        ),
        'merge_min_taxa_module_sc_abundance_across_samples': (
            ['Params_vc_project', '05_merge_tax_Module_SumCoverage', 'abund'], 0.0, float,
            '05_merge_tax_Module_SumCoverage: Minimum mean sum coverage [0-1e6] of a Taxon-Module across samples. Taxa-Modules with mean sum coverage below this threshold are excluded. Note: GPM/TPM sum of a sample is 1e6.'
        ),
        'mclink_project_min_ko_num_per_taxon': (
            ['Params_vc_project', '06_mclink_taxon', 'min_KO_num_in_one_taxon'], 100, int,
            '06_mclink_taxon: Minimum number of KOs per taxon for calculating taxon-level metabolic pathway/module completeness on the project-averaged taxon-function matrix. Taxa with fewer than this number of detected KOs in the project-averaged matrix will be filtered out. KO abundances are first averaged across all samples to generate a single project-level matrix before calculation.'
        ),
        'mclink_project_min_ko_abundance_per_taxon': (
            ['Params_vc_project', '06_mclink_taxon', 'min_mean_KO_RA_in_one_taxon'], 0, float,
            '06_mclink_taxon: Minimum mean KO abundance [0-1e6] per taxon for calculating taxon-level metabolic pathway/module completeness on the project-averaged taxon-function matrix. Taxa with mean KO abundance below this threshold in the project-averaged matrix will be filtered out. Note: GPM/TPM sum of all taxa is 1e6 per sample; values are averaged across samples prior to filtering and calculation.'
        ),
        'mclink_project_topn_taxon': (
            ['Params_vc_project', '06_mclink_taxon', 'TopN'], 100, int,
            '06_mclink_taxon: Number of top taxa (ranked by mean KO abundance in the project-averaged taxon-function matrix) retained for calculating taxon-level metabolic pathway/module completeness. Set to 0 to disable filtering.'
        ),
        'mclink_project_completeness_per_taxon': (
            ['Params_vc_project', '06_mclink_taxon', 'completeness_threshold'], 0.5, float,
            '06_mclink_taxon: Completeness threshold [0-1] of modules for calculating taxon-level metabolic pathway/module abundance on the project-averaged taxon-function matrix. Higher values retain only more complete profiles.'),
        'mclink_project_plus_scale_method_per_taxon': (
            ['Params_vc_project', '06_mclink_taxon', 'plus_scale_method'], 'mean', str,
            "06_mclink_taxon: Aggregation method for '+'-separated KOs in a module [mean/max/min] on the project-averaged taxon-function matrix. 'mean': average abundance across all required KOs; 'max': highest KO abundance; 'min': lowest KO abundance."),
        'mclink_project_comma_scale_method_per_taxon': (
            ['Params_vc_project', '06_mclink_taxon', 'comma_scale_method'], 'sum', str,
            "06_mclink_taxon: Aggregation method for ','-separated alternative KOs in a module [sum/max] on the project-averaged taxon-function matrix. 'sum': cumulative abundance of all alternative KOs; 'max': abundance of the highest single KO among alternatives. Note: this applies to abundance calculation only; module completeness is always computed with 'max'."),
    },
    'vcr': {
        'raw_stats_sample_threads': (
            ['Params_b_sample', '00_raw_stats', 'threads'], 4, int,
            '00_raw_stats: Threads for seqkit.'),
        'trim_sample_threads': (
            ['Params_b_sample', '01_trim', 'threads'], 16, int,
            '01_trim: Threads for fastp (next-generation short reads) or fastplong (third-generation long reads).'),
        'trim_length_required': (
            ['Params_b_sample', '01_trim', 'length_required'], 15, int,
            '01_trim: Reads shorter than length_required will be discarded'),
        'trim_n_base_limit': (
            ['Params_b_sample', '01_trim', 'n_base_limit'], 5, int,
            "01_trim: if one read's number of N base is >n_base_limit, then this read/pair is discarded."),
        'trim_compression': (
            ['Params_b_sample', '01_trim', 'compression'], 4, int,
            '01_trim: Compression level for gzip output [1-9]. 1 is fastest, 9 is smallest.'),
        'trim_extra_args': (
            ['Params_b_sample', '01_trim', 'extra_args'], '', str,
            '01_trim: Extra arguments forwarded to fastp or fastplong for raw reads trimming (e.g., "--reads_to_process 0 -Q").'),
        'trim_stats_sample_threads': (
            ['Params_b_sample', '01_trim_stats', 'threads'], 16, int,
            '01_trim_stats: Threads.'),
        'assembly_method': (
            ['Params_b_sample', '02_assembly', 'method'], 'megahit', str,
            '\b\n02_assembly: Underlying assembly software using one of [megahit/spades/metaflye].\n'
            'Must be compatible with --read-type (config key read_type):\n'
            'megahit: MEGAHIT (second-generation sequencing; for paired or single DNA reads).\n'
            'spades: SPAdes (second-generation sequencing; metaSPAdes for paired DNA reads, rnaSPAdes for paired-rna or single-rna RNA reads).\n'
            'metaflye: metaFlye (third-generation sequencing; required for pacbio-raw/corr/hifi or nano-raw/corr/hq long reads).'),
        'assembly_threads': (
            ['Params_b_sample', '02_assembly', 'threads'], 40, int,
            '02_assembly: Threads for [megahit/spades/metaflye].'),
        'megahit_presets': (
            ['Params_b_sample', '02_megahit', 'presets'], 'meta-sensitive', str,
            '\b\n02_megahit: MEGAHIT presets using one of [meta-sensitive/meta-large].\n'
            "meta-sensitive: '--min-count 1 --k-list 21,29,39,49,...,129,141'.\n"
            "meta-large: '--k-min 27 --k-max 127 --k-step 10'.\n"
            '(large & complex metagenomes, like soil).'),
        'megahit_k_list': (
            ['Params_b_sample', '02_megahit', 'k_list'], '', str,
            '02_megahit: MEGAHIT --k-list (optional). Explicit comma-separated k-mer sizes (e.g. "21,29,39,59,79,99,119,141"). Overrides --megahit-presets. If set, --megahit-k-min, --megahit-k-max, --megahit-k-step must NOT be set. Optional --megahit-min-count can be used together.'),
        'megahit_min_count': (
            ['Params_b_sample', '02_megahit', 'min_count'], 2, int,
            '02_megahit: MEGAHIT --min-count (optional). Minimum multiplicity for filtering (k_min+1)-mers. Can be used with --megahit-k-list or --megahit-k-min/--k-max/--k-step.'),
        'megahit_k_min': (
            ['Params_b_sample', '02_megahit', 'k_min'], '', str,
            '02_megahit: MEGAHIT --k-min (optional). Minimum kmer size (<= 255), must be odd number. If set with --k-max and --k-step, overrides --megahit-presets and uses custom k-mer mode.'),
        'megahit_k_max': (
            ['Params_b_sample', '02_megahit', 'k_max'], '', str,
            '02_megahit: MEGAHIT --k-max (optional). Maximum kmer size (<= 255), must be odd number.'),
        'megahit_k_step': (
            ['Params_b_sample', '02_megahit', 'k_step'], '', str,
            '02_megahit: MEGAHIT --k-step (optional). Increment of kmer size of each iteration (<= 28), must be even number.'),
        'megahit_extra_args': (
            ['Params_b_sample', '02_megahit', 'extra_args'], '', str,
            '02_megahit: Extra arguments forwarded to MEGAHIT (e.g. "--no-mercy --tmp-dir /path/to/your/tmp/dir"). Do NOT use this for k-mer settings (--k-list, --k-min, --presets, etc.).'),
        'spades_extra_args': (
            ['Params_b_sample', '02_spades', 'extra_args'], '', str,
            '02_spades: Extra arguments forwarded to SPAdes. Do NOT use this for input/output or thread settings.'),
        'metaflye_extra_args': (
            ['Params_b_sample', '02_metaflye', 'extra_args'], '', str,
            '02_metaflye: Extra arguments forwarded to Flye for third-generation long reads assembly (e.g., "--iterations 2"). Do NOT use this for read type or input/output settings.'),
        'mapping_threads': (
            ['Params_b_sample', '03_mapping', 'threads'], 20, int,
            '03_mapping: Threads for [minibwa/bwamem/bwamem2/strobealign].'),
        'mapping_method': (
            ['Params_b_sample', '03_mapping', 'method'], 'minibwa', str,
            '\b\n03_mapping: Underlying mapping software used [minibwa/bwamem/bwamem2/strobealign].\n'
            'minibwa: minibwa using default parameters.\n'
            'bwamem: bwa mem using default parameters.\n'
            'bwamem2: bwa-mem2 using default parameters.\n'
            'strobealign: strobealign using default parameters.'),
        'mapping_samtools_sort_mem': (
            ['Params_b_sample', '03_mapping', 'samtools_sort_mem'], '768M', str,
            '03_mapping: Set maximum memory per thread; suffix K/M/G recognized, e.g. 2G, 768M.'),
        'mapping_extra_args': (
            ['Params_b_sample', '03_mapping', 'extra_args'], '', str,
            '03_mapping: Extra arguments forwarded to the active mapping aligner (minibwa map, bwa mem, bwa-mem2 mem, or strobealign --use-index), (e.g. "-B 8" for minibwa map).'),
        'contigs_depth_sample_threads': (
            ['Params_b_sample', '07_contigs_depth', 'threads'], 8, int,
            '07_contigs_depth: Threads.'),
        'depth_min_read_percent_identity': (
            ['Params_b_sample', '07_contigs_depth', 'min_read_percent_identity'], 95, float,
            '\b\n07_contigs_depth: Exclude reads by overall percent identity [0-100] for contig depth estimation\n'
            'using jgi_summarize_bam_contig_depths --percentIdentity, e.g. 95 for 95%.\n'
            'Reads with percent identity below this threshold are excluded from depth calculation.\n'
            'Higher values increase stringency, retaining only high-quality alignments.'),
        'depth_contig_end_exclusion': (
            ['Params_b_sample', '07_contigs_depth', 'contig_end_exclusion'], 75, int,
            '\b\n07_contigs_depth: Exclude bases at the ends of contig sequences for depth estimation using jgi_summarize_bam_contig_depths --maxEdgeBases.\n'
            'Bases within this distance from either end of a contig are ignored when calculating depth,\n'
            'preventing edge effects from partial alignments.'),
        'depth_min_contig_length': (
            ['Params_b_sample', '07_contigs_depth', 'min_contig_length'], 300, int,
            '\b\n07_contigs_depth: Minimum contig length for depth estimation [bp], e.g. 300 for 300 bp.\n'
            'Contigs shorter than this threshold are excluded from depth calculation.'),
        'contigs_prodigal_extra_args': (
            ['Params_b_sample', '07_contigs_depth', 'extra_args'], '', str,
            '07_contigs_depth: Extra arguments forwarded to prodigal for contig gene prediction (e.g., "-m -n").'),
        'raw_stats_project_threads': (
            ['Params_b_project', '01_cat_raw_stats', 'threads'], 4, int,
            '01_cat_raw_stats: Threads.'),
        'trimmed_stats_project_threads': (
            ['Params_b_project', '02_cat_trimmed_stats', 'threads'], 4, int,
            '02_cat_trimmed_stats: Threads.'),
        'contigs_stats_project_threads': (
            ['Params_b_project', '03_contigs_stats', 'threads'], 4, int,
            '03_contigs_stats: Threads.'),
        'merge_stats_project_threads': (
            ['Params_b_project', '04_merge_reads_contigs_stats', 'threads'], 4, int,
            '04_merge_reads_contigs_stats: Threads.'),
        'virus_identify_method': (
            ['Params_v_sample', '01_virus_identify', 'identify_method'], 'genomad_vs_vb', str,
            '\b\n01_virus_identify: Virus identification method, should be one of the [genomad/vs_vb/genomad_vs_vb].\n'
            'genomad: Use geNomad;\n'
            'vs_vb: Use VirSorter2 and CheckV, then use VIBRANT to check the low quality viruses from VirSorter2;\n'
            'genomad_vs_vb: Combine genomad and vs_vb virus to get the union results.'
        ),
        'virus_length_limit': (
            ['Params_v_sample', '01_virus_identify', 'virus_length_limit'], 2000, int,
            '01_virus_identify: Minimum length (bp) for viral contigs.'
        ),
        'virus_identify_threads': (
            ['Params_v_sample', '01_virus_identify', 'threads'], 8, int,
            '01_virus_identify: Threads.'),
        'virus_contig_taxonomy_cat_range': (
            ['Params_v_sample', '03_taxonomy_CAT_pack', 'range'], 10, float,
            '\b\n03_taxonomy_CAT_pack: Percentage range [0-100] of DIAMOND hits relative to top bit-score to include in LCA calculation.\n'
            'Higher values yield more conservative classifications by incorporating divergent hits, while lower values enable specific assignments.'
        ),
        'virus_contig_taxonomy_cat_fraction': (
            ['Params_v_sample', '03_taxonomy_CAT_pack', 'fraction'], 0.50, float,
            '\b\n03_taxonomy_CAT_pack: Minimum bit-score support fraction [0-0.99].\n'
            'Higher values enforce stricter consensus (one classification per sequence), while lower values allow multiple lineage detections.'
        ),
        'virus_contig_taxonomy_cat_pack_threads': (
            ['Params_v_sample', '03_taxonomy_CAT_pack', 'threads'], 8, int,
            '03_taxonomy_CAT_pack: Threads.'),
        'virus_contig_taxonomy_cat_pack_extra_args': (
            ['Params_v_sample', '03_taxonomy_CAT_pack', 'extra_args'], '', str,
            '03_taxonomy_CAT_pack: Extra arguments forwarded to CAT_pack contigs (e.g., "--block_size 12.0").'),
        'virus_mapping_threads': (
            ['Params_v_sample', '03_mapping', 'threads'], 20, int,
            '03_mapping: Threads for [minibwa/bwamem/bwamem2/strobealign].'),
        'virus_mapping_method': (
            ['Params_v_sample', '03_mapping', 'method'], 'minibwa', str,
            '\b\n03_mapping: Underlying mapping software used [minibwa/bwamem/bwamem2/strobealign].\n'
            'minibwa: minibwa using default parameters.\n'
            'bwamem: bwa mem using default parameters.\n'
            'bwamem2: bwa-mem2 using default parameters.\n'
            'strobealign: strobealign using default parameters.'),
        'virus_mapping_samtools_sort_mem': (
            ['Params_v_sample', '03_mapping', 'samtools_sort_mem'], '768M', str,
            '03_mapping: Set maximum memory per thread; suffix K/M/G recognized, e.g. 2G, 768M.'),
        'virus_mapping_extra_args': (
            ['Params_v_sample', '03_mapping', 'extra_args'], '', str,
            '03_mapping: Extra arguments forwarded to the active mapping aligner (minibwa map, bwa mem, bwa-mem2 mem, or strobealign --use-index), (e.g. "-B 8" for minibwa map).'),
        'virus_binning_threads': (
            ['Params_v_sample', '04_binning', 'threads'], 8, int,
            '04_binning_vRhyme: Threads.'),
        'virus_checkv_threads': (
            ['Params_v_sample', '05_checkv', 'threads'], 8, int,
            '05_checkv: Threads.'),
        'virus_vclust_ani_sample': (
            ['Params_v_sample', '06_vclust', 'ani'], 0.95, float,
            '06_vclust: ANI threshold [0-1] of vclust for sample-level virus clustering.'
        ),
        'virus_vclust_threads_sample': (
            ['Params_v_sample', '06_vclust', 'threads'], 8, int,
            '06_vclust: Threads.'),
        'vclust_sample_prefilter_extra_args': (
            ['Params_v_sample', '06_vclust', 'prefilter_extra_args'], '', str,
            '06_vclust: Extra arguments forwarded to vclust prefilter (e.g., "-k 25").'),
        'vclust_sample_align_extra_args': (
            ['Params_v_sample', '06_vclust', 'align_extra_args'], '', str,
            '06_vclust: Extra arguments forwarded to vclust align (e.g., "--mal 11").'),
        'vclust_sample_cluster_extra_args': (
            ['Params_v_sample', '06_vclust', 'cluster_extra_args'], '', str,
            '06_vclust: Extra arguments forwarded to vclust cluster (e.g., "--len_ratio 0").'),
        'virus_iphop_threads_sample': (
            ['Params_v_sample', '07_iPHoP', 'threads'], 20, int,
            '07_iPHoP: Threads.'),
        'virus_vmags_taxonomy_bat_range': (
            ['Params_v_sample', '08_taxonomy_BAT', 'range'], 10, float,
            '\b\n08_taxonomy_BAT: Percentage range [0-100] of DIAMOND hits relative to top bit-score to include in LCA calculation.\n'
            'Higher values yield more conservative classifications by incorporating divergent hits, while lower values enable specific assignments.'
        ),
        'virus_vmags_taxonomy_bat_fraction': (
            ['Params_v_sample', '08_taxonomy_BAT', 'fraction'], 0.50, float,
            '\b\n08_taxonomy_BAT: Minimum bit-score support fraction [0-0.99].\n'
            'Higher values enforce stricter consensus (one classification per sequence), while lower values allow multiple lineage detections.'
        ),
        'virus_vmags_taxonomy_bat_threads': (
            ['Params_v_sample', '08_taxonomy_BAT', 'threads'], 8, int,
            '08_taxonomy_BAT: Threads.'),
        'virus_vmags_taxonomy_bat_extra_args': (
            ['Params_v_sample', '08_taxonomy_BAT', 'extra_args'], '', str,
            '08_taxonomy_BAT: Extra arguments forwarded to CAT_pack bins (e.g., "--block_size 12.0").'),
        'virus_coverm_method_sample': (
            ['Params_v_sample', '09_coverm_eachsample_vMAGs', 'method'], 'trimmed_mean', str,
            '\b\n09_coverm_eachsample_vMAGs: Method(s) for coverm-based coverage calculation, must be one of [trimmed_mean/mean].\n'
            'A more detailed description of the different methods is available at https://github.com/wwood/CoverM#calculation-methods'
        ),
        'virus_coverm_min_read_percent_identity_sample': (
            ['Params_v_sample', '09_coverm_eachsample_vMAGs', 'min_read_percent_identity'], 95, float,
            '09_coverm_eachsample_vMAGs: Exclude reads by overall percent identity [0-100], e.g. 95 for 95%.'
        ),
        'virus_coverm_min_read_aligned_length_sample': (
            ['Params_v_sample', '09_coverm_eachsample_vMAGs', 'min_read_aligned_length'], 45, int,
            '09_coverm_eachsample_vMAGs: Exclude reads with fewer aligned bases.'
        ),
        'virus_coverm_contig_end_exclusion_sample': (
            ['Params_v_sample', '09_coverm_eachsample_vMAGs', 'contig_end_exclusion'], 75, int,
            '09_coverm_eachsample_vMAGs: Exclude bases at the ends of contig sequences for coverm calculation.'
        ),
        'virus_coverm_threads_sample': (
            ['Params_v_sample', '09_coverm_eachsample_vMAGs', 'threads'], 8, int,
            '09_coverm_eachsample_vMAGs: Threads.'),
        'virus_coverm_sample_extra_args': (
            ['Params_v_sample', '09_coverm_eachsample_vMAGs', 'extra_args'], '', str,
            '09_coverm_eachsample_vMAGs: Extra arguments forwarded to coverm genome (e.g., "--quiet").'),
        'contigs_taxonomy_classify_method': (
            ['Params_c_sample', '01_taxonomy', 'contigs_taxonomy_classify_method'], 'kraken2', str,
            '01_taxonomy: Contigs taxonomy classifier using one of [CAT_pack/kraken2].'
        ),
        'contigs_taxonomy_threads': (
            ['Params_c_sample', '01_taxonomy', 'threads'], 60, int,
            '01_taxonomy: Threads for [CAT_pack/kraken2].'
        ),
        'contigs_taxonomy_cat_range': (
            ['Params_c_sample', '01_taxonomy_CAT_pack', 'range'], 10, float,
            '\b\n01_taxonomy_CAT_pack: Percentage range [0-100] of DIAMOND hits relative to top bit-score to include in LCA calculation.\n'
            'Higher values yield more conservative classifications by incorporating divergent hits, while lower values enable specific assignments.\n'
            'Applicable when --contigs-taxonomy-classify-method is set to CAT_pack.'
        ),
        'contigs_taxonomy_cat_fraction': (
            ['Params_c_sample', '01_taxonomy_CAT_pack', 'fraction'], 0.50, float,
            '\b\n01_taxonomy_CAT_pack: Minimum bit-score support fraction [0-0.99].\n'
            'Higher values enforce stricter consensus (one classification per sequence), while lower values allow multiple lineage detections.\n'
            'Applicable when --contigs-taxonomy-classify-method is set to CAT_pack.'
        ),
        'contigs_taxonomy_cat_pack_extra_args': (
            ['Params_c_sample', '01_taxonomy_CAT_pack', 'extra_args'], '', str,
            '01_taxonomy_CAT_pack: Extra arguments forwarded to CAT_pack contigs (e.g., "--block_size 12.0").'),
        'contigs_taxonomy_kraken2_confidence': (
            ['Params_c_sample', '01_taxonomy_kraken2', 'confidence'], 0.7, float,
            '\b\n01_taxonomy_kraken2: Confidence threshold [0-1] for taxonomic classification, applicable when --contigs-taxonomy-classify-method is set to kraken2.\n'
            'Default value was tuned on the CAMI2 marine dataset. For samples with higher complexity, a lower confidence threshold is usually recommended.'
        ),
        'contigs_taxonomy_kraken2_extra_args': (
            ['Params_c_sample', '01_taxonomy_kraken2', 'extra_args'], '', str,
            '01_taxonomy_kraken2: Extra arguments forwarded to kraken2 (e.g., "--memory-mapping").'),
        'contigs_function_classify_method': (
            ['Params_c_sample', '02_function', 'contigs_function_classify_method'], 'kraken2', str,
            '02_function: Contigs function classifier using one of [kraken2/emapper].'
        ),
        'contigs_function_threads': (
            ['Params_c_sample', '02_function', 'threads'], 60, int,
            '02_function: Threads for [kraken2/emapper].'),
        'contigs_function_emapper_evalue': (
            ['Params_c_sample', '02_function_emapper', 'seed_ortholog_evalue'], 1e-10, float,
            '02_function_emapper: Report only alignments with e-value below or equal to the threshold, applicable when --contigs-function-classify-method is set to emapper.'
        ),
        'contigs_function_emapper_method': (
            ['Params_c_sample', '02_function_emapper', 'method'], 'diamond', str,
            '\b\n02_function_emapper: search seed orthologs using one of [diamond/mmseqs], applicable when --contigs-function-classify-method is set to emapper.\n'
            'diamond: using default parameters of diamond mode of emapper.\n'
            'mmseqs: using default parameters of mmseqs mode of emapper'
        ),
        'contigs_function_emapper_query_cover': (
            ['Params_c_sample', '02_function_emapper', 'query_cover'], 30, float,
            '\b\n02_function_emapper: Report only alignments with query cover above or equal to the given percentage [0-100],.\n'
            'applicable when --contigs-function-classify-method is set to emapper'
        ),
        'contigs_function_emapper_extra_args': (
            ['Params_c_sample', '02_function_emapper', 'extra_args'], '', str,
            '02_function_emapper: Extra arguments forwarded to eggNOG-mapper (e.g., "--pident 30").'),
        'contigs_function_kraken2_confidence': (
            ['Params_c_sample', '02_function_kraken2', 'confidence'], 0.05, float,
            '02_function_kraken2: Confidence threshold [0-1] for functional classification, applicable when --contigs-function-classify-method is set to kraken2.'
        ),
        'contigs_function_kraken2_extra_args': (
            ['Params_c_sample', '02_function_kraken2', 'extra_args'], '', str,
            '02_function_kraken2: Extra arguments forwarded to kraken2 (e.g., "--memory-mapping").'),
        'vcr_unmapped_reads_merge_threads': (
            ['Params_vcr_sample', '01_unmapped_reads_trim_merge', 'threads'], 16, int,
            '01_unmapped_reads_trim_merge_fastp: Threads for fastp, samtools, and seqkit (paired-end), or samtools and seqkit only (single-end).'),
        'vcr_unmapped_reads_merge_length_required': (
            ['Params_vcr_sample', '01_unmapped_reads_trim_merge', 'length_required'], 15, int,
            '01_unmapped_reads_trim_merge_fastp: Reads shorter than length_required will be discarded. Only applicable when read_type is paired.'),
        'vcr_unmapped_reads_merge_n_base_limit': (
            ['Params_vcr_sample', '01_unmapped_reads_trim_merge', 'n_base_limit'], 5, int,
            "01_unmapped_reads_trim_merge_fastp: if one read's number of N base is >n_base_limit, then this read/pair is discarded. Only applicable when read_type is paired."),
        'vcr_merge_compression': (
            ['Params_vcr_sample', '01_unmapped_reads_trim_merge', 'compression'], 4, int,
            '01_unmapped_reads_trim_merge_fastp: Compression level for gzip output [1-9]. 1 is fastest, 9 is smallest. Only applicable when read_type is paired.'),
        'vcr_unmapped_reads_fastp_extra_args': (
            ['Params_vcr_sample', '01_unmapped_reads_trim_merge', 'extra_args'], '', str,
            '01_unmapped_reads_trim_merge: Extra arguments forwarded to fastp for unmapped reads trimming (e.g., "--reads_to_process 0 -Q"). Only applicable when read_type is paired.'),
        'vcr_fraggenescanrs_model': (
            ['Params_vcr_sample', '02_FragGeneScanRs', 'model'], 'illumina_1', str,
            '\b\n02_FragGeneScanRs: Model parameter used for short-read open reading frame prediction.\n'
            '[complete] for complete genomic sequences or short sequence reads without sequencing error\n'
            '[sanger_5] for Sanger sequencing reads with about 0.5% error rate\n'
            '[sanger_10] for Sanger sequencing reads with about 1% error rate\n'
            '[454_5] for 454 pyrosequencing reads with about 0.5% error rate\n'
            '[454_10] for 454 pyrosequencing reads with about 1% error rate\n'
            '[454_30] for 454 pyrosequencing reads with about 3% error rate\n'
            '[illumina_1] for Illumina sequencing reads with about 0.1% error rate\n'
            '[illumina_5] for Illumina sequencing reads with about 0.5% error rate\n'
            '[illumina_10] for Illumina sequencing reads with about 1% error rate'
        ),
        'vcr_fraggenescanrs_threads': (
            ['Params_vcr_sample', '02_FragGeneScanRs', 'threads'], 20, int,
            '02_FragGeneScanRs: Threads.'),
        'vcr_fraggenescanrs_stats_threads': (
            ['Params_vcr_sample', '02_FragGeneScanRs_stats', 'threads'], 10, int,
            '02_FragGeneScanRs_stats: Threads.'),
        'vcr_taxonomy_threads': (
            ['Params_vcr_sample', '03_taxonomy', 'threads'], 60, int,
            '03_taxonomy: Threads.'),
        'vcr_taxonomy_kraken2_confidence': (
            ['Params_vcr_sample', '03_taxonomy_kraken2', 'confidence'], 0.7, float,
            '\b\n03_taxonomy_kraken2: Confidence threshold [0-1] for taxonomic classification of unmapped reads with kraken2.\n'
            'Default value was tuned on the CAMI2 marine dataset. For samples with higher complexity, a lower confidence threshold is usually recommended.'
        ),
        'vcr_taxonomy_kraken2_extra_args': (
            ['Params_vcr_sample', '03_taxonomy_kraken2', 'extra_args'], '', str,
            '03_taxonomy_kraken2: Extra arguments forwarded to kraken2 (e.g., "--memory-mapping").'),
        'vcr_function_classify_method': (
            ['Params_vcr_sample', '04_function', 'reads_function_classify_method'], 'kraken2', str,
            '04_function: reads functional classifier using one of [kraken2/emapper].'
        ),
        'vcr_function_threads': (
            ['Params_vcr_sample', '04_function', 'threads'], 60, int,
            '04_function: Threads for [kraken2/emapper].'),
        'vcr_function_kraken2_confidence': (
            ['Params_vcr_sample', '04_function_kraken2', 'confidence'], 0.05, float,
            '04_function_kraken2: Confidence threshold [0-1] for functional classification for reads when --vcr-function-classify-method is set to kraken2.'
        ),
        'vcr_function_kraken2_extra_args': (
            ['Params_vcr_sample', '04_function_kraken2', 'extra_args'], '', str,
            '04_function_kraken2: Extra arguments forwarded to kraken2 (e.g., "--memory-mapping").'),
        'vcr_function_emapper_evalue': (
            ['Params_vcr_sample', '04_function_emapper', 'seed_ortholog_evalue'], 1e-10, float,
            '04_function_emapper: Report only alignments with e-value below or equal to the threshold, applicable when --vcr-function-classify-method is set to emapper.'
        ),
        'vcr_function_emapper_method': (
            ['Params_vcr_sample', '04_function_emapper', 'method'], 'diamond', str,
            '\b\n04_function_emapper: search seed orthologs using one of [diamond/mmseqs], applicable when --vcr-function-classify-method is set to emapper.\n'
            'diamond: using default parameters of diamond mode of emapper.\n'
            'mmseqs: using default parameters of mmseqs mode of emapper'
        ),
        'vcr_function_emapper_query_cover': (
            ['Params_vcr_sample', '04_function_emapper', 'query_cover'], 30, float,
            '\b\n04_function_emapper: Report only alignments with query cover above or equal to the given percentage [0-100],.\n'
            'applicable when --vcr-function-classify-method is set to emapper'
        ),
        'vcr_function_emapper_extra_args': (
            ['Params_vcr_sample', '04_function_emapper', 'extra_args'], '', str,
            '04_function_emapper: Extra arguments forwarded to eggNOG-mapper (e.g., "--pident 30").'),
        'merge_min_orf_length': (
            ['Params_vcr_sample', '05_merge_tax_fun_depth', 'minORFLen'], 50, int,
            '05_merge_tax_fun_depth: Minimum ORF length for finally used ORFs of vMAGs, contigs, and reads.'
        ),
        'merge_vmags_checkv_quality': (
            ['Params_vcr_sample', '05_merge_tax_fun_depth', 'checkv_quality'], 'High-quality', str,
            '\b\n05_merge_tax_fun_depth: CheckV quality filter for viral MAGs (lower bound filter).\n'
            'viral MAGs with quality equal to or higher than the specified threshold will be retained.\n'
            'Options: Complete, High-quality, Medium-quality, Low-quality, Not-determined.\n'
            'Complete: complete genome (highest quality)\n'
            'High-quality: genome completeness ≥90%\n'
            'Medium-quality: genome completeness 50-90%\n'
            'Low-quality: genome completeness <50%\n'
            'Not-determined: completeness could not be determined (lowest quality)'
        ),
        'merge_multiple_ko_handling_mode': (
            ['Params_vcr_sample', '05_merge_tax_fun_depth', 'KO_mode'], 'ignore_multiple_KOs', str,
            '\b\n05_merge_tax_fun_depth: Method using one of [split_multiple_KOs/keep_multiple_KOs/ignore_multiple_KOs] when an ORF has multiple KO annotations:.\n'
            'split_multiple_KOs: evenly distributing the total abundance across all assigned KOs;\n'
            'keep_multiple_KOs: retain original multiple KO annotations;\n'
            'ignore_multiple_KOs: ignore the ORFs with multiple KO annotations and only use the ORFs with single KO annotation.'
        ),
        'merge_min_ko_num_per_taxon': (
            ['Params_vcr_sample', '05_merge_tax_fun_depth', 'min_KO_num_in_one_taxon'], 100, int,
            '\b\n05_merge_tax_fun_depth: Minimum number of KOs per taxon, used to filter outputs in KO_Taxa_RelativeAbundance and KO_Taxa_SumCoverage.\n'
            'This threshold filters out taxa with insufficient KOs,\n'
            'retaining only KO-rich taxa to ensure robust metabolic reconstruction\n'
            'while preventing excessive output file sizes caused by sparse matrices.'
        ),
        'merge_tax_fun_depth_threads': (
            ['Params_vcr_sample', '05_merge_tax_fun_depth', 'threads'], 20, int,
            '05_merge_tax_fun_depth: Threads.'),
        'merge_cami_format_output': (
            ['Params_vcr_sample', '05_merge_tax_fun_depth', 'CAMI_format_output'], False, bool,
            '05_merge_tax_fun_depth: Report a taxonomic profile using the CAMI output format (flag only).'
        ),
        'merge_output_detail': (
            ['Params_vcr_sample', '05_merge_tax_fun_depth', 'output_detail'], False, bool,
            '05_merge_tax_fun_depth: Whether to output detail tables (flag only).'),
        'mclink_sample_min_ko_num_per_taxon': (
            ['Params_vcr_sample', '06_mclink_taxon', 'min_KO_num_in_one_taxon'], 100, int,
            '06_mclink_taxon: Minimum number of KOs per taxon for calculating taxon-level metabolic pathway/module completeness or abundance per sample. Taxa with fewer than this number of detected KOs in a given sample will be filtered out. Calculated independently for each sample.'
        ),
        'mclink_sample_min_ko_abundance_per_taxon': (
            ['Params_vcr_sample', '06_mclink_taxon', 'min_mean_KO_RA_in_one_taxon'], 0, float,
            '06_mclink_taxon: Minimum mean KO abundance [0-1e6] per taxon for calculating taxon-level metabolic pathway/module completeness or abundance per sample. Taxa with mean KO abundance below this threshold in a given sample will be filtered out. Note: GPM/TPM sum of all taxa in each sample is 1e6.'
        ),
        'mclink_sample_topn_taxon': (
            ['Params_vcr_sample', '06_mclink_taxon', 'TopN'], 100, int,
            '06_mclink_taxon: Number of top taxa (ranked by mean KO abundance) retained for calculating taxon-level metabolic pathway/module completeness or abundance per sample. Set to 0 to disable filtering.'
        ),
        'mclink_sample_completeness_per_taxon': (
            ['Params_vcr_sample', '06_mclink_taxon', 'completeness_threshold'], 0.5, float,
            '06_mclink_taxon: Completeness threshold [0-1] of modules for calculating taxon-level metabolic pathway/module abundance. Higher values retain only more complete profiles.'),
        'mclink_sample_plus_scale_method_per_taxon': (
            ['Params_vcr_sample', '06_mclink_taxon', 'plus_scale_method'], 'mean', str,
            "06_mclink_taxon: Aggregation method for '+'-separated KOs in a module [mean/max/min]. 'mean': average abundance across all required KOs; 'max': highest KO abundance; 'min': lowest KO abundance."),
        'mclink_sample_comma_scale_method_per_taxon': (
            ['Params_vcr_sample', '06_mclink_taxon', 'comma_scale_method'], 'sum', str,
            "06_mclink_taxon: Aggregation method for ','-separated alternative KOs in a module [sum/max]. 'sum': cumulative abundance of all alternative KOs; 'max': abundance of the highest single KO among alternatives. Note: this applies to abundance calculation only; module completeness is always computed with 'max'."),
        'barplot_topn_taxon': (
            ['Params_vcr_project', '01_Taxa_visualize', 'TopN'], 10, int,
            '01_Taxa_visualize: Number of top taxa to display in stacked barplots. Remaining categories are merged into Others. Maximum 29.'
        ),
        'barplot_topn_fun': (
            ['Params_vcr_project', '03_Fun_visualize', 'TopN'], 10, int,
            '03_Fun_visualize: Number of top functional categories to display in stacked barplots. Remaining categories are merged into Others. Maximum 29.'
        ),
        'mclink_min_ko_num_per_sample': (
            ['Params_vcr_project', '03_mclink_sample', 'min_sample_num_per_KO'], 1, int,
            '03_mclink_sample: Minimum number of samples in which a KO must be detected for calculating community-level metabolic pathway/module completeness or abundance.'
        ),
        'mclink_min_ko_prevalence_per_sample': (
            ['Params_vcr_project', '03_mclink_sample', 'preva'], 0, float,
            '03_mclink_sample: Minimum detection rate [0-1] of a KO across all samples for calculating community-level metabolic pathway/module completeness or abundance.'
        ),
        'mclink_min_ko_abundance_per_sample': (
            ['Params_vcr_project', '03_mclink_sample', 'abund'], 0, float,
            '03_mclink_sample: Minimum mean relative abundance [0-1e6] of a KO across all samples for calculating community-level metabolic pathway/module completeness or abundance. KOs with mean abundance below this threshold are excluded. Note: GPM/TPM sum of a sample is 1e6.'
        ),
        'mclink_completeness_per_sample': (
            ['Params_vcr_project', '03_mclink_sample', 'completeness_threshold'], 0.5, float,
            '03_mclink_sample: Completeness threshold [0-1] of modules for calculating community-level metabolic pathway/module abundance. Higher values retain only more complete profiles.'
        ),
        'mclink_plus_scale_method_per_sample': (
            ['Params_vcr_project', '03_mclink_sample', 'plus_scale_method'], 'mean', str,
            "03_mclink_sample: Aggregation method for '+'-separated KOs in a module [mean/max/min]. 'mean': average abundance across all required KOs; 'max': highest KO abundance; 'min': lowest KO abundance."
        ),
        'mclink_comma_scale_method_per_sample': (
            ['Params_vcr_project', '03_mclink_sample', 'comma_scale_method'], 'sum', str,
            "03_mclink_sample: Aggregation method for ','-separated alternative KOs in a module [sum/max]. 'sum': cumulative abundance of all alternative KOs; 'max': abundance of the highest single KO among alternatives. Note: this applies to abundance calculation only; module completeness is always computed with 'max'."
        ),
        'merge_min_taxa_ko_prevalence_across_samples': (
            ['Params_vcr_project', '05_merge_tax_KO_depth', 'preva'], 0.10, float,
            '05_merge_tax_KO_depth: Minimum prevalence [0-1] of a Taxon-KO across samples. Taxa-KOs detected in fewer samples than this proportion are excluded.'
        ),
        'merge_min_taxa_ko_abundance_across_samples': (
            ['Params_vcr_project', '05_merge_tax_KO_depth', 'abund'], 0.0, float,
            '05_merge_tax_KO_depth: Minimum mean relative abundance [0-1e6] of a Taxon-KO across samples. Taxa-KOs with mean abundance below this threshold are excluded. Note: GPM/TPM sum of a sample is 1e6.'
        ),
        'merge_min_taxa_module_prevalence_across_samples': (
            ['Params_vcr_project', '05_merge_tax_Module_RelativeAbundance', 'preva'], 0.10, float,
            '05_merge_tax_Module_RelativeAbundance: Minimum prevalence [0-1] of a Taxon-Module across samples. Taxa-Modules detected in fewer samples than this proportion are excluded.'
        ),
        'merge_min_taxa_module_abundance_across_samples': (
            ['Params_vcr_project', '05_merge_tax_Module_RelativeAbundance', 'abund'], 0.0, float,
            '05_merge_tax_Module_RelativeAbundance: Minimum mean relative abundance [0-1e6] of a Taxon-Module across samples. Taxa-Modules with mean abundance below this threshold are excluded. Note: GPM/TPM sum of a sample is 1e6.'
        ),
        'merge_min_taxa_module_sc_prevalence_across_samples': (
            ['Params_vcr_project', '05_merge_tax_Module_SumCoverage', 'preva'], 0.10, float,
            '05_merge_tax_Module_SumCoverage: Minimum prevalence [0-1] of a Taxon-Module across samples. Taxa-Modules detected in fewer samples than this proportion are excluded.'
        ),
        'merge_min_taxa_module_sc_abundance_across_samples': (
            ['Params_vcr_project', '05_merge_tax_Module_SumCoverage', 'abund'], 0.0, float,
            '05_merge_tax_Module_SumCoverage: Minimum mean sum coverage [0-1e6] of a Taxon-Module across samples. Taxa-Modules with mean sum coverage below this threshold are excluded. Note: GPM/TPM sum of a sample is 1e6.'
        ),
        'mclink_project_min_ko_num_per_taxon': (
            ['Params_vcr_project', '06_mclink_taxon', 'min_KO_num_in_one_taxon'], 100, int,
            '06_mclink_taxon: Minimum number of KOs per taxon for calculating taxon-level metabolic pathway/module completeness on the project-averaged taxon-function matrix. Taxa with fewer than this number of detected KOs in the project-averaged matrix will be filtered out. KO abundances are first averaged across all samples to generate a single project-level matrix before calculation.'
        ),
        'mclink_project_min_ko_abundance_per_taxon': (
            ['Params_vcr_project', '06_mclink_taxon', 'min_mean_KO_RA_in_one_taxon'], 0, float,
            '06_mclink_taxon: Minimum mean KO abundance [0-1e6] per taxon for calculating taxon-level metabolic pathway/module completeness on the project-averaged taxon-function matrix. Taxa with mean KO abundance below this threshold in the project-averaged matrix will be filtered out. Note: GPM/TPM sum of all taxa is 1e6 per sample; values are averaged across samples prior to filtering and calculation.'
        ),
        'mclink_project_topn_taxon': (
            ['Params_vcr_project', '06_mclink_taxon', 'TopN'], 100, int,
            '06_mclink_taxon: Number of top taxa (ranked by mean KO abundance in the project-averaged taxon-function matrix) retained for calculating taxon-level metabolic pathway/module completeness. Set to 0 to disable filtering.'
        ),
        'mclink_project_completeness_per_taxon': (
            ['Params_vcr_project', '06_mclink_taxon', 'completeness_threshold'], 0.5, float,
            '06_mclink_taxon: Completeness threshold [0-1] of modules for calculating taxon-level metabolic pathway/module abundance on the project-averaged taxon-function matrix. Higher values retain only more complete profiles.'),
        'mclink_project_plus_scale_method_per_taxon': (
            ['Params_vcr_project', '06_mclink_taxon', 'plus_scale_method'], 'mean', str,
            "06_mclink_taxon: Aggregation method for '+'-separated KOs in a module [mean/max/min] on the project-averaged taxon-function matrix. 'mean': average abundance across all required KOs; 'max': highest KO abundance; 'min': lowest KO abundance."),
        'mclink_project_comma_scale_method_per_taxon': (
            ['Params_vcr_project', '06_mclink_taxon', 'comma_scale_method'], 'sum', str,
            "06_mclink_taxon: Aggregation method for ','-separated alternative KOs in a module [sum/max] on the project-averaged taxon-function matrix. 'sum': cumulative abundance of all alternative KOs; 'max': abundance of the highest single KO among alternatives. Note: this applies to abundance calculation only; module completeness is always computed with 'max'."),
    },
    'vm': {
        'raw_stats_sample_threads': (
            ['Params_b_sample', '00_raw_stats', 'threads'], 4, int,
            '00_raw_stats: Threads for seqkit.'),
        'trim_sample_threads': (
            ['Params_b_sample', '01_trim', 'threads'], 16, int,
            '01_trim: Threads for fastp (next-generation short reads) or fastplong (third-generation long reads).'),
        'trim_length_required': (
            ['Params_b_sample', '01_trim', 'length_required'], 15, int,
            '01_trim: Reads shorter than length_required will be discarded'),
        'trim_n_base_limit': (
            ['Params_b_sample', '01_trim', 'n_base_limit'], 5, int,
            "01_trim: if one read's number of N base is >n_base_limit, then this read/pair is discarded."),
        'trim_compression': (
            ['Params_b_sample', '01_trim', 'compression'], 4, int,
            '01_trim: Compression level for gzip output [1-9]. 1 is fastest, 9 is smallest.'),
        'trim_extra_args': (
            ['Params_b_sample', '01_trim', 'extra_args'], '', str,
            '01_trim: Extra arguments forwarded to fastp or fastplong for raw reads trimming (e.g., "--reads_to_process 0 -Q").'),
        'trim_stats_sample_threads': (
            ['Params_b_sample', '01_trim_stats', 'threads'], 16, int,
            '01_trim_stats: Threads.'),
        'assembly_method': (
            ['Params_b_sample', '02_assembly', 'method'], 'megahit', str,
            '\b\n02_assembly: Underlying assembly software using one of [megahit/spades/metaflye].\n'
            'Must be compatible with --read-type (config key read_type):\n'
            'megahit: MEGAHIT (second-generation sequencing; for paired or single DNA reads).\n'
            'spades: SPAdes (second-generation sequencing; metaSPAdes for paired DNA reads, rnaSPAdes for paired-rna or single-rna RNA reads).\n'
            'metaflye: metaFlye (third-generation sequencing; required for pacbio-raw/corr/hifi or nano-raw/corr/hq long reads).'),
        'assembly_threads': (
            ['Params_b_sample', '02_assembly', 'threads'], 40, int,
            '02_assembly: Threads for [megahit/spades/metaflye].'),
        'megahit_presets': (
            ['Params_b_sample', '02_megahit', 'presets'], 'meta-sensitive', str,
            '\b\n02_megahit: MEGAHIT presets using one of [meta-sensitive/meta-large].\n'
            "meta-sensitive: '--min-count 1 --k-list 21,29,39,49,...,129,141'.\n"
            "meta-large: '--k-min 27 --k-max 127 --k-step 10'.\n"
            '(large & complex metagenomes, like soil).'),
        'megahit_k_list': (
            ['Params_b_sample', '02_megahit', 'k_list'], '', str,
            '02_megahit: MEGAHIT --k-list (optional). Explicit comma-separated k-mer sizes (e.g. "21,29,39,59,79,99,119,141"). Overrides --megahit-presets. If set, --megahit-k-min, --megahit-k-max, --megahit-k-step must NOT be set. Optional --megahit-min-count can be used together.'),
        'megahit_min_count': (
            ['Params_b_sample', '02_megahit', 'min_count'], 2, int,
            '02_megahit: MEGAHIT --min-count (optional). Minimum multiplicity for filtering (k_min+1)-mers. Can be used with --megahit-k-list or --megahit-k-min/--k-max/--k-step.'),
        'megahit_k_min': (
            ['Params_b_sample', '02_megahit', 'k_min'], '', str,
            '02_megahit: MEGAHIT --k-min (optional). Minimum kmer size (<= 255), must be odd number. If set with --k-max and --k-step, overrides --megahit-presets and uses custom k-mer mode.'),
        'megahit_k_max': (
            ['Params_b_sample', '02_megahit', 'k_max'], '', str,
            '02_megahit: MEGAHIT --k-max (optional). Maximum kmer size (<= 255), must be odd number.'),
        'megahit_k_step': (
            ['Params_b_sample', '02_megahit', 'k_step'], '', str,
            '02_megahit: MEGAHIT --k-step (optional). Increment of kmer size of each iteration (<= 28), must be even number.'),
        'megahit_extra_args': (
            ['Params_b_sample', '02_megahit', 'extra_args'], '', str,
            '02_megahit: Extra arguments forwarded to MEGAHIT (e.g. "--no-mercy --tmp-dir /path/to/your/tmp/dir"). Do NOT use this for k-mer settings (--k-list, --k-min, --presets, etc.).'),
        'spades_extra_args': (
            ['Params_b_sample', '02_spades', 'extra_args'], '', str,
            '02_spades: Extra arguments forwarded to SPAdes. Do NOT use this for input/output or thread settings.'),
        'metaflye_extra_args': (
            ['Params_b_sample', '02_metaflye', 'extra_args'], '', str,
            '02_metaflye: Extra arguments forwarded to Flye for third-generation long reads assembly (e.g., "--iterations 2"). Do NOT use this for read type or input/output settings.'),
        'mapping_threads': (
            ['Params_b_sample', '03_mapping', 'threads'], 20, int,
            '03_mapping: Threads for [minibwa/bwamem/bwamem2/strobealign].'),
        'mapping_method': (
            ['Params_b_sample', '03_mapping', 'method'], 'minibwa', str,
            '\b\n03_mapping: Underlying mapping software used [minibwa/bwamem/bwamem2/strobealign].\n'
            'minibwa: minibwa using default parameters.\n'
            'bwamem: bwa mem using default parameters.\n'
            'bwamem2: bwa-mem2 using default parameters.\n'
            'strobealign: strobealign using default parameters.'),
        'mapping_samtools_sort_mem': (
            ['Params_b_sample', '03_mapping', 'samtools_sort_mem'], '768M', str,
            '03_mapping: Set maximum memory per thread; suffix K/M/G recognized, e.g. 2G, 768M.'),
        'mapping_extra_args': (
            ['Params_b_sample', '03_mapping', 'extra_args'], '', str,
            '03_mapping: Extra arguments forwarded to the active mapping aligner (minibwa map, bwa mem, bwa-mem2 mem, or strobealign --use-index), (e.g. "-B 8" for minibwa map).'),
        'binning_metawrap_depth_threads': (
            ['Params_b_sample', '05_metawrap_binning_depth', 'threads'], 8, int,
            '05_metawrap_binning_depth: Threads.'),
        'binning_min_read_percent_identity': (
            ['Params_b_sample', '05_metawrap_binning_depth', 'min_read_percent_identity'], 97, float,
            '05_metawrap_binning_depth: Exclude reads by overall percent identity [0-100] during the binning process, e.g. 97 for 97%.'),
        'binning_contig_end_exclusion': (
            ['Params_b_sample', '05_metawrap_binning_depth', 'contig_end_exclusion'], 75, int,
            '05_metawrap_binning_depth: Exclude bases at the ends of contig sequences used for depth calculation and bining.'),
        'binning_min_contig_length': (
            ['Params_b_sample', '05_metawrap_binning_depth', 'min_contig_length'], 1500, int,
            '05_metawrap_binning_depth: Minimum contig length for binning.'),
        'binning_tools': (
            ['Params_b_sample', '05_metawrap_binning_depth', 'binning_tools'], 'metabat2,maxbin2,concoct', str,
            '\b\n05_metawrap_binning_depth: Comma-separated list of binning tools, choose 2 or 3 from [metabat2,maxbin2,concoct].\n'
            'At least 2 tools recommended for bin_refinement.'),
        'binning_max_memory_gb': (
            ['Params_b_sample', '05_metawrap_binning_depth', 'max_memory_gb'], 512, int,
            '05_metawrap_binning_depth: Maximum memory (GB) for MetaWRAP binning.'),
        'binning_metawrap_extra_args': (
            ['Params_b_sample', '05_metawrap_binning_depth', 'extra_args'], '', str,
            '05_metawrap_binning_depth: Extra arguments forwarded to metaWRAP binning (e.g., "-m 5").'),
        'metawrap_refinement_max_memory_gb': (
            ['Params_b_sample', '06_metawrap_bin_refinement', 'max_memory_gb'], 512, int,
            '06_metawrap_bin_refinement: Maximum memory in GB for MetaWRAP bin_refinement.'),
        'metawrap_refinement_extra_args': (
            ['Params_b_sample', '06_metawrap_bin_refinement', 'extra_args'], '', str,
            '06_metawrap_bin_refinement: Extra arguments forwarded to metaWRAP bin_refinement (e.g., "--remove-ambiguous").'),
        'bin_checkm_threads': (
            ['Params_b_sample', '06_bin_checkm', 'threads'], 8, int,
            '06_bin_checkm or 06_metawrap_bin_refinement: Threads.'),
        'bin_min_completeness': (
            ['Params_b_sample', '06_bin_checkm', 'min_completeness'], 50, float,
            '06_bin_checkm or 06_metawrap_bin_refinement: Minimum completeness (%) bin filtering, e.g. 50 for 50%.'),
        'bin_max_contamination': (
            ['Params_b_sample', '06_bin_checkm', 'max_contamination'], 10, float,
            '06_bin_checkm or 06_metawrap_bin_refinement: Maximum contamination (%) bin filtering, e.g. 10 for 10%.'),
        'bin_min_quality_score': (
            ['Params_b_sample', '06_bin_checkm', 'min_quality_score'], 50, float,
            '06_bin_checkm: Minimum quality score for bin filtering (GTDB-like), QS = completeness - 5*contamination, e.g. 50 for QS>=50.'),
        'contigs_depth_sample_threads': (
            ['Params_b_sample', '07_contigs_depth', 'threads'], 8, int,
            '07_contigs_depth: Threads.'),
        'depth_min_read_percent_identity': (
            ['Params_b_sample', '07_contigs_depth', 'min_read_percent_identity'], 95, float,
            '\b\n07_contigs_depth: Exclude reads by overall percent identity [0-100] for contig depth estimation\n'
            'using jgi_summarize_bam_contig_depths --percentIdentity, e.g. 95 for 95%.\n'
            'Reads with percent identity below this threshold are excluded from depth calculation.\n'
            'Higher values increase stringency, retaining only high-quality alignments.'),
        'depth_contig_end_exclusion': (
            ['Params_b_sample', '07_contigs_depth', 'contig_end_exclusion'], 75, int,
            '\b\n07_contigs_depth: Exclude bases at the ends of contig sequences for depth estimation using jgi_summarize_bam_contig_depths --maxEdgeBases.\n'
            'Bases within this distance from either end of a contig are ignored when calculating depth,\n'
            'preventing edge effects from partial alignments.'),
        'depth_min_contig_length': (
            ['Params_b_sample', '07_contigs_depth', 'min_contig_length'], 300, int,
            '\b\n07_contigs_depth: Minimum contig length for depth estimation [bp], e.g. 300 for 300 bp.\n'
            'Contigs shorter than this threshold are excluded from depth calculation.'),
        'contigs_prodigal_extra_args': (
            ['Params_b_sample', '07_contigs_depth', 'extra_args'], '', str,
            '07_contigs_depth: Extra arguments forwarded to prodigal for contig gene prediction (e.g., "-m -n").'),
        'bin_prodigal_threads': (
            ['Params_b_sample', '08_bin_prodigal', 'threads'], 8, int,
            '08_bin_prodigal: Threads.'),
        'bin_prodigal_extra_args': (
            ['Params_b_sample', '08_bin_prodigal', 'extra_args'], '', str,
            '08_bin_prodigal: Extra arguments forwarded to prodigal for MAG gene prediction (e.g., "-m -n").'),
        'raw_stats_project_threads': (
            ['Params_b_project', '01_cat_raw_stats', 'threads'], 4, int,
            '01_cat_raw_stats: Threads.'),
        'trimmed_stats_project_threads': (
            ['Params_b_project', '02_cat_trimmed_stats', 'threads'], 4, int,
            '02_cat_trimmed_stats: Threads.'),
        'contigs_stats_project_threads': (
            ['Params_b_project', '03_contigs_stats', 'threads'], 4, int,
            '03_contigs_stats: Threads.'),
        'merge_stats_project_threads': (
            ['Params_b_project', '04_merge_reads_contigs_stats', 'threads'], 4, int,
            '04_merge_reads_contigs_stats: Threads.'),
        'gather_checkm_project_threads': (
            ['Params_b_project', '05_gather_project_bin_checkm_stats', 'threads'], 4, int,
            '05_gather_project_bin_checkm_stats: Threads.'),
        'virus_identify_method': (
            ['Params_v_sample', '01_virus_identify', 'identify_method'], 'genomad_vs_vb', str,
            '\b\n01_virus_identify: Virus identification method, should be one of the [genomad/vs_vb/genomad_vs_vb].\n'
            'genomad: Use geNomad;\n'
            'vs_vb: Use VirSorter2 and CheckV, then use VIBRANT to check the low quality viruses from VirSorter2;\n'
            'genomad_vs_vb: Combine genomad and vs_vb virus to get the union results.'
        ),
        'virus_length_limit': (
            ['Params_v_sample', '01_virus_identify', 'virus_length_limit'], 2000, int,
            '01_virus_identify: Minimum length (bp) for viral contigs.'
        ),
        'virus_identify_threads': (
            ['Params_v_sample', '01_virus_identify', 'threads'], 8, int,
            '01_virus_identify: Threads.'),
        'virus_contig_taxonomy_cat_range': (
            ['Params_v_sample', '03_taxonomy_CAT_pack', 'range'], 10, float,
            '\b\n03_taxonomy_CAT_pack: Percentage range [0-100] of DIAMOND hits relative to top bit-score to include in LCA calculation.\n'
            'Higher values yield more conservative classifications by incorporating divergent hits, while lower values enable specific assignments.'
        ),
        'virus_contig_taxonomy_cat_fraction': (
            ['Params_v_sample', '03_taxonomy_CAT_pack', 'fraction'], 0.50, float,
            '\b\n03_taxonomy_CAT_pack: Minimum bit-score support fraction [0-0.99].\n'
            'Higher values enforce stricter consensus (one classification per sequence), while lower values allow multiple lineage detections.'
        ),
        'virus_contig_taxonomy_cat_pack_threads': (
            ['Params_v_sample', '03_taxonomy_CAT_pack', 'threads'], 8, int,
            '03_taxonomy_CAT_pack: Threads.'),
        'virus_contig_taxonomy_cat_pack_extra_args': (
            ['Params_v_sample', '03_taxonomy_CAT_pack', 'extra_args'], '', str,
            '03_taxonomy_CAT_pack: Extra arguments forwarded to CAT_pack contigs (e.g., "--block_size 12.0").'),
        'virus_mapping_threads': (
            ['Params_v_sample', '03_mapping', 'threads'], 20, int,
            '03_mapping: Threads for [minibwa/bwamem/bwamem2/strobealign].'),
        'virus_mapping_method': (
            ['Params_v_sample', '03_mapping', 'method'], 'minibwa', str,
            '\b\n03_mapping: Underlying mapping software used [minibwa/bwamem/bwamem2/strobealign].\n'
            'minibwa: minibwa using default parameters.\n'
            'bwamem: bwa mem using default parameters.\n'
            'bwamem2: bwa-mem2 using default parameters.\n'
            'strobealign: strobealign using default parameters.'),
        'virus_mapping_samtools_sort_mem': (
            ['Params_v_sample', '03_mapping', 'samtools_sort_mem'], '768M', str,
            '03_mapping: Set maximum memory per thread; suffix K/M/G recognized, e.g. 2G, 768M.'),
        'virus_mapping_extra_args': (
            ['Params_v_sample', '03_mapping', 'extra_args'], '', str,
            '03_mapping: Extra arguments forwarded to the active mapping aligner (minibwa map, bwa mem, bwa-mem2 mem, or strobealign --use-index), (e.g. "-B 8" for minibwa map).'),
        'virus_binning_threads': (
            ['Params_v_sample', '04_binning', 'threads'], 8, int,
            '04_binning_vRhyme: Threads.'),
        'virus_checkv_threads': (
            ['Params_v_sample', '05_checkv', 'threads'], 8, int,
            '05_checkv: Threads.'),
        'virus_vclust_ani_sample': (
            ['Params_v_sample', '06_vclust', 'ani'], 0.95, float,
            '06_vclust: ANI threshold [0-1] of vclust for sample-level virus clustering.'
        ),
        'virus_vclust_threads_sample': (
            ['Params_v_sample', '06_vclust', 'threads'], 8, int,
            '06_vclust: Threads.'),
        'vclust_sample_prefilter_extra_args': (
            ['Params_v_sample', '06_vclust', 'prefilter_extra_args'], '', str,
            '06_vclust: Extra arguments forwarded to vclust prefilter (e.g., "-k 25").'),
        'vclust_sample_align_extra_args': (
            ['Params_v_sample', '06_vclust', 'align_extra_args'], '', str,
            '06_vclust: Extra arguments forwarded to vclust align (e.g., "--mal 11").'),
        'vclust_sample_cluster_extra_args': (
            ['Params_v_sample', '06_vclust', 'cluster_extra_args'], '', str,
            '06_vclust: Extra arguments forwarded to vclust cluster (e.g., "--len_ratio 0").'),
        'virus_iphop_threads_sample': (
            ['Params_v_sample', '07_iPHoP', 'threads'], 20, int,
            '07_iPHoP: Threads.'),
        'virus_vmags_taxonomy_bat_range': (
            ['Params_v_sample', '08_taxonomy_BAT', 'range'], 10, float,
            '\b\n08_taxonomy_BAT: Percentage range [0-100] of DIAMOND hits relative to top bit-score to include in LCA calculation.\n'
            'Higher values yield more conservative classifications by incorporating divergent hits, while lower values enable specific assignments.'
        ),
        'virus_vmags_taxonomy_bat_fraction': (
            ['Params_v_sample', '08_taxonomy_BAT', 'fraction'], 0.50, float,
            '\b\n08_taxonomy_BAT: Minimum bit-score support fraction [0-0.99].\n'
            'Higher values enforce stricter consensus (one classification per sequence), while lower values allow multiple lineage detections.'
        ),
        'virus_vmags_taxonomy_bat_threads': (
            ['Params_v_sample', '08_taxonomy_BAT', 'threads'], 8, int,
            '08_taxonomy_BAT: Threads.'),
        'virus_vmags_taxonomy_bat_extra_args': (
            ['Params_v_sample', '08_taxonomy_BAT', 'extra_args'], '', str,
            '08_taxonomy_BAT: Extra arguments forwarded to CAT_pack bins (e.g., "--block_size 12.0").'),
        'virus_coverm_method_sample': (
            ['Params_v_sample', '09_coverm_eachsample_vMAGs', 'method'], 'trimmed_mean', str,
            '\b\n09_coverm_eachsample_vMAGs: Method(s) for coverm-based coverage calculation, must be one of [trimmed_mean/mean].\n'
            'A more detailed description of the different methods is available at https://github.com/wwood/CoverM#calculation-methods'
        ),
        'virus_coverm_min_read_percent_identity_sample': (
            ['Params_v_sample', '09_coverm_eachsample_vMAGs', 'min_read_percent_identity'], 95, float,
            '09_coverm_eachsample_vMAGs: Exclude reads by overall percent identity [0-100], e.g. 95 for 95%.'
        ),
        'virus_coverm_min_read_aligned_length_sample': (
            ['Params_v_sample', '09_coverm_eachsample_vMAGs', 'min_read_aligned_length'], 45, int,
            '09_coverm_eachsample_vMAGs: Exclude reads with fewer aligned bases.'
        ),
        'virus_coverm_contig_end_exclusion_sample': (
            ['Params_v_sample', '09_coverm_eachsample_vMAGs', 'contig_end_exclusion'], 75, int,
            '09_coverm_eachsample_vMAGs: Exclude bases at the ends of contig sequences for coverm calculation.'
        ),
        'virus_coverm_threads_sample': (
            ['Params_v_sample', '09_coverm_eachsample_vMAGs', 'threads'], 8, int,
            '09_coverm_eachsample_vMAGs: Threads.'),
        'virus_coverm_sample_extra_args': (
            ['Params_v_sample', '09_coverm_eachsample_vMAGs', 'extra_args'], '', str,
            '09_coverm_eachsample_vMAGs: Extra arguments forwarded to coverm genome (e.g., "--quiet").'),
        'virus_vclust_ani_project': (
            ['Params_v_project', '06_vclust', 'ani'], 0.95, float,
            '06_vclust: ANI threshold [0-1] of vclust for project-level virus clustering.'
        ),
        'virus_vclust_threads_project': (
            ['Params_v_project', '06_vclust', 'threads'], 8, int,
            '06_vclust: Threads.'),
        'vclust_project_prefilter_extra_args': (
            ['Params_v_project', '06_vclust', 'prefilter_extra_args'], '', str,
            '06_vclust: Extra arguments forwarded to vclust prefilter (e.g., "-k 25".'),
        'vclust_project_align_extra_args': (
            ['Params_v_project', '06_vclust', 'align_extra_args'], '', str,
            '06_vclust: Extra arguments forwarded to vclust align (e.g., "--mal 11").'),
        'vclust_project_cluster_extra_args': (
            ['Params_v_project', '06_vclust', 'cluster_extra_args'], '', str,
            '06_vclust: Extra arguments forwarded to vclust cluster (e.g., "--len_ratio 0").'),
        'virus_iphop_custom_mags_dir_project': (
            ['Params_v_project', '07_iPHoP', 'custom_MAGs_dir'], 'auto_drep95', str,
            '\b\n07_iPHoP: Custom mMAGs directory for iPHoP host prediction, using one of ["none"/user-defined-path].\n'
            'Use "none" to disable custom mMAGs.\n'
            'Use an absolute path to a custom MAGs directory (including *.fna files).\n'
            'NOTE: "auto_drep99" and "auto_drep95" (automatically using dereplicated mMAGs with metalink default path) are only effective in vm mode.\n'
            'This parameter is not applicable in vc/vcr/vmc/vmcr modes (iPHoP host prediction is skipped).'
        ),
        'virus_iphop_threads_project': (
            ['Params_v_project', '07_iPHoP', 'threads'], 60, int,
            '07_iPHoP: Threads.'),
        'mmags_checkm2_threads': (
            ['Params_m_sample', '01_checkm2', 'threads'], 60, int,
            '01_checkm2: Threads.'),
        'mmags_checkm2_extra_args': (
            ['Params_m_sample', '01_checkm2', 'extra_args'], '', str,
            '01_checkm2: Extra arguments forwarded to checkm2 predict (e.g., "--remove_intermediates").'),
        'mmags_taxonomy_classify_method': (
            ['Params_m_sample', '02_taxonomy', 'MAGs_taxonomy_classify_method'], 'gtdbtk', str,
            '02_taxonomy: mMAGs taxonomy classifier, should be one of the [gtdbtk/BAT].'
        ),
        'mmags_taxonomy_threads': (
            ['Params_m_sample', '02_taxonomy', 'threads'], 60, int,
            '02_taxonomy_gtdbtk: Threads for [gtdbtk/BAT].'),
        'mmags_taxonomy_gtdbtk_extra_args': (
            ['Params_m_sample', '02_taxonomy_gtdbtk', 'extra_args'], '', str,
            '02_taxonomy_gtdbtk: Extra arguments forwarded to gtdbtk classify_wf (e.g., "--min_af 0.5").'),
        'mmags_taxonomy_bat_range': (
            ['Params_m_sample', '02_taxonomy_BAT', 'range'], 10, float,
            '\b\n02_taxonomy_BAT: Percentage range [0-100] of DIAMOND hits relative to top bit-score to include in LCA calculation.\n'
            'Higher values yield more conservative classifications by incorporating divergent hits, while lower values enable specific assignments.\n'
            'Applicable when --mMAGs-taxonomy-classify-method is set to BAT.'
        ),
        'mmags_taxonomy_bat_fraction': (
            ['Params_m_sample', '02_taxonomy_BAT', 'fraction'], 0.50, float,
            '\b\n02_taxonomy_BAT: Minimum bit-score support fraction [0-0.99].\n'
            'Higher values enforce stricter consensus (one classification per sequence), while lower values allow multiple lineage detections.\n'
            'Applicable when --mMAGs-taxonomy-classify-method is set to BAT.'
        ),
        'mmags_taxonomy_bat_extra_args': (
            ['Params_m_sample', '02_taxonomy_BAT', 'extra_args'], '', str,
            '02_taxonomy_BAT: Extra arguments forwarded to CAT_pack bins (e.g., "--block_size 12.0").'),
        'mmags_drep99_completeness_per_sample': (
            ['Params_m_sample', '05_drep_eachsample_99', 'comp'], 50, float,
            '05_drep_eachsample_99: The first level dereplication of mMAGs in a sample using completeness threshold [0-100], e.g., 50 for 50%.'
        ),
        'mmags_drep99_contamination_per_sample': (
            ['Params_m_sample', '05_drep_eachsample_99', 'con'], 10, float,
            '05_drep_eachsample_99: The first level dereplication of mMAGs in a sample using contamination threshold [0-100], e.g., 10 for 10%.'
        ),
        'mmags_drep99_ani_per_sample': (
            ['Params_m_sample', '05_drep_eachsample_99', 'sa'], 0.99, float,
            '05_drep_eachsample_99: The first level dereplication of mMAGs in a sample using ANI threshold [0-1], e.g., 0.99 for 99%.'
        ),
        'mmags_drep99_sample_threads': (
            ['Params_m_sample', '05_drep_eachsample_99', 'threads'], 20, int,
            '05_drep_eachsample_99: Threads.'),
        'mmags_drep99_sample_extra_args': (
            ['Params_m_sample', '05_drep_eachsample_99', 'extra_args'], '', str,
            '05_drep_eachsample: Extra arguments forwarded to dRep dereplicate (e.g., "--completeness_weight 1 --contamination_weight 5").'),
        'mmags_drep95_completeness_per_sample': (
            ['Params_m_sample', '06_drep_eachsample_95', 'comp'], 50, float,
            '06_drep_eachsample_95: The second level dereplication of mMAGs in a sample using completeness threshold [0-100], e.g., 50 for 50%.'
        ),
        'mmags_drep95_contamination_per_sample': (
            ['Params_m_sample', '06_drep_eachsample_95', 'con'], 10, float,
            '06_drep_eachsample_95: The second level dereplication of mMAGs in a sample using contamination threshold [0-100], e.g., 10 for 10%.'
        ),
        'mmags_drep95_ani_per_sample': (
            ['Params_m_sample', '06_drep_eachsample_95', 'sa'], 0.95, float,
            '\b\n06_drep_eachsample_95: The second level dereplication of mMAGs in a sample using ANI threshold [0-1], e.g., 0.95 for 95%.\n'
            'This value should be smaller than --drep99-ANI-per-sample'
        ),
        'mmags_drep95_sample_threads': (
            ['Params_m_sample', '06_drep_eachsample_95', 'threads'], 20, int,
            '06_drep_eachsample_95: Threads.'),
        'mmags_drep95_sample_extra_args': (
            ['Params_m_sample', '06_drep_eachsample_95', 'extra_args'], '', str,
            '06_drep_eachsample: Extra arguments forwarded to dRep dereplicate (e.g., "--completeness_weight 1 --contamination_weight 5").'),
        'mmags_coverm_method': (
            ['Params_m_sample', '07_coverm_eachsample', 'method'], 'trimmed_mean', str,
            '\b\n07_coverm_eachsample: Method(s) for coverm-based coverage calculation, must be one of [trimmed_mean/mean].\n'
            'A more detailed description of the different methods is available at https://github.com/wwood/CoverM#calculation-methods'
        ),
        'mmags_coverm_min_read_percent_identity': (
            ['Params_m_sample', '07_coverm_eachsample', 'min_read_percent_identity'], 95, float,
            '07_coverm_eachsample: Exclude reads by overall percent identity [0-100], e.g. 95 for 95%.'
        ),
        'mmags_coverm_min_read_aligned_length': (
            ['Params_m_sample', '07_coverm_eachsample', 'min_read_aligned_length'], 45, int,
            '07_coverm_eachsample: Exclude reads with fewer aligned bases.'
        ),
        'mmags_coverm_contig_end_exclusion': (
            ['Params_m_sample', '07_coverm_eachsample', 'contig_end_exclusion'], 75, int,
            '07_coverm_eachsample: Exclude bases at the ends of contig sequences for coverm calculation.'
        ),
        'mmags_coverm_threads': (
            ['Params_m_sample', '07_coverm_eachsample', 'threads'], 20, int,
            '07_coverm_eachsample: Threads.'),
        'mmags_coverm_sample_extra_args': (
            ['Params_m_sample', '07_coverm_eachsample', 'extra_args'], '', str,
            '07_coverm_eachsample: Extra arguments forwarded to coverm genome (e.g., "--quiet").'),
        'mmags_function_classify_method': (
            ['Params_m_sample', '08_function', 'MAGs_function_classify_method'], 'emapper', str,
            '08_function: mMAGs function classifier using one of [kraken2/emapper].'
        ),
        'mmags_function_threads': (
            ['Params_m_sample', '08_function', 'threads'], 60, int,
            '08_function: Threads for [kraken2/emapper].'),
        'mmags_function_kraken2_confidence': (
            ['Params_m_sample', '08_function_kraken2', 'confidence'], 0.05, float,
            '08_function_kraken2: Confidence threshold [0-1] for functional classification for mMAGs when --mMAGs-function-classify-method is set to kraken2.'
        ),
        'mmags_function_kraken2_extra_args': (
            ['Params_m_sample', '08_function_kraken2', 'extra_args'], '', str,
            '08_function_kraken2: Extra arguments forwarded to kraken2 (e.g., "--memory-mapping").'),
        'mmags_function_emapper_evalue': (
            ['Params_m_sample', '08_function_emapper', 'seed_ortholog_evalue'], 1e-10, float,
            '08_function_emapper: Report only alignments with e-value below or equal to the threshold, applicable when --mMAGs-function-classify-method is set to emapper.'
        ),
        'mmags_function_emapper_method': (
            ['Params_m_sample', '08_function_emapper', 'method'], 'diamond', str,
            '\b\n08_function_emapper: search seed orthologs using one of [diamond/mmseqs], applicable when --mMAGs-function-classify-method is set to emapper.\n'
            'diamond: using default parameters of diamond mode of emapper.\n'
            'mmseqs: using default parameters of mmseqs mode of emapper'
        ),
        'mmags_function_emapper_query_cover': (
            ['Params_m_sample', '08_function_emapper', 'query_cover'], 30, float,
            '\b\n08_function_emapper: Report only alignments with query cover above or equal to the given percentage [0-100],.\n'
            'applicable when --mMAGs-function-classify-method is set to emapper'
        ),
        'mmags_function_emapper_extra_args': (
            ['Params_m_sample', '08_function_emapper', 'extra_args'], '', str,
            '08_function_emapper: Extra arguments forwarded to eggNOG-mapper (e.g., "--pident 30").'),
        'mmags_drep99_completeness_per_project': (
            ['Params_m_project', '01_drep_project_99', 'comp'], 50, float,
            '01_drep_project_99: Dereplicate mMAGs from all samples in a project with completeness threshold [0-100], e.g., 50 for 50%.'
        ),
        'mmags_drep99_contamination_per_project': (
            ['Params_m_project', '01_drep_project_99', 'con'], 10, float,
            '01_drep_project_99: Dereplicate mMAGs from all samples in a project with contamination threshold [0-100], e.g., 5 for 5%.'
        ),
        'mmags_drep99_project_threads': (
            ['Params_m_project', '01_drep_project_99', 'threads'], 60, int,
            '01_drep_project_99: Threads.'),
        'mmags_drep99_project_extra_args': (
            ['Params_m_project', '01_drep_project_99', 'extra_args'], '', str,
            '01_drep_project_99: Extra arguments forwarded to dRep dereplicate ("--completeness_weight 1 --contamination_weight 5").'),
        'mmags_drep95_completeness_per_project': (
            ['Params_m_project', '01_drep_project_95', 'comp'], 50, float,
            '01_drep_project_95: Dereplicate mMAGs from all samples in a project with completeness threshold [0-100], e.g., 50 for 50%.'
        ),
        'mmags_drep95_contamination_per_project': (
            ['Params_m_project', '01_drep_project_95', 'con'], 10, float,
            '01_drep_project_95: Dereplicate mMAGs from all samples in a project with contamination threshold [0-100], e.g., 5 for 5%.'
        ),
        'mmags_drep95_project_threads': (
            ['Params_m_project', '01_drep_project_95', 'threads'], 60, int,
            '01_drep_project_95: Threads.'),
        'mmags_drep95_project_extra_args': (
            ['Params_m_project', '01_drep_project_95', 'extra_args'], '', str,
            '01_drep_project_95: Extra arguments forwarded to dRep dereplicate (e.g., "--completeness_weight 1 --contamination_weight 5").'),
        'merge_mmags_min_completeness1': (
            ['Params_vm_sample', '05_merge_tax_fun_depth', 'MAGs_min_completeness1'], 50, float,
            '02_quality_gate: Minimum completeness threshold (%) for CheckM1-based mMAG quality filtering [0-100], e.g., 50 for 50%.'
        ),
        'merge_mmags_max_contamination1': (
            ['Params_vm_sample', '05_merge_tax_fun_depth', 'MAGs_max_contamination1'], 1, float,
            '02_quality_gate: Maximum contamination threshold (%) for CheckM1-based mMAG quality filtering [0-100], e.g., 5 for 5%.'
        ),
        'merge_mmags_min_completeness2': (
            ['Params_vm_sample', '05_merge_tax_fun_depth', 'MAGs_min_completeness2'], 50, float,
            '02_quality_gate: Minimum completeness threshold (%) for CheckM2-based mMAG quality filtering [0-100], e.g., 50 for 50%.'
        ),
        'merge_mmags_max_contamination2': (
            ['Params_vm_sample', '05_merge_tax_fun_depth', 'MAGs_max_contamination2'], 1, float,
            '02_quality_gate: Maximum contamination threshold (%) for CheckM2-based mMAG quality filtering [0-100], e.g., 5 for 5%.'
        ),
        'merge_min_orf_length': (
            ['Params_vm_sample', '05_merge_tax_fun_depth', 'minORFLen'], 50, int,
            '05_merge_tax_fun_depth: Minimum ORF length for finally used ORFs of vMAGs and mMAGs.'
        ),
        'merge_vmags_checkv_quality': (
            ['Params_vm_sample', '05_merge_tax_fun_depth', 'checkv_quality'], 'High-quality', str,
            '\b\n05_merge_tax_fun_depth: CheckV quality filter for viral MAGs (lower bound filter).\n'
            'viral MAGs with quality equal to or higher than the specified threshold will be retained.\n'
            'Options: Complete, High-quality, Medium-quality, Low-quality, Not-determined.\n'
            'Complete: complete genome (highest quality)\n'
            'High-quality: genome completeness ≥90%\n'
            'Medium-quality: genome completeness 50-90%\n'
            'Low-quality: genome completeness <50%\n'
            'Not-determined: completeness could not be determined (lowest quality)'
        ),
        'merge_multiple_ko_handling_mode': (
            ['Params_vm_sample', '05_merge_tax_fun_depth', 'KO_mode'], 'ignore_multiple_KOs', str,
            '\b\n05_merge_tax_fun_depth: Method using one of [split_multiple_KOs/keep_multiple_KOs/ignore_multiple_KOs] when an ORF has multiple KO annotations:.\n'
            'split_multiple_KOs: evenly distributing the total abundance across all assigned KOs;\n'
            'keep_multiple_KOs: retain original multiple KO annotations;\n'
            'ignore_multiple_KOs: ignore the ORFs with multiple KO annotations and only use the ORFs with single KO annotation.'
        ),
        'merge_min_ko_num_per_taxon': (
            ['Params_vm_sample', '05_merge_tax_fun_depth', 'min_KO_num_in_one_taxon'], 100, int,
            '\b\n05_merge_tax_fun_depth: Minimum number of KOs per taxon, used to filter outputs in KO_Taxa_RelativeAbundance and KO_Taxa_SumCoverage.\n'
            'This threshold filters out taxa with insufficient KOs,\n'
            'retaining only KO-rich taxa to ensure robust metabolic reconstruction\n'
            'while preventing excessive output file sizes caused by sparse matrices.'
        ),
        'merge_tax_fun_depth_threads': (
            ['Params_vm_sample', '05_merge_tax_fun_depth', 'threads'], 20, int,
            '05_merge_tax_fun_depth: Threads.'),
        'merge_cami_format_output': (
            ['Params_vm_sample', '05_merge_tax_fun_depth', 'CAMI_format_output'], False, bool,
            '05_merge_tax_fun_depth: Report a taxonomic profile using the CAMI output format (flag only). Disabled when mMAG mode uses gtdbtk.'
        ),
        'merge_output_detail': (
            ['Params_vm_sample', '05_merge_tax_fun_depth', 'output_detail'], False, bool,
            '05_merge_tax_fun_depth: Whether to output detail tables (flag only).'),
        'mclink_sample_min_ko_num_per_taxon': (
            ['Params_vm_sample', '06_mclink_taxon', 'min_KO_num_in_one_taxon'], 100, int,
            '06_mclink_taxon: Minimum number of KOs per taxon for calculating taxon-level metabolic pathway/module completeness or abundance per sample. Taxa with fewer than this number of detected KOs in a given sample will be filtered out. Calculated independently for each sample.'
        ),
        'mclink_sample_min_ko_abundance_per_taxon': (
            ['Params_vm_sample', '06_mclink_taxon', 'min_mean_KO_RA_in_one_taxon'], 0, float,
            '06_mclink_taxon: Minimum mean KO abundance [0-1e6] per taxon for calculating taxon-level metabolic pathway/module completeness or abundance per sample. Taxa with mean KO abundance below this threshold in a given sample will be filtered out. Note: GPM/TPM sum of all taxa in each sample is 1e6.'
        ),
        'mclink_sample_topn_taxon': (
            ['Params_vm_sample', '06_mclink_taxon', 'TopN'], 100, int,
            '06_mclink_taxon: Number of top taxa (ranked by mean KO abundance) retained for calculating taxon-level metabolic pathway/module completeness or abundance per sample. Set to 0 to disable filtering.'
        ),
        'mclink_sample_completeness_per_taxon': (
            ['Params_vm_sample', '06_mclink_taxon', 'completeness_threshold'], 0.5, float,
            '06_mclink_taxon: Completeness threshold [0-1] of modules for calculating taxon-level metabolic pathway/module abundance. Higher values retain only more complete profiles.'),
        'mclink_sample_plus_scale_method_per_taxon': (
            ['Params_vm_sample', '06_mclink_taxon', 'plus_scale_method'], 'mean', str,
            "06_mclink_taxon: Aggregation method for '+'-separated KOs in a module [mean/max/min]. 'mean': average abundance across all required KOs; 'max': highest KO abundance; 'min': lowest KO abundance."),
        'mclink_sample_comma_scale_method_per_taxon': (
            ['Params_vm_sample', '06_mclink_taxon', 'comma_scale_method'], 'sum', str,
            "06_mclink_taxon: Aggregation method for ','-separated alternative KOs in a module [sum/max]. 'sum': cumulative abundance of all alternative KOs; 'max': abundance of the highest single KO among alternatives. Note: this applies to abundance calculation only; module completeness is always computed with 'max'."),
        'barplot_topn_taxon': (
            ['Params_vm_project', '01_Taxa_visualize', 'TopN'], 10, int,
            '01_Taxa_visualize: Number of top taxa to display in stacked barplots. Remaining categories are merged into Others. Maximum 29.'
        ),
        'barplot_topn_fun': (
            ['Params_vm_project', '03_Fun_visualize', 'TopN'], 10, int,
            '03_Fun_visualize: Number of top functional categories to display in stacked barplots. Remaining categories are merged into Others. Maximum 29.'
        ),
        'mclink_min_ko_num_per_sample': (
            ['Params_vm_project', '03_mclink_sample', 'min_sample_num_per_KO'], 1, int,
            '03_mclink_sample: Minimum number of samples in which a KO must be detected for calculating community-level metabolic pathway/module completeness or abundance.'
        ),
        'mclink_min_ko_prevalence_per_sample': (
            ['Params_vm_project', '03_mclink_sample', 'preva'], 0, float,
            '03_mclink_sample: Minimum detection rate [0-1] of a KO across all samples for calculating community-level metabolic pathway/module completeness or abundance.'
        ),
        'mclink_min_ko_abundance_per_sample': (
            ['Params_vm_project', '03_mclink_sample', 'abund'], 0, float,
            '03_mclink_sample: Minimum mean relative abundance [0-1e6] of a KO across all samples for calculating community-level metabolic pathway/module completeness or abundance. KOs with mean abundance below this threshold are excluded. Note: GPM/TPM sum of a sample is 1e6.'
        ),
        'mclink_completeness_per_sample': (
            ['Params_vm_project', '03_mclink_sample', 'completeness_threshold'], 0.5, float,
            '03_mclink_sample: Completeness threshold [0-1] of modules for calculating community-level metabolic pathway/module abundance. Higher values retain only more complete profiles.'
        ),
        'mclink_plus_scale_method_per_sample': (
            ['Params_vm_project', '03_mclink_sample', 'plus_scale_method'], 'mean', str,
            "03_mclink_sample: Aggregation method for '+'-separated KOs in a module [mean/max/min]. 'mean': average abundance across all required KOs; 'max': highest KO abundance; 'min': lowest KO abundance."
        ),
        'mclink_comma_scale_method_per_sample': (
            ['Params_vm_project', '03_mclink_sample', 'comma_scale_method'], 'sum', str,
            "03_mclink_sample: Aggregation method for ','-separated alternative KOs in a module [sum/max]. 'sum': cumulative abundance of all alternative KOs; 'max': abundance of the highest single KO among alternatives. Note: this applies to abundance calculation only; module completeness is always computed with 'max'."
        ),
        'merge_min_taxa_ko_prevalence_across_samples': (
            ['Params_vm_project', '05_merge_tax_KO_depth', 'preva'], 0.10, float,
            '05_merge_tax_KO_depth: Minimum prevalence [0-1] of a Taxon-KO across samples. Taxa-KOs detected in fewer samples than this proportion are excluded.'
        ),
        'merge_min_taxa_ko_abundance_across_samples': (
            ['Params_vm_project', '05_merge_tax_KO_depth', 'abund'], 0.0, float,
            '05_merge_tax_KO_depth: Minimum mean relative abundance [0-1e6] of a Taxon-KO across samples. Taxa-KOs with mean abundance below this threshold are excluded. Note: GPM/TPM sum of a sample is 1e6.'
        ),
        'merge_min_taxa_module_prevalence_across_samples': (
            ['Params_vm_project', '05_merge_tax_Module_RelativeAbundance', 'preva'], 0.10, float,
            '05_merge_tax_Module_RelativeAbundance: Minimum prevalence [0-1] of a Taxon-Module across samples. Taxa-Modules detected in fewer samples than this proportion are excluded.'
        ),
        'merge_min_taxa_module_abundance_across_samples': (
            ['Params_vm_project', '05_merge_tax_Module_RelativeAbundance', 'abund'], 0.0, float,
            '05_merge_tax_Module_RelativeAbundance: Minimum mean relative abundance [0-1e6] of a Taxon-Module across samples. Taxa-Modules with mean abundance below this threshold are excluded. Note: GPM/TPM sum of a sample is 1e6.'
        ),
        'merge_min_taxa_module_sc_prevalence_across_samples': (
            ['Params_vm_project', '05_merge_tax_Module_SumCoverage', 'preva'], 0.10, float,
            '05_merge_tax_Module_SumCoverage: Minimum prevalence [0-1] of a Taxon-Module across samples. Taxa-Modules detected in fewer samples than this proportion are excluded.'
        ),
        'merge_min_taxa_module_sc_abundance_across_samples': (
            ['Params_vm_project', '05_merge_tax_Module_SumCoverage', 'abund'], 0.0, float,
            '05_merge_tax_Module_SumCoverage: Minimum mean sum coverage [0-1e6] of a Taxon-Module across samples. Taxa-Modules with mean sum coverage below this threshold are excluded. Note: GPM/TPM sum of a sample is 1e6.'
        ),
        'mclink_project_min_ko_num_per_taxon': (
            ['Params_vm_project', '06_mclink_taxon', 'min_KO_num_in_one_taxon'], 100, int,
            '06_mclink_taxon: Minimum number of KOs per taxon for calculating taxon-level metabolic pathway/module completeness on the project-averaged taxon-function matrix. Taxa with fewer than this number of detected KOs in the project-averaged matrix will be filtered out. KO abundances are first averaged across all samples to generate a single project-level matrix before calculation.'
        ),
        'mclink_project_min_ko_abundance_per_taxon': (
            ['Params_vm_project', '06_mclink_taxon', 'min_mean_KO_RA_in_one_taxon'], 0, float,
            '06_mclink_taxon: Minimum mean KO abundance [0-1e6] per taxon for calculating taxon-level metabolic pathway/module completeness on the project-averaged taxon-function matrix. Taxa with mean KO abundance below this threshold in the project-averaged matrix will be filtered out. Note: GPM/TPM sum of all taxa is 1e6 per sample; values are averaged across samples prior to filtering and calculation.'
        ),
        'mclink_project_topn_taxon': (
            ['Params_vm_project', '06_mclink_taxon', 'TopN'], 100, int,
            '06_mclink_taxon: Number of top taxa (ranked by mean KO abundance in the project-averaged taxon-function matrix) retained for calculating taxon-level metabolic pathway/module completeness. Set to 0 to disable filtering.'
        ),
        'mclink_project_completeness_per_taxon': (
            ['Params_vm_project', '06_mclink_taxon', 'completeness_threshold'], 0.5, float,
            '06_mclink_taxon: Completeness threshold [0-1] of modules for calculating taxon-level metabolic pathway/module abundance on the project-averaged taxon-function matrix. Higher values retain only more complete profiles.'),
        'mclink_project_plus_scale_method_per_taxon': (
            ['Params_vm_project', '06_mclink_taxon', 'plus_scale_method'], 'mean', str,
            "06_mclink_taxon: Aggregation method for '+'-separated KOs in a module [mean/max/min] on the project-averaged taxon-function matrix. 'mean': average abundance across all required KOs; 'max': highest KO abundance; 'min': lowest KO abundance."),
        'mclink_project_comma_scale_method_per_taxon': (
            ['Params_vm_project', '06_mclink_taxon', 'comma_scale_method'], 'sum', str,
            "06_mclink_taxon: Aggregation method for ','-separated alternative KOs in a module [sum/max] on the project-averaged taxon-function matrix. 'sum': cumulative abundance of all alternative KOs; 'max': abundance of the highest single KO among alternatives. Note: this applies to abundance calculation only; module completeness is always computed with 'max'."),
    },
    'vmc': {
        'raw_stats_sample_threads': (
            ['Params_b_sample', '00_raw_stats', 'threads'], 4, int,
            '00_raw_stats: Threads for seqkit.'),
        'trim_sample_threads': (
            ['Params_b_sample', '01_trim', 'threads'], 16, int,
            '01_trim: Threads for fastp (next-generation short reads) or fastplong (third-generation long reads).'),
        'trim_length_required': (
            ['Params_b_sample', '01_trim', 'length_required'], 15, int,
            '01_trim: Reads shorter than length_required will be discarded'),
        'trim_n_base_limit': (
            ['Params_b_sample', '01_trim', 'n_base_limit'], 5, int,
            "01_trim: if one read's number of N base is >n_base_limit, then this read/pair is discarded."),
        'trim_compression': (
            ['Params_b_sample', '01_trim', 'compression'], 4, int,
            '01_trim: Compression level for gzip output [1-9]. 1 is fastest, 9 is smallest.'),
        'trim_extra_args': (
            ['Params_b_sample', '01_trim', 'extra_args'], '', str,
            '01_trim: Extra arguments forwarded to fastp or fastplong for raw reads trimming (e.g., "--reads_to_process 0 -Q").'),
        'trim_stats_sample_threads': (
            ['Params_b_sample', '01_trim_stats', 'threads'], 16, int,
            '01_trim_stats: Threads.'),
        'assembly_method': (
            ['Params_b_sample', '02_assembly', 'method'], 'megahit', str,
            '\b\n02_assembly: Underlying assembly software using one of [megahit/spades/metaflye].\n'
            'Must be compatible with --read-type (config key read_type):\n'
            'megahit: MEGAHIT (second-generation sequencing; for paired or single DNA reads).\n'
            'spades: SPAdes (second-generation sequencing; metaSPAdes for paired DNA reads, rnaSPAdes for paired-rna or single-rna RNA reads).\n'
            'metaflye: metaFlye (third-generation sequencing; required for pacbio-raw/corr/hifi or nano-raw/corr/hq long reads).'),
        'assembly_threads': (
            ['Params_b_sample', '02_assembly', 'threads'], 40, int,
            '02_assembly: Threads for [megahit/spades/metaflye].'),
        'megahit_presets': (
            ['Params_b_sample', '02_megahit', 'presets'], 'meta-sensitive', str,
            '\b\n02_megahit: MEGAHIT presets using one of [meta-sensitive/meta-large].\n'
            "meta-sensitive: '--min-count 1 --k-list 21,29,39,49,...,129,141'.\n"
            "meta-large: '--k-min 27 --k-max 127 --k-step 10'.\n"
            '(large & complex metagenomes, like soil).'),
        'megahit_k_list': (
            ['Params_b_sample', '02_megahit', 'k_list'], '', str,
            '02_megahit: MEGAHIT --k-list (optional). Explicit comma-separated k-mer sizes (e.g. "21,29,39,59,79,99,119,141"). Overrides --megahit-presets. If set, --megahit-k-min, --megahit-k-max, --megahit-k-step must NOT be set. Optional --megahit-min-count can be used together.'),
        'megahit_min_count': (
            ['Params_b_sample', '02_megahit', 'min_count'], 2, int,
            '02_megahit: MEGAHIT --min-count (optional). Minimum multiplicity for filtering (k_min+1)-mers. Can be used with --megahit-k-list or --megahit-k-min/--k-max/--k-step.'),
        'megahit_k_min': (
            ['Params_b_sample', '02_megahit', 'k_min'], '', str,
            '02_megahit: MEGAHIT --k-min (optional). Minimum kmer size (<= 255), must be odd number. If set with --k-max and --k-step, overrides --megahit-presets and uses custom k-mer mode.'),
        'megahit_k_max': (
            ['Params_b_sample', '02_megahit', 'k_max'], '', str,
            '02_megahit: MEGAHIT --k-max (optional). Maximum kmer size (<= 255), must be odd number.'),
        'megahit_k_step': (
            ['Params_b_sample', '02_megahit', 'k_step'], '', str,
            '02_megahit: MEGAHIT --k-step (optional). Increment of kmer size of each iteration (<= 28), must be even number.'),
        'megahit_extra_args': (
            ['Params_b_sample', '02_megahit', 'extra_args'], '', str,
            '02_megahit: Extra arguments forwarded to MEGAHIT (e.g. "--no-mercy --tmp-dir /path/to/your/tmp/dir"). Do NOT use this for k-mer settings (--k-list, --k-min, --presets, etc.).'),
        'spades_extra_args': (
            ['Params_b_sample', '02_spades', 'extra_args'], '', str,
            '02_spades: Extra arguments forwarded to SPAdes. Do NOT use this for input/output or thread settings.'),
        'metaflye_extra_args': (
            ['Params_b_sample', '02_metaflye', 'extra_args'], '', str,
            '02_metaflye: Extra arguments forwarded to Flye for third-generation long reads assembly (e.g., "--iterations 2"). Do NOT use this for read type or input/output settings.'),
        'mapping_threads': (
            ['Params_b_sample', '03_mapping', 'threads'], 20, int,
            '03_mapping: Threads for [minibwa/bwamem/bwamem2/strobealign].'),
        'mapping_method': (
            ['Params_b_sample', '03_mapping', 'method'], 'minibwa', str,
            '\b\n03_mapping: Underlying mapping software used [minibwa/bwamem/bwamem2/strobealign].\n'
            'minibwa: minibwa using default parameters.\n'
            'bwamem: bwa mem using default parameters.\n'
            'bwamem2: bwa-mem2 using default parameters.\n'
            'strobealign: strobealign using default parameters.'),
        'mapping_samtools_sort_mem': (
            ['Params_b_sample', '03_mapping', 'samtools_sort_mem'], '768M', str,
            '03_mapping: Set maximum memory per thread; suffix K/M/G recognized, e.g. 2G, 768M.'),
        'mapping_extra_args': (
            ['Params_b_sample', '03_mapping', 'extra_args'], '', str,
            '03_mapping: Extra arguments forwarded to the active mapping aligner (minibwa map, bwa mem, bwa-mem2 mem, or strobealign --use-index), (e.g. "-B 8" for minibwa map).'),
        'binning_metawrap_depth_threads': (
            ['Params_b_sample', '05_metawrap_binning_depth', 'threads'], 8, int,
            '05_metawrap_binning_depth: Threads.'),
        'binning_min_read_percent_identity': (
            ['Params_b_sample', '05_metawrap_binning_depth', 'min_read_percent_identity'], 97, float,
            '05_metawrap_binning_depth: Exclude reads by overall percent identity [0-100] during the binning process, e.g. 97 for 97%.'),
        'binning_contig_end_exclusion': (
            ['Params_b_sample', '05_metawrap_binning_depth', 'contig_end_exclusion'], 75, int,
            '05_metawrap_binning_depth: Exclude bases at the ends of contig sequences used for depth calculation and bining.'),
        'binning_min_contig_length': (
            ['Params_b_sample', '05_metawrap_binning_depth', 'min_contig_length'], 1500, int,
            '05_metawrap_binning_depth: Minimum contig length for binning.'),
        'binning_tools': (
            ['Params_b_sample', '05_metawrap_binning_depth', 'binning_tools'], 'metabat2,maxbin2,concoct', str,
            '\b\n05_metawrap_binning_depth: Comma-separated list of binning tools, choose 2 or 3 from [metabat2,maxbin2,concoct].\n'
            'At least 2 tools recommended for bin_refinement.'),
        'binning_max_memory_gb': (
            ['Params_b_sample', '05_metawrap_binning_depth', 'max_memory_gb'], 512, int,
            '05_metawrap_binning_depth: Maximum memory (GB) for MetaWRAP binning.'),
        'binning_metawrap_extra_args': (
            ['Params_b_sample', '05_metawrap_binning_depth', 'extra_args'], '', str,
            '05_metawrap_binning_depth: Extra arguments forwarded to metaWRAP binning (e.g., "-m 5").'),
        'metawrap_refinement_max_memory_gb': (
            ['Params_b_sample', '06_metawrap_bin_refinement', 'max_memory_gb'], 512, int,
            '06_metawrap_bin_refinement: Maximum memory in GB for MetaWRAP bin_refinement.'),
        'metawrap_refinement_extra_args': (
            ['Params_b_sample', '06_metawrap_bin_refinement', 'extra_args'], '', str,
            '06_metawrap_bin_refinement: Extra arguments forwarded to metaWRAP bin_refinement (e.g., "--remove-ambiguous").'),
        'bin_checkm_threads': (
            ['Params_b_sample', '06_bin_checkm', 'threads'], 8, int,
            '06_bin_checkm or 06_metawrap_bin_refinement: Threads.'),
        'bin_min_completeness': (
            ['Params_b_sample', '06_bin_checkm', 'min_completeness'], 50, float,
            '06_bin_checkm or 06_metawrap_bin_refinement: Minimum completeness (%) bin filtering, e.g. 50 for 50%.'),
        'bin_max_contamination': (
            ['Params_b_sample', '06_bin_checkm', 'max_contamination'], 10, float,
            '06_bin_checkm or 06_metawrap_bin_refinement: Maximum contamination (%) bin filtering, e.g. 10 for 10%.'),
        'bin_min_quality_score': (
            ['Params_b_sample', '06_bin_checkm', 'min_quality_score'], 50, float,
            '06_bin_checkm: Minimum quality score for bin filtering (GTDB-like), QS = completeness - 5*contamination, e.g. 50 for QS>=50.'),
        'contigs_depth_sample_threads': (
            ['Params_b_sample', '07_contigs_depth', 'threads'], 8, int,
            '07_contigs_depth: Threads.'),
        'depth_min_read_percent_identity': (
            ['Params_b_sample', '07_contigs_depth', 'min_read_percent_identity'], 95, float,
            '\b\n07_contigs_depth: Exclude reads by overall percent identity [0-100] for contig depth estimation\n'
            'using jgi_summarize_bam_contig_depths --percentIdentity, e.g. 95 for 95%.\n'
            'Reads with percent identity below this threshold are excluded from depth calculation.\n'
            'Higher values increase stringency, retaining only high-quality alignments.'),
        'depth_contig_end_exclusion': (
            ['Params_b_sample', '07_contigs_depth', 'contig_end_exclusion'], 75, int,
            '\b\n07_contigs_depth: Exclude bases at the ends of contig sequences for depth estimation using jgi_summarize_bam_contig_depths --maxEdgeBases.\n'
            'Bases within this distance from either end of a contig are ignored when calculating depth,\n'
            'preventing edge effects from partial alignments.'),
        'depth_min_contig_length': (
            ['Params_b_sample', '07_contigs_depth', 'min_contig_length'], 300, int,
            '\b\n07_contigs_depth: Minimum contig length for depth estimation [bp], e.g. 300 for 300 bp.\n'
            'Contigs shorter than this threshold are excluded from depth calculation.'),
        'contigs_prodigal_extra_args': (
            ['Params_b_sample', '07_contigs_depth', 'extra_args'], '', str,
            '07_contigs_depth: Extra arguments forwarded to prodigal for contig gene prediction (e.g., "-m -n").'),
        'bin_prodigal_threads': (
            ['Params_b_sample', '08_bin_prodigal', 'threads'], 8, int,
            '08_bin_prodigal: Threads.'),
        'bin_prodigal_extra_args': (
            ['Params_b_sample', '08_bin_prodigal', 'extra_args'], '', str,
            '08_bin_prodigal: Extra arguments forwarded to prodigal for MAG gene prediction (e.g., "-m -n").'),
        'raw_stats_project_threads': (
            ['Params_b_project', '01_cat_raw_stats', 'threads'], 4, int,
            '01_cat_raw_stats: Threads.'),
        'trimmed_stats_project_threads': (
            ['Params_b_project', '02_cat_trimmed_stats', 'threads'], 4, int,
            '02_cat_trimmed_stats: Threads.'),
        'contigs_stats_project_threads': (
            ['Params_b_project', '03_contigs_stats', 'threads'], 4, int,
            '03_contigs_stats: Threads.'),
        'merge_stats_project_threads': (
            ['Params_b_project', '04_merge_reads_contigs_stats', 'threads'], 4, int,
            '04_merge_reads_contigs_stats: Threads.'),
        'gather_checkm_project_threads': (
            ['Params_b_project', '05_gather_project_bin_checkm_stats', 'threads'], 4, int,
            '05_gather_project_bin_checkm_stats: Threads.'),
        'virus_identify_method': (
            ['Params_v_sample', '01_virus_identify', 'identify_method'], 'genomad_vs_vb', str,
            '\b\n01_virus_identify: Virus identification method, should be one of the [genomad/vs_vb/genomad_vs_vb].\n'
            'genomad: Use geNomad;\n'
            'vs_vb: Use VirSorter2 and CheckV, then use VIBRANT to check the low quality viruses from VirSorter2;\n'
            'genomad_vs_vb: Combine genomad and vs_vb virus to get the union results.'
        ),
        'virus_length_limit': (
            ['Params_v_sample', '01_virus_identify', 'virus_length_limit'], 2000, int,
            '01_virus_identify: Minimum length (bp) for viral contigs.'
        ),
        'virus_identify_threads': (
            ['Params_v_sample', '01_virus_identify', 'threads'], 8, int,
            '01_virus_identify: Threads.'),
        'virus_contig_taxonomy_cat_range': (
            ['Params_v_sample', '03_taxonomy_CAT_pack', 'range'], 10, float,
            '\b\n03_taxonomy_CAT_pack: Percentage range [0-100] of DIAMOND hits relative to top bit-score to include in LCA calculation.\n'
            'Higher values yield more conservative classifications by incorporating divergent hits, while lower values enable specific assignments.'
        ),
        'virus_contig_taxonomy_cat_fraction': (
            ['Params_v_sample', '03_taxonomy_CAT_pack', 'fraction'], 0.50, float,
            '\b\n03_taxonomy_CAT_pack: Minimum bit-score support fraction [0-0.99].\n'
            'Higher values enforce stricter consensus (one classification per sequence), while lower values allow multiple lineage detections.'
        ),
        'virus_contig_taxonomy_cat_pack_threads': (
            ['Params_v_sample', '03_taxonomy_CAT_pack', 'threads'], 8, int,
            '03_taxonomy_CAT_pack: Threads.'),
        'virus_contig_taxonomy_cat_pack_extra_args': (
            ['Params_v_sample', '03_taxonomy_CAT_pack', 'extra_args'], '', str,
            '03_taxonomy_CAT_pack: Extra arguments forwarded to CAT_pack contigs (e.g., "--block_size 12.0").'),
        'virus_mapping_threads': (
            ['Params_v_sample', '03_mapping', 'threads'], 20, int,
            '03_mapping: Threads for [minibwa/bwamem/bwamem2/strobealign].'),
        'virus_mapping_method': (
            ['Params_v_sample', '03_mapping', 'method'], 'minibwa', str,
            '\b\n03_mapping: Underlying mapping software used [minibwa/bwamem/bwamem2/strobealign].\n'
            'minibwa: minibwa using default parameters.\n'
            'bwamem: bwa mem using default parameters.\n'
            'bwamem2: bwa-mem2 using default parameters.\n'
            'strobealign: strobealign using default parameters.'),
        'virus_mapping_samtools_sort_mem': (
            ['Params_v_sample', '03_mapping', 'samtools_sort_mem'], '768M', str,
            '03_mapping: Set maximum memory per thread; suffix K/M/G recognized, e.g. 2G, 768M.'),
        'virus_mapping_extra_args': (
            ['Params_v_sample', '03_mapping', 'extra_args'], '', str,
            '03_mapping: Extra arguments forwarded to the active mapping aligner (minibwa map, bwa mem, bwa-mem2 mem, or strobealign --use-index), (e.g. "-B 8" for minibwa map).'),
        'virus_binning_threads': (
            ['Params_v_sample', '04_binning', 'threads'], 8, int,
            '04_binning_vRhyme: Threads.'),
        'virus_checkv_threads': (
            ['Params_v_sample', '05_checkv', 'threads'], 8, int,
            '05_checkv: Threads.'),
        'virus_vclust_ani_sample': (
            ['Params_v_sample', '06_vclust', 'ani'], 0.95, float,
            '06_vclust: ANI threshold [0-1] of vclust for sample-level virus clustering.'
        ),
        'virus_vclust_threads_sample': (
            ['Params_v_sample', '06_vclust', 'threads'], 8, int,
            '06_vclust: Threads.'),
        'vclust_sample_prefilter_extra_args': (
            ['Params_v_sample', '06_vclust', 'prefilter_extra_args'], '', str,
            '06_vclust: Extra arguments forwarded to vclust prefilter (e.g., "-k 25").'),
        'vclust_sample_align_extra_args': (
            ['Params_v_sample', '06_vclust', 'align_extra_args'], '', str,
            '06_vclust: Extra arguments forwarded to vclust align (e.g., "--mal 11").'),
        'vclust_sample_cluster_extra_args': (
            ['Params_v_sample', '06_vclust', 'cluster_extra_args'], '', str,
            '06_vclust: Extra arguments forwarded to vclust cluster (e.g., "--len_ratio 0").'),
        'virus_iphop_threads_sample': (
            ['Params_v_sample', '07_iPHoP', 'threads'], 20, int,
            '07_iPHoP: Threads.'),
        'virus_vmags_taxonomy_bat_range': (
            ['Params_v_sample', '08_taxonomy_BAT', 'range'], 10, float,
            '\b\n08_taxonomy_BAT: Percentage range [0-100] of DIAMOND hits relative to top bit-score to include in LCA calculation.\n'
            'Higher values yield more conservative classifications by incorporating divergent hits, while lower values enable specific assignments.'
        ),
        'virus_vmags_taxonomy_bat_fraction': (
            ['Params_v_sample', '08_taxonomy_BAT', 'fraction'], 0.50, float,
            '\b\n08_taxonomy_BAT: Minimum bit-score support fraction [0-0.99].\n'
            'Higher values enforce stricter consensus (one classification per sequence), while lower values allow multiple lineage detections.'
        ),
        'virus_vmags_taxonomy_bat_threads': (
            ['Params_v_sample', '08_taxonomy_BAT', 'threads'], 8, int,
            '08_taxonomy_BAT: Threads.'),
        'virus_vmags_taxonomy_bat_extra_args': (
            ['Params_v_sample', '08_taxonomy_BAT', 'extra_args'], '', str,
            '08_taxonomy_BAT: Extra arguments forwarded to CAT_pack bins (e.g., "--block_size 12.0").'),
        'virus_coverm_method_sample': (
            ['Params_v_sample', '09_coverm_eachsample_vMAGs', 'method'], 'trimmed_mean', str,
            '\b\n09_coverm_eachsample_vMAGs: Method(s) for coverm-based coverage calculation, must be one of [trimmed_mean/mean].\n'
            'A more detailed description of the different methods is available at https://github.com/wwood/CoverM#calculation-methods'
        ),
        'virus_coverm_min_read_percent_identity_sample': (
            ['Params_v_sample', '09_coverm_eachsample_vMAGs', 'min_read_percent_identity'], 95, float,
            '09_coverm_eachsample_vMAGs: Exclude reads by overall percent identity [0-100], e.g. 95 for 95%.'
        ),
        'virus_coverm_min_read_aligned_length_sample': (
            ['Params_v_sample', '09_coverm_eachsample_vMAGs', 'min_read_aligned_length'], 45, int,
            '09_coverm_eachsample_vMAGs: Exclude reads with fewer aligned bases.'
        ),
        'virus_coverm_contig_end_exclusion_sample': (
            ['Params_v_sample', '09_coverm_eachsample_vMAGs', 'contig_end_exclusion'], 75, int,
            '09_coverm_eachsample_vMAGs: Exclude bases at the ends of contig sequences for coverm calculation.'
        ),
        'virus_coverm_threads_sample': (
            ['Params_v_sample', '09_coverm_eachsample_vMAGs', 'threads'], 8, int,
            '09_coverm_eachsample_vMAGs: Threads.'),
        'virus_coverm_sample_extra_args': (
            ['Params_v_sample', '09_coverm_eachsample_vMAGs', 'extra_args'], '', str,
            '09_coverm_eachsample_vMAGs: Extra arguments forwarded to coverm genome (e.g., "--quiet").'),
        'mmags_checkm2_threads': (
            ['Params_m_sample', '01_checkm2', 'threads'], 60, int,
            '01_checkm2: Threads.'),
        'mmags_checkm2_extra_args': (
            ['Params_m_sample', '01_checkm2', 'extra_args'], '', str,
            '01_checkm2: Extra arguments forwarded to checkm2 predict (e.g., "--remove_intermediates").'),
        'mmags_taxonomy_classify_method': (
            ['Params_m_sample', '02_taxonomy', 'MAGs_taxonomy_classify_method'], 'gtdbtk', str,
            '02_taxonomy: mMAGs taxonomy classifier, should be one of the [gtdbtk/BAT].'
        ),
        'mmags_taxonomy_threads': (
            ['Params_m_sample', '02_taxonomy', 'threads'], 60, int,
            '02_taxonomy_gtdbtk: Threads for [gtdbtk/BAT].'),
        'mmags_taxonomy_gtdbtk_extra_args': (
            ['Params_m_sample', '02_taxonomy_gtdbtk', 'extra_args'], '', str,
            '02_taxonomy_gtdbtk: Extra arguments forwarded to gtdbtk classify_wf (e.g., "--min_af 0.5").'),
        'mmags_taxonomy_bat_range': (
            ['Params_m_sample', '02_taxonomy_BAT', 'range'], 10, float,
            '\b\n02_taxonomy_BAT: Percentage range [0-100] of DIAMOND hits relative to top bit-score to include in LCA calculation.\n'
            'Higher values yield more conservative classifications by incorporating divergent hits, while lower values enable specific assignments.\n'
            'Applicable when --mMAGs-taxonomy-classify-method is set to BAT.'
        ),
        'mmags_taxonomy_bat_fraction': (
            ['Params_m_sample', '02_taxonomy_BAT', 'fraction'], 0.50, float,
            '\b\n02_taxonomy_BAT: Minimum bit-score support fraction [0-0.99].\n'
            'Higher values enforce stricter consensus (one classification per sequence), while lower values allow multiple lineage detections.\n'
            'Applicable when --mMAGs-taxonomy-classify-method is set to BAT.'
        ),
        'mmags_taxonomy_bat_extra_args': (
            ['Params_m_sample', '02_taxonomy_BAT', 'extra_args'], '', str,
            '02_taxonomy_BAT: Extra arguments forwarded to CAT_pack bins (e.g., "--block_size 12.0").'),
        'mmags_drep99_completeness_per_sample': (
            ['Params_m_sample', '05_drep_eachsample_99', 'comp'], 50, float,
            '05_drep_eachsample_99: The first level dereplication of mMAGs in a sample using completeness threshold [0-100], e.g., 50 for 50%.'
        ),
        'mmags_drep99_contamination_per_sample': (
            ['Params_m_sample', '05_drep_eachsample_99', 'con'], 10, float,
            '05_drep_eachsample_99: The first level dereplication of mMAGs in a sample using contamination threshold [0-100], e.g., 10 for 10%.'
        ),
        'mmags_drep99_ani_per_sample': (
            ['Params_m_sample', '05_drep_eachsample_99', 'sa'], 0.99, float,
            '05_drep_eachsample_99: The first level dereplication of mMAGs in a sample using ANI threshold [0-1], e.g., 0.99 for 99%.'
        ),
        'mmags_drep99_sample_threads': (
            ['Params_m_sample', '05_drep_eachsample_99', 'threads'], 20, int,
            '05_drep_eachsample_99: Threads.'),
        'mmags_drep99_sample_extra_args': (
            ['Params_m_sample', '05_drep_eachsample_99', 'extra_args'], '', str,
            '05_drep_eachsample: Extra arguments forwarded to dRep dereplicate (e.g., "--completeness_weight 1 --contamination_weight 5").'),
        'mmags_drep95_completeness_per_sample': (
            ['Params_m_sample', '06_drep_eachsample_95', 'comp'], 50, float,
            '06_drep_eachsample_95: The second level dereplication of mMAGs in a sample using completeness threshold [0-100], e.g., 50 for 50%.'
        ),
        'mmags_drep95_contamination_per_sample': (
            ['Params_m_sample', '06_drep_eachsample_95', 'con'], 10, float,
            '06_drep_eachsample_95: The second level dereplication of mMAGs in a sample using contamination threshold [0-100], e.g., 10 for 10%.'
        ),
        'mmags_drep95_ani_per_sample': (
            ['Params_m_sample', '06_drep_eachsample_95', 'sa'], 0.95, float,
            '\b\n06_drep_eachsample_95: The second level dereplication of mMAGs in a sample using ANI threshold [0-1], e.g., 0.95 for 95%.\n'
            'This value should be smaller than --drep99-ANI-per-sample'
        ),
        'mmags_drep95_sample_threads': (
            ['Params_m_sample', '06_drep_eachsample_95', 'threads'], 20, int,
            '06_drep_eachsample_95: Threads.'),
        'mmags_drep95_sample_extra_args': (
            ['Params_m_sample', '06_drep_eachsample_95', 'extra_args'], '', str,
            '06_drep_eachsample: Extra arguments forwarded to dRep dereplicate (e.g., "--completeness_weight 1 --contamination_weight 5").'),
        'mmags_coverm_method': (
            ['Params_m_sample', '07_coverm_eachsample', 'method'], 'trimmed_mean', str,
            '\b\n07_coverm_eachsample: Method(s) for coverm-based coverage calculation, must be one of [trimmed_mean/mean].\n'
            'A more detailed description of the different methods is available at https://github.com/wwood/CoverM#calculation-methods'
        ),
        'mmags_coverm_min_read_percent_identity': (
            ['Params_m_sample', '07_coverm_eachsample', 'min_read_percent_identity'], 95, float,
            '07_coverm_eachsample: Exclude reads by overall percent identity [0-100], e.g. 95 for 95%.'
        ),
        'mmags_coverm_min_read_aligned_length': (
            ['Params_m_sample', '07_coverm_eachsample', 'min_read_aligned_length'], 45, int,
            '07_coverm_eachsample: Exclude reads with fewer aligned bases.'
        ),
        'mmags_coverm_contig_end_exclusion': (
            ['Params_m_sample', '07_coverm_eachsample', 'contig_end_exclusion'], 75, int,
            '07_coverm_eachsample: Exclude bases at the ends of contig sequences for coverm calculation.'
        ),
        'mmags_coverm_threads': (
            ['Params_m_sample', '07_coverm_eachsample', 'threads'], 20, int,
            '07_coverm_eachsample: Threads.'),
        'mmags_coverm_sample_extra_args': (
            ['Params_m_sample', '07_coverm_eachsample', 'extra_args'], '', str,
            '07_coverm_eachsample: Extra arguments forwarded to coverm genome (e.g., "--quiet").'),
        'mmags_function_classify_method': (
            ['Params_m_sample', '08_function', 'MAGs_function_classify_method'], 'emapper', str,
            '08_function: mMAGs function classifier using one of [kraken2/emapper].'
        ),
        'mmags_function_threads': (
            ['Params_m_sample', '08_function', 'threads'], 60, int,
            '08_function: Threads for [kraken2/emapper].'),
        'mmags_function_kraken2_confidence': (
            ['Params_m_sample', '08_function_kraken2', 'confidence'], 0.05, float,
            '08_function_kraken2: Confidence threshold [0-1] for functional classification for mMAGs when --mMAGs-function-classify-method is set to kraken2.'
        ),
        'mmags_function_kraken2_extra_args': (
            ['Params_m_sample', '08_function_kraken2', 'extra_args'], '', str,
            '08_function_kraken2: Extra arguments forwarded to kraken2 (e.g., "--memory-mapping").'),
        'mmags_function_emapper_evalue': (
            ['Params_m_sample', '08_function_emapper', 'seed_ortholog_evalue'], 1e-10, float,
            '08_function_emapper: Report only alignments with e-value below or equal to the threshold, applicable when --mMAGs-function-classify-method is set to emapper.'
        ),
        'mmags_function_emapper_method': (
            ['Params_m_sample', '08_function_emapper', 'method'], 'diamond', str,
            '\b\n08_function_emapper: search seed orthologs using one of [diamond/mmseqs], applicable when --mMAGs-function-classify-method is set to emapper.\n'
            'diamond: using default parameters of diamond mode of emapper.\n'
            'mmseqs: using default parameters of mmseqs mode of emapper'
        ),
        'mmags_function_emapper_query_cover': (
            ['Params_m_sample', '08_function_emapper', 'query_cover'], 30, float,
            '\b\n08_function_emapper: Report only alignments with query cover above or equal to the given percentage [0-100],.\n'
            'applicable when --mMAGs-function-classify-method is set to emapper'
        ),
        'mmags_function_emapper_extra_args': (
            ['Params_m_sample', '08_function_emapper', 'extra_args'], '', str,
            '08_function_emapper: Extra arguments forwarded to eggNOG-mapper (e.g., "--pident 30").'),
        'contigs_taxonomy_classify_method': (
            ['Params_c_sample', '01_taxonomy', 'contigs_taxonomy_classify_method'], 'kraken2', str,
            '01_taxonomy: Contigs taxonomy classifier using one of [CAT_pack/kraken2].'
        ),
        'contigs_taxonomy_threads': (
            ['Params_c_sample', '01_taxonomy', 'threads'], 60, int,
            '01_taxonomy: Threads for [CAT_pack/kraken2].'
        ),
        'contigs_taxonomy_cat_range': (
            ['Params_c_sample', '01_taxonomy_CAT_pack', 'range'], 10, float,
            '\b\n01_taxonomy_CAT_pack: Percentage range [0-100] of DIAMOND hits relative to top bit-score to include in LCA calculation.\n'
            'Higher values yield more conservative classifications by incorporating divergent hits, while lower values enable specific assignments.\n'
            'Applicable when --contigs-taxonomy-classify-method is set to CAT_pack.'
        ),
        'contigs_taxonomy_cat_fraction': (
            ['Params_c_sample', '01_taxonomy_CAT_pack', 'fraction'], 0.50, float,
            '\b\n01_taxonomy_CAT_pack: Minimum bit-score support fraction [0-0.99].\n'
            'Higher values enforce stricter consensus (one classification per sequence), while lower values allow multiple lineage detections.\n'
            'Applicable when --contigs-taxonomy-classify-method is set to CAT_pack.'
        ),
        'contigs_taxonomy_cat_pack_extra_args': (
            ['Params_c_sample', '01_taxonomy_CAT_pack', 'extra_args'], '', str,
            '01_taxonomy_CAT_pack: Extra arguments forwarded to CAT_pack contigs (e.g., "--block_size 12.0").'),
        'contigs_taxonomy_kraken2_confidence': (
            ['Params_c_sample', '01_taxonomy_kraken2', 'confidence'], 0.7, float,
            '\b\n01_taxonomy_kraken2: Confidence threshold [0-1] for taxonomic classification, applicable when --contigs-taxonomy-classify-method is set to kraken2.\n'
            'Default value was tuned on the CAMI2 marine dataset. For samples with higher complexity, a lower confidence threshold is usually recommended.'
        ),
        'contigs_taxonomy_kraken2_extra_args': (
            ['Params_c_sample', '01_taxonomy_kraken2', 'extra_args'], '', str,
            '01_taxonomy_kraken2: Extra arguments forwarded to kraken2 (e.g., "--memory-mapping").'),
        'contigs_function_classify_method': (
            ['Params_c_sample', '02_function', 'contigs_function_classify_method'], 'kraken2', str,
            '02_function: Contigs function classifier using one of [kraken2/emapper].'
        ),
        'contigs_function_threads': (
            ['Params_c_sample', '02_function', 'threads'], 60, int,
            '02_function: Threads for [kraken2/emapper].'),
        'contigs_function_emapper_evalue': (
            ['Params_c_sample', '02_function_emapper', 'seed_ortholog_evalue'], 1e-10, float,
            '02_function_emapper: Report only alignments with e-value below or equal to the threshold, applicable when --contigs-function-classify-method is set to emapper.'
        ),
        'contigs_function_emapper_method': (
            ['Params_c_sample', '02_function_emapper', 'method'], 'diamond', str,
            '\b\n02_function_emapper: search seed orthologs using one of [diamond/mmseqs], applicable when --contigs-function-classify-method is set to emapper.\n'
            'diamond: using default parameters of diamond mode of emapper.\n'
            'mmseqs: using default parameters of mmseqs mode of emapper'
        ),
        'contigs_function_emapper_query_cover': (
            ['Params_c_sample', '02_function_emapper', 'query_cover'], 30, float,
            '\b\n02_function_emapper: Report only alignments with query cover above or equal to the given percentage [0-100],.\n'
            'applicable when --contigs-function-classify-method is set to emapper'
        ),
        'contigs_function_emapper_extra_args': (
            ['Params_c_sample', '02_function_emapper', 'extra_args'], '', str,
            '02_function_emapper: Extra arguments forwarded to eggNOG-mapper (e.g., "--pident 30").'),
        'contigs_function_kraken2_confidence': (
            ['Params_c_sample', '02_function_kraken2', 'confidence'], 0.05, float,
            '02_function_kraken2: Confidence threshold [0-1] for functional classification, applicable when --contigs-function-classify-method is set to kraken2.'
        ),
        'contigs_function_kraken2_extra_args': (
            ['Params_c_sample', '02_function_kraken2', 'extra_args'], '', str,
            '02_function_kraken2: Extra arguments forwarded to kraken2 (e.g., "--memory-mapping").'),
        'merge_mmags_min_completeness1': (
            ['Params_vmc_sample', '05_merge_tax_fun_depth', 'MAGs_min_completeness1'], 50, float,
            '02_quality_gate: Minimum completeness threshold (%) for CheckM1-based mMAG quality filtering [0-100], e.g., 50 for 50%.'
        ),
        'merge_mmags_max_contamination1': (
            ['Params_vmc_sample', '05_merge_tax_fun_depth', 'MAGs_max_contamination1'], 1, float,
            '02_quality_gate: Maximum contamination threshold (%) for CheckM1-based mMAG quality filtering [0-100], e.g., 5 for 5%.'
        ),
        'merge_mmags_min_completeness2': (
            ['Params_vmc_sample', '05_merge_tax_fun_depth', 'MAGs_min_completeness2'], 50, float,
            '02_quality_gate: Minimum completeness threshold (%) for CheckM2-based mMAG quality filtering [0-100], e.g., 50 for 50%.'
        ),
        'merge_mmags_max_contamination2': (
            ['Params_vmc_sample', '05_merge_tax_fun_depth', 'MAGs_max_contamination2'], 1, float,
            '02_quality_gate: Maximum contamination threshold (%) for CheckM2-based mMAG quality filtering [0-100], e.g., 5 for 5%.'
        ),
        'merge_min_orf_length': (
            ['Params_vmc_sample', '05_merge_tax_fun_depth', 'minORFLen'], 50, int,
            '05_merge_tax_fun_depth: Minimum ORF length for finally used ORFs of vMAGs, mMAGs, and contigs.'
        ),
        'merge_vmags_checkv_quality': (
            ['Params_vmc_sample', '05_merge_tax_fun_depth', 'checkv_quality'], 'High-quality', str,
            '\b\n05_merge_tax_fun_depth: CheckV quality filter for viral MAGs (lower bound filter).\n'
            'viral MAGs with quality equal to or higher than the specified threshold will be retained.\n'
            'Options: Complete, High-quality, Medium-quality, Low-quality, Not-determined.\n'
            'Complete: complete genome (highest quality)\n'
            'High-quality: genome completeness ≥90%\n'
            'Medium-quality: genome completeness 50-90%\n'
            'Low-quality: genome completeness <50%\n'
            'Not-determined: completeness could not be determined (lowest quality)'
        ),
        'merge_multiple_ko_handling_mode': (
            ['Params_vmc_sample', '05_merge_tax_fun_depth', 'KO_mode'], 'ignore_multiple_KOs', str,
            '\b\n05_merge_tax_fun_depth: Method using one of [split_multiple_KOs/keep_multiple_KOs/ignore_multiple_KOs] when an ORF has multiple KO annotations:.\n'
            'split_multiple_KOs: evenly distributing the total abundance across all assigned KOs;\n'
            'keep_multiple_KOs: retain original multiple KO annotations;\n'
            'ignore_multiple_KOs: ignore the ORFs with multiple KO annotations and only use the ORFs with single KO annotation.'
        ),
        'merge_min_ko_num_per_taxon': (
            ['Params_vmc_sample', '05_merge_tax_fun_depth', 'min_KO_num_in_one_taxon'], 100, int,
            '\b\n05_merge_tax_fun_depth: Minimum number of KOs per taxon, used to filter outputs in KO_Taxa_RelativeAbundance and KO_Taxa_SumCoverage.\n'
            'This threshold filters out taxa with insufficient KOs,\n'
            'retaining only KO-rich taxa to ensure robust metabolic reconstruction\n'
            'while preventing excessive output file sizes caused by sparse matrices.'
        ),
        'merge_tax_fun_depth_threads': (
            ['Params_vmc_sample', '05_merge_tax_fun_depth', 'threads'], 20, int,
            '05_merge_tax_fun_depth: Threads.'),
        'merge_cami_format_output': (
            ['Params_vmc_sample', '05_merge_tax_fun_depth', 'CAMI_format_output'], False, bool,
            '05_merge_tax_fun_depth: Report a taxonomic profile using the CAMI output format (flag only). Disabled when mMAG mode uses gtdbtk.'
        ),
        'merge_output_detail': (
            ['Params_vmc_sample', '05_merge_tax_fun_depth', 'output_detail'], False, bool,
            '05_merge_tax_fun_depth: Whether to output detail tables (flag only).'),
        'mclink_sample_min_ko_num_per_taxon': (
            ['Params_vmc_sample', '06_mclink_taxon', 'min_KO_num_in_one_taxon'], 100, int,
            '06_mclink_taxon: Minimum number of KOs per taxon for calculating taxon-level metabolic pathway/module completeness or abundance per sample. Taxa with fewer than this number of detected KOs in a given sample will be filtered out. Calculated independently for each sample.'
        ),
        'mclink_sample_min_ko_abundance_per_taxon': (
            ['Params_vmc_sample', '06_mclink_taxon', 'min_mean_KO_RA_in_one_taxon'], 0, float,
            '06_mclink_taxon: Minimum mean KO abundance [0-1e6] per taxon for calculating taxon-level metabolic pathway/module completeness or abundance per sample. Taxa with mean KO abundance below this threshold in a given sample will be filtered out. Note: GPM/TPM sum of all taxa in each sample is 1e6.'
        ),
        'mclink_sample_topn_taxon': (
            ['Params_vmc_sample', '06_mclink_taxon', 'TopN'], 100, int,
            '06_mclink_taxon: Number of top taxa (ranked by mean KO abundance) retained for calculating taxon-level metabolic pathway/module completeness or abundance per sample. Set to 0 to disable filtering.'
        ),
        'mclink_sample_completeness_per_taxon': (
            ['Params_vmc_sample', '06_mclink_taxon', 'completeness_threshold'], 0.5, float,
            '06_mclink_taxon: Completeness threshold [0-1] of modules for calculating taxon-level metabolic pathway/module abundance. Higher values retain only more complete profiles.'),
        'mclink_sample_plus_scale_method_per_taxon': (
            ['Params_vmc_sample', '06_mclink_taxon', 'plus_scale_method'], 'mean', str,
            "06_mclink_taxon: Aggregation method for '+'-separated KOs in a module [mean/max/min]. 'mean': average abundance across all required KOs; 'max': highest KO abundance; 'min': lowest KO abundance."),
        'mclink_sample_comma_scale_method_per_taxon': (
            ['Params_vmc_sample', '06_mclink_taxon', 'comma_scale_method'], 'sum', str,
            "06_mclink_taxon: Aggregation method for ','-separated alternative KOs in a module [sum/max]. 'sum': cumulative abundance of all alternative KOs; 'max': abundance of the highest single KO among alternatives. Note: this applies to abundance calculation only; module completeness is always computed with 'max'."),
        'barplot_topn_taxon': (
            ['Params_vmc_project', '01_Taxa_visualize', 'TopN'], 10, int,
            '01_Taxa_visualize: Number of top taxa to display in stacked barplots. Remaining categories are merged into Others. Maximum 29.'
        ),
        'barplot_topn_fun': (
            ['Params_vmc_project', '03_Fun_visualize', 'TopN'], 10, int,
            '03_Fun_visualize: Number of top functional categories to display in stacked barplots. Remaining categories are merged into Others. Maximum 29.'
        ),
        'mclink_min_ko_num_per_sample': (
            ['Params_vmc_project', '03_mclink_sample', 'min_sample_num_per_KO'], 1, int,
            '03_mclink_sample: Minimum number of samples in which a KO must be detected for calculating community-level metabolic pathway/module completeness or abundance.'
        ),
        'mclink_min_ko_prevalence_per_sample': (
            ['Params_vmc_project', '03_mclink_sample', 'preva'], 0, float,
            '03_mclink_sample: Minimum detection rate [0-1] of a KO across all samples for calculating community-level metabolic pathway/module completeness or abundance.'
        ),
        'mclink_min_ko_abundance_per_sample': (
            ['Params_vmc_project', '03_mclink_sample', 'abund'], 0, float,
            '03_mclink_sample: Minimum mean relative abundance [0-1e6] of a KO across all samples for calculating community-level metabolic pathway/module completeness or abundance. KOs with mean abundance below this threshold are excluded. Note: GPM/TPM sum of a sample is 1e6.'
        ),
        'mclink_completeness_per_sample': (
            ['Params_vmc_project', '03_mclink_sample', 'completeness_threshold'], 0.5, float,
            '03_mclink_sample: Completeness threshold [0-1] of modules for calculating community-level metabolic pathway/module abundance. Higher values retain only more complete profiles.'
        ),
        'mclink_plus_scale_method_per_sample': (
            ['Params_vmc_project', '03_mclink_sample', 'plus_scale_method'], 'mean', str,
            "03_mclink_sample: Aggregation method for '+'-separated KOs in a module [mean/max/min]. 'mean': average abundance across all required KOs; 'max': highest KO abundance; 'min': lowest KO abundance."
        ),
        'mclink_comma_scale_method_per_sample': (
            ['Params_vmc_project', '03_mclink_sample', 'comma_scale_method'], 'sum', str,
            "03_mclink_sample: Aggregation method for ','-separated alternative KOs in a module [sum/max]. 'sum': cumulative abundance of all alternative KOs; 'max': abundance of the highest single KO among alternatives. Note: this applies to abundance calculation only; module completeness is always computed with 'max'."
        ),
        'merge_min_taxa_ko_prevalence_across_samples': (
            ['Params_vmc_project', '05_merge_tax_KO_depth', 'preva'], 0.10, float,
            '05_merge_tax_KO_depth: Minimum prevalence [0-1] of a Taxon-KO across samples. Taxa-KOs detected in fewer samples than this proportion are excluded.'
        ),
        'merge_min_taxa_ko_abundance_across_samples': (
            ['Params_vmc_project', '05_merge_tax_KO_depth', 'abund'], 0.0, float,
            '05_merge_tax_KO_depth: Minimum mean relative abundance [0-1e6] of a Taxon-KO across samples. Taxa-KOs with mean abundance below this threshold are excluded. Note: GPM/TPM sum of a sample is 1e6.'
        ),
        'merge_min_taxa_module_prevalence_across_samples': (
            ['Params_vmc_project', '05_merge_tax_Module_RelativeAbundance', 'preva'], 0.10, float,
            '05_merge_tax_Module_RelativeAbundance: Minimum prevalence [0-1] of a Taxon-Module across samples. Taxa-Modules detected in fewer samples than this proportion are excluded.'
        ),
        'merge_min_taxa_module_abundance_across_samples': (
            ['Params_vmc_project', '05_merge_tax_Module_RelativeAbundance', 'abund'], 0.0, float,
            '05_merge_tax_Module_RelativeAbundance: Minimum mean relative abundance [0-1e6] of a Taxon-Module across samples. Taxa-Modules with mean abundance below this threshold are excluded. Note: GPM/TPM sum of a sample is 1e6.'
        ),
        'merge_min_taxa_module_sc_prevalence_across_samples': (
            ['Params_vmc_project', '05_merge_tax_Module_SumCoverage', 'preva'], 0.10, float,
            '05_merge_tax_Module_SumCoverage: Minimum prevalence [0-1] of a Taxon-Module across samples. Taxa-Modules detected in fewer samples than this proportion are excluded.'
        ),
        'merge_min_taxa_module_sc_abundance_across_samples': (
            ['Params_vmc_project', '05_merge_tax_Module_SumCoverage', 'abund'], 0.0, float,
            '05_merge_tax_Module_SumCoverage: Minimum mean sum coverage [0-1e6] of a Taxon-Module across samples. Taxa-Modules with mean sum coverage below this threshold are excluded. Note: GPM/TPM sum of a sample is 1e6.'
        ),
        'mclink_project_min_ko_num_per_taxon': (
            ['Params_vmc_project', '06_mclink_taxon', 'min_KO_num_in_one_taxon'], 100, int,
            '06_mclink_taxon: Minimum number of KOs per taxon for calculating taxon-level metabolic pathway/module completeness on the project-averaged taxon-function matrix. Taxa with fewer than this number of detected KOs in the project-averaged matrix will be filtered out. KO abundances are first averaged across all samples to generate a single project-level matrix before calculation.'
        ),
        'mclink_project_min_ko_abundance_per_taxon': (
            ['Params_vmc_project', '06_mclink_taxon', 'min_mean_KO_RA_in_one_taxon'], 0, float,
            '06_mclink_taxon: Minimum mean KO abundance [0-1e6] per taxon for calculating taxon-level metabolic pathway/module completeness on the project-averaged taxon-function matrix. Taxa with mean KO abundance below this threshold in the project-averaged matrix will be filtered out. Note: GPM/TPM sum of all taxa is 1e6 per sample; values are averaged across samples prior to filtering and calculation.'
        ),
        'mclink_project_topn_taxon': (
            ['Params_vmc_project', '06_mclink_taxon', 'TopN'], 100, int,
            '06_mclink_taxon: Number of top taxa (ranked by mean KO abundance in the project-averaged taxon-function matrix) retained for calculating taxon-level metabolic pathway/module completeness. Set to 0 to disable filtering.'
        ),
        'mclink_project_completeness_per_taxon': (
            ['Params_vmc_project', '06_mclink_taxon', 'completeness_threshold'], 0.5, float,
            '06_mclink_taxon: Completeness threshold [0-1] of modules for calculating taxon-level metabolic pathway/module abundance on the project-averaged taxon-function matrix. Higher values retain only more complete profiles.'),
        'mclink_project_plus_scale_method_per_taxon': (
            ['Params_vmc_project', '06_mclink_taxon', 'plus_scale_method'], 'mean', str,
            "06_mclink_taxon: Aggregation method for '+'-separated KOs in a module [mean/max/min] on the project-averaged taxon-function matrix. 'mean': average abundance across all required KOs; 'max': highest KO abundance; 'min': lowest KO abundance."),
        'mclink_project_comma_scale_method_per_taxon': (
            ['Params_vmc_project', '06_mclink_taxon', 'comma_scale_method'], 'sum', str,
            "06_mclink_taxon: Aggregation method for ','-separated alternative KOs in a module [sum/max] on the project-averaged taxon-function matrix. 'sum': cumulative abundance of all alternative KOs; 'max': abundance of the highest single KO among alternatives. Note: this applies to abundance calculation only; module completeness is always computed with 'max'."),
    }
}

def _build_snakemake_cmd(config, jobs, verbose, silent, dryrun,
                        snakemake_args):
    if verbose and not silent:
        quiet_flag = ''
    else:
        quiet_flag = '--quiet'
    env_prefix = 'METALINK_SILENT=1 ' if silent else ''
    return (
        '{env_prefix}snakemake --snakefile {snakefile} --directory {work_dir} '
        '--jobs {jobs} --configfile {config_f} '
        '--use-conda {dryrun} {quiet} {args}'
    ).format(
        env_prefix=env_prefix,
        snakefile=shlex.quote(config['workflow_snakefile']),
        work_dir=shlex.quote(config.get('Directory', {}).get('work', '.')),
        jobs=jobs,
        config_f=shlex.quote(config['_path']),
        dryrun='--dry-run' if dryrun else '',
        quiet=quiet_flag,
        args=' '.join(shlex.quote(a) for a in snakemake_args)
    )

def exec_workflow(mode, input_tsv, working_dir, jobs,
                  verbose, silent, dryrun,
                  db_dir, project_name, snakemake_args,
                  cli_provided=None, **kwargs):

    working_dir = os.path.abspath(working_dir)
    input_tsv = os.path.abspath(input_tsv)
    os.makedirs(working_dir, exist_ok=True)
    os.chdir(working_dir)
    logging.info(f"Working directory is set to: {os.getcwd()}")

    config_dir = os.path.join(working_dir, 'config')
    os.makedirs(config_dir, exist_ok=True)
    config_f = os.path.join(config_dir, _generated_config_name(mode, kwargs))
    make_config(mode=mode,
                input_tsv=input_tsv,
                threads=jobs,
                cli_overrides=kwargs,
                cli_provided=cli_provided,
                config_f=config_f,
                db_dir=db_dir,
                working_dir=working_dir,
                project_name=project_name)
    if not silent:
        logging.info(f'Config file written to {config_f}')
        logging.info(f'Using config file: {config_f}')
    config = load_configfile(config_f)
    config['_path'] = config_f
    cmd = _build_snakemake_cmd(config, jobs, verbose, silent,
                              dryrun, snakemake_args)
    if not silent:
        logging.info('Executing: %s' % cmd)
    _t0 = time.time()
    subprocess.run(cmd, shell=True, check=True)
    _warn_zero_kraken2_classification(
        working_dir, _t0,
        _get_nested(config, ['Params_r_sample', '03_taxonomy', 'bracken_level'])
        if mode == 'r' else None,
        _get_nested(config, ['Params_r_sample', '03_taxonomy', 'bracken_threshold'])
        if mode == 'r' else 0)

def exec_workflow_with_existing_config(mode, config_f, jobs, verbose, silent,
                             dryrun, snakemake_args, **kwargs):
    config = load_configfile(config_f)
    if not silent:
        logging.info(f'Using config file: {config_f}')
    if config.get('mode') != mode:
        raise ValueError(
            f"Config mode '{config.get('mode')}' != CLI mode '{mode}'"
        )

    # CLI mode-specific parameters override config values when explicitly provided
    _mags = None
    if mode == 'pmr':
        # An explicit --mags on the command line wins over the config value;
        # it also decides which Params_pmr_{m,v,vm}_* sections the {mags}
        # placeholder yaml_paths below resolve to.
        _mags = kwargs.get('mags') or config.get('mags') or 'm'
        if _mags not in _PMR_MAGS_VALUES:
            raise ValueError(
                "--mags must be one of {}, got: {}".format(list(_PMR_MAGS_VALUES), _mags)
            )
        # The config must carry the sections of the resolved --mags chain;
        # otherwise it came from an older version or a different --mags value
        # and rules would fail later with opaque KeyErrors.
        _required = MODE_SPECIFIC.get('pmr_{}'.format(_mags), [])
        _missing = [s for s in _required if s not in config]
        if _missing:
            raise ValueError(
                "Config {} lacks section(s) {} required for --mags {}; it was likely "
                "generated by an older Metalink version or a different --mags value. "
                "Please regenerate it (metalink pmr --mags {} --generate-config ...)."
                .format(config_f, _missing, _mags, _mags)
            )
    overridden = False
    for opt_name, (yaml_path, default, typ, help_txt) in PARAMS_HELP[mode].items():
        user_val = kwargs.get(opt_name)
        if user_val is None:
            continue
        if typ == bool and not user_val:
            continue  # bool flag=False treated as not provided
        option_name = f"--{opt_name.replace('_', '-')}"
        if mode == 'pmr':
            yaml_path = _resolve_mags_placeholder(yaml_path, _mags)
            if not _pmr_option_applicable(opt_name, yaml_path, config, _mags):
                continue
        if not silent:
            logging.info(f'Setting config {" -> ".join(yaml_path)} to {user_val!r} (from CLI {option_name})')
        set_nested(config, yaml_path, user_val)
        overridden = True

    # result-neutral or additive options that may also override an existing
    # config (they are not PARAMS_HELP entries): taxonomy ranks, project-level
    # switch, intermediate-file cleanup, kraken2 taxonomy concurrency
    _tr = kwargs.get('taxonomy_rank')
    if _tr is not None:
        config['taxonomy_rank'] = [r.strip() for r in _tr.split(',')] if isinstance(_tr, str) else _tr
        overridden = True
        if not silent:
            logging.info(f"Setting config taxonomy_rank to {config['taxonomy_rank']!r} (from CLI --taxonomy-rank)")
    if kwargs.get('run_project'):
        config['run_project'] = True
        overridden = True
        if not silent:
            logging.info("Setting config run_project to True (from CLI --run-project)")
    _pn = kwargs.get('project_name')
    if _pn is not None:
        config['name_project'] = [_pn] if isinstance(_pn, str) else _pn
        overridden = True
        if not silent:
            logging.info(f"Setting config name_project to {config['name_project']!r} (from CLI --project-name)")
    _ki = config.get('keep_intermediate')
    config['keep_intermediate'] = dict(_ki) if isinstance(_ki, dict) else \
        {'bam': True, 'drep_data': False, 'vs2_iter0': False, 'assembly': False}
    for _flag, _cat in (('delete_bam', 'bam'),
                        ('keep_drep_data', 'drep_data'),
                        ('keep_vs2_iter0', 'vs2_iter0'),
                        ('keep_assembly_intermediate', 'assembly')):
        if kwargs.get(_flag, False):
            config['keep_intermediate'][_cat] = _flag.startswith('keep_')
            overridden = True
            if not silent:
                logging.info(f"Setting config keep_intermediate.{_cat} to {config['keep_intermediate'][_cat]!r} (from CLI --{_flag.replace('_','-')})")
    _ktj = kwargs.get('kraken2_taxonomy_jobs')
    if _ktj is not None:
        config['kraken2_taxonomy_jobs'] = int(_ktj)
        overridden = True
        if not silent:
            logging.info(f"Setting config kraken2_taxonomy_jobs to {int(_ktj)!r} (from CLI --kraken2-taxonomy-jobs)")

    # Enforce parameter compatibility after applying CLI overrides.
    if mode == 'pmr':
        if _migrate_pmr_quality_gate(config, _mags):
            overridden = True
        _validate_pmr_iphop_custom_mags_dir(config, _mags)
    if _validate_cami_gtdbtk_compatibility(config, mode):
        overridden = True

    # skip_level is recorded in the config by make_config; for configs written
    # by older versions it is inferred from workflow_includes. skip_level is
    # authoritative: workflow_includes are rebuilt from it so that editing
    # skip_level in the config file takes effect.
    skip_level = _infer_skip_level(config)
    if skip_level not in WORKFLOW_MAP:
        raise ValueError(f"skip_level in config must be one of [0/1/2/3], got: {skip_level}")
    config['skip_level'] = skip_level
    _rules_dir = _rules_dir_for(config.get('read_type'))
    _workflow_includes = ["{}/{}".format(_rules_dir, f) for f in _workflow_smk_files(mode, skip_level, _mags)]
    _common_includes = ["{}/init.smk".format(_rules_dir)]
    if (config.get('workflow_includes') != _workflow_includes
            or config.get('workflow_common_includes') != _common_includes):
        config['workflow_includes'] = _workflow_includes
        config['workflow_common_includes'] = _common_includes
        overridden = True
        if not silent:
            logging.info(f'workflow_includes rebuilt from skip_level={skip_level}')

    _validate_read_type_assembly_method(config, mode, skip_level)
    _validate_bracken_level(config, mode)

    if overridden:
        tmp_dir = os.path.dirname(os.path.abspath(config_f))
        fd, tmp_path = tempfile.mkstemp(
            suffix='.yaml',
            prefix=f'config.{mode}.runtime.',
            dir=tmp_dir
        )
        with os.fdopen(fd, 'w') as f:
            yaml.dump(config, f)
        config_f = tmp_path
        if not silent:
            logging.info(f'CLI overrides merged to temporary config: {tmp_path}')

    config['_path'] = config_f
    cmd = _build_snakemake_cmd(config, jobs, verbose, silent,
                              dryrun, snakemake_args)
    if not silent:
        logging.info('Executing: %s' % cmd)
    _t0 = time.time()
    subprocess.run(cmd, shell=True, check=True)
    _warn_zero_kraken2_classification(
        config.get('Directory', {}).get('work', '.'), _t0,
        _get_nested(config, ['Params_r_sample', '03_taxonomy', 'bracken_level'])
        if mode == 'r' else None,
        _get_nested(config, ['Params_r_sample', '03_taxonomy', 'bracken_threshold'])
        if mode == 'r' else 0)

def exec_workflow_generate_config_only(mode, input_tsv, jobs, working_dir,
                                       db_dir, project_name=None,
                                       cli_provided=None, **kwargs):
    config_dir = os.path.join(os.path.abspath(working_dir), 'config')
    os.makedirs(config_dir, exist_ok=True)
    config_f = os.path.join(config_dir, _generated_config_name(mode, kwargs))
    make_config(mode=mode,
                input_tsv=input_tsv,
                threads=jobs,
                cli_overrides=kwargs,
                cli_provided=cli_provided,
                config_f=config_f,
                db_dir=db_dir,
                working_dir=working_dir,
                project_name=project_name)
    logging.info('Config file written to %s', config_f)

def add_mode_options(mode):
    def decorator(f):
        from collections import OrderedDict
        grouped = OrderedDict()
        for opt_name, (yaml_path, default, typ, help_txt) in PARAMS_HELP[mode].items():
            prefix = yaml_path[0] if yaml_path else ''
            if prefix not in grouped:
                grouped[prefix] = []
            grouped[prefix].append((opt_name, yaml_path, default, typ, help_txt))

        section_items = list(grouped.items())
        _SECTION_ORDER[mode] = [get_section_name(prefix) for prefix, _ in section_items]

        # click decorators are applied bottom-to-top; reversed() ensures display order matches source order
        for prefix, params in reversed(section_items):
            section_name = get_section_name(prefix)
            for opt_name, yaml_path, default, typ, help_txt in reversed(params):
                if yaml_path and yaml_path[-1] == 'confidence' and 'kraken2' in str(yaml_path[-2]):
                    help_txt += ('\nWith a reduced kraken2 database (e.g. kraken2_16GB), generated '
                                 'configs default this to 0; the default shown applies to the full '
                                 'database (kraken2_raw).')
                if typ == bool:
                    f = click.option(
                        f'--{opt_name.replace("_", "-")}',
                        is_flag=True,
                        default=default,
                        help=f'{help_txt} (default: {default})'
                    )(f)
                else:
                    f = click.option(
                        f'--{opt_name.replace("_", "-")}',
                        type=typ,
                        default=None,
                        help=f'{help_txt} (default: {default})'
                    )(f)
                if hasattr(f, '__click_params__') and f.__click_params__:
                    last_opt = f.__click_params__[-1]
                    last_opt._section_name = section_name
        return f
    return decorator

common_options = [
    click.option('--read-type',
                 type=str, default="paired", help=
                '\b\nInput read type, must be one of [paired/paired-rna/single/single-rna/pacbio-raw/pacbio-corr/pacbio-hifi/nano-raw/nano-corr/nano-hq].\n'
                'paired: paired-end DNA short reads;\n'
                'paired-rna: paired-end RNA short reads;\n'
                'single: single-end DNA short reads;\n'
                'single-rna: single-end RNA short reads;\n'
                'pacbio-raw/pacbio-corr/pacbio-hifi: PacBio long reads (raw/corrected/HiFi);\n'
                'nano-raw/nano-corr/nano-hq: Oxford Nanopore long reads (raw/corrected/high-quality).\n'
                'Mutually exclusive with -c/--configfile; edit read_type in the config file instead.'
                 ),
    click.option('--skip-level',
                 type=click.IntRange(0, 3), default=0, help=
                 '\b\nSkip level for input data preparation, using one of [0/1/2/3].\n'
                 '   0: Run full pipeline including trim, assembly, binning, and other processes (default).\n'
                 '   1: Clean reads provided; skip trim, start from assembly.\n'
                 '   2: Clean reads and contigs provided; skip trim and assembly.\n'
                 '   3: Clean reads, contigs, and mMAGs provided; skip trim, assembly, and binning.\n'
                 '   Mutually exclusive with -c/--configfile; use the config file to change the skip level.'
                 ),
    click.option('-i', '--input-tsv', required=False,
                type=click.Path(exists=True), help=
                '\b\nA tab-separated file specifying sample paths. Metalink will auto-generate a config.<mode>.yaml from this file for the current mode.\n'
                'Format: header required. The first column must be "sample"; remaining columns are file paths required by the selected mode.\n'
                '   --skip-level 0: sample\tread_1\tread_2\n'
                '   --skip-level 1: sample\tread_1\tread_2\n'
                '   --skip-level 2: sample\tread_1\tread_2\tcontig_file\n'
                '   --skip-level 3: sample\tread_1\tread_2\tcontig_file\tmMAG_dir\n'
                'Column descriptions:\n'
                '   read_1, read_2: FASTQ (*.gz) files. When read type is single, use read_s.\n'
                '   contig_file:   contigs FASTA (*.fa) file.\n'
                '   mMAG_dir:      mMAGs FASTA (*.fa) files.\n'
                'Constraints:\n'
                '   Each mMAG file must end with .fa suffix.\n'
                '   Internal sequence headers of mMAGs must be original contig IDs exactly matching those in contig_file.\n'
                '   All files must not contain three consecutive dashes (---) as they are reserved as internal delimiters.\n'
                'Mutually exclusive with -c/--configfile.'
                ),
    click.option('-d', '--db-dir',
                type=click.Path(dir_okay=True, resolve_path=True),
                default=None,
                help=
                '\b\nPath to the reference database directory.\n'
                'Must contain all required databases for the selected mode (e.g., Kraken2\n'
                'databases for taxonomy, EggNOG databases for function).\n'
                'This value is written to Directory.database.DBDIR and is only used to resolve\n'
                'database paths that are relative. Absolute paths in the config are preserved\n'
                'unchanged, so per-tool database overrides remain effective.\n'
                'If not specified, the DBDIR value from the config template will be used.\n'
                'Mutually exclusive with -c/--configfile; edit Directory.database.DBDIR in the config file instead.'
                 ),
    click.option('-w', '--working-dir', default='.',
                 type=click.Path(dir_okay=True, writable=True),
                 help=
                '\b\nWorking directory for the workflow execution.\n'
                'All output files, logs, and intermediate results will be written here.\n'
                'Default is the current directory (.).\n'
                'Mutually exclusive with -c/--configfile; edit Directory.work in the config file instead.'
                 ),
    click.option('-c', '--configfile', metavar='PATH',
                 type=click.Path(exists=True),
                 help=
                 '\b\nUse an existing config.yaml instead of auto-generating one from the built-in template.\n'
                 'Config-source options (--read-type, --skip-level, -i, -d, -w) cannot be used together with -c; edit them directly in the existing config file.\n'
                 'Other options (e.g., --project-name, --run-project, --taxonomy-rank, --kraken2-taxonomy-jobs,\n'
                 'and --trim-sample-threads) override the corresponding config values for that run.'
                 ),
    click.option('-g', '--generate-config', is_flag=True,
                 help=
                '\b\nGenerate a config.yaml template to <working_dir>/config/config.<mode>.yaml and exit (flag only).\n'
                'You can then manually edit all workflow parameters,\n'
                'e.g., rule-level threads, tool-specific flags and run with:\n'
                'metalink <mode> --configfile <working_dir>/config/config.<mode>.yaml.\n'
                'If used with dry-run (-pn), config is written to the working directory.'
                 ),
    click.option('-j', '--jobs', default=multiprocessing.cpu_count(),
                 type=int,
                 help=
                '\b\nMaximum number of parallel jobs to run simultaneously.\n'
                'In local mode, this limits the number of concurrent threads/processes.\n'
                'In cluster mode (with Snakemake --profile, passed through as an extra\n'
                'argument), this limits the number of concurrent submitted jobs.\n'
                'Default is the number of CPU cores on the current machine.'
                 ),
    click.option('-n', '--dry-run', 'dryrun', is_flag=True,
                 help=
                '\b\nPerform a dry-run without executing any jobs (flag only).\n'
                'The execution plan will be printed, showing which rules would be run\n'
                'and their dependencies, but no actual computation is performed.\n'
                'Useful for validating workflow logic and checking input/output expectations.'
                 ),
    click.option('--run-project', is_flag=True,
                help='Enable project-level analysis (flag only). '
                    'If not set, only sample-level analysis is performed. '
                    'With -c/--configfile, enables it for that run (runtime override).',
                default=False),
    click.option('--delete-bam', is_flag=True,
                help='Delete mapping BAMs (b/v-sample/v-project) after last use '
                    '(flag only; default: kept; config key keep_intermediate.bam).',
                default=False),
    click.option('--keep-drep-data', is_flag=True,
                help='Keep dRep data/ working dirs '
                    '(flag only; default: deleted; config key keep_intermediate.drep_data).',
                default=False),
    click.option('--keep-vs2-iter0', is_flag=True,
                help='Keep VirSorter2 iter-0 dirs '
                    '(flag only; default: deleted; config key keep_intermediate.vs2_iter0).',
                default=False),
    click.option('--keep-assembly-intermediate', is_flag=True,
                help='Keep assembler/metaWRAP intermediate files '
                    '(flag only; default: deleted; config key keep_intermediate.assembly).',
                default=False),
    click.option('--project-name',
                 help='Project name for the samples. '
                      'With -c/--configfile, overrides name_project for that run (runtime override).',
                 default="test.metalink"),
    click.option('--taxonomy-rank',
                 type=str, default="Phylum,Class,Order,Family,Genus,Species", help=
                '\b\nComma-separated taxonomy ranks controlling aLL rank-resolved outputs, using any of the [Phylum/Class/Order/Family/Genus/Species].\n'
                'Overrides taxonomy_rank in the config template.\n'
                'With -c/--configfile, overrides taxonomy_rank for that run (runtime override).'
                 ),
    click.option('--kraken2-taxonomy-jobs',
                 type=int, default=None,
                 help=
                '\b\nMaximum number of kraken2 taxonomy jobs to run concurrently. kraken2 taxonomy rules are resource-tagged (kraken2_taxonomy=1). \n'
                'The value is written to config key kraken2_taxonomy_jobs and applied by the workflow as a Snakemake global resource limit.\n'
                'With -c/--configfile, overrides it for that run (runtime override). Default: unlimited.'),
    click.option('--verbose', is_flag=True,
                 help=
                '\b\nEnable verbose output during workflow execution (flag only).\n'
                'When enabled, detailed execution logs are printed, including rule activation,\n'
                'job submission, and resource usage. Useful for debugging workflow failures.\n'
                'When disabled (default), the workflow runs in quiet mode.'
                 ),
    click.option('--silent', is_flag=True,
                 help=
                '\b\nSilence all Metalink messages (flag only).\n'
                'When enabled, Metalink initialization and execution logs are suppressed.'
                 ),
    click.argument('snakemake_args', nargs=-1, type=click.UNPROCESSED)
]

def _validate_passthrough_args(mode, snakemake_args):
    """Reject Metalink-looking options that were silently swallowed into
    snakemake_args (ignore_unknown_options=True passes them through).

    Catches three mistake classes:
      1. single-dash long options: -merge-output-detail (meant: --merge-output-detail)
      2. valid Metalink options used with the wrong mode (silently ignored)
      3. close misspellings of Metalink long options (--merge-output-detai)
    Genuine Snakemake options (e.g. --ri, -k, --rerun-incomplete) pass through.
    """
    import difflib
    # option name in dash form -> set of modes that accept it
    opt_to_modes = {}
    for _m, _specs in PARAMS_HELP.items():
        for _opt in _specs:
            opt_to_modes.setdefault(_opt.replace('_', '-'), set()).add(_m)
    # mode-agnostic options defined outside PARAMS_HELP (consumed by click
    # when spelled correctly; listed here only for diagnostics)
    for _opt in ('read-type', 'skip-level', 'input-tsv', 'db-dir',
                 'working-dir', 'jobs', 'verbose', 'silent', 'dry-run',
                 'generate-config', 'configfile', 'project-name',
                 'run-project', 'taxonomy-rank', 'kraken2-taxonomy-jobs',
                 'mags', 'databases', 'kraken2-taxonomy-db-size',
                 'kraken2-function-db-size', 'delete-archives'):
        opt_to_modes.setdefault(_opt, set()).add('*')

    # Single-dash long-option typos (e.g. -merge-output-detail) are mangled by
    # click's short-option clustering before they reach snakemake_args, so they
    # must be detected on the raw command line. A single leading dash followed
    # by letters containing another dash can never be a valid short-option
    # cluster (-pn), so this is unambiguous.
    for tok in sys.argv[2:]:
        if (tok.startswith('-') and not tok.startswith('--')
                and len(tok) > 2 and tok[1].isalpha() and '-' in tok[1:]):
            name = tok[1:].split('=', 1)[0]
            hint = (" Did you mean '--{}'?".format(name)
                    if name in opt_to_modes else '')
            raise click.BadParameter(
                "'{}' looks like a long option written with a single dash; "
                "it was NOT applied.{}".format(tok, hint))

    for tok in snakemake_args:
        if not tok.startswith('--'):
            continue
        name = tok.lstrip('-').split('=', 1)[0]
        if not name or not name[0].isalpha():
            continue
        if name in opt_to_modes and mode not in opt_to_modes[name] \
                and '*' not in opt_to_modes[name]:
            raise click.BadParameter(
                "'{}' is not a Metalink option for mode '{}' (valid for: {}); "
                "it was NOT applied.".format(tok, mode,
                                             ', '.join(sorted(opt_to_modes[name]))))
        if name not in opt_to_modes:
            close = difflib.get_close_matches(name, opt_to_modes, n=1, cutoff=0.85)
            if close:
                raise click.BadParameter(
                    "unknown option '{}'. Did you mean '--{}'?".format(tok, close[0]))


def build_command(mode):
    """Generate a top-level command: metalink <mode>"""
    def cmd_func(input_tsv, working_dir, jobs, verbose, silent, dryrun,
                 db_dir, project_name, snakemake_args,
                 generate_config, configfile,
                 **kw):

        set_logger(silent=silent)
        if not silent:
            logging.info(f'Metalink {__version__}')
            logging.info(' '.join(map(shlex.quote, sys.argv)))

        # Unknown options land in snakemake_args (ignore_unknown_options=True);
        # --mags is pmr-only, so reject it explicitly for every other mode.
        if mode != 'pmr' and any(a == '--mags' or a.startswith('--mags=') for a in snakemake_args):
            raise click.BadParameter('--mags is only valid for the pmr mode.')

        _validate_passthrough_args(mode, snakemake_args)

        if generate_config and configfile:
            raise click.BadParameter(
                '--generate-config and --configfile are mutually exclusive.'
            )

        if not configfile and not input_tsv:
            raise click.BadParameter(
                'Either -i/--input-tsv or -c/--configfile must be given.'
            )

        # Parameters that determine the content of the config file.
        # They are mutually exclusive with -c/--configfile because the
        # existing config already encodes these values.
        config_source_params = {
            'input_tsv': '-i/--input-tsv',
            'db_dir': '-d/--db-dir',
            'working_dir': '-w/--working-dir',
            'read_type': '--read-type',
            'skip_level': '--skip-level',
        }

        ctx = click.get_current_context()
        provided_config_sources = {
            name for name in config_source_params
            if ctx.get_parameter_source(name) == click.core.ParameterSource.COMMANDLINE
        }
        provided_options = [config_source_params[name] for name in provided_config_sources]

        if configfile and provided_options:
            raise click.BadParameter(
                'The following parameters cannot be used with -c/--configfile: '
                f"{', '.join(provided_options)}. "
                'Either remove -c/--configfile and use these options to generate a new config, '
                'or use the config file directly.'
            )

        # 1. --generate-config，must spicify -i
        if generate_config:
            if not input_tsv:
                raise click.BadParameter(
                    '-i/--input-tsv is required when using --generate-config.'
                )
            exec_workflow_generate_config_only(
                mode, input_tsv, jobs, working_dir, db_dir, project_name,
                cli_provided=provided_config_sources, **kw
            )
            return

        # 2. --configfile or use template config
        if configfile:
            # project_name is an explicit cmd_func parameter (not in **kw);
            # thread it through only when the user actually passed it
            if ctx.get_parameter_source('project_name') == click.core.ParameterSource.COMMANDLINE:
                kw['project_name'] = project_name
            exec_workflow_with_existing_config(
                mode=mode, config_f=configfile, jobs=jobs, verbose=verbose,
                silent=silent, dryrun=dryrun,
                snakemake_args=snakemake_args, **kw
            )
        else:
            if not input_tsv:
                raise click.BadParameter('-i/--input-tsv or -c/--configfile is required.')
            exec_workflow(mode, input_tsv, working_dir, jobs,
                        verbose, silent, dryrun,
                        db_dir, project_name, snakemake_args,
                        cli_provided=provided_config_sources, **kw)

    # Apply common options in reverse so that the first item in the list appears first in help
    for opt in reversed(common_options):
        cmd_func = opt(cmd_func)
    cmd_func = add_mode_options(mode)(cmd_func)
    for _p in getattr(cmd_func, '__click_params__', []):
        if _p.name in ('delete_bam', 'keep_drep_data', 'keep_vs2_iter0',
                       'keep_assembly_intermediate', 'kraken2_taxonomy_jobs'):
            _p._section_name = 'Intermediate files & resources'

    doc = {
        'b': """
Base pipeline: assembly → mapping → binning.
""",
        'm': """
mMAG-level analysis; only evidence from mMAGs is used, several mMAG-related annotations will also be processed.
""",
        'c': """
Contig-level analysis; only evidence from contigs is used.
""",
        'r': """
Read-level analysis; only evidence from reads is used.
""",
        'mc': """
mMAG-contig-level analysis; use mMAG taxonomy to refine contig annotation.
""",
        'mcr': """
mMAG-contig-read-level analysis; use taxonomy/function of unmapped read to refine mc results.
""",
        'cr': """
Contig-read-level analysis; use taxonomy/function of unmapped read to refine c results.
""",
        'v': """
Virus analysis; vMAG identification from contigs, several virus-related annotations will also be processed.
""",
        'vc': """
vMAG-contig-level analysis; use vMAG taxonomy to refine contig annotation.
""",
        'vcr': """
vMAG-contig-read-level analysis; use taxonomy/function of unmapped read to refine vc results.
""",
        'vm': """
vMAG-mMAG-level analysis; use vMAG taxonomy to refine mMAG annotation.
""",
        'vmc': """
vMAG-mMAG-contig-level analysis; use vMAG taxonomy to refine mc annotation.
""",
        'vmcr': """
vMAG-mMAG-contig-read-level analysis; use taxonomy/function of unmapped read to refine vmc results.
""",
        'pmr': """
Project-MAG-read-level analysis: map reads to the project-level MAG reference and integrate taxonomy, function and depth into abundance matrices. --mags selects the MAG source: m (mMAGs, default; p+mMAG+r), v (vMAGs; p+vMAG+r), or vm (vMAGs+mMAGs combined into one reference with a single mapping pass; p+vMAG+mMAG+r).
""",

    }[mode]

    cmd_func = cli.command(
        mode,
        cls=OrderedCommand,
        context_settings=dict(ignore_unknown_options=True),
        short_help=doc.strip()
    )(cmd_func)
    return cmd_func

import shutil
def get_width():
    try:
        return shutil.get_terminal_size().columns
    except Exception:
        return 200
class OrderedGroup(click.Group):
    """Group that lists setup commands before analysis modes, in sections."""

    setup_commands = ('get-envs', 'download-database')

    def list_commands(self, ctx):
        names = [n for n in self.setup_commands if n in self.commands]
        names += [n for n in self.commands if n not in self.setup_commands]
        return names

    def format_commands(self, ctx, formatter):
        sections = [
            ('Setup (run once before analysis)', self.setup_commands),
            ('Analysis modes',
             [n for n in self.commands if n not in self.setup_commands]),
        ]
        for title, names in sections:
            rows = [(n, self.commands[n].get_short_help_str())
                    for n in names if n in self.commands]
            if rows:
                with formatter.section(title):
                    formatter.write_dl(rows)

@click.group(cls=OrderedGroup,
    context_settings=dict(
        help_option_names=["-h", "--help"],
        max_content_width=get_width(),
        terminal_width=get_width(),
        show_default=True,
    )
)
@click.version_option(__version__, '-v', '--version')
@click.pass_context
def cli(obj):
    """Metalink: Joint taxonomic and functional profiling of metagenomes and metatranscriptomes across resolution scales. Metalink can be run in different modes to obtain joint outputs of functional and taxonomic profiles at multiple levels."""


build_command('b')
build_command('m')
build_command('c')
build_command('r')
build_command('mc')
build_command('mcr')
build_command('cr')
build_command('v')
build_command('vc')
build_command('vcr')
build_command('vm')
build_command('vmc')
build_command('vmcr')
build_command('pmr')


def get_snakefile(filename='Snakefile'):
    """Locate a snakemake workflow file bundled with this package."""
    sf = os.path.join(os.path.dirname(os.path.abspath(__file__)), filename)
    if not os.path.exists(sf):
        sys.exit('Unable to locate the Snakemake workflow file; tried %s' % sf)
    return sf


class CommaSeparatedChoice(click.ParamType):
    """Click type that validates a comma-separated list against allowed choices."""

    name = 'comma-separated-choice'

    def __init__(self, choices):
        self.choices = tuple(choices)

    def convert(self, value, param, ctx):
        if value is None:
            return None
        items = [item.strip() for item in str(value).split(',') if item.strip()]
        for item in items:
            if item not in self.choices:
                self.fail(
                    f"'{item}' is not a valid database name. "
                    f"Choose from: {', '.join(self.choices)}",
                    param,
                    ctx,
                )
        return ','.join(items)


DOWNLOAD_DATABASE_CHOICES = [
    'kraken2_taxonomy',
    'kraken2_function',
    'cat_pack_db',
    'eggnog',
    'gtdbtk',
    'checkm',
    'checkm2',
    'checkv',
    'genomad',
    'vibrant',
    'virsorter',
    'iphop',
]

# Database -> required conda environment key(s).  Keys are resolved to actual
# conda environment names via the conda_envs mapping in the download config
# (generated from config/template.config.general.yaml defaults).
DB_REQUIRED_CONDA_ENVS = {
    'kraken2_taxonomy': ['metalink-main'],
    'kraken2_function': ['metalink-main'],
    'cat_pack_db': ['metalink-annotation'],
    'eggnog': ['metalink-annotation'],
    'genomad': ['metalink-genomad'],
    'checkv': ['metalink-annotation'],
    'vibrant': ['metalink-vibrant'],
    'iphop': ['metalink-iphop'],
    'gtdbtk': ['metalink-gtdbtk'],
    'checkm': ['metalink-binning'],
    'virsorter': ['metalink-virsorter'],
    'checkm2': ['metalink-checkm2'],
}


def _load_conda_env_map():
    """Return the {env key: conda env name} mapping from the general template."""
    template = Path(__file__).with_name('config') / 'template.config.general.yaml'
    try:
        with open(template) as f:
            data = yaml.load(f)
        mapping = (data or {}).get('conda_envs', {}) or {}
        return {str(k): str(v) for k, v in mapping.items() if v}
    except Exception:
        return {}


def _env_yaml_name(key):
    """Return the name: field of metalink/envs/<key>.yaml, or None."""
    env_yaml = Path(__file__).with_name('envs') / '{}.yaml'.format(key)
    try:
        with open(env_yaml) as f:
            for line in f:
                if line.startswith('name:'):
                    return line.split(':', 1)[1].strip()
    except OSError:
        pass
    return None


def _check_required_conda_envs(selected_dbs, env_map, warn_drift=True):
    """Verify that all conda environments required by selected databases exist."""
    required_keys = {'metalink-main'}
    for db in selected_dbs:
        required_keys.update(DB_REQUIRED_CONDA_ENVS.get(db, []))

    # key -> actual conda env name (fall back to the key itself)
    required = {key: env_map.get(key) or key for key in required_keys}

    # Catch drift between the env yaml name: (used by `mamba env create`) and
    # the conda_envs value (used by the workflow): a mismatch means the created
    # environment will never be found.  Only checked for template defaults --
    # when the user supplies a config file they own the mapping, and a
    # deliberate rename is not drift.
    if warn_drift:
        for key, name in sorted(required.items()):
            yaml_name = _env_yaml_name(key)
            if yaml_name and yaml_name != name:
                logging.warning(
                    'conda env name mismatch for %s: metalink/envs/%s.yaml defines '
                    '"%s" but the config template maps it to "%s". '
                    'Keep the two in sync.', key, key, yaml_name, name)

    try:
        result = subprocess.run(
            ['conda', 'env', 'list'],
            capture_output=True, text=True, check=True
        )
    except FileNotFoundError:
        logging.error('conda is not available in PATH; cannot verify conda environments.')
        sys.exit(1)
    except subprocess.CalledProcessError as e:
        logging.error('Failed to list conda environments: %s', e.stderr.strip())
        sys.exit(1)

    existing = set()
    for line in result.stdout.splitlines():
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        parts = line.split()
        # The active env is prefixed with '*'; skip that marker.
        if parts and parts[0] == '*':
            parts = parts[1:]
        if len(parts) >= 2:
            existing.add(parts[0])

    missing = sorted(name for name in required.values() if name not in existing)
    if missing:
        logging.error('The following required conda environments are missing:')
        for env in missing:
            logging.error('  - %s', env)
        logging.error('Please create them first, e.g.:')
        logging.error('  conda env create -f metalink/envs/<env>.yaml')
        sys.exit(1)


def _generated_download_config_path(db_dir):
    return os.path.join(os.path.abspath(db_dir), 'config', 'config.download.yaml')


def make_download_config(db_dir, databases, kraken2_taxonomy_db_size,
                         kraken2_function_db_size, delete_archives, config_f):
    """Build config.download.yaml from the download template plus CLI options."""
    config_dir = Path(__file__).with_name('config')
    download_template = config_dir / 'template.config.download.yaml'
    if not download_template.exists():
        raise FileNotFoundError('missing download template {}'.format(download_template))
    with open(download_template) as f:
        cfg = yaml.load(f)

    cfg['mode'] = 'download'
    cfg['Directory']['database']['DBDIR'] = str(Path(db_dir).resolve())

    selected = [x.strip() for x in databases.split(',') if x.strip()]
    invalid = [d for d in selected if d not in DOWNLOAD_DATABASE_CHOICES]
    if invalid:
        raise ValueError(
            'unknown database(s): {}; choose from {}'.format(
                ', '.join(invalid), ', '.join(DOWNLOAD_DATABASE_CHOICES)))
    cfg['Databases'] = selected

    cfg['Kraken2_taxonomy_db_size'] = kraken2_taxonomy_db_size
    cfg['Kraken2_function_db_size'] = kraken2_function_db_size
    cfg['Delete_archives'] = bool(delete_archives)
    # Re-assigning cfg['Databases'] above drops the template comment that
    # followed the list; attach the size-variant documentation explicitly.
    cfg.yaml_set_comment_before_after_key(
        'Kraken2_taxonomy_db_size',
        before='\nSize variants of the Metalink-built kraken2 databases:\n'
        '  Kraken2_taxonomy_db_size: kraken2_16GB / kraken2_64GB / kraken2_256GB / kraken2_raw\n'
        '  Kraken2_function_db_size: kraken2_16GB / kraken2_raw')

    # The conda_envs mapping is shared with the analysis modes; pull it from
    # the general template (comments included) so it is maintained in one
    # place and users rename environments the same way in both configs.
    general_template = config_dir / 'template.config.general.yaml'
    with open(general_template) as f:
        general_cfg = yaml.load(f)
    if general_cfg and general_cfg.get('conda_envs'):
        cfg['conda_envs'] = general_cfg['conda_envs']
        cfg.yaml_set_comment_before_after_key(
            'conda_envs',
            before='\nMapping from Metalink-internal env keys to the actual conda environment\n'
            'names. You may point a key at a differently-named environment of your own\n'
            'here; both the pre-download environment check and the workflow use them.')

    # Record the absolute path of the bundle manifest so setup.smk does not
    # have to infer it from __file__ (unreliable inside Snakemake).
    cfg['DB_bundle'] = str(config_dir / 'db_bundle.yaml')
    cfg.yaml_set_comment_before_after_key(
        'DB_bundle',
        before='Database bundle manifest shipped with this Metalink installation.')

    os.makedirs(os.path.dirname(config_f), exist_ok=True)
    with open(config_f, 'w') as f:
        f.write('# Metalink v{}\n'.format(__version__))
        f.write('# Note: when using this config with -c/--configfile, the config-source options\n')
        f.write('#       (-d, --databases, --kraken2-taxonomy-db-size, --kraken2-function-db-size)\n')
        f.write('#       can still be given on the command line to override the values here for\n')
        f.write('#       that run; this file is not modified. Edit it to change the defaults.\n')
        # Keep long database paths on a single line instead of folding them.
        _old_width = getattr(yaml, 'width', None)
        yaml.width = 4096
        yaml.dump(cfg, f)
        yaml.width = _old_width
    _compact_yaml(config_f)


@cli.command(
    'download-database',
    context_settings=dict(ignore_unknown_options=True),
    short_help='download MetaLink reference databases',
)
@click.option('-d',
    '--db-dir',
    help='directory path for databases',
    type=click.Path(dir_okay=True, writable=True, resolve_path=True),
    default=None
)
@click.option(
    '-j',
    '--jobs',
    default=multiprocessing.cpu_count(),
    type=int,
    show_default=True,
    help='number of simultaneous downloads',
)
@click.option(
    '--databases',
    default=','.join(DOWNLOAD_DATABASE_CHOICES),
    type=CommaSeparatedChoice(DOWNLOAD_DATABASE_CHOICES),
    show_default=True,
    help=(
        'comma-separated list of database names to download, using any of the '
        '[kraken2_taxonomy/kraken2_function/cat_pack_db/eggnog/gtdbtk/checkm/checkm2/checkv/genomad/vibrant/virsorter/iphop]'
    ),
)
@click.option(
    '--kraken2-taxonomy-db-size',
    default='kraken2_raw',
    type=str,
    show_default=True,
    help='size variant of the MetaLink kraken2 taxonomy database, using one of [kraken2_16GB/kraken2_64GB/kraken2_256GB/kraken2_raw]',
)
@click.option(
    '--kraken2-function-db-size',
    default='kraken2_raw',
    type=str,
    show_default=True,
    help='size variant of the MetaLink kraken2 function database, using one of [kraken2_16GB/kraken2_raw]',
)
@click.option(
    '--delete-archives', 'delete_archives', is_flag=True,
    help=
    '\b\nDelete the downloaded archives after successful extraction (flag only).\n'
    'By default the .tar.gz archives and their .md5 files are kept in DBDIR/downloads/\n'
    'so databases can be re-extracted or copied to another machine without\n'
    're-downloading. With this flag they are removed to save disk space.'
)
@click.option(
    '--generate-config', 'generate_config', is_flag=True,
    help=
    '\b\nOnly generate the config file (flag only).\n'
    'Writes <db_dir>/config/config.download.yaml merged from the package template\n'
    'and the CLI options, then exits without running the workflow.\n'
    'Mutually exclusive with -c/--configfile.'
)
@click.option(
    '-c',
    '--configfile',
    metavar='PATH',
    type=click.Path(exists=True, dir_okay=False, resolve_path=True),
    help=(
        '\b\nrun from an existing config.download.yaml (see --generate-config); config-source options\n'
        '(-d, --databases, --kraken2-taxonomy-db-size, --kraken2-function-db-size)\n'
        'given on the command line override the file values for this run without modifying the file.'
    ),
)
@click.option('--verbose', is_flag=True,
              help=
              '\b\nEnable verbose output during workflow execution (flag only).\n'
              'When enabled, detailed execution logs are printed, including rule activation,\n'
              'job submission, and resource usage. Useful for debugging workflow failures.\n'
              'When disabled (default), the workflow runs in quiet mode.'
              )
@click.option('--silent', is_flag=True,
              help=
              '\b\nSilence all Metalink messages (flag only).\n'
              'When enabled, Metalink initialization and execution logs are suppressed.'
              )
@click.option('-n', '--dry-run', 'dryrun', is_flag=True,
              help=
              '\b\nPerform a dry run without executing any jobs (flag only).\n'
              'Useful for validating workflow logic and checking input/output expectations.'
              )
@click.argument('snakemake_args', nargs=-1, type=click.UNPROCESSED)
def download_database(db_dir, jobs, databases,
                      kraken2_taxonomy_db_size, kraken2_function_db_size,
                      delete_archives, generate_config, configfile,
                      verbose, silent, dryrun, snakemake_args):
    '''Download MetaLink reference databases.

    Executes a Snakemake workflow to download the MetaLink-built kraken2 and
    CAT database archives and third-party tool databases, extract them, and
    verify expected marker files. Required conda environments are checked
    before any download starts; missing environments are reported and the
    command stops.

    Like the analysis modes, there are three ways to run:

    \b
      1. CLI options:  metalink download-database -d DBDIR [--databases ...]
         (a config is written to DBDIR/config/config.download.yaml)
      2. Config only:  metalink download-database --generate-config -d DBDIR
      3. From config:  metalink download-database -c CONFIG [--verbose ...]
    '''
    set_logger(silent=silent)

    _validate_passthrough_args('download-database', snakemake_args)

    if generate_config and configfile:
        raise click.BadParameter(
            '--generate-config and --configfile are mutually exclusive.'
        )
    if not configfile and not db_dir:
        raise click.BadParameter(
            'Either -d/--db-dir or -c/--configfile must be given.'
        )

    # Config-source parameters explicitly given on the command line override
    # the values in a -c/--configfile for this run (merged into a temporary
    # runtime config; the file itself is left unchanged).
    config_source_params = {
        'db_dir': '-d/--db-dir',
        'databases': '--databases',
        'kraken2_taxonomy_db_size': '--kraken2-taxonomy-db-size',
        'kraken2_function_db_size': '--kraken2-function-db-size',
        'delete_archives': '--delete-archives',
    }
    ctx = click.get_current_context()
    provided_config_sources = {
        name for name in config_source_params
        if ctx.get_parameter_source(name) == click.core.ParameterSource.COMMANDLINE
    }

    # 1. --generate-config: write the config and stop.
    if generate_config:
        config_f = _generated_download_config_path(db_dir)
        make_download_config(db_dir, databases, kraken2_taxonomy_db_size,
                             kraken2_function_db_size, delete_archives, config_f)
        logging.info('Config file written to %s', config_f)
        return

    # 2. Resolve the config to run from: an existing one (-c), or a freshly
    #    generated one at <db_dir>/config/config.download.yaml.
    if configfile:
        config_f = configfile
        cfg = load_configfile(config_f)
        if not silent:
            logging.info('Using config file: %s', config_f)
        if cfg.get('mode') != 'download':
            raise ValueError(
                "Config mode '{}' != 'download'; generate a download config with "
                "'metalink download-database --generate-config -d <DBDIR>'."
                .format(cfg.get('mode'))
            )
        # CLI config-source options win over the file values; the overrides go
        # into a temporary runtime config so the user's file stays untouched.
        overridden = False
        if 'db_dir' in provided_config_sources:
            cfg['Directory']['database']['DBDIR'] = str(Path(db_dir).resolve())
            overridden = True
            if not silent:
                logging.info("Setting config Directory.database.DBDIR to %r (from CLI -d/--db-dir)",
                             cfg['Directory']['database']['DBDIR'])
        if 'databases' in provided_config_sources:
            cfg['Databases'] = [x.strip() for x in databases.split(',') if x.strip()]
            overridden = True
            if not silent:
                logging.info("Setting config Databases to %r (from CLI --databases)",
                             cfg['Databases'])
        for _name, _key, _val in (
                ('kraken2_taxonomy_db_size', 'Kraken2_taxonomy_db_size', kraken2_taxonomy_db_size),
                ('kraken2_function_db_size', 'Kraken2_function_db_size', kraken2_function_db_size)):
            if _name in provided_config_sources:
                cfg[_key] = _val
                overridden = True
                if not silent:
                    logging.info("Setting config %s to %r (from CLI --%s)",
                                 _key, _val, _name.replace('_', '-'))
        if 'delete_archives' in provided_config_sources:
            cfg['Delete_archives'] = delete_archives
            overridden = True
            if not silent:
                logging.info("Setting config Delete_archives to %r (from CLI --delete-archives)",
                             delete_archives)
        if overridden:
            tmp_dir = os.path.dirname(os.path.abspath(config_f))
            fd, tmp_path = tempfile.mkstemp(
                suffix='.yaml',
                prefix='config.download.runtime.',
                dir=tmp_dir
            )
            with os.fdopen(fd, 'w') as f:
                yaml.dump(cfg, f)
            config_f = tmp_path
            if not silent:
                logging.info('CLI overrides merged to temporary config: %s', tmp_path)
    else:
        db_dir = os.path.abspath(db_dir)
        os.makedirs(db_dir, exist_ok=True)
        config_f = _generated_download_config_path(db_dir)
        make_download_config(db_dir, databases, kraken2_taxonomy_db_size,
                             kraken2_function_db_size, delete_archives, config_f)
        if not silent:
            logging.info('Config file written to %s', config_f)
            logging.info('Using config file: %s', config_f)
        cfg = load_configfile(config_f)

    db_dir = str(cfg.get('Directory', {}).get('database', {}).get('DBDIR') or '')
    if not db_dir:
        raise ValueError(
            'Config {} does not set Directory.database.DBDIR.'.format(config_f))
    if not os.path.isabs(db_dir):
        db_dir = os.path.abspath(os.path.join(os.path.dirname(config_f), db_dir))
    os.makedirs(db_dir, exist_ok=True)

    selected_dbs = cfg.get('Databases') or list(DOWNLOAD_DATABASE_CHOICES)
    if isinstance(selected_dbs, str):
        selected_dbs = [x.strip() for x in selected_dbs.split(',') if x.strip()]
    invalid = [d for d in selected_dbs if d not in DOWNLOAD_DATABASE_CHOICES]
    if invalid:
        raise ValueError(
            'unknown database(s) in {}: {}; choose from {}'.format(
                config_f, ', '.join(invalid), ', '.join(DOWNLOAD_DATABASE_CHOICES)))

    # Conda environment names come from the config's conda_envs section
    # (copied from the general template at generation time); fall back to the
    # template for hand-written configs that lack the section.
    env_map = {str(k): str(v) for k, v in (cfg.get('conda_envs') or {}).items() if v}
    if not env_map:
        env_map = _load_conda_env_map()
    _check_required_conda_envs(selected_dbs, env_map,
                               warn_drift=(configfile is None))

    quiet_flag = '' if (verbose and not silent) else '--quiet'
    dryrun_flag = '--dry-run' if dryrun else ''
    env_prefix = 'METALINK_SILENT=1 ' if silent else ''

    cmd = (
        '{env_prefix}snakemake --snakefile {snakefile} '
        '--directory {db_dir} '
        '--configfile {config_f} '
        '--jobs {jobs} --rerun-incomplete --latency-wait 600 '
        '--nolock --use-conda '
        '{dryrun} {quiet} {args}'
    )
    cmd_str = cmd.format(
        env_prefix=env_prefix,
        snakefile=get_snakefile('rules/setup.smk'),
        db_dir=shlex.quote(db_dir),
        config_f=shlex.quote(config_f),
        jobs=jobs,
        dryrun=dryrun_flag,
        quiet=quiet_flag,
        args=' '.join(shlex.quote(a) for a in snakemake_args),
    )

    if not silent:
        logging.info('Setting up MetaLink databases; this only needs to be done once.')
        logging.info('Executing: %s' % cmd_str)
    try:
        subprocess.run(cmd_str, check=True, shell=True)
    except subprocess.CalledProcessError:
        sys.exit(1)


@cli.command(
    'get-envs',
    short_help='locate or export the bundled conda environment YAML files',
)
@click.option('-o',
    '--output-dir',
    help='copy the environment YAML files into this directory instead of just printing their location',
    type=click.Path(dir_okay=True, writable=True, resolve_path=True),
    default=None
)
def get_envs(output_dir):
    '''Locate the conda environment YAML files bundled with this package.

    When Metalink is installed from PyPI or Bioconda, the YAML files live
    inside the installed package rather than in a source checkout, so the
    relative paths shown in the documentation (e.g. metalink/envs/metalink-main.yaml)
    do not exist on disk. Without options this prints the absolute path of the
    bundled envs directory; with -o/--output-dir the YAML files are copied to
    a directory of your choice.

    \b
    Examples:
      metalink get-envs                       # print the bundled envs directory
      mamba env create -f "$(metalink get-envs)/metalink-main.yaml"
      metalink get-envs -o ./metalink_envs    # copy the YAML files out
      mamba env create -f ./metalink_envs/metalink-main.yaml
    '''
    set_logger(silent=False)
    envs_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'envs')
    if not os.path.isdir(envs_dir):
        sys.exit('Unable to locate the bundled envs directory; tried %s' % envs_dir)

    if output_dir is None:
        click.echo(envs_dir)
        return

    os.makedirs(output_dir, exist_ok=True)
    copied = []
    for name in sorted(os.listdir(envs_dir)):
        if name.endswith('.yaml'):
            shutil.copy2(os.path.join(envs_dir, name),
                         os.path.join(output_dir, name))
            copied.append(name)
    click.echo(output_dir)
    logging.info('Copied %d environment YAML file(s) to %s', len(copied), output_dir)
    for name in copied:
        logging.info('  mamba env create -f %s', os.path.join(output_dir, name))


if __name__ == '__main__':
    cli()
