#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
import os
import sys
import argparse


def _fmt_num(x):
    """Format a number for file/directory names; drop trailing .0 from floats."""
    if isinstance(x, float) and x.is_integer():
        return str(int(x))
    return str(x)


def build_refinement_cmd(bin_dirs, out_path, threads, memory_gb,
                         min_completeness, max_contamination,
                         dir_script_py, extra_args=None):
    """
    Build metawrap bin_refinement command based on available binning results.
    Compatible with Python 2.7+ and Python 3.x
    """
    flag_map = {'metabat2': '-A', 'maxbin2': '-B', 'concoct': '-C'}
    available = []

    for tool, path in bin_dirs.items():
        if path and os.path.isdir(path):
            available.append((flag_map[tool], path))

    n_tools = len(available)

    # 0 tools: create placeholder so Snakemake output target exists
    if n_tools == 0:
        placeholder = os.path.join(out_path, 'this_is_not_the_bin_you_are_looking_for')
        return "mkdir -p '{0}'".format(placeholder)

    # 1 tool: copy bins directly, no refinement needed
    if n_tools == 1:
        flag, src_dir = available[0]
        out_bins = os.path.join(
            out_path,
            'metawrap_{0}_{1}_bins'.format(_fmt_num(min_completeness), _fmt_num(max_contamination))
        )
        return (
            "mkdir -p '{0}' && "
            "cp '{1}'/bin.* '{0}'/ 2>/dev/null || true".format(out_bins, src_dir)
        )

    # 2 or 3 tools: run metawrap refinement
    cmd_parts = ['bash', os.path.join(dir_script_py, 'bin_refinement.sh')]
    cmd_parts += ['-o', out_path]
    cmd_parts += ['-t', threads]
    cmd_parts += ['-m', str(memory_gb)]

    for flag, path in available:
        cmd_parts += [flag, path]

    cmd_parts += ['-c', _fmt_num(min_completeness)]
    cmd_parts += ['-x', _fmt_num(max_contamination)]

    # extra_args: explicit split on whitespace, filter empty strings
    if extra_args:
        extra_list = [arg for arg in extra_args.split() if arg]
        cmd_parts += extra_list

    return ' '.join("'{0}'".format(p) if ' ' in p else p for p in cmd_parts)


def main():
    parser = argparse.ArgumentParser(
        description='MetaWRAP bin refinement wrapper with flexible tool combinations'
    )
    parser.add_argument('metabat2',
                        help='Path to metabat2 bins directory, or "None"')
    parser.add_argument('maxbin2',
                        help='Path to maxbin2 bins directory, or "None"')
    parser.add_argument('concoct',
                        help='Path to concoct bins directory, or "None"')
    parser.add_argument('out_path',
                        help='Output directory path')
    parser.add_argument('threads',
                        help='Number of threads')
    parser.add_argument('dir_script_py',
                        help='Path to script directory containing bin_refinement.sh')
    parser.add_argument('--memory-gb', type=float, default=64,
                        help='Max memory in GB (default: 64)')
    parser.add_argument('--min-completeness', type=float, default=50,
                        help='Minimum completeness %% (default: 50)')
    parser.add_argument('--max-contamination', type=float, default=10,
                        help='Maximum contamination %% (default: 10)')
    parser.add_argument('--extra-args', dest='extra_args', default='',
                        help='Extra arguments forwarded to metaWRAP bin_refinement '
                             '(space-separated, e.g. "--remove-ambiguous --quick")')

    args = parser.parse_args()

    # Explicit None handling: convert string "None" to Python None
    bin_dirs = {
        'metabat2': args.metabat2 if args.metabat2 != "None" else None,
        'maxbin2':  args.maxbin2  if args.maxbin2  != "None" else None,
        'concoct':  args.concoct  if args.concoct  != "None" else None,
    }

    cmd = build_refinement_cmd(
        bin_dirs=bin_dirs,
        out_path=args.out_path,
        threads=args.threads,
        memory_gb=args.memory_gb,
        min_completeness=args.min_completeness,
        max_contamination=args.max_contamination,
        dir_script_py=args.dir_script_py,
        extra_args=args.extra_args
    )

    print("[CMD] {0}".format(cmd), file=sys.stderr)
    ret = os.system(cmd)
    if ret != 0:
        print("[ERROR] Command failed with exit code {0}".format(ret), file=sys.stderr)
        sys.exit(1)
    else:
        print("[OK] Command completed successfully", file=sys.stderr)


if __name__ == '__main__':
    main()
