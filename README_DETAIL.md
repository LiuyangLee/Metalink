# Metalink — Detailed Documentation

[![License](https://img.shields.io/badge/license-GPL%20v3-blue.svg)](LICENSE)

> **Metalink**: Joint taxonomic and functional profiling of metagenomes and metatranscriptomes across resolution scales.

🔙 **Back to the [project overview](README.md).**

## Overview

Metalink implements an integrated strategy that functionally and taxonomically profiles metagenomic and metatranscriptomic sequencing. It provides 14 processing modes organized by sequence type inclusion and viral analysis scope. All modes except read-only (`r`) depend on the base workflow (`b`) for contig assembly, and binning-related modes (`m`, `mc`, `mcr`, `pmr --mags m`/`vm`, `vm`, `vmc`, `vmcr`) additionally run the base mMAG binning workflow. All virus-related modes additionally depend on the virus mode (`v`) for viral prediction and viral binning workflows.

1. **Comprehensive & Scalable Workflow**:  14 analysis modes (b to vmcr) covering assembly, binning, and joint taxonomic and functional profiling. Supports 2nd-generation short reads (single-end / paired-end), 3rd-generation long reads (PacBio / Nanopore), and metatranscriptomes (read and assembly-level profiling). Built on Snakemake with robust dependency tracking, dry-run (`-n`), resume (`--ri`/`--rerun-incomplete`), and cluster execution.
2. **Accurate Taxonomic & Functional Profiling**: Fast k-mer-based classification (kraken2) with a unified taxonomy integrating GTDB (prokaryotes), ICTV (viruses), and NCBI (microeukaryotes). The LCFH (Last Common Function-based functional Hierarchy) reorganizes ~424 million KO annotations into a seven-level tree (lv1–lv7), enabling LCA-like functional assignment.
3. **Joint Taxonomic-Functional Output**: Per-sequence co-annotation directly answers which species contribute to which functions, moving analysis beyond cataloging species and functions toward mechanistic inference of species–function–environment triadic associations.

---

## Table of Contents

- [Overview](#overview)
- [Installation & Setup](#installation--setup)
  - [Metalink CLI Installation](#metalink-cli-installation)
    - [Option 1: PyPI Installation](#option-1-pypi-installation-recommended-for-most-users)
    - [Option 2: Manual Installation](#option-2-manual-installation-from-source)
    - [Option 3: Bioconda Installation](#option-3-bioconda-installation)
    - [Verify installation](#verify-installation)
  - [Workflow Environment Setup](#workflow-environment-setup)
    - [Environment overview](#environment-overview)
    - [Minimal setup (read-only mode `r`)](#minimal-setup-read-only-mode-r)
    - [Full setup (all analysis modes)](#full-setup-all-analysis-modes)
  - [Database Setup](#database-setup)
    - [Default path (recommended)](#default-path-recommended)
    - [Custom paths](#custom-paths)
    - [Default path layout](#default-path-layout)
    - [Required databases](#required-databases)
    - [Quick setup](#quick-setup)
    - [Which databases do I need?](#which-databases-do-i-need)
- [Input File Formats](#input-file-formats)
  - [Paired-end Format](#paired-end-format)
  - [Single-end Format (short single-end or long read)](#single-end-format-short-single-end-or-long-read)
  - [Advanced Format (pre-computed contigs and/or MAGs)](#advanced-format-pre-computed-contigs-andor-mags)
- [Analysis Modes](#analysis-modes)
  - [Mode Hierarchy](#mode-hierarchy)
  - [Mode Details](#mode-details)
    - [`b` — Base Analysis](#b--base-analysis)
    - [`c` — Contig Analysis](#c--contig-analysis)
    - [`cr` — Contig-Read Analysis](#cr--contig-read-analysis)
    - [`r` — Read Analysis](#r--read-analysis)
    - [`m` — mMAG Analysis](#m--mmag-analysis)
    - [`mc` / `mcr` — Integrated MAG Modes](#mc--mcr--integrated-mag-modes)
    - [`v` — vMAG Analysis](#v--vmag-analysis)
    - [`vm` — vMAG + mMAG Integrated Analysis](#vm--vmag--mmag-integrated-analysis)
    - [`vc` / `vcr` / `vmc` / `vmcr` — Integrated vMAG Modes](#vc--vcr--vmc--vmcr--integrated-vmag-modes)
    - [`pmr` — Project-Level MAG + Read Analysis](#pmr--project-level-mag--read-analysis)
- [Usage](#usage)
  - [Choosing a Mode](#choosing-a-mode)
    - [Progressive Workflow (Optional)](#progressive-workflow-optional)
  - [Mode-Specific Differences and Special Cases](#mode-specific-differences-and-special-cases)
    - [`r`: Read-Only Direct Profiling](#r-read-only-direct-profiling)
    - [`v` vs `vm`: iPHoP Host Prediction Strategies](#v-vs-vm-iphop-host-prediction-strategies)
    - [`pmr`: Project-Level Pooled MAGs + Reads](#pmr-project-level-pooled-mags--reads)
  - [Skip-Level: Reusing Pre-Computed Results](#skip-level-reusing-pre-computed-results)
    - [Advanced skip-level use cases](#advanced-skip-level-use-cases)
  - [Read Type](#read-type)
    - [Read Type Behavior](#read-type-behavior)
    - [Paired-End Data](#paired-end-data)
    - [Single-End Data](#single-end-data)
    - [Long-Read Data](#long-read-data)
    - [Metatranscriptomic Data](#metatranscriptomic-data)
  - [Configuration & Parameters](#configuration--parameters)
    - [How the Config Is Built](#how-the-config-is-built)
    - [Reusing an Existing Config (`-c/--configfile`)](#reusing-an-existing-config--c--configfile)
    - [Key Global Parameters](#key-global-parameters)
    - [Parameter Help](#parameter-help)
    - [Customizing Tool Parameters with --extra-args](#customizing-tool-parameters-with--extra-args)
  - [Workflow Management](#workflow-management)
    - [Common Snakemake flags](#common-snakemake-flags)
    - [Generating a DAG](#generating-a-dag)
    - [Resuming, unlocking, and cleaning up](#resuming-unlocking-and-cleaning-up)
    - [Deleting large intermediate files](#deleting-large-intermediate-files)
    - [Limiting kraken2 taxonomy concurrency](#limiting-kraken2-taxonomy-concurrency)
  - [How to Tell if a Run Completed Successfully](#how-to-tell-if-a-run-completed-successfully)
    - [1. Exit Code](#1-exit-code)
    - [2. Target File Presence](#2-target-file-presence)
    - [3. Snakemake Summary](#3-snakemake-summary)
    - [4. Log Inspection](#4-log-inspection)
- [Output File Formats](#output-file-formats)
  - [Working Directory Structure (Generated at Runtime)](#working-directory-structure-generated-at-runtime)
  - [Results by Mode and Scope](#results-by-mode-and-scope)
  - [Sample-level Outputs](#sample-level-outputs)
    - [Base mode (`b`) sample layout](#base-mode-b-sample-layout)
    - [Contig mode (`c`) sample layout (extends `b`)](#contig-mode-c-sample-layout-extends-b)
    - [Read mode (`r`) sample layout](#read-mode-r-sample-layout)
    - [MAG mode (`m`) sample layout](#mag-mode-m-sample-layout)
    - [Sample-level `05_visualization/` files](#sample-level-05_visualization-files)
    - [Virus mode (`v`) sample layout](#virus-mode-v-sample-layout)
    - [Integrated modes: merge-level fallback and the `mode` flag](#integrated-modes-merge-level-fallback-and-the-mode-flag)
  - [Project-level Outputs (when `--run-project`)](#project-level-outputs-when--run-project)
    - [Base mode (`b`) project summaries](#base-mode-b-project-summaries)
    - [MAG modes (`m`, `vm`) project summaries](#mag-modes-m-vm-project-summaries)
    - [Virus modes (`v` / `vm`) project summaries](#virus-modes-v--vm-project-summaries)
    - [Taxonomy and Function Outputs](#taxonomy-and-function-outputs)
      - [Long-format table columns](#long-format-table-columns)
      - [ReadCount interpretation for 2nd-generation vs. 3rd-generation reads](#readcount-interpretation-for-2nd-generation-vs-3rd-generation-reads)
      - [Wide-format tables](#wide-format-tables)
      - [CAMI-format taxonomic profile](#cami-format-taxonomic-profile)
  - [Joint Taxonomy–Function Outputs](#joint-taxonomyfunction-outputs)
    - [Merge taxonomy and function tables](#merge-taxonomy-and-function-tables)
      - [Columns of `detail_for_each_seq.tsv.gz`](#columns-of-detail_for_each_seqtsvgz)
      - [Key columns of the other merge tables](#key-columns-of-the-other-merge-tables)
    - [Metabolic pathway/module abundance tables (mclink outputs)](#metabolic-pathwaymodule-abundance-tables-mclink-outputs)
  - [Other Summary Tables (mMAGs and vMAGs)](#other-summary-tables-mmags-and-vmags)
    - [`mMAG_summary_info.tsv` — mMAG summary](#mmag_summary_infotsv--mmag-summary)
    - [`vMAG_summary_info.tsv` — vMAG summary](#vmag_summary_infotsv--vmag-summary)
    - [Cross-sample species abundance matrices](#cross-sample-species-abundance-matrices)
  - [Visualization Outputs](#visualization-outputs)
- [Optimizations & Patches to Third-Party Tools](#optimizations--patches-to-third-party-tools)
  - [MetaWRAP Binning Module](#metawrap-binning-module)
    - [1. FASTQ.GZ File Recognition](#1-fastqgz-file-recognition)
    - [2. CheckM / dRep Result Consistency & Avoidance of Redundant CheckM Runs](#2-checkm--drep-result-consistency--avoidance-of-redundant-checkm-runs)
- [FAQ & Troubleshooting](#faq--troubleshooting)
  - [Q: How do I choose the kraken2 confidence threshold for taxonomic classification?](#q-how-do-i-choose-the-kraken2-confidence-threshold-for-taxonomic-classification)
  - [Q: How should I interpret abundance and relative abundance — are they cell counts or molecule counts?](#q-how-should-i-interpret-abundance-and-relative-abundance--are-they-cell-counts-or-molecule-counts)
  - [Q: Can I run Metalink on a single sample?](#q-can-i-run-metalink-on-a-single-sample)
  - [Q: Can I mix paired-end and single-end samples in one run?](#q-can-i-mix-paired-end-and-single-end-samples-in-one-run)
  - [Q: How much memory does Metalink need?](#q-how-much-memory-does-metalink-need)
  - [Q: Can I use pre-trimmed reads, pre-assembled contigs, or pre-binned MAGs?](#q-can-i-use-pre-trimmed-reads-pre-assembled-contigs-or-pre-binned-mags)
  - [Q: What if my run fails on the cluster?](#q-what-if-my-run-fails-on-the-cluster)
  - [Q: What read types and assemblers are compatible?](#q-what-read-types-and-assemblers-are-compatible)
  - [Q: How often is the Metalink database updated?](#q-how-often-is-the-metalink-database-updated)
  - [Q: Where can I get help?](#q-where-can-i-get-help)
- [Citation](#citation)
  - [Workflow Engine](#workflow-engine)
  - [Read Preprocessing & Assembly](#read-preprocessing--assembly)
  - [Read Mapping](#read-mapping)
  - [MAG Binning & Quality Assessment](#mag-binning--quality-assessment)
  - [Taxonomic Classification](#taxonomic-classification)
  - [Functional Annotation](#functional-annotation)
  - [Viral Analysis](#viral-analysis)
  - [MAG Dereplication & Coverage](#mag-dereplication--coverage)
  - [ORF Prediction](#orf-prediction)
  - [Metalink Custom Tools](#metalink-custom-tools)
- [License](#license)

---

## Installation & Setup

Setting up Metalink involves three steps: installing the Metalink CLI, creating the workflow tool environments, and downloading the databases.

### Metalink CLI Installation

Metalink can be installed via three methods. The core Metalink CLI (Python package + Snakemake wrapper) is separate from the workflow-specific tool environments, which must be created independently.

---

#### Option 1: PyPI Installation

The fastest way to get the Metalink CLI:

```bash
# 1. Create a Python environment
conda create -n metalink-0.3.7
conda activate metalink-0.3.7
mamba install -c conda-forge -c bioconda python=3.11 snakemake=7.32.4 click ruamel.yaml graphviz

# 2. Install Metalink from PyPI
pip install metalink
```

> [!NOTE]
> The PyPI package contains **only** the Metalink CLI and its core Python dependencies (`snakemake>=7.32.4`, `click`, `ruamel.yaml`).
> The Snakemake workflow-specific tool environments (e.g., `metalink-main-0.3` with `kraken2`, `mclink`) are **not** bundled in the PyPI package and must be created separately via the provided YAML files. See [Workflow Environment Setup](#workflow-environment-setup) below.

---

#### Option 2: Manual Installation

Best for developers or users who need the latest unreleased features.

Metalink is compatible with a range of Snakemake versions. Tested combinations:

  | Snakemake | Python |
  |-----------|--------|
  | 7.32.4 | 3.11 |
  | 8.24.1 | 3.12 |
  | 9.25.2 | 3.14 |

```bash
# 1. Clone the repository
git clone https://github.com/LiuyangLee/Metalink.git
cd Metalink

# 2. Create a Python environment for the Metalink CLI
conda create -n metalink-0.3.7
conda activate metalink-0.3.7
mamba install -c conda-forge -c bioconda python=3.11 snakemake=7.32.4 click ruamel.yaml graphviz

# 3. Install Metalink in editable mode
pip install -e .
```

---

#### Option 3: Bioconda Installation

```bash
conda create -n metalink-0.3.7
conda activate metalink-0.3.7
mamba install -c bioconda metalink
```

> [!NOTE]
> Bioconda availability is pending approval. Until the package is accepted into the Bioconda channel, use Option 1 (PyPI) or Option 2 (manual).


#### Verify installation

```bash
conda activate metalink-0.3.7
metalink --help
metalink b --help
```

---

### Workflow Environment Setup

Metalink is a Snakemake workflow that orchestrates dozens of third-party bioinformatics tools. To avoid dependency conflicts, these tools are partitioned into **separate conda environments**. After installing the Metalink CLI (via any of the three methods above), you **must** create the workflow environments before running any analysis.

The environment YAML files are bundled inside the installed `metalink` package. Export them once with the built-in `get-envs` command (this works for all three installation methods; with a source checkout the same files also live at `metalink/envs/`):

```bash
metalink get-envs -o ./metalink_envs
```

(`metalink get-envs` without `-o` just prints the bundled envs directory path, e.g. for `mamba env create -f "$(metalink get-envs)/metalink-main.yaml"`.)

#### Environment overview

| Environment | Purpose | Required for modes |
|-------------|---------|-------------------|
| `metalink-main-0.3` | QC (`fastp`/`fastplong`), assembly (MEGAHIT/SPAdes/Flye), mapping (bwa/bwa-mem2/minibwa/strobealign), coverage (CoverM), read-level taxonomy/function (kraken2/bracken/FragGeneScanRs/prodigal), R + `mclink` reporting | **All modes** (read-only `r` needs only this, unless `reads_function_classify_method` is set to `emapper`) |
| `metalink-binning-0.3` | metaWRAP binning (metabat2/MaxBin2/CONCOCT), CheckM | All binning-related modes: `b`, `m`, `mc`, `mcr`, `vm`, `vmc`, `vmcr`, `pmr --mags m`/`vm` |
| `metalink-gtdbtk-2.4.1` | GTDB-Tk taxonomy classification | `m`, `mc`, `mcr`, `vm`, `vmc`, `vmcr`, `pmr --mags m`/`vm` when `MAGs_taxonomy_classify_method = gtdbtk` (default; not needed with `BAT`); also `v`/`vm` when custom-host iPHoP (`07_iPHoP.custom_MAGs_dir`) is configured |
| `metalink-checkm2-1.1.0` | CheckM2 MAG quality assessment | `m`, `mc`, `mcr`, `vm`, `vmc`, `vmcr`, `pmr --mags m`/`vm` |
| `metalink-annotation-0.3` | eggNOG-mapper, dRep, CheckV, CAT/BAT | `m`, `mc`, `mcr`, `pmr`, `vm`, `vmc`, `vmcr` (MAG dRep always runs; MAG function defaults to `emapper`); `v`, `vc`, `vcr`, and `pmr --mags v`/`vm` (CheckV and viral CAT_pack/BAT always run); `c`, `cr` only with `CAT_pack` taxonomy or `emapper` function; `r` only with `emapper` function |
| `metalink-genomad-1.12.0` | geNomad viral identification | `v`, `vc`, `vcr`, `vm`, `vmc`, `vmcr`, and `pmr --mags v`/`vm` with `identify_method = genomad` or `genomad_vs_vb` (default); also used by custom-host iPHoP in `v`/`vm` |
| `metalink-vclust-1.3.1` | vClust viral clustering | `v`, `vc`, `vcr`, `vm`, `vmc`, `vmcr`, and `pmr --mags v`/`vm` |
| `metalink-vibrant-1.2.1` | VIBRANT viral identification/annotation | `v`, `vc`, `vcr`, `vm`, `vmc`, `vmcr`, and `pmr --mags v`/`vm` (used under every `identify_method`) |
| `metalink-virsorter-2.2.4` | VirSorter2 viral identification | `v`, `vc`, `vcr`, `vm`, `vmc`, `vmcr`, and `pmr --mags v`/`vm` with `identify_method = vs_vb` or `genomad_vs_vb` (default) |
| `metalink-vrhyme-1.1.0` | vRhyme viral binning | `v`, `vc`, `vcr`, `vm`, `vmc`, `vmcr`, and `pmr --mags v`/`vm` |
| `metalink-iphop-1.4.1` | iPHoP host prediction | All virus modes: `v`, `vc`, `vcr`, `vm`, `vmc`, `vmcr`, and `pmr --mags v`/`vm` (sample-level iPHoP runs in every v-containing mode; the project-level `custom_MAGs_dir` host strategy is available in `v`, `vm`, `pmr --mags v`/`vm`, with `auto_drep99`/`auto_drep95` only in `vm` and `pmr --mags vm`) |

> [!IMPORTANT]
> Keep the environment **names** exactly as defined in the YAML files (e.g. `metalink-main-0.3`, `metalink-gtdbtk-2.4.1`, ...). These names are referenced by the workflow rules and scripts via the `conda_envs` mapping in the config. If you use custom names, you must manually edit the generated **config.\<mode\>.yaml** file (under `out_dir/config/`) to update each tool's entry under the `conda_envs` section.
> Naming convention: multi-tool bundle environments carry the Metalink major.minor version (e.g. `metalink-main-0.3`; the name only changes on big releases), while single-tool environments carry the tool version (e.g. `metalink-gtdbtk-2.4.1`). This lets multiple Metalink versions coexist on the same machine without conda environment name clashes, while still sharing environments whose contents have not changed.

#### Minimal setup (read-only mode `r`)

If you only need read-level taxonomic and functional profiling (`metalink r`),
create only the `metalink-main-0.3` environment:

```bash
mamba env create -f ./metalink_envs/metalink-main.yaml
```

This covers the default `r` configuration (`reads_function_classify_method = kraken2`).
If you set `reads_function_classify_method = emapper`, the `metalink-annotation-0.3` environment is additionally required.


> [!IMPORTANT]
> Bioconda availability of `r-mclink` is pending approval (20260825). Until the package is accepted into the Bioconda channel, install it manually after creating the environment:
> ```bash
> conda activate metalink-main-0.3
> R
> install.packages("mclink")
> ```

#### Full setup (all analysis modes)

For assembly-based, MAG-level, or virus-level analysis, create all environments:

```bash
# Base environment (required for all modes)
mamba env create -f ./metalink_envs/metalink-main.yaml

# binning
mamba env create -f ./metalink_envs/metalink-binning.yaml

# MAG taxonomy & quality
mamba env create -f ./metalink_envs/metalink-gtdbtk.yaml
mamba env create -f ./metalink_envs/metalink-checkm2.yaml
mamba env create -f ./metalink_envs/metalink-annotation.yaml

# Viral analysis
mamba env create -f ./metalink_envs/metalink-genomad.yaml
mamba env create -f ./metalink_envs/metalink-vclust.yaml
mamba env create -f ./metalink_envs/metalink-vibrant.yaml
mamba env create -f ./metalink_envs/metalink-virsorter.yaml
mamba env create -f ./metalink_envs/metalink-vrhyme.yaml
mamba env create -f ./metalink_envs/metalink-iphop.yaml
```

YAML definitions for all environments are bundled in the `metalink` package (exported above with `metalink get-envs`); in a source checkout they are located directly at `metalink/envs/`.

### Database Setup

By default, Metalink requires two categories of databases under a common **DBDIR**:

1. **Metalink-built databases** — provided by the Metalink authors (Zenodo or institutional mirror):
   - Metalink kraken2 taxonomy DB
   - Metalink kraken2 function DB (LCFH)
   - Metalink-customized CAT_pack DB
2. **Third-party tool databases** — downloaded from official sources (optional):
   - eggNOG-mapper, GTDB-Tk, CheckM2, geNomad, VirSorter2, VIBRANT, CheckV, iPHoP

> [!NOTE]
> **Unified taxonomy ranks**: when the Metalink-built databases merge the three source taxonomies (GTDB for prokaryotes, ICTV for viruses, NCBI for microeukaryotes), only the seven canonical ranks are kept — Domain, Phylum, Class, Order, Family, Genus, Species (`D/P/C/O/F/G/S`). Here Domain (`D`) covers Bacteria, Archaea, Eukaryota, and Viruses as the top level. Intermediate ranks present in NCBI and ICTV (e.g. subclade) are omitted so that the three taxonomic systems remain mutually compatible within a single hierarchy.

#### Default path (recommended)

By default, all database paths in `out_dir/config/config.<mode>.yaml` point to a unified **DBDIR**. To set this up, run:

```bash
metalink download-database -d db_dir
```

> [!IMPORTANT]
> `metalink download-database` verifies the required conda environments **before** downloading anything: the `metalink-main-0.3` environment (plus the environment of each selected database's tool) must already exist, otherwise the command exits with an error listing the missing environments. Create the workflow environments first (see [Workflow Environment Setup](#workflow-environment-setup)).

This downloads all databases into a fixed layout under **DBDIR**. No manual path editing is needed unless you want to override specific databases.

#### Custom paths

If you already have a compatible version of any database stored elsewhere (e.g., a team-shared location or a previous download), you can override its path by setting the **absolute path** directly in `config.<mode>.yaml`. **Once an absolute path is explicitly specified in the config, Metalink will use that path directly and ignore the default DBDIR location for that database.** You do not need to run `metalink download-database` for databases with custom paths.

> [!NOTE]
> If you override only some databases, the remaining ones will still be resolved relative to the default **DBDIR**. Ensure those databases exist at the default location, or override their paths as well.

#### Default path layout
The default DBDIR layout (used by `out_dir/config/config.<mode>.yaml`) is shown below.
You can use a different layout, but then you must update the absolute or relative paths in the config accordingly.

During and after `metalink download-database`, the following directories are created under **DBDIR**:

```text
DBDIR/
├── kraken2_taxonomy_db/          # Metalink kraken2 taxonomy DB
│   └── v20250930/                # Build version
│       ├── kraken2_raw/          # Default
│       ├── kraken2_16GB/         # Optional
│       ├── kraken2_64GB/         # Optional
│       └── kraken2_256GB/        # Optional
├── kraken2_function_db/          # Metalink kraken2 function DB (LCFH)
│   └── v20250930/                # Build version
│       ├── kraken2_raw/          # Default
│       └── kraken2_16GB/         # Optional
├── cat_pack_db/                  # Metalink-customized CAT DB
│   └── v20250930/                # Build version
│       ├── cat_pack_gtdb_ncbi_ictv/
│       │   ├── db/
│       │   └── tax/
│       └── cat_pack_ictv/
│           ├── db/
│           └── tax/
├── eggnog_db/                    # eggNOG-mapper data
│   └── v5.0.2/
│       ├── eggnog_proteins.dmnd
│       └── mmseqs/
├── gtdbtk_db/                    # GTDB-Tk reference data
│   └── R226/
├── checkm2_db/                   # CheckM2 database
│   └── v1.1.0/
│       └── uniref100.KO.1.dmnd
├── checkm_db/                    # CheckM database (used by MetaWRAP binning)
│   └── checkm_data_2015_01_16/
├── checkv_db/                    # CheckV database
│   └── checkv-db-v1.5/
├── genomad_db/                   # geNomad database
│   └── genomad_db_v1.9/
├── virsorter2_db/                # VirSorter2 database
│   └── v2.2.4/
├── vibrant_db/                   # VIBRANT database
│   └── v1.2.1/
│       ├── databases/
│       └── files/
├── iphop_db/                     # iPHoP database
│   └── Jun_2025_pub_rw/
├── downloads/                    # Temporary download archives
├── logs/                         # Download and setup logs
├── markers/                      # Completion flags
├── config/                       # config.download.yaml used for the run
└── .metalink_db_bundle.yaml      # Installed bundle version
```
#### Required databases

| Config key | Database | Version / note | Expected path / files |
|------------|----------|----------------|-----------------------|
| `kraken2_taxonomy_db` | Metalink kraken2 taxonomy DB | v20250930 | `kraken2_taxonomy_db/v20250930/kraken2_raw/` |
| `kraken2_function_db` | Metalink kraken2 function DB (LCFH) | v20250930 | `kraken2_function_db/v20250930/kraken2_raw/` |
| `cat_pack_db` | Metalink-customized CAT/BAT databases (GTDB-NCBI-ICTV + ICTV) | v20250930 | `cat_pack_db/v20250930/` |
| `eggnog_db` | eggNOG-mapper | v5.0.2 | `eggnog_db/v5.0.2/` |
| `genomad_db` | geNomad | v1.9 | `genomad_db/genomad_db_v1.9/` |
| `checkv_db` | CheckV | v1.5 | `checkv_db/checkv-db-v1.5/` |
| `virsorter2_db` | VirSorter2 | v2.2.4 | `virsorter2_db/v2.2.4/` |
| `vibrant_db` | VIBRANT | v1.2.1 | `vibrant_db/v1.2.1/` |
| `iphop_db` | iPHoP | Jun_2025_pub_rw | `iphop_db/Jun_2025_pub_rw/` |
| `gtdbtk_db` | GTDB-Tk | R226 | `gtdbtk_db/R226/` |
| `checkm2_db` | CheckM2 | v1.1.0 | `checkm2_db/v1.1.0/uniref100.KO.1.dmnd` |
| `checkm_db`— (set via `checkm data setRoot`) | CheckM | 2015_01_16 | `checkm_db/checkm_data_2015_01_16/` |

#### Quick setup

1. Download all databases (Metalink-built + third-party):

    ```bash
    # -n for dry run
    metalink download-database -d db_dir
    ```

    Or download only selected databases by listing comma-separated names:

    ```bash
    metalink download-database -d db_dir --databases kraken2_taxonomy,kraken2_function,cat_pack_db
    metalink download-database -d db_dir --databases gtdbtk
    ```

    The Metalink kraken2 databases come in size variants; default is `kraken2_raw`:

    ```bash
    metalink download-database -d db_dir --kraken2-taxonomy-db-size kraken2_16GB
    metalink download-database -d db_dir --kraken2-function-db-size kraken2_16GB
    ```

    Available sizes: taxonomy `kraken2_16GB` / `kraken2_64GB` / `kraken2_256GB` / `kraken2_raw`; function `kraken2_16GB` / `kraken2_raw`.

    Like the analysis modes, `download-database` is config-driven: the options above are written to `db_dir/config/config.download.yaml` before the workflow starts. You can also generate the config without downloading, edit it (database selection, size variants, custom conda environment names in `conda_envs`), and run from it:

    ```bash
    metalink download-database --generate-config -d db_dir
    # edit db_dir/config/config.download.yaml, then:
    metalink download-database -c db_dir/config/config.download.yaml
    ```

    With `-c/--configfile` for `download-database`, the config-source options (`-d`, `--databases`, `--kraken2-taxonomy-db-size`, `--kraken2-function-db-size`) given on the command line override the file values for that run (merged into a temporary runtime config; the file itself is not modified). Edit the file to change the persistent defaults. If you created the conda environments under custom names, update the `conda_envs` mapping in this config; both the pre-download environment check and the workflow itself will use your names.

    Approximate disk and memory requirements of the databases (Metalink-built: v20250930; third-party: versions as shipped by `metalink download-database`):

    | Database | Variant | Disk space | kraken2 RAM |
    |----------|---------|-----------|-------------|
    | kraken2 + bracken taxonomy DB | `kraken2_16GB` | ~35 GB | ~15 GB |
    | kraken2 + bracken taxonomy DB | `kraken2_64GB` | ~79 GB | ~60 GB |
    | kraken2 + bracken taxonomy DB | `kraken2_256GB` | ~256 GB | ~240 GB |
    | kraken2 + bracken taxonomy DB | `kraken2_raw` (default) | ~652 GB | ~640 GB |
    | kraken2 function DB (LCFH) | `kraken2_16GB` | ~29 GB | ~15 GB |
    | kraken2 function DB (LCFH) | `kraken2_raw` (default) | ~74 GB | ~60 GB |
    | CAT_pack DB | — | ~198 GB | — |
    | eggNOG-mapper DB | v5.0.2 | ~59 GB | — |
    | GTDB-Tk DB | R226 | ~139 GB | — |
    | CheckM2 DB | v1.1.0 | ~2.9 GB | — |
    | CheckM DB | 2015_01_16 | ~1.4 GB | — |
    | CheckV DB | v1.5 | ~6.4 GB | — |
    | geNomad DB | v1.9 | ~1.4 GB | — |
    | VirSorter2 DB | v2.2.4 | ~12 GB | — |
    | VIBRANT DB | v1.2.1 | ~8.0 GB | — |
    | iPHoP DB | Jun_2025_pub_rw | ~294 GB | — |

2. You can manually edit `out_dir/config/config.<mode>.yaml` so that all `Directory.database.*` paths point to your DBDIR.

#### Which databases do I need?

- **All modes** need the Metalink kraken2 taxonomy DB, kraken2 function DB, and the corresponding lineage files.
- **Modes involving contig taxonomy** (`c`, `cr`) use kraken2 by default for both contig taxonomy and function and need no extra database; the CAT databases are only required when `contigs_taxonomy_classify_method = CAT_pack`, and the eggNOG-mapper database only when `*_function_classify_method = emapper`.
- **Modes involving MAG taxonomy** (`m`, `mc`, `mcr`, `pmr --mags m`/`vm`) additionally need the CheckM2 database, plus the GTDB-Tk database under the default `MAGs_taxonomy_classify_method = gtdbtk` (not needed with `BAT`, which uses the CAT databases instead). The eggNOG-mapper database is required under the default `MAGs_function_classify_method = emapper` (not needed with `kraken2`).
- **Virus modes** (`v`, `vc`, `vcr`, `vm`, `vmc`, `vmcr`, and `pmr --mags v`/`vm`) additionally need geNomad, VirSorter2, VIBRANT, CheckV, and iPHoP databases.

If you do not plan to run a tool, you can skip its database, but you must also avoid the modes that trigger that tool.

---


## Input File Formats

Metalink uses a single TSV (tab-separated values) file to define all samples. The format is flexible and supports both **paired-end** and **single-end** reads.

### Paired-end Format

| sample | read_1 | read_2 |
|--------|--------|--------|
| sample1 | /path/to/sample1_r1.fastq.gz | /path/to/sample1_r2.fastq.gz |
| sample2 | /path/to/sample2_R1.fastq.gz | /path/to/sample2_R2.fastq.gz |
| sample3 | /path/to/sample3_read1.fq.gz | /path/to/sample3_read2.fq.gz |

### Single-end Format (short single-end or long read)

| sample | read_s |
|--------|--------|
| sample1 | /path/to/sample1.fastq.gz |
| sample2 | /path/to/sample2.fastq.gz |

### Advanced Format (pre-computed contigs and/or MAGs)
For skip-level ≥ 2, you can provide pre-computed contigs and/or MAGs. See [Skip-Level: Reusing Pre-Computed Results](#skip-level-reusing-pre-computed-results)

| sample | read_1 | read_2 | contig_file | mMAG_dir |
|--------|--------|--------|-------------|----------|
| sample1 | /path/to/sample1_r1.fastq.gz | /path/to/sample1_r2.fastq.gz | /path/to/sample1.contigs.fa | /path/to/sample1_mags_dir/ |
| sample2 | /path/to/sample2_r1.fastq.gz | /path/to/sample2_r2.fastq.gz | /path/to/sample2.contigs.fa | /path/to/sample2_mags_dir/ |

For single-end or long-read samples, use a `read_s` column instead of `read_1`/`read_2` (same `contig_file` and `mMAG_dir` columns apply).

> [!IMPORTANT]
>- Supports `.fastq.gz` and `.fq.gz` extensions.
>- The input TSV is validated when the workflow starts: it must be tab-separated, contain a `sample` column plus the columns required by your read type and skip-level (`read_1`/`read_2` for paired, `read_s` for single-end/long reads, `contig_file` for skip-level ≥ 2, `mMAG_dir` for skip-level ≥ 3), and required columns must not be empty. Malformed files are rejected with the exact line number and reason.
>- All FASTA files (contigs, MAGs) must **not** contain `---` in headers (used as internal delimiter).
>- MAG directories should contain `.fa` files.
>- The contig IDs must match between MAG files and the contigs file.

---


## Analysis Modes

Metalink supports 14 analysis modes, specified as the first positional argument. The hierarchy reflects analytical level inclusion and viral scope:

### Mode Hierarchy

```
b (base)
├── r (read)            → read-based direct analysis
├── c (contig)          → contig taxonomy/function
│   └── cr              → c + unmapped read analysis
├── m (mMAG)            → mMAG quality/taxonomy/function/drep/coverage
│   ├── mc              → m + contig analysis
│   └── mcr             → m + contig + read analysis
├── v (vMAG)            → viral identification/clustering/host/coverage
│   ├── vc              → v + contig taxonomy/function
│   ├── vcr             → v + contig + read analysis
│   ├── vm              → v + m + vMAG–mMAG integration
│   ├── vmc             → vm + contig analysis
│   └── vmcr            → vm + contig + read analysis (full integration)
└── pmr                 → project-level pooled MAG reference + read analysis
    ├── --mags m        → pooled mMAGs (dRep ANI 95%) + reads
    ├── --mags v        → pooled vMAGs (vClust ANI 95%) + reads
    └── --mags vm       → pooled mMAGs + vMAGs (deduplicated) + reads
```

### Mode Details

#### `b` — Base Analysis

Steps:
1. Raw read statistics (seqkit)
2. Trimming and QC (automatically selected by `--read-type`):
   - **Paired-end / single-end DNA short reads** (`paired`, `single`): fastp
   - **Paired-end / single-end RNA short reads** (`paired-rna`, `single-rna`): fastp
   - **Long reads** (`pacbio-raw`, `pacbio-corr`, `pacbio-hifi`, `nano-raw`, `nano-corr`, `nano-hq`): fastplong
3. Assembly (automatically selected by `--read-type`; validated against `--assembly-method`):
   - **Paired-end short reads**: MEGAHIT (default runs `meta-sensitive`), spades (default runs `metaSPAdes` for metagenomes and `rnaSPAdes` for metatranscriptomes)
   - **Single-end short reads**: MEGAHIT
   - **Long reads (PacBio / Nanopore)**: metaFlye (supports `pacbio-raw`, `pacbio-corr`, `pacbio-hifi`, `nano-raw`, `nano-corr`, `nano-hq`)
4. Contig length filtering
5. Read mapping (minibwa / bwa-mem2 / bwa mem / strobealign)
6. BAM summary (SAMtools)
7. Binning (MetaWRAP: CONCOCT, MaxBin2, MetaBAT2)
8. Bin refinement (MetaWRAP)
9. CheckM quality assessment
10. Contig depth estimation (jgi_summarize_bam_contig_depths)
11. Prodigal gene prediction on contigs and bins

#### `c` — Contig Analysis

Additional steps:
- Contig taxonomy classification (CAT_pack or kraken2)
- Contig function classification (eggNOG-Mapper or kraken2 against LCFH)
- Merge taxonomy + function + depth
- Taxonomic visualization (stacked barplots)
- Functional visualization

#### `cr` — Contig-Read Analysis

Additional steps:
- Extract unmapped reads from BAM
- Trim and merge unmapped reads (fastp with `--merge`)
- ORF prediction on unmapped reads (FragGeneScanRs)
- Taxonomy classification (kraken2)
- Function classification (kraken2 against LCFH or eggNOG-Mapper)
- Merge contig + read results

#### `r` — Read Analysis

Steps:
- Optional read downsampling (seqkit sample)
- Trim and merge reads (if paired-end) (fastp)
- ORF prediction (FragGeneScanRs)
- Taxonomy classification (kraken2 + Bracken)
- Function classification (kraken2 against LCFH or eggNOG-Mapper)
- Merge results

#### `m` — mMAG Analysis

Additional steps:
- CheckM2 quality assessment
- Taxonomy (GTDB-Tk or CAT_pack BAT mode)
- Sample-level dereplication (dRep at 99% and 95% ANI)
- Coverage estimation (CoverM: `trimmed_mean` or `mean`)
- Function annotation (eggNOG-Mapper or kraken2 against LCFH)
- Merge taxonomy + function + coverage

#### `mc` / `mcr` — Integrated MAG Modes

These extend `m` with contig and/or read analysis at sample level. The mcr mode implements the pipeline: clean reads → contigs → mMAGs → read mapping → unassembled reads → overlap-based deduplication → unified taxonomic-functional profiles.

#### `v` — vMAG Analysis

Additional steps:
- Viral identification (geNomad, or VirSorter2 + CheckV + VIBRANT, or union)
- Prophage prediction (Modified ViWrap scripts)
- Viral contig taxonomy (CAT_pack)
- Viral read mapping and coverage (minibwa / bwa-mem2 / bwa mem / strobealign + SAMtools)
- Viral binning (vRhyme)
- CheckV quality assessment
- Sample-level viral clustering (vClust)
- Host prediction (iPHoP)
- vMAG taxonomy (CAT_pack BAT)
- vMAG coverage (CoverM)

#### `vm` — vMAG + mMAG Integrated Analysis

`vm` runs the **complete** `v` and `m` workflows and adds a dedicated vMAG–mMAG integration layer — it is not a simple stacking of the two modes:

- All `m` steps: CheckM2 quality assessment, GTDB-Tk/BAT taxonomy, within-sample dRep (99% and 95% ANI), CoverM coverage, mMAG function annotation.
- All `v` steps: viral identification, vRhyme binning, CheckV, vClust, iPHoP, vMAG taxonomy and coverage.
- **Joint vMAG–mMAG profiling**: mMAGs passing the CheckM completeness/contamination filter and vMAGs passing the CheckV quality filter are merged into unified taxonomic–functional profiles (`05_merge_tax_fun_depth/`), including a vMAG/mMAG–contig taxonomy consistency check (`MAG_contig_consistency/`).
- **Workflow-derived host prediction**: unlike `v`, `vm` can supply iPHoP with the dereplicated mMAGs produced by the workflow itself (`--virus-iphop-custom-mags-dir-project auto_drep95`, the default in `vm`; see [`v` vs `vm`: iPHoP Host Prediction Strategies](#v-vs-vm-iphop-host-prediction-strategies)).
- **Robust to missing MAGs**: samples without vMAGs fall back to m-mode merging (recorded in `{sample}.mode.tsv`); samples without mMAGs cannot be merged and exit with an error.


#### `vc` / `vcr` / `vmc` / `vmcr` — Integrated vMAG Modes

These extend `v` with contig, read, mMAG, or full integration analysis. The `vmcr` mode implements the complete pipeline: clean reads → contigs → mMAGs + vMAGs → read mapping → unassembled reads → overlap-based deduplication → unified taxonomic-functional profiles.


#### `pmr` — Project-Level MAG + Read Analysis

- Constructs a **project-level non-redundant MAG reference** pooled from ALL samples (structurally different from the per-sample MAG modes), maps every sample's reads against it in a single pass, and integrates taxonomy × function × depth into project-level abundance matrices.
- The `--mags m|v|vm` option selects the MAG source: mMAGs (default), vMAGs, or both combined into one deduplicated reference.
- Variant details, outputs, and usage: [`pmr`: Project-Level Pooled MAGs + Reads](#pmr-project-level-pooled-mags--reads).

> [!NOTE]
> Metabolic completeness analysis (mclink R package) is implemented in 12 modes (`r`, `c`, `cr`, `m`, `mc`, `mcr`, `vc`, `vcr`, `vm`, `vmc`, `vmcr`, `pmr`).
---


## Usage

Metalink is designed around **analysis goals** rather than individual parameters. The CLI provides extensive `--help` output for every mode; this section focuses on **which mode to choose for which analysis goal** and how to construct practical command lines.


### Choosing a Mode

Metalink modes form a progressive hierarchy. Because the workflow is built on Snakemake, **you do not need to run the base mode first**: starting from any mode (e.g. `mcr` or `vmcr`) will automatically build the full dependency graph and execute all required upstream steps (`b`, assembly, mapping, binning, etc.) in the same working directory (`-w`). Completed steps are cached, so you can also run modes sequentially into the same directory if you prefer to inspect intermediate results.

#### Progressive Workflow (Optional)

**Microbial analysis paths**

```bash
# Path A: contig-focused (no MAGs)
metalink b  -i sample.tsv -d db_dir -w out_dir -j 4 --run-project  # optional
metalink c  -i sample.tsv -d db_dir -w out_dir -j 4 --run-project
metalink cr -i sample.tsv -d db_dir -w out_dir -j 4 --run-project

# Path B: MAG-focused
metalink m   -i sample.tsv -d db_dir -w out_dir -j 4 --run-project
metalink mc  -i sample.tsv -d db_dir -w out_dir -j 4 --run-project
metalink mcr -i sample.tsv -d db_dir -w out_dir -j 4 --run-project

# Path C: MAG + reads (no contig integration)
metalink pmr -i sample.tsv -d db_dir -w out_dir -j 4 --run-project
# Need a virus abundance matrix?            → metalink pmr --mags v  ...
# Joint mMAG + vMAG quantification in one mapping pass → metalink pmr --mags vm ...
```

**Virus-related analysis paths**

```bash
# Path D: viral contig-focused (no MAGs)
metalink v   -i sample.tsv -d db_dir -w out_dir -j 4 --run-project  # optional
metalink vc  -i sample.tsv -d db_dir -w out_dir -j 4 --run-project
metalink vcr -i sample.tsv -d db_dir -w out_dir -j 4 --run-project

# Path E: viral MAG-focused
metalink vm   -i sample.tsv -d db_dir -w out_dir -j 4 --run-project
metalink vmc  -i sample.tsv -d db_dir -w out_dir -j 4 --run-project
metalink vmcr -i sample.tsv -d db_dir -w out_dir -j 4 --run-project
```
> [!IMPORTANT]
> Snakemake resolves the entire DAG for the requested mode. If upstream outputs (e.g. trimmed reads, contigs, MAGs) already exist in the working directory, they are reused; otherwise they are generated automatically.

### Mode-Specific Differences and Special Cases


#### `r`: Read-Only Direct Profiling

The `r` mode is independent of assembly and binning. It performs direct taxonomic and functional profiling on reads without any assembly- or binning-level steps. This makes it ideal for:

- **Metatranscriptomes** (after rRNA removal).
- **Quick profiling** when assembly is not needed or not feasible
- **Low-biomass samples** where assembly yields few contigs

```bash
# r mode: direct read profiling (no assembly, no binning)
metalink r -i sample.tsv -d db_dir -w out_dir -j 4 --run-project

# r mode with target stage control (taxonomy only, function only, or both)
metalink r -i sample.tsv -d db_dir -w out_dir -j 4 --run-project --target-stage taxonomy               # Taxonomy only
metalink r -i sample.tsv -d db_dir -w out_dir -j 4 --run-project --target-stage function               # Function only
metalink r -i sample.tsv -d db_dir -w out_dir -j 4 --run-project --target-stage taxonomy_and_function  # Taxonomy and function (not merge)
metalink r -i sample.tsv -d db_dir -w out_dir -j 4 --run-project --target-stage all                    # Taxonomy and function (merge)
```

#### `v` vs `vm`: iPHoP Host Prediction Strategies

The most important difference between `v` (viral-only) and `vm` (viral + microbial) modes is **iPHoP host prediction configuration**:

| Mode | iPHoP custom_MAGs_dir Options | Behavior |
|------|------------------------------|----------|
| `v`, `pmr --mags v` | `none` (default) or **user-defined path** only | Host prediction with iphop default database or user-provided genomes |
| `vm`,`pmr --mags vm` | `none`, `auto_drep99`, `auto_drep95` (default), or **user-defined path** | Automatically uses dereplicated microbial MAGs from the workflow for host prediction; `auto_drep99` uses 99% ANI dereplicated mMAGs, `auto_drep95` uses 95% ANI |


```bash
# v mode: host prediction with iphop default database
metalink v -i sample.tsv -d db_dir -w out_dir -j 4 --virus-iphop-custom-mags-dir-project none

# v mode: host prediction with external MAGs (user provides path)
metalink v -i sample.tsv -d db_dir -w out_dir -j 4 --virus-iphop-custom-mags-dir-project /path/to/external/mags/

# vm mode: host prediction with auto-derived microbial MAGs
metalink vm -i sample.tsv -d db_dir -w out_dir -j 4 --virus-iphop-custom-mags-dir-project auto_drep95
```

#### `pmr`: Project-Level Pooled MAGs + Reads

The `pmr` mode is structurally different from all other modes:

- **Other modes** (`m`, `mc`, `mcr`, `vm`, `vmc`, `vmcr`): Perform **sample-level** MAG binning and dereplication. Each sample is analyzed independently; cross-sample comparison occurs only at the project-level aggregation step.
- **pmr mode**: Constructs a **project-level non-redundant MAG reference** pooled from ALL samples, maps every sample's reads against this reference, and integrates taxonomy × function × depth into project-level abundance matrices.

**`--mags m|v|vm` selects the MAG source** (default `m`):

| `--mags` | MAG reference | Chains | Output directories | Cross-sample abundance matrix |
|----------|---------------|--------|--------------------|-------------------------------|
| `m` (default) | Project-level dRep ANI 95% mMAGs | b + m-s + m-p + pmr_m | `results/pmr_m-sample/`, `results/pmr_m-project/` | `mMAG_species_abundance_matrix.tsv` |
| `v` | Project-level vClust ANI 95% vMAGs | b + v-s + v-p + pmr_v | `results/pmr_v-sample/`, `results/pmr_v-project/` | `vMAG_species_abundance_matrix.tsv` |
| `vm` | mMAGs and vMAGs combined into one reference | b + v-s + v-p + m-s + m-p + pmr_vm | `results/pmr_vm-sample/`, `results/pmr_vm-project/` | `mMAG_vMAG_species_abundance_matrix.tsv` |

Variant-specific notes:

- `--mags m` (p+mMAG+r, default): pools the project-level dRep ANI 95% mMAGs into one reference and maps every sample's reads against it, producing the cross-sample abundance matrix `mMAG_species_abundance_matrix.tsv`, the full taxonomy × function × depth downstream, mclink module completeness, and project-level visualization.
- `--mags v` (p+vMAG+r): pools the project-level vClust ANI 95% vMAGs into one reference and maps every sample's reads against it, producing the cross-sample abundance matrix `vMAG_species_abundance_matrix.tsv`, the full taxonomy × function × depth downstream, mclink module completeness, and project-level visualization.
- `--mags vm` (p+vMAG+mMAG+r) quantifies mMAGs and vMAGs jointly in a single pass: one combined reference, one index build, one mapping, one CoverM run, and one integrated downstream. The combined reference is deduplicated before indexing: proviral regions on mMAG contigs are N-masked and fully-viral contigs in mMAGs are dropped, and host ORFs overlapping viral regions are excluded from the host KO tables. Audit outputs are written under `results/pmr_vm-project/02_quality_gate/` (`{project}.viral_intervals.bed`, `{project}.excluded_ORFs.txt`, `{project}.fully_viral_contigs.txt`).

**Quality gate**: before building the mapping reference, `pmr` runs `02_quality_gate` (config section `Params_pmr_{m,v,vm}_project → 02_quality_gate`). It filters the mMAG and vMAG sets with the thresholds below. For `--mags m`/`v` the passing genomes are linked into a gated directory (`results/pmr_m-project/02_quality_gate/{project}.mMAGs_gated/`, `results/pmr_v-project/02_quality_gate/{project}.vMAGs_gated/`); for `--mags vm` the mMAG gating is applied inside `{project}.mMAGs_masked_genomes` and vMAGs are linked into `{project}.vMAGs_gated/`. The mapping reference, CoverM output, and downstream taxonomy/function tables are all built from the gated set.

| Option | Description | Default |
|--------|-------------|---------|
| `--mmags-min-completeness1/2` | Minimum CheckM1/CheckM2 completeness (%) for mMAG gating | 0 |
| `--mmags-max-contamination1/2` | Maximum CheckM1/CheckM2 contamination (%) for mMAG gating | 100 |
| `--vmags-checkv-quality` | Minimum CheckV tier for vMAG gating, one of `Complete` / `High-quality` / `Medium-quality` / `Low-quality` / `Not-determined` | `Not-determined` (keep all vMAGs) |
| `--merge-vmags-checkv-quality` | CheckV tier filter applied when merging taxonomy × function × depth in the `vc`/`vcr`/`vm`/`vmc`/`vmcr` modes | `High-quality` |

```bash
# pmr mode: project-level pooled mMAGs + read mapping analysis (default, --mags m)
metalink pmr --mags m -i sample.tsv -d db_dir -w out_dir -j 4 --run-project

# virus abundance matrix: project-level vMAGs + reads
metalink pmr --mags v -i sample.tsv -d db_dir -w out_dir -j 4 --run-project

# joint mMAG + vMAG quantification with a single mapping pass
metalink pmr --mags vm -i sample.tsv -d db_dir -w out_dir -j 4 --run-project
```

### Skip-Level: Reusing Pre-Computed Results

Use `--skip-level` (0–3) to skip upstream steps if you already have pre-computed results:

| Skip-level | Skips | Required Inputs | Use Case |
|------------|-------|-----------------|----------|
| `0` (default) | Nothing | Raw FASTQ | Start from scratch |
| `1` | Trimming | Clean FASTQ | Pre-trimmed reads available |
| `2` | Trimming + Assembly | Clean FASTQ + Contigs | Pre-assembled contigs available |
| `3` | Trimming + Assembly + Binning | Clean FASTQ + Contigs + mMAGs | Pre-binned mMAGs available |

> [!NOTE]
> - The choice of `--read-type` determines both the trimming tool and the allowed assembler at skip-level `0` and `1`.
> - When using externally assembled contigs, the FASTA requirements in [Input File Formats](#input-file-formats) apply: contig IDs must not contain `---` (reserved as internal delimiter), and MAG contig IDs must match the original contigs file when skip-level `3` is used.
> - For `v`/`vc`/`vcr` and `pmr --mags v`, skip-level `3` is equivalent to skip-level `2` (these modes do not use mMAGs).

#### Advanced skip-level use cases

Metalink's built-in assembler and binner combinations cover the most common scenarios (MEGAHIT, SPAdes, metaFlye; MetaBAT2, MaxBin2, CONCOCT). For workflows that require tools or parameters not directly exposed by Metalink, use `--skip-level` to inject externally generated results and continue with Metalink's downstream taxonomy, function, and integration analyses.

| Use case | Skip-level | What you provide | Metalink mode | What Metalink does |
|----------|-----------|------------------|---------------|-------------------|
| **Hybrid assembly** (2nd-gen + 3rd-gen mixed) | `2` | Pre-assembled contigs (from external hybrid assembler) + clean reads (paired or long reads) | Any (`c`, `m`, `mc`, `mcr`, etc.) | Skips trimming and assembly; runs mapping, binning, annotation, and MAGs/contigs/reads taxonomy/function annotation and project-level aggregation |
| **Long-read binning** | `3` | Pre-binned MAGs (from external binner) + contigs (from external assembler) + clean reads | Any (`m`, `mc`, `mcr`, etc.) | Skips trimming, assembly, and binning; runs MAGs/contigs/reads taxonomy/function annotation and project-level aggregation |
| **Metatranscriptome read mapping to MAGs** | `3` | Pre-binned MAGs (from external binner) + metatranscriptome reads | `pmr` | Skips trimming, assembly, and binning; maps reads to provided MAGs and performs project-level MAG + read analysis |

> [!TIP]
> - For hybrid assembly, generate contigs with your preferred tool, then supply them in the TSV's `contig_file` column and run with `--skip-level 2`.
> - For metatranscriptomes, rRNA should be removed.
> - For `--skip-level 3`, place each sample's MAG directory (with existing MAGs) in the `mMAG_dir` column of the TSV.


**Skip-level 1 example** (pre-trimmed reads):

```bash
# Provide clean FASTQ paths in the TSV (or reuse from previous Metalink run)
metalink c -i samples_clean.tsv -d db_dir -w out_dir --skip-level 1 -j 4 --run-project
```

**Skip-level 2 example** (pre-assembled contigs):

```bash
# TSV format: sample  read_1  read_2  contig_file
# Column `contig_file` must contain contig paths
metalink c -i samples_with_contigs.tsv -d db_dir -w out_dir --skip-level 2 -j 4 --run-project
```

**Skip-level 2 example** (metatranscriptome contigs assembled externally):

```bash
# Assemble metatranscriptomes externally (e.g., Trinity) with rRNA removed,
# then supply contigs to Metalink for taxonomy/function annotation.
# TSV: sample  read_1  read_2  contig_file
metalink c -i mt_samples_with_contigs.tsv -d db_dir -w out_dir --skip-level 2 --read-type paired -j 4 --run-project
```

> [!NOTE]
> For metatranscriptome datasets, virus identification with the `vs_vb` or `genomad_vs_vb` method (VirSorter2 + VIBRANT, optionally combined with geNomad) runs VirSorter2 with the `--include-groups dsDNAphage,NCLDV,RNA,ssDNA,lavidaviridae` parameter, so RNA viruses are included in the search space.

**Skip-level 2 example** (hybrid 2nd-gen + 3rd-gen assembly):

```bash
# Assemble mixed short and long reads externally,
# then supply the hybrid contigs to Metalink for downstream analysis.
# TSV: sample  read_1  read_2  contig_file
metalink mc -i hybrid_samples.tsv -d db_dir -w out_dir --skip-level 2 --read-type paired -j 4 --run-project
```

**Skip-level 3 example** (pre-binned MAGs):

```bash
# TSV format: sample  read_1  read_2  contig_file  mMAG_dir/
# Column `mMAG_dir` must contain paths of MAGs of each sample
metalink vmcr -i samples_with_mags.tsv -d db_dir -w out_dir --skip-level 3 -j 4 --run-project
```

**Skip-level 3 example** (metatranscriptome reads mapped to pre-existing MAGs):

```bash
# Map metatranscriptome reads to MAGs recovered from a matching metagenomic sample.
# TSV: sample  read_1  read_2  contig_file  mMAG_dir/
metalink pmr -i mt_samples_with_mags.tsv -d db_dir -w out_dir --skip-level 3 --read-type paired -j 4 --run-project
```
### Read Type
#### Read Type Behavior

| `--read-type` | Trimming tool | Assembler (default) | Behavior |
|---------------|---------------|---------------------|------------------------|
| `paired` | fastp | MEGAHIT / metaSPAdes | Uses both R1 and R2. Paired-end specific rules are active. |
| `single` | fastp | MEGAHIT | Uses the single-end reads (`read_s` column of the sample sheet). Single-end specific rules are active. |
| `paired-rna` | fastp | rnaSPAdes | Paired-end RNA short reads (rRNA must be removed before input). |
| `single-rna` | fastp | rnaSPAdes | Single-end RNA short reads (rRNA must be removed before input). |
| `pacbio-raw` | fastplong | metaFlye | PacBio raw CLR long reads (`read_s` column of the sample sheet). |
| `pacbio-corr` | fastplong | metaFlye | PacBio corrected CLR long reads. |
| `pacbio-hifi` | fastplong | metaFlye | PacBio HiFi long reads. |
| `nano-raw` | fastplong | metaFlye | Oxford Nanopore raw long reads. |
| `nano-corr` | fastplong | metaFlye | Oxford Nanopore corrected long reads. |
| `nano-hq` | fastplong | metaFlye | Oxford Nanopore high-quality (Q20+) long reads. |

> [!IMPORTANT]
> `paired`, `single`, `paired-rna`, and `single-rna` represent 2nd-generation short-read data, while `pacbio-raw`, `pacbio-corr`, `pacbio-hifi`, `nano-raw`, `nano-corr`, and `nano-hq` represent 3rd-generation long-read data. Each long-read `read-type` maps to the corresponding Flye mode:
>  - `pacbio-raw` → `--pacbio-raw`
>  - `pacbio-corr` → `--pacbio-corr`
>  - `pacbio-hifi` → `--pacbio-hifi`
>  - `nano-raw` → `--nano-raw`
>  - `nano-corr` → `--nano-corr`
>  - `nano-hq` → `--nano-hq`

#### Paired-End Data

All modes support paired-end short reads via `--read-type paired`.

```bash
# Paired-end with custom contig length threshold
metalink c -i pe_sample.tsv -d db_dir -w out_dir --read-type paired --depth-min-contig-length 1000 -j 4 --run-project

# Paired-end read-only profiling
metalink r -i pe_sample.tsv -d db_dir -w out_dir --read-type paired -j 4 --run-project
```

#### Single-End Data

All modes support single-end short reads via `--read-type single`.

```bash
# Single-end with custom contig length threshold
metalink c -i se_sample.tsv -d db_dir -w out_dir --read-type single --depth-min-contig-length 1000 -j 4 --run-project

# Single-end read-only profiling
metalink r -i se_sample.tsv -d db_dir -w out_dir --read-type single -j 4 --run-project
```

#### Long-Read Data

The 3rd-generation long reads (PacBio and Oxford Nanopore) are supported via dedicated `--read-type` values (`pacbio-hifi`, `nano-hq`, etc.):

```bash
# PacBio HiFi assembly + contig analysis
metalink c -i pacbio_sample.tsv -d db_dir -w out_dir --read-type pacbio-hifi -j 4 --run-project
```

> [!NOTE]
>- Long-read preprocessing uses `fastplong` instead of `fastp`.
>- `--assembly-method` is automatically validated against `--read-type`: long-read types require `metaflye`.

> [!WARNING]
> Metalink does **not** currently ship a gene-calling tool specifically designed for 3rd-generation (PacBio / Oxford Nanopore) raw or corrected long reads. Read-level ORF prediction is performed by **FragGeneScanRs**, which was designed for short-read sequences. It may complete on high-quality 3rd-generation reads such as `pacbio-hifi` or `nano-hq`, but results should be interpreted with caution.


#### Metatranscriptomic Data

Metalink supports metatranscriptomic analysis via the `c`, `cr`, `r`, `v`, `vc`, and `vcr` modes. **rRNA and other non-coding RNA sequences must be removed prior to input** (e.g., with [RiboDetector](https://github.com/hzi-bifo/RiboDetector)). Metalink assumes the input reads are already protein-coding-enriched.

```bash
# Paired-end metatranscriptome (provide rRNA-depleted reads)
metalink r -i mt_samples_rrnaremoved.tsv -d db_dir -w out_dir --read-type paired-rna -j 4

# Single-end metatranscriptome (provide rRNA-depleted reads)
metalink r -i mt_samples_rrnaremoved.tsv -d db_dir -w out_dir --read-type single-rna -j 4
```

> [!IMPORTANT]
> Metalink does not include a metatranscriptome binner. For MAG-level analysis, you can run `metalink pmr --skip-level 3` with pre-binned MAGs (e.g., metagenomic bins with metatranscriptomic reads).
>
> When analyzing metatranscriptomes with a virus mode, the `vs_vb`/`genomad_vs_vb` virus identification method runs VirSorter2 with the `--include-groups dsDNAphage,NCLDV,RNA,ssDNA,lavidaviridae` parameter.


### Configuration & Parameters

#### How the Config Is Built

When you run Metalink without `-c/--configfile`, it builds `{work_dir}/config/config.<mode>.yaml` (`config.pmr_{m,v,vm}.yaml` for `pmr`, keyed by `--mags`) from layered sources. Later layers override earlier ones:

1. `metalink/config/template.config.<mode>.yaml` — defaults (`template.config.pmr_{m,v,vm}.yaml` for `pmr`)
2. CLI overrides — values passed on the command line

#### Reusing an Existing Config (`-c/--configfile`)

To re-run or continue an analysis from an existing config file:

```bash
cp config/config.c.yaml config/config.c.manually.yaml
metalink c -c config/config.c.manually.yaml
```

When `-c` is used, Metalink loads the existing config directly and does **not** regenerate it from the templates.

**Which parameters can be used with `-c`?**

- **Config-source parameters** (`--read-type`, `--skip-level`, `-i`, `-d`, `-w`) define the run's identity (read type, reuse strategy, input, databases, working directory). They are **mutually exclusive with `-c`**; to change them, edit the config file directly or regenerate it with `-i`.
- **Additive parameters and mode-specific parameters** (e.g. `--taxonomy-rank`, `--run-project` / `--project-name`, `--kraken2-taxonomy-jobs`, `--delete-bam`, `--depth-min-contig-length`, `--megahit-presets`) override values inside the existing config and **can be used with `-c`**. In this case Metalink writes the overrides to a temporary `config/config.{mode}.runtime.*.yaml` and passes that file to Snakemake; the original `config.{mode}.yaml` remains unchanged:

```bash
metalink c -c config/config.c.manually.yaml --depth-min-contig-length 1000
```
> [!WARNING]
> With project-level chains, a new `--project-name` places the whole chain under the new project directory, i.e. that part reruns.


#### Key Global Parameters

| Parameter | Description | Typical Value |
|-----------|-------------|---------------|
| `--read-type` | `paired`, `single`, `paired-rna`, `single-rna`, `pacbio-raw`, `pacbio-corr`, `pacbio-hifi`, `nano-raw`, `nano-corr`, `nano-hq` | `paired` (default) |
| `--skip-level` | 0, 1, 2, or 3 | 0 (default) |
| `-i, --input-tsv` | Sample sheet TSV | Required unless `-c` is used |
| `-d, --db-dir` | Database root directory | Optional; defaults to the `DBDIR` value from the config template |
| `-w, --working-dir` | Working/output directory | Optional; defaults to `.` (current directory) |
| `--jobs, -j` | Max parallel jobs/cores (local: total concurrent threads; cluster: submitted jobs) | Number of CPU cores |
| `--run-project` | Enable project-level aggregation | Flag |
| `--project-name` | Project name for outputs | Optional |
| `-g, --generate-config` | Generate `config/config.<mode>.yaml` and exit (edit it, then re-run with `-c`) | Flag |
| `--taxonomy-rank` | Comma-separated ranks | Phylum,Class,Order,Family,Genus,Species |
| `--kraken2-taxonomy-jobs` | Cap concurrent kraken2 taxonomy jobs (see [Limiting kraken2 taxonomy concurrency](#limiting-kraken2-taxonomy-concurrency)) | Optional; default unlimited |
| `--silent` | Suppress all Metalink log messages | Flag |

#### Parameter Help

Run `metalink <mode> --help` for the complete list of mode-specific parameters. Parameters are grouped by analytical section (e.g. `base sample-level options`, `vMAG project-level options`, `vMAG-mMAG sample-level options`) for easy navigation. The cleanup flags and `--kraken2-taxonomy-jobs` live in their own `Intermediate files & resources` group.

#### Customizing Tool Parameters with --extra-args

Metalink exposes `extra_args` parameters for key third-party tools, allowing fine-grained control without modifying the workflow source code. These are passed directly to the tool's command line.

> [!IMPORTANT]
> **Always quote `extra_args` values**. Whether you pass a single flag/value or multiple flags, wrap the entire argument in quotes so the shell passes it as one token to Metalink:
> ```bash
> # Single value
> --spades-extra-args "--rna"
> # Multiple values
> --metaflye-extra-args "--iterations 2 --meta"
> --mmags-drep99-sample-extra-args "--completeness_weight 1 --contamination_weight 5"
> ```

**Common use cases**:

```bash
# fastp: limit reads for testing (process all reads = 0)
metalink b -i sample.tsv -d db_dir -w out_dir -j 4 --trim-extra-args "--reads_to_process 0"

# fastplong (long-read trimming): adjust quality threshold
metalink b -i sample.tsv -d db_dir -w out_dir -j 4 --read-type nano-hq --trim-extra-args "--qualified_quality_phred 10"

# MEGAHIT: disable mercy k-mer extension for faster assembly
metalink b -i sample.tsv -d db_dir -w out_dir -j 4 --megahit-extra-args "--no-mercy"

# metaFlye: add iterations
metalink b -i sample.tsv -d db_dir -w out_dir -j 4 --read-type pacbio-hifi --metaflye-extra-args "--iterations 2"

# Mapping: adjust alignment scoring (minibwa / bwa-mem2)
metalink b -i sample.tsv -d db_dir -w out_dir -j 4 --mapping-extra-args "-B 8"

# MetaWRAP binning: adjust minimum contig count per bin
metalink b -i sample.tsv -d db_dir -w out_dir -j 4 --binning-metawrap-extra-args "-m 5"

# MetaWRAP refinement: remove ambiguously assigned contigs
metalink b -i sample.tsv -d db_dir -w out_dir -j 4 --metawrap-refinement-extra-args "--remove-ambiguous"

# Prodigal: mask runs of N (-m) and bypass the Shine-Dalgarno trainer (-n)
# (by default contigs are predicted with `-p meta` and MAGs with `-p single`)
metalink b -i sample.tsv -d db_dir -w out_dir -j 4 --contigs-prodigal-extra-args "-m -n" --bin-prodigal-extra-args "-m -n"

# dRep: custom weighting for completeness/contamination
metalink m -i sample.tsv -d db_dir -w out_dir -j 4 --mmags-drep99-sample-extra-args "--completeness_weight 1 --contamination_weight 5" --mmags-drep95-sample-extra-args "--completeness_weight 1 --contamination_weight 5"

# GTDB-Tk: minimum alignment fraction
metalink m -i sample.tsv -d db_dir -w out_dir -j 4 --mmags-taxonomy-gtdbtk-extra-args "--min_af 0.5"

# CAT_pack: block size for DIAMOND
metalink c -i sample.tsv -d db_dir -w out_dir -j 4 --contigs-taxonomy-cat-pack-extra-args "--block_size 12.0"

# eggNOG-Mapper: minimum percent identity
metalink c -i sample.tsv -d db_dir -w out_dir -j 4 --contigs-function-emapper-extra-args "--pident 30"

# vClust: custom prefilter k-mer size
metalink v -i sample.tsv -d db_dir -w out_dir -j 4 --vclust-sample-prefilter-extra-args "-k 25" --vclust-sample-align-extra-args "--mal 11" --vclust-sample-cluster-extra-args "--len_ratio 0"

# CoverM: suppress output for cleaner logs
metalink m -i sample.tsv -d db_dir -w out_dir -j 4 --mmags-coverm-sample-extra-args "--quiet"
```


### Workflow Management

Metalink is fundamentally a Snakemake workflow. Any valid Snakemake CLI argument can be appended after the Metalink options and is passed through to Snakemake unchanged.

- **Dependency tracking**: Only re-runs jobs when inputs or parameters change.
- **Parallel execution**: Scales from single-core laptops to HPC clusters (Slurm, PBS, etc.).
- **Resume**: Interrupted workflows can be safely resumed without losing progress.

> **Snakemake docs**: [https://snakemake.readthedocs.io](https://snakemake.readthedocs.io)

#### Common Snakemake flags

| Flag | Purpose | Example |
| --- | --- | --- |
| `-n`, `--dry-run` | Preview the execution plan without running anything. Highly recommended before a new run. | `metalink c -i sample.tsv -d db_dir -w out_dir -n` |
| `--summary` | Show status of all output files. | `metalink c --configfile config/config.{mode}.yaml --summary` |
| `--forcerun <rule>` | Force re-execution of a rule and all its dependents. | `metalink c ... --forcerun b_sample_02_assembly_megahit` |
| `--ri`/`--rerun-incomplete` | Re-run jobs that produced incomplete output files. | `metalink c ... --rerun-incomplete` |
| `--verbose` | Show Snakemake progress. By default Metalink passes `--quiet` to Snakemake; `--verbose` disables this. | `metalink c ... --verbose` |
| `--dag` / `--rulegraph` | Print the workflow DAG / rule graph in DOT format (pipe to Graphviz). | `metalink c ... --dag \| dot -Tpng > workflow.png` |
| `--profile <dir>` | Submit jobs via a Snakemake profile **directory** (containing a `config.yaml` with scheduler settings such as `cluster: "sbatch ..."`). It is a path, not a scheduler name. | `metalink c ... --profile /path/to/cluster/config --jobs 100` |

For how to write a profile directory (e.g., a `config.yaml` with `cluster: "sbatch ..."` or `qsub` settings), see the [Snakemake profiles documentation](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles).

> [!NOTE]
> `--verbose` only changes whether Metalink passes `--quiet` to Snakemake (it does by default); it does not change how much individual tools log. Tool output is still written to `logs/{mode}-{scope}/{step_dir}/{sample}.{rule}.log` (e.g. `logs/b-sample/02_assembly/megahit/sample1.megahit.log`).

#### Generating a DAG

In graph mode (`--dag` / `--rulegraph` / `--filegraph`), Metalink suppresses log output so that stdout is pure DOT and can be piped directly to Graphviz:

```bash
metalink b -i sample.tsv -w out_dir --rulegraph | dot -Tpdf > Metalink.b.sample.dag.pdf

metalink m -i sample.tsv -w out_dir --run-project --rulegraph | dot -Tpdf > Metalink.m.project.dag.pdf

metalink vmcr -i sample.tsv -w out_dir --rulegraph | dot -Tpdf > Metalink.vmcr.sample.dag.pdf
```

#### Resuming, unlocking, and cleaning up

Snakemake tracks completion by output file existence and timestamps, so re-running the same command usually resumes safely:

```bash
metalink c -i sample.tsv -d db_dir -w out_dir -j 4
```

For partial outputs (e.g. truncated BAMs or logs), add `--ri`/`--rerun-incomplete`:

```bash
metalink c ... --ri
```

If a run was killed abruptly, Snakemake may leave a lock in `.snakemake/locks/` and subsequent runs fail with `Error: Directory cannot be locked`. Release the lock with:

```bash
metalink <mode> -i sample.tsv -d db_dir -w out_dir --unlock
```

To force a fresh start, remove Snakemake metadata or specific outputs (use with caution):

```bash
# Remove Snakemake metadata (forces full re-run)
rm -r {work_dir}/.snakemake/

# Or remove a specific output
rm {work_dir}/results/b-sample/03_mapping/minibwa/sample1.bam
rm {work_dir}/markers/b-sample/03_mapping/minibwa/sample1.mk
```

> [!NOTE]
> **Why delete both?** Metalink rules do not all declare their outputs the same way: some rules track completion only through a marker file (`markers/.../{sample}.mk`), some only through the result file(s) in `results/`, and some declare both. Snakemake considers a step finished when **all** of its declared outputs exist, so delete both the marker and the result files for the step, or check which files a step actually tracks with `metalink <mode> --configfile {work_dir}/config/config.{mode}.yaml --summary`.


#### Deleting large intermediate files

Large intermediate files are deleted per category. The four CLI flags are the quick way to control this (`--delete-bam` deletes mapping BAMs; `--keep-assembly-intermediate`, `--keep-drep-data`, `--keep-vs2-iter0`  keep the respective category). They work both when generating a new config and with `-c` (runtime override for that run). **Running several modes sequentially into the same working directory is safe by default**: the BAMs are kept, so later modes reuse the existing mapping instead of re-running it. If disk space matters, simply add `--delete-bam` to the *last* mode you run into that directory.

Under the hood the flags write the `keep_intermediate` mapping in `config/config.{mode}.yaml` (`true` keeps everything, `false` deletes everything), which you can also edit directly for persistent defaults. Each deletion happens only after the last consuming step, so results are never removed before they are needed:

> | Category key | What is deleted when `false` | Default | Why |
> |---|---|---|---|
> | `bam` | Mapping BAMs — the `b-sample/03_mapping` BAMs (after binning/CoverM), the `cr`/`vcr` unmapped-read BAMs (after extraction), the virus sample-level BAMs under `v-sample` (after CoverM), and the pmr chains' project-reference BAMs under `pmr_m/pmr_v/pmr_vm-sample/03_mapping` (after CoverM) | `true` (keep) | Deleting them forces a costly re-mapping when additional modes are later run into the same working directory |
> | `drep_data` | The `data/` working subdirectories (MASH/prodigal/checkM files) inside the sample- and project-level dRep result folders (`05_drep_eachsample_99/95`, `01_drep_project_99/95`); `dereplicated_genomes/` and `data_tables/` are always kept | `false` (delete) | Pure dRep working files, never read downstream |
> | `vs2_iter0` | The VirSorter2 `iter-0` working directories under `v-sample/01_identification/{sample}/vs_vb_result/` (right after identification) | `false` (delete) | Multi-GB per sample; only the small `final_*` files are results |
> | `assembly` | Assembler-internal intermediates under `b-sample/02_assembly/` (MEGAHIT `intermediate_contigs/`, SPAdes `K*`/`corrected`/`tmp`/`misc`, Flye `00-assembly`/`10-consensus`/`20-repeat`/`30-contigger`/`40-polishing`, removed once `03_assembly_filter` completes and only for `--skip-level` 0/1) plus the metaWRAP binning `work_files/assembly*` copy (after binning/refinement) | `false` (delete) | Assembler/binning internals, never read downstream |

#### Limiting kraken2 taxonomy concurrency

kraken2 taxonomy rules are resource-tagged (`kraken2_taxonomy=1`): each running job holds a full copy of the kraken2 taxonomy DB in memory, so on memory-tight machines too many concurrent kraken2 jobs can exhaust RAM. `--kraken2-taxonomy-jobs N` (config key `kraken2_taxonomy_jobs`) caps how many of them run at once — e.g. `-j 330 --kraken2-taxonomy-jobs 2` lets everything else scale freely while at most two kraken2 taxonomy jobs run concurrently. The default (absent) is unlimited.

### How to Tell if a Run Completed Successfully

Metalink uses Snakemake's target-based execution model. A successful run is confirmed by:

#### 1. Exit Code

A successful run returns exit code `0`.

#### 2. Target File Presence

All expected final outputs exist in the working directory (e.g., final summary tables).

#### 3. Snakemake Summary

Run the summary command to verify all jobs are `ok`:

```bash
metalink b --configfile {work_dir}/config/config.{mode}.yaml --summary
```

Expected output (all rows should show `ok`):

```
markers/b-sample/02_assembly/megahit/sample1.mk	ok	b_sample_02_assembly_megahit	...
results/b-sample/03_mapping/minibwa/sample1.bam	ok	 b_sample_03_mapping_minibwa	...
...
```

#### 4. Log Inspection

Check the main log and rule logs for errors:

```bash
# snakemake logs
find {work_dir}/.snakemake/logs/ -name "*.log" -exec grep -l "ERROR in rule\|Error in rule" {} \;

# Rule-specific logs (should be empty or contain only completion messages)
find {work_dir}/logs/ -name "*.log" -exec grep -l "ERROR\|Error\|error" {} \;

# check log
tail {work_dir}/logs/b-sample/02_assembly/megahit/sample1.megahit.log
```

---

## Output File Formats


### Working Directory Structure (Generated at Runtime)

When you run Metalink with `-w /path/to/output`, the following structure is created. Results are separated by **mode** (`b`, `c`, `cr`, `m`, `mc`, `mcr`, `pmr`, `r`, `v`, `vc`, `vcr`, `vm`, `vmc`, `vmcr`) and by **scope** (`sample` vs `project`):

```
output/
├── config/                                  # Generated configuration files
│   ├── config.{mode}.yaml                   # e.g. config.c.yaml, config.vmcr.yaml
│   └── config.{mode}.runtime.*.yaml         # Temporary runtime configs (CLI overrides)
├── logs/                                    # Rule-level logs
│   └── {mode}-{scope}/                      # e.g. logs/b-sample/, logs/r-project/
│       └── {step_dir}/{sample}.{rule}.log   # e.g. logs/b-sample/02_assembly/megahit/sample1.megahit.log
├── benchmarks/                              # Snakemake benchmark files (runtime/resources per job)
├── markers/                                 # Completion flags (empty .mk files) per mode
│   ├── {mode}-sample/
│   └── {mode}-project/
├── results/                                 # Actual analysis outputs
│   ├── {mode}-sample/                       # Per-sample results (e.g. b-sample/, vm-sample/)
│   └── {mode}-project/                      # Project-level aggregation (with --run-project)
└── .snakemake/                              # Snakemake internal metadata
    ├── metadata/
    ├── locks/
    └── ...
```

Key points:

- `markers/{mode}-{scope}/` contains empty `.mk` flag files; their presence means Snakemake considers the corresponding target complete.
- `results/{mode}-sample/` holds per-sample intermediate and final files.
- `results/{mode}-project/` holds project-level aggregated tables (created when `--run-project` is used). Project outputs are typically grouped under a `{project}/` sub-directory (default or user-defined via `--project-name`).
- Each mode only creates the directories it needs.
- The `pmr` mode writes to `pmr_m-*`, `pmr_v-*`, or `pmr_vm-*` directories (results/logs/markers/benchmarks alike) depending on `--mags` (default `m`).
- Metalink's own messages are printed to stderr; Snakemake's internal logs live under `.snakemake/log/`.

---

### Results by Mode and Scope

All final and intermediate outputs live under `results/{mode}-{scope}/`.

- `b-sample/` and `b-project/` — base mode (trimming, assembly, binning, QC).
- `c-sample/` and `c-project/` — contig taxonomy and function.
- `cr-sample/` and `cr-project/` — contigs + unmapped reads integration.
- `m-sample/` and `m-project/` — mMAG quality, taxonomy, function, dRep.
- `mc-sample/` and `mc-project/` — mMAG + contig integration.
- `mcr-sample/` and `mcr-project/` — mMAG + contig + unmapped reads integration.
- `pmr_m-sample/` and `pmr_m-project/` — project-level mMAG + read analysis (`pmr --mags m`, default).
- `pmr_v-sample/` and `pmr_v-project/` — project-level vMAG + read analysis (`pmr --mags v`).
- `pmr_vm-sample/` and `pmr_vm-project/` — project-level combined vMAG + mMAG + read analysis (`pmr --mags vm`).
- `r-sample/` and `r-project/` — read-level taxonomy and function.
- `v-sample/` and `v-project/` — viral identification, vRhyme binning, CheckV, iPHoP host prediction, etc.
- `vc-sample/` and `vc-project/` — vMAG + contig taxonomy and function.
- `vcr-sample/` and `vcr-project/` — vMAG + contigs + unmapped reads integration.
- `vm-sample/` and `vm-project/` — vMAG + mMAG quality, taxonomy, function, iPHoP host prediction, etc.
- `vmc-sample/` and `vmc-project/` — vMAG + mMAG + contig integration.
- `vmcr-sample/` and `vmcr-project/` — vMAG + mMAG + contig + unmapped reads integration.

### Sample-level Outputs

Below are representative subdirectories for the main modes (actual contents vary with the mode and skip-level).

#### Base mode (`b`) sample layout

```
results/b-sample/
├── 00_raw_stats/                  # seqkit stats of raw reads ({sample}_raw_{1,2,s}.fastq.stats)
├── 01_trim/                       # Trimmed FASTQ files ({sample}_trimmed_{1,2,s}.fastq.gz)
├── 01_trim_html/                  # fastp HTML reports ({sample}_fastp.html)
├── 01_trim_json/                  # fastp JSON reports ({sample}_fastp.json)
├── 01_trim_stats/                 # seqkit stats of trimmed reads ({sample}_trimmed_{1,2,s}.fastq.stats)
├── 02_assembly/                   # Assembler outputs (method-specific subdir)
│   └── megahit/{sample}/
│   │   └── final.contigs.fa.gz    # MEGAHIT assembled contigs
│   └── spades/{sample}/
│   │   └── contigs.fasta.gz       # metaSPAdes assembled contigs
│   │   └── (transcripts.fasta.gz) # rnaSPAdes assembled transcripts
│   └── metaflye_<read_type>/{sample}
│       └── assembly.fasta.gz      # metaFlye assembled contigs (long reads)
├── 03_assembly_filter/
│   └── {sample}/
│       └── {sample}.assembly_filter.fa  # Length-filtered contigs
├── 03_index/minibwa/              # minibwa index of filtered contigs
├── 03_mapping/minibwa/            # Sorted BAM files ({sample}.bam)
│                                  # (deleted after the last consuming step — binning/CoverM in b/m-style
│                                  #  modes or the unmapped-read extraction in cr/vcr-style modes — unless
│                                  #  keep_intermediate: true)
├── 04_bam_summary/                # Mapping statistics ({sample}.contigs_bam_stats.tsv)
├── 05_metawrap_binning/
│   └── {sample}/                  # Initial binning with three binners
│       ├── concoct_bins/          # CONCOCT bins (bin.N.fa)
│       ├── maxbin2_bins/          # MaxBin2 bins (bin.N.fa)
│       ├── metabat2_bins/         # metaBAT2 bins (bin.N.fa)
│       └── work_files/            # metaWRAP intermediate files
├── 06_metawrap_bin_refinement/
│   └── {sample}/                  # Consolidation/refinement of the three binner sets
│       ├── concoct_bins/ maxbin2_bins/ metabat2_bins/   # plus *.contigs and *.stats sidecars
│       ├── metawrap_50_10_bins/   # Refined bins passing the thresholds (bin.N.fa)
│       ├── metawrap_50_10_bins.stats/.contigs/.Chdb.tsv  # Refinement reports
│       ├── figures/               # Refinement summary figures
│       └── work_files/            # (dir name reflects config, e.g. 50% completeness / 10% contamination)
├── 07_contigs/                    # Prodigal gene prediction + depth on filtered contigs
│   ├── contigs_fna/               # Contig nucleotide sequences ({sample}.fna)
│   ├── contigs_protein_ffn/       # Predicted gene sequences ({sample}.ffn)
│   ├── contigs_protein_faa/       # Predicted protein sequences ({sample}.faa)
│   ├── contigs_ffn_length/        # Per-gene lengths (extensionless file named {sample})
│   ├── contigs_filter_depth/      # Per-contig depth (extensionless file named {sample})
│   └── contigs_ORF_loci/          # ORF loci on contigs (extensionless file named {sample})
└── 08_gather_bin_high_quality/    # High-quality MAGs selected by CheckM
    ├── bin_genome_fna_raw/        # HQ bins, file renamed only, contig headers untouched
    │   └── {sample}/{sample}_bin.N.fa
    ├── bin_genome_fna_renamed/    # HQ bins, file renamed AND contig headers renamed
    │   └── {sample}/{sample}_bin.N.fna
    ├── bin_genome_fa_renamed/     # Copy of the renamed genomes, input for gene calling
    │   └── {sample}/{sample}_bin.N.fa
    ├── bin_protein_faa_renamed/   # MAG protein sequences ({sample}/{sample}_bin.N.faa)
    ├── bin_protein_ffn_renamed/   # MAG gene sequences ({sample}/{sample}_bin.N.ffn)
    ├── bin_ffn_length/            # Per-gene lengths of MAGs
    ├── bin_ORF_loci/              # ORF loci on MAGs
    ├── checkm_stats/              # CheckM quality reports ({sample}.csv)
```

#### Contig mode (`c`) sample layout (extends `b`)

```
results/c-sample/
├── 01_taxonomy/kraken2/               # Contig taxonomy by kraken2 (or CAT_pack)
├── 02_function/kraken2/               # Contig function by kraken2 (or emapper)
├── 05_merge_tax_fun_depth/            # Files are named {sample}.xxx directly inside each subdirectory
│   ├── detail_for_each_seq/           # Per-sequence and per-taxon taxonomic-functional annotations ({sample}.detail_for_each_seq.tsv.gz)
│   ├── summary/                       # Sample summary TSVs ({sample}.summary_Tax.tsv, ...)
│   ├── KO_Taxa_RelativeAbundance/     # KO × taxon relative abundance matrices
│   ├── KO_Taxa_SumCoverage/           # KO × taxon depth matrices
│   ├── KO_depth/                      # Community KO depth ({sample}.KO_depth.tsv)
│   ├── Taxa_abundance/                    # Taxon depth ({sample}.Taxa_abundance.tsv)
│   └── Taxa_abundance/                    # Taxon depth ({sample}.Taxa_abundance.tsv)
├── 05_visualization/                  # Sample-level summary barplots (see below)
└── 06_module_completeness/            # Sample-level module completeness
```

Modes combining multiple resolution levels additionally produce under `05_merge_tax_fun_depth/`: `mode/` (records which merge level actually ran, `{sample}.mode.tsv`; in `cr`, `mc`, `mcr`, `vc`, `vcr`, `vm`, `vmc`, `vmcr`) and `MAG_contig_consistency/` (mMAG/vMAG–contig taxonomy consistency check; in `mc`, `mcr`, `vc`, `vcr`, `vm`, `vmc`, `vmcr`).

#### Read mode (`r`) sample layout

```
results/r-sample/
├── 00_sampled_trim/               # Subsampled/trimmed FASTQ
├── 01_sampled_trim_merge/         # Merged reads
├── 01_sampled_trim_merge_html/    # fastp reports
├── 01_sampled_trim_merge_json/    # fastp JSON summaries
├── 02_FragGeneScanRs/             # Gene prediction on reads
├── 02_FragGeneScanRs_stats/       # Gene prediction statistics
├── 03_taxonomy/{kraken2,bracken}/ # Read taxonomy
├── 04_function/kraken2/           # Read function (or emapper)
├── 05_merge_tax_fun_depth/        # Joint taxonomy-function matrices
├── 05_visualization/              # Sample-level summary barplots (see below)
└── 06_module_completeness/        # Sample-level module completeness
```

#### MAG mode (`m`) sample layout

```
results/m-sample/
├── 01_checkm2/{sample}/           # CheckM2 quality assessment
├── 02_taxonomy/gtdbtk/            # GTDB-Tk taxonomy (or BAT)
├── 03_genome_stat/                # MAG genome statistics
├── 04_protein_stat/               # MAG protein statistics
├── 05_drep_eachsample_99/         # Within-sample dereplication (ANI 99%)
├── 06_drep_eachsample_95/         # Within-sample dereplication (ANI 95%)
├── 07_coverm_eachsample/          # Per-sample MAG abundance
├── 08_function/emapper/           # eggNOG-mapper functional annotation (or kraken2)
├── 09_summary/{sample}/           # MAG summary tables: mMAG_summary_info.tsv (per genome)
│   │                              # and mMAG_species_summary_info.tsv (abundance aggregated by dRep_95 cluster)
├── 05_merge_tax_fun_depth/        # Joint taxonomy-function matrices
├── 05_visualization/              # Sample-level summary barplots (see below)
└── 06_module_completeness/        # Sample-level module completeness
```

#### Sample-level `05_visualization/` files

Generated for `c`, `r`, `m`, and integrated modes:

```
{sample}_contig_count.pdf    # Contig taxonomy/function counts
{sample}_contig_prop.pdf     # Contig taxonomy/function proportions
{sample}_ORFKO_count.pdf     # ORF-KO counts
{sample}_ORFKO_prop.pdf      # ORF-KO proportions
{sample}_ORFTaxKO_count.pdf  # ORF-Taxon-KO counts
{sample}_ORFTaxKO_prop.pdf   # ORF-Taxon-KO proportions
```

#### Virus mode (`v`) sample layout

```
results/v-sample/
├── 01_identification/{sample}/
│   ├── genomad_result/                 # geNomad viral identification
│   ├── vs_vb_result/                   # VirSorter2 & VIBRANT viral identification
│   └── genomad_vs_vb_result/           # geNomad | (VirSorter2 & VIBRANT) results
├── 03_index/minibwa/                   # Index of viral scaffolds ({sample}.final_virus.fna; subdir reflects the configured mapping method)
├── 03_mapping/minibwa/{sample}/        # Sorted BAMs of reads mapped to viral scaffolds
│                                       # (deleted after 09_coverm unless keep_intermediate: true)
├── 03_taxonomy/                        # CAT_pack viral contig taxonomy
├── 04_vRhyme/{sample}/                 # vRhyme viral binning
├── 05_CheckV/{sample}/                 # CheckV quality
├── 06_vclust/{sample}/                 # vClust clustering (clusters.tsv, ani.tsv, ...)
├── 07_iPHoP/{sample}/                  # Host prediction
├── 08_taxonomy/                        # BAT vMAG taxonomy
├── 09_coverm/{sample}/                 # vMAG abundance
├── 10_summary/{sample}/
│   ├── Virus_annotation_results.tsv    # Annotation summary
│   ├── vMAG_summary_info.tsv          # vMAG quality/taxonomy summary
│   ├── Host_prediction_to_genome_m90.csv
│   ├── Host_prediction_to_genus_m90.csv
│   └── Virus_genomes_files/            # Final vMAG .fna / .faa / .ffn
└── 11_visualization/{sample}/
    ├── Result_visualization_inputs/    # Tables used for plotting
    └── Result_visualization_outputs/   # PNG/PDF plots (virus stats, taxonomy, KO metabolism)
```


#### Integrated modes: merge-level fallback and the `mode` flag

Integrated modes (`mc`, `mcr`, `vc`, `vcr`, `vm`, `vmc`, `vmcr`) merge up to four evidence groups — reads (`r`), contigs (`c`), mMAGs (`m`), vMAGs (`v`) — into unified taxonomic–functional profiles (`05_merge_tax_fun_depth/`). The merge is **per-sample adaptive**: before merging, the workflow checks which groups actually have non-empty inputs for that sample, and drops any group that is missing or empty (e.g. a sample with no high-quality MAGs, or no detected viruses). The merge then falls back to the matching lower-level — so a `vmcr` run may execute a `vcr`, `mcr`, `cr`, or even `c`-level merge for a given sample.


### Project-level Outputs (when `--run-project`)

Project-level aggregation is written under `results/{mode}-project/`, with the project name as a sub-directory inside each output section (e.g. `results/{mode}-project/01_Taxa_Sample/{project}/`). The default project name is `test.metalink` unless `--project-name` is set.

#### Base mode (`b`) project summaries

```
results/b-project/
├── 01_cat_raw_stats/{project}/01_raw_reads_stats.tsv
├── 02_cat_trimmed_stats/{project}/02_clean_reads_stats.tsv
├── 03_contigs_stats/{project}/03_contigs_stats.tsv
├── 04_merge_reads_contigs_stats/{project}/04_reads_contigs_ffn_stats.tsv
└── 05_gather_project_bin_checkm_stats/{project}/05_checkm_high_quality_bins.csv
```

#### MAG modes (`m`, `vm`) project summaries

Project-level mMAG dereplication and summary:

```
results/m-project/
├── 01_drep_project_99/{project}/    # Cross-sample dRep at 99% ANI (dereplicated_genomes/, data_tables/, ...)
├── 01_drep_project_95/{project}/    # Cross-sample dRep at 95% ANI
└── 04_summary/{project}/mMAG_summary_info.tsv   # Project-wide mMAG summary table
```

#### Virus modes (`v` / `vm`) project summaries

Project-level vMAG dereplication, summary, and host prediction:

```
results/v-project/
├── 05_CheckV/{project}/                # Project-level CheckV summary of pooled viral genomes
├── 06_vclust/{project}/clusters.tsv   # Project-level vClust clusters (all samples pooled)
├── 07_iPHoP/                          # Project-level host prediction (custom/auto_drep MAGs strategies; `v`, `vm`, and `pmr --mags v`/`vm`)
├── 10_summary/                        # Project-level vMAG summary
└── 11_visualization/                  # Project-level plots
```


#### Taxonomy and Function Outputs

For contig (`c`), read (`r`), MAG (`m`), and integrated modes (`cr`, `mc`, `mcr`, `pmr`, `vm`, `vc`, `vmc`, `vcr`, `vmcr`), the project outputs include the following standard tables:

```
results/{mode}-project/
├── 01_Taxa_Sample/{project}/
│   ├── Taxa_Sample.long.tsv                       # Long-format taxonomy table
│   ├── Taxa_Sample.Taxa_ReadCount.wide.tsv        # Wide-format read counts (m/pmr modes: Taxa_Sample.Taxa_SumCoverage.wide.tsv, summed CoverM depth)
│   └── Taxa_Sample.Taxa_RelativeAbundance.wide.tsv
├── 01_Taxa_visualization/{project}/
│   ├── Barplot.Taxa.{D,P,C,O,F,G,S}.tsv
│   └── Barplot.Taxa.{D,P,C,O,F,G,S}.pdf
├── 02_KO_Sample/{project}/
│   ├── KO_Sample.long.tsv                         # Long-format KO table
│   └── KO_Sample.KO_SumCoverage.wide.tsv          # Wide-format KO coverages
├── 03_Module_Sample/{project}/                    # Community module tables (details see below)
│   ├── Completeness/                              # Shared completeness (computed once for both scales)
│   │   ├── mc_completeness.tsv
│   │   ├── mc_detected_KOs.tsv
│   │   ├── mc_by_pathway/                         # Per-pathway completeness tables ...
│   │   └── mc_by_pathway_heatmap/                 # ... each with a PDF/HTML heatmap directory
│   ├── RelativeAbundance/                         # RA-scale (× 1e6) module abundances
│   │   ├── mc_abundance.tsv
│   │   ├── mc_abundance_filtered_by_completeness.tsv
│   │   ├── log_scale_by_pathway/
│   │   ├── log_scale_by_pathway_filtered/
│   │   ├── log_scale_by_pathway_filtered_heatmap/
│   │   ├── avg_scale_by_pathway/                  # avg_scale_* only when ≥ 2 samples
│   │   ├── avg_scale_by_pathway_filtered/
│   │   └── avg_scale_by_pathway_filtered_heatmap/
│   └── SumCoverage/                               # Coverage-scale module abundances, same structure
├── 03_Fun_visualization/{project}/
│   ├── Barplot.Module.Level_2.{tsv,pdf}
│   ├── Barplot.Module.Level_3.{tsv,pdf}
│   ├── Barplot.Module.Module_Name.{tsv,pdf}
│   └── Barplot.Module.Energy.{Level_3}.{tsv,pdf}  # Per-Level_3 breakdown within "Energy metabolism"
├── 05_merge_tax_KO_depth/{project}/               # KO × taxon matrices
├── 05_merge_tax_Module_RelativeAbundance/{project}/   # Module × taxon matrices (RA scale)
├── 05_merge_tax_Module_SumCoverage/{project}/     # Module × taxon matrices (coverage scale)
└── 06_module_completeness/{project}/{rank}/       # Taxon-level module tables (see below)
```

##### Long-format table columns

**`01_Taxa_Sample/{project}/Taxa_Sample.long.tsv`**:

In fully integrated modes (e.g. `vmcr`), the long table contains both source-specific and final integrated values:

```
sample  Level  Taxonomy  Taxa_ReadCount_vmc  Taxa_ReadCount_read  Taxa_ReadCount_vmcr
        Taxa_RelativeAbundance_vmc  Taxa_RelativeAbundance_read  Taxa_RelativeAbundance_vmcr
        Taxa_ReadCount  Taxa_RelativeAbundance
```

- `Level`: taxonomic rank code (`D`, `P`, `C`, `O`, `F`, `G`, `S` = Domain/Phylum/Class/Order/Family/Genus/Species). (Strain-level / genome-level (`T`) outputs are also written, but are experimental.)
- `Taxa_ReadCount_*`: read counts from each analytical source (`_vmc` = vMAG + mMAG + contig, `_read` = reads only, `_vmcr` = integrated vmc + read).
- `Taxa_RelativeAbundance_*`: relative abundance percentage. The formula and denominator scope depend on the mode:

  | Modes | Column(s) | Formula | Denominator scope |
  |---|---|---|---|
  | `c`, `vc`, `vm`, `vmc` | `Taxa_RelativeAbundance` | `100 * taxon_bases / total_bases` (bases = `SeqLen * Avg_fold`; MAGs are expanded to their contigs first) | all contigs, including the `unclassified` category |
  | `r` | `Taxa_RelativeAbundance` | `100 * taxon_read_count / total_read_count` | all reads, including `unclassified` |
  | `m`, `pmr` | `Taxa_RelativeAbundance` | `100 * taxon_depth / total_depth` (CoverM `trimmed_mean` depth, fold units — not bases) | quality-filtered MAGs only (`unmapped` removed first) |
  | `cr`, `mcr`, `vcr`, `vmcr` | `Taxa_RelativeAbundance_<src>`, `Taxa_RelativeAbundance_read` (`<src>` = `_contig`/`_mc`/`_vc`/`_vmc`) | within-source percentage of that source's own total (bases; `_read` covers the reads that did not map to the assembly/MAG set) | each source separately — each column sums to 100 on its own |
  | `cr`, `mcr`, `vcr`, `vmcr` | `Taxa_RelativeAbundance_<mode>`, copied to the final `Taxa_RelativeAbundance` | `100 * (bases_src + bases_read) / (total_src + total_read)` — bases-weighted combination of the two sources, **not** the average of the two within-source values | both sources combined |
- `Taxa_ReadCount` / `Taxa_RelativeAbundance`: final integrated values used for downstream summaries.

In simpler modes (e.g. `c`), the table only contains the source-appropriate columns. For example, contig mode (`c`) uses:

```
sample  Level  Taxonomy  Taxa_ContigNum  Taxa_ReadCount  Taxa_RelativeAbundance
```

where `Taxa_ContigNum` is the number of contigs assigned to the taxon, and `Taxa_ReadCount` is calculated based on taxon_bases/mean_read_length.

##### ReadCount interpretation for 2nd-generation vs. 3rd-generation reads

`Taxa_ReadCount` is reported consistently across Metalink modes (except `m`/`pmr`, where the corresponding column is named `Taxa_SumCoverage`), but the way it is estimated differs between 2nd-generation (short) and 3rd-generation (long) reads:

- **2nd-generation short reads** (`paired`, `single`, `paired-rna`, `single-rna`): all reads have approximately the same length. Metalink therefore uses a single average read length per sample to convert sequence bases into estimated read counts: for paired data this is the average length of the R1 reads (from the trimmed R1 fastq stats), and `cr`-type modes assume R2 has the same length (paired reads count `2 * ReadLength` bases). As a result, `Taxa_ReadCount` values are directly comparable across modes (e.g. `c` vs. `cr`) for short-read data.

- **3rd-generation long reads** (`pacbio-raw`, `pacbio-corr`, `pacbio-hifi`, `nano-raw`, `nano-corr`, `nano-hq`): read lengths are highly variable. Metalink uses per-read lengths where they are available:
  - In read-only modes (`r`), each classified read contributes one count.
  - In contig-only modes (`c`), contig coverage is converted to an estimated read count using the **average read length of the sample** (for paired short reads, the average R1 length).
  - In integrated contig + read modes (`cr`, `mcr`, `vcr`, `vmcr`), the contig component is converted using the **average length of reads that mapped to contigs**, while the unmapped-read component is converted using the **average length of unmapped reads**. For long reads (when per-read lengths are passed via `--in_read_fna_length`), the bases in the numerator are accumulated from per-read lengths of individual unmapped reads, but the denominator is still the single average unmapped-read length. The final integrated `Taxa_ReadCount` is the sum of these two components.

- **MAG modes (`m`, `pmr`)**: this column is named **`Taxa_SumCoverage`** (not `Taxa_ReadCount`) in the per-sample `Taxa_abundance/{sample}.Taxa_abundance.tsv` and the project-wide `Taxa_Sample.Taxa_SumCoverage.wide.tsv` — it is **not** a read count at all, but the sum of CoverM `trimmed_mean` depths (`coverm_abundance`, in fold units) of the MAGs assigned to the taxon, with no bases/read-length conversion. Only the `c`/`cr`/`mc`/`mcr`/`vc`/`vm`/`vmc`/`vcr`/`vmcr`/`r` modes compute `round(Taxa_Bases / ReadLength)` (or per-read counts in `r`) and keep the `Taxa_ReadCount` column name. Consequently `Taxa_SumCoverage` from `m`/`pmr` is on a different scale and is **not directly comparable** with the `Taxa_ReadCount` estimates of the other modes.

- **`vm` mode**: the merged `Taxa_abundance/{sample}.Taxa_abundance.tsv` stays read-count-based (`Taxa_ReadCount`), but a coverage-based companion table `Taxa_abundance/{sample}.Taxa_SumCoverage.tsv` is written next to it — same layout as the `m`/`pmr` table (`sample, Level, Taxonomy, Taxa_ContigNum, Taxa_SumCoverage, Taxa_RelativeAbundance`), computed as the per-taxon sum of the input mMAG/vMAG `coverm_abundance` values. (Samples that fall back to the m-mode merge produce only the m-style `Taxa_abundance.tsv`, which is already coverage-based.)

Because of this difference, **absolute `Taxa_ReadCount` values for contigs/taxa are not directly comparable between `c` mode and `cr` mode when using 3rd-generation data**: `c` mode uses the whole-sample average length for contig conversion, whereas `cr` mode uses the mapped-read average length for the contig component.

**`02_KO_Sample/{project}/KO_Sample.long.tsv`**:

```
sample  KO  KO_ORFNum  KO_SumCoverage  KO_RelativeAbundance  ORFLen_median
```

- `KO`: KEGG Orthology identifier.
- `KO_ORFNum`: number of ORFs annotated with the KO.
- `KO_SumCoverage`: total coverage depth of the KO: sum of per-ORF `Scaled_Avg_fold`.
- `KO_RelativeAbundance`: GPM/TPM-style normalized relative abundance: the sum of per-ORF `ORFRA = 1e6 * Scaled_Avg_fold / sum(Scaled_Avg_fold)` per sample. Because per-base coverage depth already divides by gene length (`depth = bases / length`), TPM's gene-length normalization is embedded in the depth share itself — for uniform-length short reads this measure is equivalent to TPM. The bases-share (CPM-style) variant is exported separately as per-ORF `ORFRA_Bases` in `detail_for_each_seq`.
- `ORFLen_median`: `floor(median(ORFLen))` of ORFs annotated with the KO.

> [!NOTE]
> Two pre-filters are applied to ORFs **before** any KO statistics or normalization denominators are computed (all `05_merge_tax_fun_depth` rules): ORFs shorter than `minORFLen` (default 50 bp) are discarded, and with the default `KO_mode = ignore_multiple_KOs` any ORF carrying multiple KO annotations is discarded entirely (`keep_multiple_KOs`/`split_multiple_KOs` behave differently; see `KO` in the `detail_for_each_seq` column table below). Consequently, the "each sample sums to 1e6" property of `KO_RelativeAbundance` holds only over the ORFs that pass both filters.

##### Wide-format tables

Wide tables pivot the long table so that each column is a sample:

- `Taxa_Sample.Taxa_ReadCount.wide.tsv`: `Taxonomy  Level  <sample1>  <sample2>  ...` (in `m`/`pmr` modes this file is `Taxa_Sample.Taxa_SumCoverage.wide.tsv`, holding summed CoverM depth instead of read counts)
- `Taxa_Sample.Taxa_RelativeAbundance.wide.tsv`: same layout, relative abundance.
- `KO_Sample.KO_SumCoverage.wide.tsv`: `KO  <sample1>  <sample2>  ...`

##### CAMI-format taxonomic profile

When `--merge-cami-format-output` is set (or `CAMI_format_output: true` in the mode-specific config), Metalink writes an additional CAMI-style taxonomic profiling file under `Taxa_abundance/`:

```
results/{mode}-sample/05_merge_tax_fun_depth/Taxa_abundance/{sample}.Taxa_abundance.CAMI_format.tsv
```

In `vm` mode, a coverage-based companion profile is written alongside it — `Taxa_abundance/{sample}.Taxa_SumCoverage.CAMI_format.tsv` — whose `PERCENTAGE` is computed from the summed mMAG/vMAG genome depths instead of mapped bases (mirroring the two `Taxa_abundance` tables).

The file follows the CAMI taxonomic profiling format with a metadata header and a `@@TAXID RANK TAXPATH TAXPATHSN PERCENTAGE` body:

```
# Taxonomic Profiling Output
@SampleID:sample1
@Version:0.3.7
@Ranks:superkingdom|phylum|class|order|family|genus|species|strain
@TaxonomyID:/path/to/taxid_tax_lineage.tsv.gz
@__program__:metalink c
@@TAXID	RANK	TAXPATH	TAXPATHSN	PERCENTAGE
0	root	0	unclassified	13.9186156137398
1	root	1	root	86.0813843862602
793830449	superkingdom	793830449	Archaea	8.25497903276371
...
```

- `TAXID`: Metalink internal taxon identifier.
- `RANK`: taxonomic rank.
- `TAXPATH`: pipe-delimited ancestor taxids.
- `TAXPATHSN`: pipe-delimited scientific names.
- `PERCENTAGE`: relative abundance percentage.

**Availability**: available in `c`, `cr`, `r`, `vc`, `vcr`, `m`, `mc`, `mcr`, `pmr`, `vm`, `vmc`, and `vmcr` modes. It is disabled when the mMAG branch uses GTDB-Tk for taxonomy. Metalink automatically detects this conflict before execution: if `--merge-cami-format-output` is set together with `--mmags-taxonomy-classify-method gtdbtk` (the default in mMAG modes), CAMI output is forced to `False` and a warning is logged. Use `--mmags-taxonomy-classify-method BAT` to enable CAMI output in these modes.

### Joint Taxonomy–Function Outputs

#### Merge taxonomy and function tables
`05_merge_tax_fun_depth/` (sample level) and `05_merge_tax_KO_depth/` / `05_merge_tax_Module_RelativeAbundance/` / `05_merge_tax_Module_SumCoverage/`(project level) contain the joint profiles:

| File / Directory | Description |
|------------------|-------------|
| `detail_for_each_seq/{sample}.detail_for_each_seq.tsv.gz` | Per-ORF/sequence annotation table (columns documented below). **Only exported when `--merge-output-detail` is set** (or `output_detail: true` in the config); otherwise an empty file is created |
| `summary/{sample}.summary_Tax.tsv` | Taxon counts/abundance summary per source |
| `summary/{sample}.summary_ORF_KO.tsv` | KO count summary per source |
| `summary/{sample}.summary_ORF_Tax_KO.tsv` | Combined taxon–KO source summary |
| `KO_Taxa_RelativeAbundance/{sample}.KO_Taxa_RelativeAbundance.Taxa.{rank}.tsv` | KO relative abundance per taxon at a given rank |
| `KO_Taxa_SumCoverage/{sample}.KO_Taxa_SumCoverage.Taxa.{rank}.tsv` | KO depth per taxon at a given rank |
| `Taxa_abundance/{sample}.Taxa_abundance.tsv` | Taxon depth table |
| `Taxa_abundance/{sample}.Taxa_abundance.CAMI_format.tsv` | **CAMI-format taxonomic profile** (optional; see "CAMI-format taxonomic profile" above) |
| `KO_depth/{sample}.KO_depth.tsv` | Community KO depth table |

##### Columns of `detail_for_each_seq.tsv.gz`

The per-ORF detail table contains one row per ORF (per KO in `split_multiple_KOs` mode):

| Column | Meaning |
|--------|---------|
| `qaccver` | ORF identifier |
| `KO` | KEGG Orthology identifier; `Unclassified` if no KO was assigned. With `KO_mode=keep_multiple_KOs`, multiple KOs are joined by `_`; with `split_multiple_KOs`, the ORF is expanded into one row per KO |
| `Taxonomy` | Full GTDB-style taxonomy string of the sequence (`d__;p__;c__;...`) |
| `NewTax` | Novelty status of the taxon: `Known species`, `New species`/`New genus`/`New family`/... (novel at that rank), `Unclassified`, or `Root` |
| `taxid`, `Lineage` | Taxon id and its taxid lineage path; present only when CAMI output is enabled (`--merge-cami-format-output`) |
| `Scaled_Avg_fold` | Per-KO depth attributed to this ORF row: the sequence's average fold (contig/MAG coverm depth, or 1 per read) divided by the number of KO annotations in `split_multiple_KOs` mode |
| `ORFLen` | ORF length in bp |
| `ORFRA` | ORF relative abundance by coverage-depth share: `1e6 * Scaled_Avg_fold / sum(Scaled_Avg_fold)` — GPM/TPM-style (depth = bases/length embeds gene-length normalization); this is the measure aggregated downstream, e.g. `KO_RelativeAbundance = sum(ORFRA)` per KO |
| `ORFRA_Bases` | ORF relative abundance by covered-base share: `1e6 * Scaled_Avg_fold * ORFLen / sum(Scaled_Avg_fold * ORFLen)` — CPM-style (reads/bases share) |
| `Source` | Analytical source of the sequence: `Read`, `Contig`, `mMAG`, `vMAG`, etc. |

##### Key columns of the other merge tables

**`Taxa_abundance/{sample}.Taxa_abundance.tsv`** (per-taxon abundance, one row per taxon × rank):

| Column | Meaning |
|--------|---------|
| `Level` | Taxonomic rank code: `D/P/C/O/F/G/S` = Domain/Phylum/Class/Order/Family/Genus/Species |
| `Taxonomy` | Full taxonomy path of the taxon up to that rank |
| `Taxa_ContigNum` | Number of contigs/scaffolds assigned to the taxon (counting contigs/scaffolds inside the genomes of that taxon — e.g. tens of thousands per domain, not the genome count). Used consistently across all abundance tables: `c`-family modes, every mode's `Taxa_abundance.tsv` (including `m`/`pmr` and the read-count table of `vm`), and `vm`'s coverage table `Taxa_SumCoverage.tsv` (whose per-genome depth sums are still computed at genome level; only the count column is contig-based). Not present in integrated-mode (`cr`/`mcr`/`vcr`/`vmcr`) tables |
| `Taxa_ReadCount` | Read count assigned to the taxon. For contig/read sources this is a coverage-derived estimated read count (`round(Taxa_Bases / ReadLength)`); for MAG sources in `m`/`pmr` the column is instead named `Taxa_SumCoverage` and holds the summed CoverM `trimmed_mean` depth (fold), not a read count. See "ReadCount interpretation" above for the exact conversion per mode and read type |
| `Taxa_ReadCount_<src>`, `Taxa_ReadCount_read`, `Taxa_ReadCount_<mode>` | Integrated modes (`cr`/`mcr`/`vcr`/`vmcr`) only: per-source read-count estimates (`<src>` = `_contig`/`_mc`/`_vc`/`_vmc`; `_read` = reads not mapped to the assembly/MAG set) and their sum. The final `Taxa_ReadCount` equals `Taxa_ReadCount_<mode>` |
| `Taxa_RelativeAbundance` | Relative abundance (%) of the taxon within the sample at that rank; sums to 100 per level. For MAG sources in `m`/`pmr` it is a share of CoverM `trimmed_mean` depth (not of bases), and the denominator covers only the quality-filtered MAG subset. In integrated modes (`cr`/`mcr`/`vcr`/`vmcr`) it is the integrated, bases-weighted combination of the two sources (see the per-mode table in "Long-format table columns") |
| `Taxa_RelativeAbundance_<src>`, `Taxa_RelativeAbundance_read`, `Taxa_RelativeAbundance_<mode>` | Integrated modes only: within-source percentages (each source normalized to its own total, summing to 100 separately) and the integrated value. The final `Taxa_RelativeAbundance` equals `Taxa_RelativeAbundance_<mode>`; see the per-mode table in "Long-format table columns" for formulas |

**`KO_depth/{sample}.KO_depth.tsv`** (per-KO summary, one row per KO): same columns as `KO_Sample.long.tsv` above (`KO`, `KO_ORFNum`, `KO_SumCoverage`, `KO_RelativeAbundance`, `ORFLen_median`); see there for definitions. `KO_RelativeAbundance` sums to 1e6 per sample (over ORFs passing the `minORFLen`/`KO_mode` pre-filters, see the note above).

**`detail_for_each_seq/{sample}.detail_for_each_{rank}.tsv`** (KO × taxon long table, one row per KO × taxon at the given rank):

| Column | Meaning |
|--------|---------|
| `sample`, `KO`, `{rank}` | Grouping keys (taxon name is the full path up to that rank) |
| `ORF_Num` | Number of ORFs annotated with the KO in that taxon |
| `Taxa_KO_SumCoverage` | Sum of per-ORF `Scaled_Avg_fold` for the KO in that taxon |
| `Taxa_KO_RelativeAbundance` | Sum of per-ORF `ORFRA` for the KO in that taxon |

**`KO_Taxa_SumCoverage/{sample}.KO_Taxa_SumCoverage.Taxa.{rank}.tsv`** and **`KO_Taxa_RelativeAbundance/{sample}.KO_Taxa_RelativeAbundance.Taxa.{rank}.tsv`** (wide matrices): rows are KOs, columns are `{sample}---{taxon}`; cell values are `Taxa_KO_SumCoverage` (depth) or `Taxa_KO_RelativeAbundance` (relative abundance) respectively. Taxa with fewer than `min_KO_num_in_one_taxon` KOs are dropped from the wide matrices. Note that the `detail_for_each_seq/` long tables are exported **before** the `min_KO_num_in_one_taxon` filter and therefore retain all taxa; the filter only affects the wide matrices.

#### Metabolic pathway/module abundance tables (mclink outputs)

Metalink uses the `mclink` R package to compute metabolic pathway/module abundance and completeness from KO profiles. The results are written as two kinds of module-level tables that differ in what the columns represent.

**How module scores are calculated**

For each KEGG module, mclink aggregates KO-level values according to the module definition string:

- Spaces and plus signs (`+`) between KOs are treated as **AND** operators (required components).
- Commas (`,`) are treated as **OR** operators (alternative components).
- Default aggregation:
  - `mean` for plus-separated KOs,
  - `max` for comma-separated KOs in completeness mode (hard-coded; not affected by the `comma_scale_method` configuration),
  - `sum` for comma-separated KOs in abundance mode.

KO relative abundances (`KO_RelativeAbundance`, `Taxa_KO_RelativeAbundance`) are GPM/TPM-like normalized (coverage-depth share × 1e6; since depth = bases/length, gene-length normalization is embedded, so for uniform-length short reads this equals TPM), so each sample sums to 1e6 (over ORFs passing the `minORFLen`/`KO_mode` pre-filters). Taxon relative abundances (`Taxa_RelativeAbundance_*`) are per-sample percentages that sum to 100; the exact formula and denominator scope depend on the mode — see the per-mode table in "Long-format table columns".

**Community module tables** (project level, with one column per sample: `03_Module_Sample/{project}/`, with three parallel subdirectories: `Completeness/` and `RelativeAbundance/` + `SumCoverage/` for the two abundance scales):

- Rows are KEGG modules.
- Columns are `Module_Entry`, `Level_2`, `Level_3`, `Module_Name`, `Definition`, and one column per sample.
- Files:
  - `Completeness/mc_completeness.tsv`: module completeness scores (0–1); computed once and shared by both abundance scales. KO filtering (minimum sample number / prevalence / abundance) is always derived from the relative-abundance matrix and then applied to both scales, so the two scales contain exactly the same KOs and samples and differ only in the base values.
  - `RelativeAbundance/mc_abundance.tsv` / `SumCoverage/mc_abundance.tsv`: module abundance aggregated from `KO_RelativeAbundance` / `KO_SumCoverage` respectively.
  - `RelativeAbundance/mc_abundance_filtered_by_completeness.tsv` / `SumCoverage/mc_abundance_filtered_by_completeness.tsv`: completeness-filtered abundance (see below); the SumCoverage version is the input of the module barplots.
  - `Completeness/mc_detected_KOs.tsv`: space-separated KOs detected for each module in each sample.
- These are community-level profiles: each module score is calculated from the KO profile of the whole sample.
- Scaled/transformed per-pathway outputs (`mc_by_pathway*` under `Completeness/`, from the completeness run; `log_scale_*` and `avg_scale_*` inside each scale subdirectory, e.g. `RelativeAbundance/log_scale_by_pathway/`):
  - `Completeness/mc_by_pathway/`: per-pathway module completeness tables, with `mc_by_pathway_heatmap/` (PDF + interactive HTML heatmaps).
  - `log_scale_by_pathway/`: per-pathway module abundance transformed by `log10(value + 1)` to reduce the dynamic range.
  - `avg_scale_by_pathway/`: for each module, the abundance is converted to `log2(sample_value / row_mean)`; zeros are preserved as zero. This highlights whether a sample has above-average (`> 0`) or below-average (`< 0`) abundance for each module. Generated only when there are ≥ 2 samples (≥ 2 taxa for the taxon-level tables).
  - `log_scale_by_pathway_filtered/` and `avg_scale_by_pathway_filtered/` (with their `*_heatmap/` directories): completeness-filtered versions of the above. Filtering is **per cell** (module × sample): any cell whose module completeness is below the threshold (default 0.5) is set to zero — the module row itself is kept. This removes low-confidence values that are unlikely to be truly present in the sample, so the filtered tables are the recommended input for heatmaps and downstream pathway-abundance analyses. The `avg_scale_by_pathway_filtered/` tables and their heatmaps are likewise only generated with ≥ 2 samples.
- Controlled by the project-level `03_mclink_sample` parameters, e.g. `--mclink-min-ko-num-per-sample`, `--mclink-min-ko-prevalence-per-sample`, `--mclink-min-ko-abundance-per-sample`, `--mclink-completeness-per-sample`. The KO filter thresholds (minimum sample number / prevalence / abundance) are evaluated on relative abundances re-normalized to 1e6 per sample.

**Taxon-level module tables** (`06_module_completeness/{sample}/{rank}/` at sample level and `06_module_completeness/{project}/{rank}/` at project level; each contains three parallel subdirectories: `Completeness/`, `RelativeAbundance/` and `SumCoverage/`):

- Rows are KEGG modules.
- Columns are `Module_Entry`, `Level_2`, `Level_3`, `Module_Name`, `Definition`, and one column per taxon at the chosen rank (e.g. `Phylum`).
- Values are module completeness or abundance for each taxon, calculated from the KO relative abundance (or KO coverage, in the `SumCoverage/` tables) attributed to that taxon. The project-level tables are computed from the project-level mean Taxa-KO relative abundance / coverage across samples; taxon filtering is always derived from the relative-abundance matrix and applied to both scales.
- These tables answer: "which metabolic pathways are complete/abundant in which taxa?"
- Controlled by `06_mclink_taxon` parameters, e.g. `--mclink-sample-min-ko-num-per-taxon`, `--mclink-sample-min-ko-abundance-per-taxon`, `--mclink-sample-topn-taxon`, `--mclink-sample-completeness-per-taxon`.


### Other Summary Tables (mMAGs and vMAGs)

The final per-genome summary tables are the main entry point for downstream analysis. Their columns are described below.

#### `mMAG_summary_info.tsv` — mMAG summary

One row per prokaryotic MAG (bin). Found at `results/m-sample/09_summary/{sample}/` (per sample) and `results/m-project/04_summary/{project}/` (project-wide, first column named `ID`; the two dRep columns are recomputed from cross-sample dereplication).

| Column | Meaning |
|--------|---------|
| `Prokaryotes` / `ID` | MAG (bin) ID, e.g. `{sample}_bin.1`. |
| `dRep_99_representative` | Representative genome of the dRep ANI 99% cluster (strain-level dereplication) this MAG belongs to |
| `dRep_95_representative` | Representative genome of the dRep ANI 95% cluster (species-level dereplication) |
| `coverm_abundance` | CoverM trimmed-mean depth of this MAG in the sample |
| `relative_abundance` | `coverm_abundance` divided by the sum over all MAGs in the sample (0–1) |
| `closest_genome_reference` | Closest GTDB reference genome accession (GTDB-Tk only; absent with `CAT_pack` taxonomy) |
| `Raw_Taxonomy` | Raw classification string from GTDB-Tk (default) or BAT |
| `Domain`–`Species` | Seven canonical ranks parsed from `Raw_Taxonomy`; unclassified ranks are filled as e.g. `p__unclassified <Domain> phylum` |
| `NewTax` | Novelty flag: `Known species`, `Unclassified species/genus/family/...` (first unclassified rank), `Root`, or `Unclassified` |
| `Completeness2` / `Contamination2` | CheckM2 completeness / contamination (%) |
| `Completeness_Model_Used` / `Translation_Table_Used` | CheckM2 model used for estimation / translation table used for gene prediction |
| `Completeness1` / `Contamination1` / `strain_heterogeneity` | CheckM completeness / contamination / strain heterogeneity (%) |
| `GC (%)` … `N90`, `Total scaffolds/length (>= N bp)` | MAG scaffold statistics: GC content, scaffold count, total/N bases, average/longest/shortest scaffold length, L50/N50/L75/N75/L90/N90, and counts/lengths of scaffolds ≥ 1 Mb / 100 kb / 10 kb / 1 kb |
| `Protein count` / `Coding bases` / `Coding density (%)` / `Average gene length` | Statistics of prodigal-predicted proteins: protein count, total coding bases (bp), coding density (100 × coding bases / total bases), and mean gene length (bp) |

Variants of this table:

| File | Description |
|------|-------------|
| `mMAG_summary_info.filter.tsv` | Rows of `mMAG_summary_info.tsv` that pass the configured completeness/contamination thresholds; first column renamed to `genome`. Written by the `05_merge_tax_fun_depth` step of each integrated mode at `results/<mode>-sample/05_merge_tax_fun_depth/{sample}.mMAG_summary_info.filter.tsv` (`<mode>` = `m`/`mc`/`mcr`/`vm`/`vmc`/`vmcr`); in `pmr` mode the project-level equivalent is `results/pmr_m-project/02_quality_gate/{project}.mMAG_summary_info.filter.tsv` (`pmr --mags m`) or `results/pmr_vm-project/02_quality_gate/{project}.mMAG_summary_info.filter.tsv` (`pmr --mags vm`) |
| `mMAG_species_summary_info.tsv` (same directory) | Aggregated by `dRep_95_representative`: one row per species-level cluster, taxonomy/quality/statistic columns taken from the representative genome, with `sum_coverm_abundance` and `sum_relative_abundance` summed over all cluster members |

#### `vMAG_summary_info.tsv` — vMAG summary

One row per vMAG (a vRhyme bin or an unbinned viral scaffold). Found at `results/v-sample/10_summary/{sample}/` and `results/v-project/10_summary/{suffix}/{project}/`.

| Column | Meaning |
|--------|---------|
| `Virus` | vMAG ID: vRhyme bin name (`{sample}_vRhyme_bin_N`) or unbinned viral scaffold name (`{sample}_vRhyme_unbinned_N`) |
| `cluster` | Representative genome of the vClust cluster (default ANI ≥ 95%, approximately a viral species) |
| `scaffold` | Comma-separated IDs of the contigs/scaffolds composing this vMAG. Fragment scaffolds (sub-regions of a longer host contig) carry the parent-contig coordinates in their IDs: `{parent}\|provirus_S_E` (geNomad) or `{parent}\|vs2region_S_E` (VirSorter2/CheckV) — S/E are 1-based inclusive positions on the parent contig |
| `lytic_state` | Lytic/lysogenic state inferred from geNomad (provirus topology + integrase) or VirSorter2 evidence: `lytic_virus`, `lysogenic_virus`, or scaffold-level `lytic_scaffold` / `lysogenic_scaffold` / `integrated_prophage` |
| `checkv_quality` | CheckV quality tier: `Complete`, `High-quality`, `Medium-quality`, `Low-quality`, or `Not-determined` |
| `miuvig_quality` | CheckV quality tier per the MIUViG standard |
| `genome_size` / `scaffold_num` / `protein_count` | Total length (bp), scaffold count, and predicted protein (CDS) count of the vMAG |
| `AMG_KOs` | Auxiliary metabolic genes as `KO(copy_number)` pairs, e.g. `K00018(3);K01953(4)`; filtered to remove likely host contamination |
| `completeness` / `completeness_method` | CheckV-estimated genome completeness (%) and the estimation method |
| `coverm_abundance` | CoverM trimmed-mean depth of this vMAG in the sample (× coverage) |
| `relative_abundance` | `coverm_abundance` divided by the sum over all vMAGs in the sample (0–1) |
| `lineage` | Full CAT/BAT classification lineage (with Metalink taxids) |
| `superkingdom`–`species` | Classification ranks from BAT |

Variants of this table:

| File | Description |
|------|-------------|
| `vMAG_summary_info.filter.tsv` | Rows of `vMAG_summary_info.tsv` whose `checkv_quality` passes the configured threshold; first column renamed to `genome`. Written by the `05_merge_tax_fun_depth` step of each integrated mode at `results/<mode>-sample/05_merge_tax_fun_depth/{sample}.vMAG_summary_info.filter.tsv` (`<mode>` = `vc`/`vcr`/`vm`/`vmc`/`vmcr`); in `pmr` mode the project-level equivalent is `results/pmr_v-project/02_quality_gate/{project}.vMAG_summary_info.filter.tsv` (`pmr --mags v`) or `results/pmr_vm-project/02_quality_gate/{project}.vMAG_summary_info.filter.tsv` (`pmr --mags vm`) |
| `vMAG_species_summary_info.tsv` | Aggregated by vClust `cluster` (one row per viral species-level cluster); `coverm_abundance` summed over cluster members, other columns from a cluster member. Produced by the `v` chain at `results/v-sample/10_summary/{sample}/` and `results/v-project/10_summary/{suffix}/{project}/`. The `pmr`-mode equivalent is the `04_summary` abundance matrix (`vMAG_species_abundance_matrix.tsv` for `--mags v`, `mMAG_vMAG_species_abundance_matrix.tsv` for `--mags vm`) |

#### Cross-sample species abundance matrices

`mMAG_species_abundance_matrix.tsv` (`pmr --mags m`), `vMAG_species_abundance_matrix.tsv` (`pmr --mags v`), and `mMAG_vMAG_species_abundance_matrix.tsv` (`pmr --mags vm`) are built by merging the per-sample CoverM tables and annotating them from the corresponding `02_quality_gate` filter table (`{project}.mMAG_summary_info.filter.tsv` and/or `{project}.vMAG_summary_info.filter.tsv`). The annotation depends on `--mags`: for `m`/`v` (single genome type) the matrix is a self-contained species-level summary carrying **all** filter columns (with canonical `Domain`–`Species` prepended when the table uses CAT/BAT rank names, scores stripped); for `vm` (combined) only the seven rank columns `Domain`–`Species` are kept so the mMAG/vMAG rows align on the same columns (no sparse union). In both cases the rank columns are normalized with the same `rename_GTDB_taxonomy_dataframe` + `rename_strain_dataframe` pipeline used by the combine scripts (unclassified ranks are filled as e.g. `p__unclassified d__Viruses phylum`, and vMAGs with no CAT classification are reported as `d__Viruses`). The last columns are one `{sample} Trimmed Mean` column per sample holding the CoverM trimmed-mean depth of that representative genome in that sample (0 if unmapped). Rows are sorted by taxonomy (Domain → … → Species), with genomes naturally ordered by their trailing number within a group (e.g. `bin.2` before `bin.10`):

| File | Location | Rows (representative genomes) |
|------|----------|-------------------------------|
| `mMAG_species_abundance_matrix.tsv` | `results/pmr_m-project/04_summary/{project}/` (`pmr --mags m`) | Project-level dRep ANI 95% dereplicated mMAGs (`results/m-project/01_drep_project_95/{project}/dereplicated_genomes/`), against which all samples' reads are mapped |
| `vMAG_species_abundance_matrix.tsv` | `results/pmr_v-project/04_summary/{project}/` (`pmr --mags v`) | Project-level vClust cluster representative genomes |
| `mMAG_vMAG_species_abundance_matrix.tsv` | `results/pmr_vm-project/04_summary/{project}/` (`pmr --mags vm`) | Combined project-level dRep ANI 95% mMAGs and vClust vMAG cluster representatives, quantified in a single mapping pass |

### Visualization Outputs

| File pattern | Description |
|--------------|-------------|
| `Barplot.Taxa.{D,P,C,O,F,G,S}.tsv` / `.pdf` | Stacked barplot data and plot per taxonomic rank (7 ranks, Domain–Species); input is `Taxa_Sample.Taxa_ReadCount.wide.tsv` (`Taxa_Sample.Taxa_SumCoverage.wide.tsv` for `m`/`pmr` modes) |
| `Barplot.Module.Level_2.{tsv,pdf}` | Level-2 module barplot |
| `Barplot.Module.Level_3.{tsv,pdf}` | Level-3 module barplot |
| `Barplot.Module.Module_Name.{tsv,pdf}` | Module-level barplot |
| `Barplot.Module.Energy.{Level_3}.{tsv,pdf}` | Module barplots broken down by Level_3 category, generated only within Level_2 "Energy metabolism"; module barplots use `mc_abundance_filtered_by_completeness.tsv` as input |
| `*_heatmap/` directories | PDF and interactive HTML heatmaps of module abundance/completeness |

All `Barplot.*.tsv` files contain within-sample relative abundance percentages (each sample column sums to 100 over the plotted categories). Categories are ranked by their summed relative abundance across samples; the top `TopN` categories are kept, the remainder are merged into a grey `Others` category, categories that are zero in all samples are dropped, and `TopN` is capped at 29.

---

## Optimizations & Patches to Third-Party Tools

Metalink incorporates several optimizations and patches to third-party tools to improve robustness, compatibility, and computational efficiency. These modifications are applied automatically during workflow execution.

### MetaWRAP Binning Module

Metalink patches the MetaWRAP binning module to address two limitations in the original implementation:

#### 1. FASTQ.GZ File Recognition

The original MetaWRAP binning scripts only recognize uncompressed `.fastq` files and fail to properly handle `.fastq.gz` compressed inputs. Metalink patches the binning wrapper to:
- Auto-detect `.fastq.gz` and `.fq.gz` extensions in addition to `.fastq` and `.fq`
- Pass compressed inputs directly without requiring intermediate decompression

This eliminates the need for manual decompression of large FASTQ files, reducing disk I/O and temporary storage requirements.

#### 2. CheckM / dRep Result Consistency & Avoidance of Redundant CheckM Runs

In the original MetaWRAP + dRep workflow, dereplicated MAGs from dRep are often re-evaluated by CheckM, leading to redundant quality assessment. Metalink optimizes this by:
- **Caching CheckM results**: CheckM quality scores computed during the initial binning phase are cached
- **Representative selection & result propagation**: For each dRep cluster, Metalink takes the representative genome from the dRep winner table (`Wdb.csv`), i.e. the highest-scoring genome per cluster; if `Wdb.csv` is missing or invalid, it falls back to the lexicographically largest genome ID. CheckM/CheckM2 statistics (completeness, contamination, strain heterogeneity) are merged directly onto each genome's own row by genome ID, and dRep reuses the pre-computed CheckM results via `--genomeInfo`, avoiding re-running CheckM on the dereplicated set

This optimization reduces the computational cost of the binning-to-dereplication pipeline by **avoiding redundant CheckM evaluations** on dereplicated MAG sets, particularly impactful for large-scale projects with hundreds to thousands of MAGs.

---


## FAQ & Troubleshooting

### Q: How do I choose the kraken2 confidence threshold for taxonomic classification?

Metalink applies kraken2 confidence thresholds for **taxonomic classification** at three analytical levels:

| Parameter | Default | Applies to |
|-----------|---------|------------|
| `--reads-taxonomy-kraken2-confidence` | `0.6` | Read taxonomy (`r` mode) |
| `--contigs-taxonomy-kraken2-confidence` | `0.7` | Contig taxonomy (`c`, `cr`, `mc`, `mcr`, etc.) |
| `--<mode>-taxonomy-kraken2-confidence` | `0.7` | Unmapped read taxonomy in integrated modes |

**General principle**: Lower confidence thresholds increase sensitivity (recall) at the cost of precision, which is typically preferable for complex communities where strict thresholds cause substantial under-classification. The optimal choice depends on sample complexity, database composition, and the scientific question (e.g., precision vs. recall priority).

For systematic benchmarking of kraken2 confidence thresholds across simulated communities, see:

> Wright RJ, Comeau AM, Langille MGI. From defaults to databases: parameter and database choice dramatically impact the performance of metagenomic taxonomic classification tools. *Microbial Genomics*. 2023;9(3):000949. https://doi.org/10.1099/mgen.0.000949

### Q: How should I interpret abundance and relative abundance — are they cell counts or molecule counts?

It depends on the metric. Sequencing signal is proportional to copy number × sequence length, so whether an abundance value can be read as a "count" depends on whether the length has been normalized out:

- **Read counts / read bases** — the `Taxa_ReadCount` taxonomy abundance in the merged `Taxa_abundance` tables (produced by `r`, `c`, `cr`, `mc`, `mcr`, `vc`, `vcr`, `vm`, `vmc`, `vmcr`): in `r` this is the actual number of classified reads; elsewhere it is an *equivalent* read count computed as mapped bases ÷ read length, where mapped bases = sequence length × depth. Either way the value is proportional to cell number × genome length ÷ read length: a *sequence amount*, **not** a cell count — species with larger genomes get more reads at the same cell number. Valid for comparing the same taxon across samples (after normalizing for sequencing depth); biased when comparing different taxa within a sample.
- **mMAG / vMAG coverage** (CoverM trimmed-mean depth, `coverm_abundance`): per-genome depth computed as mapped bases ÷ genome length. Because genome length is normalized out, depth ∝ cellular abundance, and `relative_abundance` ≈ the cell fraction among the recovered genomes. It appears in:
  - `vMAG_summary_info.tsv` (`v-sample/`, all virus modes) — per-vMAG `coverm_abundance`;
  - `mMAG_summary_info.tsv` (`m-project/`, `m` and `pmr` modes) — per-mMAG `coverm_abundance`;
  - the three species tables (`mMAG_species_abundance_matrix.tsv`, `vMAG_species_abundance_matrix.tsv`, `mMAG_vMAG_species_abundance_matrix.tsv`) — one `{sample} Trimmed Mean` column per sample;
  - the `Taxa_SumCoverage` columns of the `m`/`pmr` taxonomy abundance tables — per-taxon sums of MAG depths.

  The caveat: the denominator only includes recovered genomes (the unassembled/unbinned part of the community is not counted), and the interpretation assumes ~1 genome copy per cell (a good approximation for bacteria/archaea; replicating viruses may carry multiple copies).
- **Function (KO) depth** (`KO_depth`, `KO_Taxa_SumCoverage`, etc.): gene-length-normalized, so ∝ the **total copy number** of that function in the sample, and the relative abundance ≈ copy-number share. This also holds for read-mode function: each read-derived ORF contributes `ORFLen / (3 × len_median)` (normalized by the KO reference median length) instead of 1, putting it on the same coverage scale. These are community-wide copy numbers — for per-cell copy numbers, divide by the cellular abundance (depth) of the carrier taxa.

Rule of thumb: any Metalink abundance given as depth/coverage is length-normalized and interpretable as a count (cells for genomes, copies for functions); raw read counts/bases are sequence amounts only.

### Q: Can I run Metalink on a single sample?

Yes. Metalink works with any number of samples (≥1). The sample sheet TSV only needs one line.

### Q: Can I mix paired-end and single-end samples in one run?

No. All samples in a single run must use the same read type.

### Q: How much memory does Metalink need?

Memory requirements vary by mode and dataset size:
- **Assembly**: ~100–800 GB RAM depending on dataset size
- **Kraken2 taxonomic database**: ~16/64/256 GB RAM for the `kraken2_16GB`/`kraken2_64GB`/`kraken2_256GB` variants (kraken2 loads the database into RAM; the bundled bracken files only add disk usage, not kraken2 RAM); `kraken2_raw` needs hundreds of GB
- **Kraken2 functional database**: ~16 GB RAM for the `kraken2_16GB` variant; `kraken2_raw` needs RAM comparable to its on-disk size

See the database size table in [Quick setup](#quick-setup) for the exact disk space each database variant occupies.

### Q: Can I use pre-trimmed reads, pre-assembled contigs, or pre-binned MAGs?

Yes. Use `--skip-level 1` (pre-trimmed reads), `--skip-level 2` (pre-assembled contigs in the `contig_file` TSV column), or `--skip-level 3` (pre-binned MAGs in the `mMAG_dir` TSV column). See [Skip-Level: Reusing Pre-Computed Results](#skip-level-reusing-pre-computed-results).

### Q: What if my run fails on the cluster?

1. Check `{work_dir}/logs/` for error messages
2. Check cluster-specific logs (e.g., Slurm `.err` files)
3. Re-run with `--ri`/`--rerun-incomplete` to resume
4. If a specific tool failed, adjust its parameters (e.g., increase `--binning-max-memory-gb`)

### Q: What read types and assemblers are compatible?

Metalink validates `--read-type` against `--assembly-method` at runtime (skip-level 0 and 1):

| `--read-type` | Allowed `--assembly-method` | Error if violated |
|---------------|----------------------------|-------------------|
| `paired` | `MEGAHIT`, `metaSPAdes` | `ValueError` |
| `single` | `MEGAHIT` | `ValueError` |
| `paired-rna`, `single-rna` | `rnaSPAdes` | `ValueError` |
| `pacbio-raw`, `pacbio-corr`, `pacbio-hifi`, `nano-raw`, `nano-corr`, `nano-hq` | `metaFlye` | `ValueError` |

At skip-level `2` or `3` (assembly is skipped), this validation is bypassed.

### Q: How often is the Metalink database updated?

Metalink databases are updated on a **two-year cycle**. Each update tracks the taxonomy-system changes of **GTDB** (prokaryotes) and **ICTV** (viruses) and rebuilds the Metalink taxonomy database from them. The **LCFH** function database is rebuilt from the functions carried by the genomes in the taxonomy database.

### Q: Where can I get help?

- **Bug reports**: [GitHub Issues](https://github.com/LiuyangLee/Metalink/issues)
- **Snakemake help**: [Stack Overflow (snakemake tag)](https://stackoverflow.com/questions/tagged/snakemake)

---

## Citation

If you use Metalink in your research, please cite:

> Li et al. (202X). Metalink: Joint taxonomic and functional profiling of metagenomes and metatranscriptomes across resolution scales. *(under review)*.

Additionally, please cite the underlying tools and databases used in your analysis:

### Workflow Engine

- Köster, J. & Rahmann, S. Snakemake—a scalable bioinformatics workflow engine. *Bioinformatics* 28, 2520–2522 (2012).

### Read Preprocessing & Assembly

- Chen, S. Ultrafast one-pass FASTQ data preprocessing, quality control, and deduplication using fastp. *iMeta* 2, e107, doi:10.1002/imt2.107 (2023).
- Li, D. et al. MEGAHIT: an ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph. *Bioinformatics* 31, 1674–1676, doi:10.1093/bioinformatics/btv033 (2015).
- Nurk, S., Meleshko, D., Korobeynikov, A. & Pevzner, P. A. metaSPAdes: a new versatile metagenomic assembler. *Genome Research* 27, 824–834, doi:10.1101/gr.213959.116 (2017).
- Bushmanova, E., Antipov, D., Lapidus, A. & Prjibelski, A. D. rnaSPAdes: a de novo transcriptome assembler and its application to RNA-Seq data. *GigaScience* 8, giz100, doi:10.1093/gigascience/giz100 (2019).
- Kolmogorov, M. et al. metaFlye: scalable long-read metagenome assembly using repeat graphs. *Nature Methods* 17, 1103–1110, doi:10.1038/s41592-020-00971-x (2020).

### Read Mapping

- Vasimuddin, M., Misra, S., Li, H. & Aluru, S. Efficient architecture-aware acceleration of BWA-MEM for multicore systems. *2019 IEEE International Parallel and Distributed Processing Symposium (IPDPS)*, 314–324.
- Li, H. Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. *arXiv* (2013).
- Li, H. & Homer, N. Fast genomic read alignment with minibwa. *arXiv* (2026).
- Sahlin, K. Strobealign: flexible seed size enables ultra-fast and accurate read alignment. *Genome Biology* 23, 260, doi:10.1186/s13059-022-02831-7 (2022).
- Danecek, P. et al. Twelve years of SAMtools and BCFtools. *GigaScience* 10, giab008, doi:10.1093/gigascience/giab008 (2021).

### MAG Binning & Quality Assessment

- Uritskiy, G. V., DiRuggiero, J. & Taylor, J. MetaWRAP—a flexible pipeline for genome-resolved metagenomic data analysis. *Microbiome* 6, 158, doi:10.1186/s40168-018-0541-1 (2018).
- Kang, D. D. et al. MetaBAT 2: an adaptive binning algorithm for robust and efficient genome reconstruction from metagenome assemblies. *PeerJ* 7, e7359, doi:10.7717/peerj.7359 (2019).
- Wu, Y. W., Simmons, B. A. & Singer, S. W. MaxBin 2.0: an automated binning algorithm to recover genomes from multiple metagenomic datasets. *Bioinformatics* 32, 605–607, doi:10.1093/bioinformatics/btv638 (2016).
- Alneberg, J. et al. Binning metagenomic contigs by coverage and composition. *Nature Methods* 11, 1144–1146, doi:10.1038/nmeth.3103 (2014).
- Parks, D. H. et al. CheckM: assessing the quality of microbial genomes recovered from isolates, single cells, and metagenomes. *Genome Research* 25, 1043–1055, doi:10.1101/gr.186072.114 (2015).
- Chklovski, A., Parks, D. H., Woodcroft, B. J. & Tyson, G. W. CheckM2: a rapid, scalable and accurate tool for assessing microbial genome quality using machine learning. *Nature Methods* 20, 1203–1212, doi:10.1038/s41592-023-01940-w (2023).

### Taxonomic Classification

- Chaumeil, P.-A., Mussig, A. J., Hugenholtz, P. & Parks, D. H. GTDB-Tk v2: memory friendly classification with the Genome Taxonomy Database. *Bioinformatics* 38, 5315–5316, doi:10.1093/bioinformatics/btac672 (2022).
- von Meijenfeldt, F. A. B. et al. Robust taxonomic classification of uncharted microbial sequences and bins with CAT and BAT. *Genome Biology* 20, 217, doi:10.1186/s13059-019-1817-x (2019).
- Wood, D. E., Lu, J. & Langmead, B. Improved metagenomic analysis with Kraken 2. *Genome Biology* 20, 257, doi:10.1186/s13059-019-1891-0 (2019).

### Functional Annotation

- Kanehisa, M., Furumichi, M., Tanabe, M., Sato, Y. & Morishima, K. KEGG: new perspectives on genomes, pathways, diseases and drugs. *Nucleic Acids Research* 45, D353–D361, doi:10.1093/nar/gkw1092 (2017).
- Huerta-Cepas, J. et al. eggNOG 5.0: a hierarchical, functionally and phylogenetically annotated orthology resource based on 5,090 organisms and 2,502 viruses. *Nucleic Acids Research* 47, D309–D314, doi:10.1093/nar/gky1085 (2019).
- Buchfink, B., Xie, C. & Huson, D. H. Fast and sensitive protein alignment using DIAMOND. *Nature Methods* 12, 59–60, doi:10.1038/nmeth.3176 (2015).
- Steinegger, M. & Söding, J. MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. *Nature Biotechnology* 35, 1026–1028, doi:10.1038/nbt.3988 (2017).

### Viral Analysis

- Camargo, A. P. et al. Identification of mobile genetic elements with geNomad. *Nature Biotechnology* 42, 1303–1312, doi:10.1038/s41587-023-01953-y (2023).
- Guo, J. et al. VirSorter2: a multi-classifier, expert-guided approach to detect diverse DNA and RNA viruses. *Microbiome* 9, 37, doi:10.1186/s40168-020-00990-y (2021).
- Nayfach, S. et al. CheckV assesses the quality and completeness of metagenome-assembled viral genomes. *Nature Biotechnology* 39, 578–585, doi:10.1038/s41587-020-00774-7 (2021).
- Kieft, K., Zhou, Z. & Anantharaman, K. VIBRANT: automated recovery, annotation and curation of microbial viruses, and evaluation of viral community function from genomic sequences. *Microbiome* 8, 90, doi:10.1186/s40168-020-00867-0 (2020).
- Kieft, K., Adams, A., Salamzade, R., Kalan, L. & Anantharaman, K. vRhyme enables binning of viral genomes from metagenomes. *Nucleic Acids Research* 50, e83, doi:10.1093/nar/gkac341 (2022).
- Roux, S. et al. iPHoP: An integrated machine learning framework to maximize host prediction for metagenome-derived viruses of archaea and bacteria. *PLOS Biology* 21, e3002083, doi:10.1371/journal.pbio.3002083 (2023).
- Zielezinski, A. et al. Ultrafast and accurate sequence alignment and clustering of viral genomes. *Nature Methods* 22, 1191–1194, doi:10.1038/s41592-025-02701-7 (2025).

### MAG Dereplication & Coverage

- Olm, M. R., Brown, C. T., Brooks, B. & Banfield, J. F. dRep: a tool for fast and accurate genomic comparisons that enables improved genome recovery from metagenomes through de-replication. *The ISME Journal* 11, 2864–2868, doi:10.1038/ismej.2017.126 (2017).
- Aroney, S. T. N. et al. CoverM: read alignment statistics for metagenomics. *Bioinformatics* 41, doi:10.1093/bioinformatics/btaf147 (2025).

### ORF Prediction

- Hyatt, D. et al. Prodigal: prokaryotic gene recognition and translation initiation site identification. *BMC Bioinformatics* 11, 119, doi:10.1186/1471-2105-11-119 (2010).
- Van der Jeugt, F., Dawyndt, P. & Mesuere, B. FragGeneScanRs: faster gene prediction for short reads. *BMC Bioinformatics* 23, doi:10.1186/s12859-022-04736-5 (2022).

### Metalink Custom Tools

- **mclink**: Li et al. (2025). mclink R package. https://github.com/LiuyangLee/mclink

---

## License

Metalink is released under the GNU General Public License v3.0 (GPL-3.0). See [LICENSE](LICENSE) for details.
