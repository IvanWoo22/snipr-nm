# SNIPR-Nm — Small RNA Nm Ratio Analysis Pipeline

SNIPR-Nm is a high-throughput analysis workflow accompanying the experimental method for studying small-RNA abundance and 3′-terminal 2′-O-methylation (Nm). It connects the biological question—how much small RNA is present and what fraction carries terminal Nm—to reproducible sequence counts and quality-control evidence.

SNIPR-Nm trims paired small-RNA reads with Cutadapt, merges pairs with PEAR, and reports separate **exact** and **polyA** counts for every supplied reference, including external controls and zero-count entries. The built-in offline QC report covers raw R1/R2, trimmed R1/R2 and merged reads.

[中文说明](README_zh.md) includes setup, code-update instructions, parameter definitions and scientific interpretation.

## Contents

- [Background: What is SNIPR-Nm?](#background-what-is-snipr-nm)
- [Directory structure](#directory-structure)
- [Pipeline architecture](#pipeline-architecture)
- [Installation and setup](#installation-and-setup)
- [Updating an existing installation](#updating-an-existing-installation)
- [Inputs](#inputs)
- [Parameters](#parameters)
- [Count definitions and outputs](#count-definitions-and-outputs)
- [QC, storage and limitations](#qc-storage-and-limitations)
- [Recount from portable checkpoints](#recount-from-portable-checkpoints)
- [FAQ](#faq)
- [References](#references)

## Background: What is SNIPR-Nm?

3′-terminal 2′-O-methylation (Nm) is an important modification affecting small-RNA stability and function. It is widespread in plant small RNAs and is well established in animal piRNAs; Nm has also been reported in mammalian miRNAs. Measuring abundance together with terminal modification helps connect these molecules to their biological roles. This is the question that motivated SNIPR-Nm. See [Tang et al., 2025](https://doi.org/10.64898/2025.12.11.693694).

**Why combine abundance and modification measurements?** Conventional small RNA-seq can underestimate Nm-modified molecules because of adapter-ligation bias. Standard stem-loop RT-qPCR measures selected targets but does not, by itself, identify their terminal methylation status. Modification-sensitive assays based on β-elimination or poly(A) tailing provide complementary information, while combining abundance and modification measurements across many targets remains the motivation for a sequencing-based approach.

**The SNIPR-Nm experimental principle:** The method combines poly(A) polymerase with a thermostable ligase. Terminal Nm affects poly(A) addition, allowing tailing differences to carry information about modification into the sequencing library. After adapter removal and paired-read merging, the computational workflow distinguishes reads equal to a reference from reads containing that full reference followed by a qualifying poly(A) tail. The method paper reports an application to mouse testis, including piRNAs and miRNAs.

The workflow supports reference-based analysis of plant or animal small RNAs and external controls when the library design and adapters are compatible. For Nm interpretation, use libraries prepared with the corresponding experimental design and its controls. The output `no_tail_fraction` is the observed `E / (E + P)` among classified reads; relating it to a biological Nm fraction requires the experimental calibration. Detailed matching rules are given under [Count definitions and outputs](#count-definitions-and-outputs).

## Directory structure

```text
snipr-nm/
├── SNIPR-Nm.nf             # Nextflow workflow
├── nextflow.config         # Parameters and resource profiles
├── samplesheet.csv         # Sample-sheet template; replace placeholder paths
├── src/fastq_metrics.cpp   # Native full-file QC and sequence counting
├── bin/                    # Input checks, preprocessing, counting and reports
│   ├── stage_metrics.py    # Per-stage QC and merged-read counting
│   ├── mirna_count.pl      # Reference exact/polyA matching
│   ├── checkpoint_tools.py # Seal and verify portable counts checkpoints
│   ├── snipr_qc.py         # Standalone Python QC compatibility entry point
│   ├── fastq2count.pl      # Standalone Perl counting compatibility entry point
│   └── aggregate_qc.py     # Batch QC report
├── scripts/nextflow-local.sh
├── README.md
└── README_zh.md
```

The tree highlights the main entry points; `bin/` also contains helpers for trimming, PEAR, compact-mode cleanup, matrix export and HTML rendering. Supply your own reference FASTA with `--mirna_fa`; the repository does not bundle a species reference.

## Pipeline architecture

The original processing sequence—adapter removal, paired-read merging, sequence counting, reference matching and sample aggregation—remains the basis of the workflow. The current implementation collects QC at each stage and combines merged-read counting and QC in one scan:

```text
Paired FASTQ (R1/R2)
    ↓
Input validation and checksum audit → raw R1/R2 QC
    ↓
Cutadapt: adapter trimming → clean R1/R2 QC
    ↓
PEAR: paired-read merging
    ↓
Native scanner: sequence counts + merged QC
    ↓
Validated counts checkpoint
    ↓
Reference matching: exact / polyA
    ↓
Sample-labelled tables → multi-sample matrices

raw + clean + merged QC → sample reports → batch SNIPR QC report
```

This diagram shows the data-processing logic. `checkpoint` runs preprocessing as separately cached steps; `compact` groups preprocessing and clean/merged QC within a sample task to reclaim its temporary FASTQ files after validation; raw QC remains independently cached. Both use the same count definitions. `--counts_sheet` enters at a verified counts checkpoint for reference recounting. The Perl sequence counter remains available for standalone use; Nextflow uses the native scanner. Reports are generated by the built-in SNIPR reporting scripts.

## Installation and setup

### 1. Install the dependencies

On Linux/macOS, use Conda; if needed, install the matching architecture from the [official Miniforge instructions](https://github.com/conda-forge/miniforge#install). Skip environment creation if you already have the dependencies. The command below pins the main tool versions validated with this workflow:

```bash
conda create -n snipr --override-channels --strict-channel-priority \
  -c conda-forge -c bioconda \
  nextflow=25.10.4 openjdk=21 python=3.12 cutadapt=5.2 pear=0.9.11 \
  perl gzip git
conda activate snipr
```

The native QC scanner also requires a C++17 compiler and zlib development files. Nextflow builds and caches it automatically. Install these once for your operating system:

```bash
# Ubuntu / Debian
sudo apt-get update
sudo apt-get install build-essential zlib1g-dev
```

```bash
# macOS: install Xcode Command Line Tools and complete the installer dialog
xcode-select --install
```

On other Linux distributions, install the equivalent C++ compiler and zlib development packages. Verify the tools in the active environment:

```bash
command -v nextflow java cutadapt pear python3 perl gzip c++
nextflow -version
java -version
python3 --version
perl -MIO::Uncompress::Gunzip -e 'print "Perl gzip support OK\n"'
```

Python 3.11+ uses only its standard library; Perl needs `IO::Uncompress::Gunzip`. Biopython, NumPy and tsv-join are no longer required. See the [Bioconda](https://bioconda.github.io/) and [Nextflow installation documentation](https://docs.seqera.io/nextflow/install) for installation details.

### 2. Get the code and prepare inputs

```bash
git clone https://github.com/IvanWoo22/snipr-nm.git
cd snipr-nm
cp samplesheet.csv my_samples.csv
nextflow run SNIPR-Nm.nf --help
```

Edit `my_samples.csv` with real sample IDs and paired FASTQ paths in the `sample_id,R1,R2` columns. Prepare a reference FASTA for your species, target RNAs and external controls. See [Inputs](#inputs) for the format and constraints. Run `conda activate snipr` before using the workflow in each new terminal.

### 3. Save the configuration and run

Create your own `run.config` in the project directory, replace the reference path, and check that the adapters match your library:

```groovy
params {
    samplesheet = 'my_samples.csv'
    mirna_fa = '/path/to/reference.fa'
    outdir = 'results'
    adaptor_r1 = 'AGATCGGAAGAGCACA'
    adaptor_r2 = 'GATCGTCGGACTGTAGAA'
}
```

Tools available on PATH need no absolute paths. For a tool installed elsewhere, set, for example, `params.pear_bin = '/full/path/to/pear'` in `run.config`; other overrides are `cutadapt_bin`, `perl_bin`, `py_bin` and `cxx_bin`. These values accept executable paths only. Keeping personal sample sheets and configuration in separate files avoids editing the bundled files when updating the repository.

```bash
nextflow run SNIPR-Nm.nf -c run.config \
  -work-dir work -with-trace trace.tsv
```

After completion, inspect `results/mirna.merged.tsv` and `results/snipr_qc_report.html`. To resume, run the same command from the same launch directory with `-resume` and a new trace filename, keeping the original inputs, options, `work/` and `.nextflow/`. See [Parameters](#parameters) for resource settings and command-line overrides.

<details>
<summary>Existing environments or offline Nextflow</summary>

An existing environment only needs the dependencies above and executables on PATH. If the Nextflow 25.10.4 JAR is already under `.runtime/nextflow/` or `NXF_HOME`, set `JAVA_HOME` to Java 21 and use `scripts/nextflow-local.sh` in place of `nextflow`; the launcher does not install dependencies. Validate tool upgrades in a new work directory because version changes to tools on PATH do not automatically invalidate the Nextflow cache.

</details>

## Updating an existing installation

Wait for active runs to finish, save local source edits and activate the original runtime environment, then run:

```bash
cd /path/to/snipr-nm
git status --short
git switch main
git pull --ff-only origin main
nextflow run SNIPR-Nm.nf --help
```

Resolve unsaved source edits or diverged Git history before updating; personal sample sheets, configuration and analysis data do not need to be committed. This updates pipeline code only; see [Installation and setup](#installation-and-setup) for dependency changes. Resume using the original command plus `-resume`, retaining the inputs, configuration, launch directory and work directory. Tasks affected by code changes may rerun.

## Inputs

The repository's `samplesheet.csv` is a format template: replace its placeholder paths with your input files before running.

CSV columns: `sample_id,R1,R2`; optional `R1_md5,R2_md5` contain MD5 checksums of the input files as delivered, including gzip compression. Relative paths resolve against the CSV directory. Quote paths containing commas.

```csv
sample_id,R1,R2
S1,/data/S1_R1.fastq.gz,/data/S1_R2.fastq.gz
S2,/data/S2_R1.fastq.gz,/data/S2_R2.fastq.gz
```

Sample IDs must begin with a letter or digit and contain only letters, digits, `_`, `-`, `.`. `Name` and IDs ending in `_polyA` are reserved. Empty sheets, duplicate IDs and reused mate files are rejected. FASTQ files can be plain or gzip, use four-line records, ACGTN sequences and Phred+33 quality. Inputs must have finished transferring before the run. SHA256/MD5 identity is checked around trimming; supplied MD5 values are compared explicitly. Without a supplied expected checksum, hashing does not establish transfer completeness.

A reference FASTA is required: no species is assumed. Multiline sequences, lowercase and RNA U are normalized to uppercase DNA. References must contain at least four ACGTN bases. Full headers and reference order are preserved; duplicate names fail. Every reference is counted: supply a miRNA/control subset if the original FASTA also contains other RNA classes. Identical-sequence aliases remain separate rows, so summing reference rows is not a deduplicated read-union count.

## Parameters

| Parameter | Default | Meaning |
|---|---|---|
| `samplesheet`, `mirna_fa`, `outdir` | `samplesheet.csv`, required, `./results` | Inputs and destination |
| `adaptor_r1` | `AGATCGGAAGAGCACA` | R1 3′ adapter |
| `adaptor_r2` | `GATCGTCGGACTGTAGAA` | R2 3′ adapter |
| `cutadapt_opts` | `-O 5 -m 17 -e 0.18 --discard-untrimmed` | Historical trimming defaults; no added quality trimming |
| `pear_min_ovlp` | `17` | Minimum overlap, PEAR **`-v`** |
| `pear_min_length` | `17` | Minimum assembled length, PEAR **`-n`** |
| `pear_memory` | `1G` | PEAR `-y` argument |
| `cutadapt_threads`, `pear_threads` | `4`, `4` | Thread counts |
| `qc_sample_reads` | `100000` | First N reads per stage for detailed QC |
| `local_cpus`, `local_memory`, `local_queue_size` | `6`, `10 GB`, `6` | Desktop scheduling budget |
| `intermediates` | `checkpoint` | Step-level cache, or task-local FASTQ retirement with `compact` |
| `counts_sheet` | unset | Validated `sample_id,counts,manifest` CSV; skips FASTQ preprocessing/QC |
| `count_max_unique` | `1000000` | Spill sorted exact-count chunks at this unique-sequence threshold |
| `qc_max_forks`, `count_max_forks`, `preprocess_max_forks` | `4`, `2`, `2` | Per-process concurrency caps |
| `compact_threads` | `4` | Threads for sequential tools in a compact task |
| `execution_reports` | `false` | Enable trace, execution report and timeline under pipeline_info |

Use one resource profile per run. Command-line `--local_cpus`, `--local_memory`, `--local_queue_size` and tool-thread options override its defaults. When overriding resources in a separate `-c` file, set `executor.cpus`, `executor.memory`, `executor.queueSize` and the appropriate `process` selectors directly; changing only a late `params.*` assignment does not recompute previously derived config values. See the [Nextflow configuration rules](https://docs.seqera.io/nextflow/config).

`cutadapt_opts` is parsed into an argument list in both modes. Quoted values are supported; shell substitutions, variable expansion, pipes and redirects are not evaluated. The resource budget applies to one Nextflow invocation and is a scheduling reservation, not an OS memory limit. A cluster requires a separate executor configuration; no LSF configuration is bundled.

## Count definitions and outputs

For each reference, **E** counts reads equal to the reference. **P** counts reads beginning with the complete reference followed by `AA`, `AAA`, or a tail longer than three bases containing at least 75% A. A single terminal A is excluded. These classes are disjoint for a given reference. Aliases and overlapping reference prefixes can share reads, preserving historical behavior.

`classified_total = E + P`. `no_tail_fraction = E / (E + P)`, or `NA` for a zero denominator. This is an untailed fraction among classified reads, not a calibrated absolute Nm measurement. The previous `(E-P)/E` interpretation was incorrect. For E=120 and P=3 the fraction is 120/123≈0.97561.

Published outputs:

- `mirna.merged.tsv`: `Name,S1,S1_polyA,S2,S2_polyA,...` (tab separated); an unsuffixed sample column contains **exact**, not E+P.
- `mirna.exact.tsv`, `mirna.polyA.tsv`, `mirna.classified_total.tsv`, `mirna.no_tail_fraction.tsv`: `Name,sequence` and one column per sample.
- `reference_annotations.tsv`: every reference name and normalized sequence, including controls.
- `<sample>/mirna.count.tsv`: `Name/exact/polyA`; `<sample>/mirna.<sample>.annotated.tsv`: sample-labelled columns.
- `<sample>/reads.count.tsv.gz`: sorted sequence/count rows without a header; `counts.manifest.json` seals content hashes, inputs, parameters and QC for independent reuse.
- `<sample>/input_manifest.json`, `cutadapt.log`, `cutadapt.json`, `pear.log`, `pear.json`: input identity, trimming statistics, merging log and validated intermediate compression records (`pear.json` in checkpoint mode; `intermediates.json` in compact mode).
- `<sample>/<sample>_qc_metrics.json`, `_quality_by_pos.json`, `_length_dist.json`, `_qc_report.html`: stage-specific QC and an offline sample report.
- `snipr_qc_report.html`, `snipr_qc_data/`: batch report, sample summary TSV and copied QC files. These are custom SNIPR outputs, not MultiQC outputs.

Matrices follow samplesheet order and FASTA order. Aggregation consumes staged task outputs, rejects missing samples or mismatched references, and does not read asynchronously published directories. The standalone Perl counter retains headerless output by default: `perl bin/mirna_count.pl [--header] counts.tsv[.gz] reference.fa`. `fastq2count.pl` accepts stdin or a single FASTQ filename.

## QC, storage and limitations

Read counts, base counts, full length histograms, global mean quality, Q30 and GC use all reads in each of the five stages. Reads longer than 50 nt remain in the 18–25 nt percentage denominator. Position quality, sequence redundancy, dinucleotide entropy and exact adapter motifs use the first N reads, explicitly labelled as a nonrandom sample. Redundancy is not PCR duplication, and exact motif prevalence is not partial-adapter contamination. Empty denominators are null/NA. QC does not filter reads.

The workflow reconciles raw mate counts, Cutadapt input/output, PEAR assembled/unassembled/discarded counts and the sum of sequence counts. A library with zero reads after trimming skips PEAR, which otherwise aborts on empty input, and retains zero reference counts.

The C++17/zlib scanner is compiled and cached by BUILD_SCANNER. Metrics are generated per stage; merged counts and QC share one FASTQ parse. Numerical assembly and HTML rendering consume small files. Only actual script dependencies are staged, so rendering changes do not invalidate numerical tasks, reference changes reuse preprocessing, and adding/changing one sample preserves other sample caches.

`checkpoint` retains compressed PEAR outputs for step-level recovery. `compact` owns clean/merged temporary files inside one sample task, scans uncompressed merged reads once, verifies counts/QC/checkpoint integrity, and retires only six explicitly generated FASTQ files. Failures retain intermediate evidence. Raw QC is independently cached; changing compact-task QC parameters or recovering a failed compact task can rerun sample preprocessing. Temporary disk use scales with concurrent compact samples, but PEAR still needs its uncompressed outputs and there is no machine-wide disk quota guarantee.

To start a new FASTQ run with compact intermediates and execution reports:

```bash
nextflow run SNIPR-Nm.nf --samplesheet samplesheet.csv \
  --mirna_fa reference.fa --intermediates compact --outdir compact-results \
  -work-dir compact-work --execution_reports true
```

Exact counts spill sorted chunks when `count_max_unique` is reached and use at most 32-way merges. This is a unique-entry threshold, not an RSS byte limit: read length, parsing buffers and bounded detailed-QC tables also consume memory. Spilling adds temporary disk I/O. Count compression and SHA256 verification also stream data rather than loading the entire count file into memory. Counts are never approximated or truncated. The legacy Perl/Python commands remain available as compatibility entry points and semantic oracles.

Profiles: `desktop` (default 6 CPU/10 GB/6 task slots), `workstation` (8 CPU/12 GB/8 slots), and `test` (2 CPU/6 GB, single-thread tools). Each process cap is separate; total resource budgets prevent oversubscription. Local reservations are not enforced memory limits, and independent workflow runs do not share budgets.

Deep caching preserves input identity checks. Published task outputs explicitly overwrite matching destinations on resume, but obsolete sample directories are not automatically removed: prefer a fresh outdir for a different sample set. Never clean active work outputs. Tool upgrades require revalidation because PATH tool versions are not automatically included in cache keys. macOS without containers may lack CPU/RSS trace metrics; missing measurements are not zero usage.

`--execution_reports true` writes `pipeline_info/trace.tsv`, `report.html`, and `timeline.html`, overwriting per run. Use new outdirs or explicit `-with-trace/-with-report/-with-timeline` paths for a history.

## Recount from portable checkpoints

```csv
sample_id,counts,manifest
S1,/results/S1/reads.count.tsv.gz,/results/S1/counts.manifest.json
```

```bash
nextflow run SNIPR-Nm.nf --counts_sheet checkpoints.csv \
  --mirna_fa updated_reference.fa --outdir recount -work-dir recount-work
```

Import validates the `snipr-counts-v1` schema, sample identity, compressed/content SHA256, unique sorted rows, read totals and QC evidence. It runs reference counting/matrix generation without raw FASTQ, scanner compilation or a new QC report, and needs no old work cache. Legacy counts without a valid manifest are not automatically trusted. Counts do not retain read IDs, qualities or pairing; preserve raw inputs.

The native scanner combines full-file validation with GC, quality-sum and Q30 accumulation in one pass and uses buffered line scanning. Full-file and first-N definitions remain unchanged. The standalone Python QC and Perl counting entry points remain available for compatibility.

## FAQ

### An input file cannot be found

Check the sample sheet and reference path, and confirm that the files exist and have finished transferring. Relative FASTQ paths are resolved against the sample-sheet directory. Software paths must identify an executable; verify its availability in the environment used to launch Nextflow.

### A sample or reference has zero counts

Check `cutadapt.json`, `pear.log` and the sample QC report to locate where reads were retained or lost. A sample with no reads after trimming is supported and keeps zero-count reference rows. For a reference with zero matches, check its sequence and the exact/polyA matching rules. Both relative and absolute `--outdir` paths are supported; changing the path alone does not resolve a failed task or a missing report. A failed workflow can leave partial published outputs.

### Adding more samples

Append unique sample IDs and mate paths to the sample sheet, then rerun the original command with `-resume` from the same launch directory and work directory. Unchanged samples can reuse their cached tasks. If removing samples or changing the analysis set, prefer a new output directory because obsolete published sample directories are not removed automatically.

### Investigating a failed step

Start with `.nextflow.log`, then inspect `.command.err`, `.command.out` and `.command.sh` in the failed task's work directory. Correct the reported input, dependency or resource problem before resuming with the original command plus `-resume`. Keep the work directory and `.nextflow/` cache while investigating; they contain the evidence and cached outputs needed for recovery.

## References

Tang Y, et al. **A high-throughput methods SNIPR-Nm reveals small RNA methylation ratio.** bioRxiv, 2025 (preprint). [doi:10.64898/2025.12.11.693694](https://doi.org/10.64898/2025.12.11.693694). Consult the method paper for experimental design, controls and interpretation of Nm ratios.
