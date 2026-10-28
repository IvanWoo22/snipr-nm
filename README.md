# SNIPR-Nm

SNIPR-Nm trims paired small-RNA reads with Cutadapt, merges pairs with PEAR, and reports separate **exact** and **polyA** counts for every supplied reference, including external controls and zero-count entries. The built-in offline QC report covers raw R1/R2, trimmed R1/R2 and merged reads.

[中文说明](README_zh.md) includes setup, code-update instructions, parameter definitions and scientific interpretation.

## Requirements and quick start

Validated locally with Nextflow 25.10.4, Java 21, Cutadapt 5.2, PEAR 0.9.11, Python 3.12, Perl, gzip, a C++17 compiler and zlib development headers. Python 3.11+ uses only its standard library; Perl needs `IO::Uncompress::Gunzip`. Biopython, NumPy and tsv-join are no longer required.

```bash
# Activate an environment with the dependencies listed above on PATH.
nextflow run SNIPR-Nm.nf \
  --samplesheet samplesheet.csv \
  --mirna_fa /path/to/hsa_miRNA_external.fa \
  --outdir results -work-dir work -with-trace trace.tsv
```

If you already have the Nextflow 25.10.4 JAR under `.runtime/nextflow/` or `NXF_HOME`, the optional `scripts/nextflow-local.sh` launcher includes the required Java 21 module options. Set `JAVA_HOME` for that installation. The launcher does not install dependencies. Use a different trace filename for each invocation.

Use `--help` for the CLI entry point. Resume with the same command and `-resume`. On other machines, install the dependencies and put executables on PATH, or pass `--cutadapt_bin`, `--pear_bin`, `--perl_bin`, `--py_bin`, `--cxx_bin` as executable paths. Changes to installed tool versions are not automatically detected by the Nextflow cache; validate upgrades in a fresh work directory.

## Updating an existing installation

Wait for active pipeline runs to finish, then update the existing checkout. `origin/main` must contain the version you intend to install. Check local changes first; commit needed source edits to a separate branch or save them elsewhere before continuing. Local data and run directories do not need to be added to Git.

```bash
cd /path/to/snipr-nm
# Review tracked source edits before proceeding.
git status --short
git diff --stat
git switch main

# Keep a reference to the currently installed commit, then fetch and review updates.
git branch "codex/before-update-$(date +%Y%m%d-%H%M%S)"
git fetch origin
git log --oneline HEAD..origin/main

# Only advance main when it can be fast-forwarded.
git merge --ff-only origin/main
git log -1 --oneline
nextflow run SNIPR-Nm.nf --help
```

If Git reports conflicting local edits or diverged history, resolve those changes before updating. The commands above do not overwrite local history or delete analysis files. Updating the repository changes pipeline code; it does not install or upgrade Nextflow, Python, Cutadapt or PEAR. This version requires Python 3.11+, a C++17 compiler and zlib headers. Nextflow automatically builds and caches the scanner.

To continue with an existing work directory after updating:

```bash
nextflow run SNIPR-Nm.nf \
  --samplesheet samplesheet.csv \
  --mirna_fa /path/to/reference.fa \
  --outdir results -work-dir work \
  -resume -with-trace trace.updated.tsv
```

Use the original inputs, output directory and work directory when resuming. Task code and dependency changes can cause affected tasks to rerun; an update does not guarantee a fully cached run. Keep `work/` and `.nextflow/` for recovery. For an independent new run, choose new output/work directories and omit `-resume`.

## Inputs

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

Exact counts spill sorted chunks when `count_max_unique` is reached and use at most 32-way merges. This is a unique-entry threshold, not an RSS byte limit: read length, parsing buffers and bounded detailed-QC tables also consume memory. Spilling adds temporary disk I/O. Counts are never approximated or truncated. The legacy Perl/Python commands remain available as compatibility entry points and semantic oracles.

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

nextflow run SNIPR-Nm.nf --samplesheet samplesheet.csv \
  --mirna_fa reference.fa --intermediates compact --outdir compact-results \
  -work-dir compact-work --execution_reports true
```

Import validates the `snipr-counts-v1` schema, sample identity, compressed/content SHA256, unique sorted rows, read totals and QC evidence. It runs reference counting/matrix generation without raw FASTQ, scanner compilation or a new QC report, and needs no old work cache. Legacy counts without a valid manifest are not automatically trusted. Counts do not retain read IDs, qualities or pairing; preserve raw inputs.

The native scanner combines full-file validation with GC, quality-sum and Q30 accumulation in one pass and uses buffered line scanning. Full-file and first-N definitions remain unchanged. The standalone Python QC and Perl counting entry points remain available for compatibility.
