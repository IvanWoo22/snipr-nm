# SNIPR-Nm

SNIPR-Nm trims paired small-RNA reads with Cutadapt, merges pairs with PEAR, and reports separate **exact** and **polyA** counts for every supplied reference, including external controls and zero-count entries. The built-in offline QC report covers raw R1/R2, trimmed R1/R2 and merged reads.

[中文说明](README_zh.md) includes parameter definitions, scientific interpretation and validation of the 20260930 batch.

## Requirements and quick start

Validated locally with Nextflow 25.10.4, Java 21, Cutadapt 5.2, PEAR 0.9.11, Python 3.12, Perl and gzip. Python uses only its standard library; Perl needs `IO::Uncompress::Gunzip`. Biopython, NumPy and tsv-join are no longer required.

```bash
# Existing local installation, when available:
source analysis/20260930/env.sh

nextflow run SNIPR-Nm.nf \
  --samplesheet samplesheet.csv \
  --mirna_fa /path/to/hsa_miRNA_external.fa \
  --outdir results -work-dir work -with-trace trace.tsv
```

The prepared local environment uses `scripts/nextflow-local.sh` with the Java 21 module options required by the official launcher; a bare `java -jar` launcher failed to persist task caches in local validation. Use a different trace filename for each invocation.

Use `--help` for the CLI entry point. Resume with the same command and `-resume`. On other machines, install the dependencies and put executables on PATH, or pass `--cutadapt_bin`, `--pear_bin`, `--perl_bin`, `--py_bin` as executable paths. Changes to installed tool versions are not automatically detected by the Nextflow cache; validate upgrades in a fresh work directory.

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
| `local_cpus`, `local_memory`, `local_queue_size` | `8`, `12 GB`, `2` | Local scheduling budget |

`cutadapt_opts` is a trusted shell option string, not an untrusted data field. The resource budget applies to one Nextflow invocation and is a scheduling reservation, not an OS memory limit. A cluster requires a separate executor configuration; no LSF configuration is bundled.

## Count definitions and outputs

For each reference, **E** counts reads equal to the reference. **P** counts reads beginning with the complete reference followed by `AA`, `AAA`, or a tail longer than three bases containing at least 75% A. A single terminal A is excluded. These classes are disjoint for a given reference. Aliases and overlapping reference prefixes can share reads, preserving historical behavior.

`classified_total = E + P`. `no_tail_fraction = E / (E + P)`, or `NA` for a zero denominator. This is an untailed fraction among classified reads, not a calibrated absolute Nm measurement. The previous `(E-P)/E` interpretation was incorrect. For E=120 and P=3 the fraction is 120/123≈0.97561.

Published outputs:

- `mirna.merged.tsv`: `Name,S1,S1_polyA,S2,S2_polyA,...` (tab separated); an unsuffixed sample column contains **exact**, not E+P.
- `mirna.exact.tsv`, `mirna.polyA.tsv`, `mirna.classified_total.tsv`, `mirna.no_tail_fraction.tsv`: `Name,sequence` and one column per sample.
- `reference_annotations.tsv`: every reference name and normalized sequence, including controls.
- `<sample>/mirna.count.tsv`: `Name/exact/polyA`; `<sample>/mirna.<sample>.annotated.tsv`: sample-labelled columns.
- `<sample>/reads.count.tsv.gz`: sorted sequence/count rows without a header.
- `<sample>/input_manifest.json`, `cutadapt.log`, `cutadapt.json`, `pear.log`, `pear.json`: input identity, trimming statistics, merging log and validated intermediate compression records.
- `<sample>/<sample>_qc_metrics.json`, `_quality_by_pos.json`, `_length_dist.json`, `_qc_report.html`: stage-specific QC and an offline sample report.
- `snipr_qc_report.html`, `snipr_qc_data/`: batch report, sample summary TSV and copied QC files. These are custom SNIPR outputs, not MultiQC outputs.

Matrices follow samplesheet order and FASTA order. Aggregation consumes staged task outputs, rejects missing samples or mismatched references, and does not read asynchronously published directories. The standalone Perl counter retains headerless output by default: `perl bin/mirna_count.pl [--header] counts.tsv[.gz] reference.fa`. `fastq2count.pl` accepts stdin or a single FASTQ filename.

## QC, storage and limitations

Read counts, base counts, full length histograms, global mean quality, Q30 and GC use all reads in each of the five stages. Reads longer than 50 nt remain in the 18–25 nt percentage denominator. Position quality, sequence redundancy, dinucleotide entropy and exact adapter motifs use the first N reads, explicitly labelled as a nonrandom sample. Redundancy is not PCR duplication, and exact motif prevalence is not partial-adapter contamination. Empty denominators are null/NA. QC does not filter reads.

The workflow reconciles raw mate counts, Cutadapt input/output, PEAR assembled/unassembled/discarded counts and the sum of sequence counts. A library with zero reads after trimming skips PEAR, which otherwise aborts on empty input, and retains zero reference counts.

PEAR's four task-local FASTQs are compressed and decompression SHA256 is checked before removing each uncompressed copy. Compressed intermediates remain in work for resume. PEAR still needs space for one sample's uncompressed outputs; work accumulates across samples. There is no automatic work cleanup or free-space gate. Confirm deliverables before explicitly cleaning work; removed tasks cannot be reused from cache.

Deep caching hashes input content, and scripts are staged as task inputs so code changes invalidate dependent tasks. This adds hashing I/O. Failed runs may leave partial published results and are explicitly reported as failures. Prefer fresh output directories when changing parameters/references. Only a successful Nextflow exit after QC and aggregation establishes batch completion.

Unique-sequence counting still uses memory proportional to the number of distinct reads. Detailed QC memory scales with its sample limit and read length, rather than full-library depth. Five-stage full-file QC adds CPU/I/O time. macOS without containers does not provide reliable Nextflow trace CPU/RSS metrics.

## Tests and historical results

```bash
python3 -m unittest discover -s tests -v
python3 tests/run_nextflow_smoke.py --evidence /path/to/smoke-evidence
```

Tests cover exhaustive short-tail rules, aliases/controls, zero counts, malformed inputs, QC denominators, safe compression, strict aggregation, real Nextflow execution/resume, empty libraries and checksum failures.

The existing `analysis/20260930/hsa_primary/` contains the earlier independent executor's **21 full-library results**. This code update does not relabel those artifacts as regenerated by the revised Nextflow workflow. Original sources are archived against their original checksum manifest in `analysis/20260930/provenance/baseline_source/`. Validation of the revised counter against all 21 full-library count tables and end-to-end pilot runs is recorded under `analysis/20260930/validation/optimized_pipeline/`.
