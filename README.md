# SNIPR-Nm — Small RNA Nm Ratio Analysis Pipeline

> **中文说明**: 本文档为英文版。如需查看中文版说明，请访问 [README_zh.md](README_zh.md)。

SNIPR-Nm (Small RNA Integration-based Polyadenylation and Ligation for Nm detection) is a high-throughput analysis pipeline designed for detecting 3' end 2'-O-methylation (Nm) in small RNAs. It simultaneously outputs **expression levels** and **Nm modification ratios** for each small RNA, suitable for plant, animal tissues, or any small RNA libraries.

> Related Literature: Tang et al. "SNIPR-Nm reveals small RNA methylation ratio" (bioRxiv, 2025)

---

## Table of Contents

- [Directory Structure](#directory-structure)
- [Background: What is SNIPR-Nm?](#background-what-is-snipr-nm)
- [Step 1: Install Dependencies](#step-1-install-dependencies)
- [Step 2: Prepare Sample Sheet](#step-2-prepare-sample-sheet)
- [Step 3: Verify Software Paths](#step-3-verify-software-paths)
- [Pipeline Architecture](#pipeline-architecture)
- [Step 4: Run the Pipeline](#step-4-run-the-pipeline)
- [Parameter Reference](#parameter-reference)
- [Output Results](#output-results)
- [QC Metrics](#qc-metrics)
- [FAQ](#faq)
- [References](#references)

---

## Directory Structure

```
snipr-nm/
├── SNIPR-Nm.nf          # Main pipeline file
├── nextflow.config      # Parameter configuration file
├── samplesheet.csv      # Sample list (needs to be filled)
├── ath_mirna_neo.fa     # miRNA reference sequences (replace with your own)
├── bin/
│   ├── fastq2count.pl   # FASTQ → sequence count table
│   ├── mirna_count.pl    # miRNA alignment counting + poly-A detection
│   ├── snipr_qc.py      # Small RNA-seq QC analysis
│   └── aggregate_qc.py   # Aggregate QC reports for all samples
├── README.md            # English version
└── README_zh.md         # Chinese version
```

---

## Background: What is SNIPR-Nm?

3' end 2'-O-methylation (Nm) is an important modification that regulates small RNA stability and function. In plants, HEN1-mediated methylation covers almost all small RNAs; in animals, Nm primarily appears in piRNAs, and recent studies have also found partial methylation modifications in mammalian miRNAs.

**Limitations of Existing Technologies:**
- Traditional small RNA-seq: Relies on adapter ligation, Nm-modified RNAs are **underestimated** due to low ligation efficiency
- Stem-loop RT-qPCR: Can accurately quantify specific small RNAs, but **cannot report methylation status**
- β-elimination or poly(A) tailing qPCR: Can detect 3' Nm, but only for specific targets, requires large sample input, and cannot provide absolute quantification

**SNIPR-Nm's Solution:**
Utilizing a combination of poly(A) polymerase and a special ligase to simultaneously detect small RNA **expression levels** and **Nm modification ratios** in the same reaction. The specific principle:
1. Poly(A) polymerase performs tailing on RNA
2. Nm-modified RNAs are **not tailed** because their 3' end is protected; unmodified RNAs are poly-A tailed
3. Differential amplification through ligase distinguishes modified from unmodified molecules
4. Combined with high-throughput sequencing, enabling transcriptome-wide Nm ratio quantification

This method has been successfully applied to mouse testicular tissue, yielding reliable piRNA and miRNA Nm modification ratio data.

---

## Step 1: Install Dependencies

### Using conda (Recommended)

```bash
# Create environment
conda create -n snipr -c bioconda -c conda-forge \
    nextflow cutadapt pear tsv-utils "python>=3.9"

# Activate and install Python dependencies
conda activate snipr
pip install biopython

# pear is recommended to be installed from KantorLab channel (newer version):
# conda install -n snipr -c kantorlab pear
# If the channel is not available, use bioconda version:
# conda install -n snipr -c bioconda pear

# Verify all tools are installed
which nextflow cutadapt pear python3 perl tsv-join
```

**Software Functions:**

| Software | Purpose |
|----------|---------|
| `nextflow` | Pipeline execution engine |
| `cutadapt` | Remove adapter sequences |
| `pear` | Merge PE reads |
| `perl` | Run miRNA counting scripts |
| `python3` + biopython | QC scripts |
| `tsv-utils` | Merge count tables from multiple samples |

### Using Existing System Software

If your system already has these tools and biopython installed, you can skip conda. Use `which <command>` to verify they are accessible.

---

## Step 2: Prepare Sample Sheet

Edit `samplesheet.csv` with your own data paths:

```csv
sample_id,R1,R2
U1,/home/username/data/U1_1.fq.gz,/home/username/data/U1_2.fq.gz
U2,/home/username/data/U2_1.fq.gz,/home/username/data/U2_2.fq.gz
```

- `sample_id`: Sample name (avoid Chinese characters and special symbols)
- `R1` / `R2`: **Complete absolute paths** to FASTQ files

**How to Get Absolute Path?**

```bash
cd /your/data/directory
pwd    # Print absolute path
ls *.fq.gz  # Confirm file names
```

---

## Step 3: Verify Software Paths

```bash
conda activate snipr
which cutadapt   # Output example: /home/xxx/miniconda3/envs/snipr/bin/cutadapt
which pear
which python3
which perl
```

Fill these paths into the corresponding lines in `nextflow.config` (e.g., `params.cutadapt_bin = "/full/path/to/cutadapt"`). If commands can be found directly in PATH, you can just use the command name (e.g., `"cutadapt"`), no need for full path.

---

## Pipeline Architecture

Before running the pipeline, it is recommended to understand the overall processing steps:

```
Input FASTQ (PE small RNA library)
    ↓
[Cutadapt]  Remove adapter sequences
    ↓
[PEAR]       Merge PE reads
    ↓
[fastq2count] Convert to sequence count table
    ↓
[mirna_count] Align with miRNA reference
    ↓           Detect poly-A tails simultaneously (= unmodified molecules)
    ↓
[annotate]    Generate unique filenames for each sample (avoid name conflicts)
    ↓
[merge]       Merge all samples horizontally into one table
    ↓
[QC]          Quality control (length distribution, quality, complexity)
    ↓
[multiqc]     Generate aggregated QC report
```

**Step Descriptions:**

| Step | Function | Output Files |
|------|----------|--------------|
| Cutadapt | Remove adapters, filter low-quality reads | clean_R1/R2.fq.gz |
| PEAR | Merge paired-end reads | assembled.assembled.fastq |
| fastq2count | Deduplicate and count sequences | reads.count.tsv |
| mirna_count | Align with miRNA reference, detect poly-A tails | mirna.count.tsv |
| annotate | Add sample name prefix | mirna.{sample}.annotated.tsv |
| merge | Merge all sample results | mirna.merged.tsv |
| QC | Quality control analysis | *_qc_metrics.json, *_qc_report.html |
| multiqc | Generate summary report | multiqc_report.html |

---

## Step 4: Run the Pipeline

```bash
cd /path/to/snipr-nm

nextflow run SNIPR-Nm.nf \
    --samplesheet samplesheet.csv \
    --mirna_fa ath_mirna_neo.fa \
    --outdir ./results
```

**Resume after interruption** (skip completed steps):

```bash
nextflow run SNIPR-Nm.nf \
    --samplesheet samplesheet.csv \
    --mirna_fa ath_mirna_neo.fa \
    --outdir ./results \
    -resume
```

---

## Parameter Reference

### Required Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--samplesheet` | `samplesheet.csv` | Sample list file path (CSV format with sample_id, R1, R2 columns) |
| `--mirna_fa` | `ath_mirna_neo.fa` | miRNA reference sequence FASTA file path |
| `--outdir` | `./results` | Output directory path (**using absolute path is recommended** to avoid file not found issues) |

### Optional Parameters

#### Adapter Sequence Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--adaptor_r1` | `AGATCGGAAGAGCACA` | R1 adapter sequence (common Illumina Small RNA TruSeq adapter) |
| `--adaptor_r2` | `GATCGTCGGACTGTAGAA` | R2 adapter sequence |

#### Thread and Performance Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--cutadapt_threads` | `8` | Number of CPU threads for Cutadapt step |
| `--pear_threads` | `4` | Number of CPU threads for PEAR step |

#### PEAR Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--pear_min_ovlp` | `17` | PEAR minimum overlap length (bp). Reads are merged only when R1 and R2 overlap ≥ this length. Recommended range: 15-25, typically 17-20 for small RNA-seq |

### Advanced Parameters

#### Cutadapt Advanced Options

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--cutadapt_opts` | `-O 5 -m 17 -e 0.18 --discard-untrimmed` | Cutadapt advanced options string. See details below: |

**`--cutadapt_opts` Detailed Explanation:**

| Option | Value | Description |
|--------|-------|-------------|
| `-O 5` | 5 | Minimum adapter overlap length. Only considered valid adapter when overlap ≥ 5bp |
| `-m 17` | 17 | Minimum read length (bp). Reads shorter than this are discarded. miRNAs are typically 18-25nt, recommend setting 15-17 |
| `-e 0.18` | 0.18 | Maximum error rate (allowed mismatches/overlap length). 0.18 means ~18% error tolerance, balancing sensitivity and specificity |
| `--discard-untrimmed` | - | Discard reads without detected adapters. These may be low-quality or off-target fragments |

**Other Available Cutadapt Options (require manual modification of `nextflow.config`):**

| Option | Example Value | Description |
|--------|---------------|-------------|
| `-q` | `20,20` | Quality trim threshold, format "front,back". Remove bases with quality below this value |
| `--max-n` | `0.5` | Maximum N proportion (0-1). Allowed proportion of N bases |
| `-M` | `30` | Maximum read length. Reads longer than this are discarded |

### Common Usage Examples

```bash
# Basic run (using absolute path for output)
nextflow run SNIPR-Nm.nf \
    --samplesheet samplesheet.csv \
    --mirna_fa ath_mirna_neo.fa \
    --outdir /home/username/srna_results

# Change adapter sequences (if your library uses different adapters)
nextflow run SNIPR-Nm.nf \
    --samplesheet samplesheet.csv \
    --mirna_fa ath_mirna_neo.fa \
    --outdir ./results \
    --adaptor_r1 your_R1_adapter_sequence \
    --adaptor_r2 your_R2_adapter_sequence

# Adjust thread count (based on server configuration)
nextflow run SNIPR-Nm.nf \
    --samplesheet samplesheet.csv \
    --mirna_fa ath_mirna_neo.fa \
    --outdir ./results \
    --cutadapt_threads 16 \
    --pear_threads 8

# Adjust PEAR overlap parameter (for special libraries)
nextflow run SNIPR-Nm.nf \
    --samplesheet samplesheet.csv \
    --mirna_fa ath_mirna_neo.fa \
    --outdir ./results \
    --pear_min_ovlp 20

# Complete parameter example
nextflow run SNIPR-Nm.nf \
    --samplesheet samplesheet.csv \
    --mirna_fa ath_mirna_neo.fa \
    --outdir /home/username/srna_results \
    --adaptor_r1 AGATCGGAAGAGCACA \
    --adaptor_r2 GATCGTCGGACTGTAGAA \
    --cutadapt_threads 16 \
    --pear_threads 8 \
    --pear_min_ovlp 17

# Resume after interruption
nextflow run SNIPR-Nm.nf \
    --samplesheet samplesheet.csv \
    --mirna_fa ath_mirna_neo.fa \
    --outdir ./results \
    -resume
```

### Configuration File

For frequently used parameter combinations, it is recommended to modify the `nextflow.config` file:

```groovy
// nextflow.config example configuration
params {
    samplesheet    = "samplesheet.csv"
    mirna_fa       = "ath_mirna_neo.fa"
    outdir         = "/home/username/srna_results"  // Use absolute path
    adaptor_r1     = "AGATCGGAAGAGCACA"
    adaptor_r2     = "GATCGTCGGACTGTAGAA"
    cutadapt_threads = 16
    pear_threads     = 8
    pear_min_ovlp    = 17
}
```

After modifying the configuration file, you only need to specify required parameters when running:

```bash
nextflow run SNIPR-Nm.nf --samplesheet samplesheet.csv --mirna_fa ath_mirna_neo.fa
```

---

## Output Results

### Directory Structure

```
results/
├── U1/
│   ├── cutadapt.log              # Adapter removal log
│   ├── reads.count.tsv           # Sequence count table
│   ├── mirna.count.tsv          # miRNA alignment results
│   ├── mirna.U1.annotated.tsv   # Annotated file with sample name
│   ├── U1_qc_metrics.json        # QC metrics (JSON)
│   ├── U1_qc_report.html         # Detailed QC HTML report
│   ├── U1_quality_by_pos.json    # Per-position quality data
│   └── U1_length_dist.json        # Length distribution data
├── U2/ ... (one subdirectory per sample)
├── mirna.merged.tsv              # Horizontally merged table for all samples
├── multiqc_report.html           # Aggregated QC report
└── multiqc_data/                 # QC raw data
```

### Detailed Output File Descriptions

#### 1. cutadapt.log (Adapter Removal Log)

Records Cutadapt processing statistics:

```log
This is cutadapt 4.9 with Python 3.11.0
Command line parameters: -a AGATCGGAAGAGCACA -A GATCGTCGGACTGTAGAA ...
Command line parameters: -a AGATCGGAAGAGCACA -A GATCGTCGGACTGTAGAA -O 5 -m 17 -e 0.18 ...
Finished 4 processes (0.4s wall time, 14.6 MB peak memory)

=== Adapter 1 (first read) ===
Sequence: AGATCGGAAGAGCACA; Type: regular 3'; Trimmed: 1250 reads

=== Adapter 2 (second read) ===
Sequence: GATCGTCGGACTGTAGAA; Type: regular 3'; Trimmed: 980 reads

=== Summary ===
Reads written (passing filters): 85000
```

**Key Metric Interpretation:**

| Metric | Description | Reference Value |
|--------|-------------|-----------------|
| Adapter 1/2 Trimmed | Number of reads with successful adapter removal | - |
| Reads written | Final reads passing quality filters | - |
| Untrimmed reads ratio | Reads discarded by `--discard-untrimmed` | < 20% is normal |

#### 2. reads.count.tsv (Sequence Count Table)

Converts FASTQ to sequence count table, one unique sequence per line with its count:

```tsv
sequencount
UUUAAAGCUAUAGUAGUAGUACCCGUACAUCUACUAUAUACUAUAUACUAUAUACU  125
AGCUAGCUAGCUAGCUAGCUA  98
UCGUACGUCAGCUAGCUAGCU  67
...
```

#### 3. mirna.count.tsv (miRNA Alignment Results)

Exact match and poly-A tail detection counts for each reference miRNA:

```tsv
Name    exact  polyA
ath-miR156a-5p  120  3
ath-miR166b     45   1
ath-miR167a     89   2
ath-miR172a     67   0
...
```

**Column Descriptions:**

| Column | Description |
|--------|-------------|
| Name | miRNA name (from reference FASTA) |
| exact | Exactly matched sequence count |
| polyA | Detected poly-A tail sequence count (unmodified) |

#### 4. mirna.{sample}.annotated.tsv (Annotated File with Sample Name)

Adds sample name prefix on top of mirna.count.tsv for merging:

```tsv
Name    U1  U1_polyA
ath-miR156a-5p  120  3
ath-miR166b     45   1
ath-miR167a     89   2
```

#### 5. mirna.merged.tsv (Merged Summary Table)

Horizontally merged results from all samples:

```tsv
Name            U1    U1_polyA  U2    U2_polyA  U3    U3_polyA
ath-miR156a-5p  120   3         98    5         145   12
ath-miR166b     45    1         67    2         78    3
ath-miR167a     89    2         102   1         156   8
ath-miR172a     67    0         89    0         123   2
ath-miR393a     23    5         34    8         45    15
...
```

#### 6. QC Output Files

- `*_qc_metrics.json`: QC metrics summary (JSON format)
- `*_qc_report.html`: Detailed QC HTML report (viewable in browser)
- `*_quality_by_pos.json`: Per-position quality distribution data
- `*_length_dist.json`: Read length distribution data
- `multiqc_report.html`: Aggregated QC report for all samples
- `multiqc_data/`: QC raw data directory

### Nm Modification Ratio Calculation Example

Based on data in `mirna.merged.tsv`, you can calculate Nm modification ratios:

**Example Data:**

```tsv
Name            U1    U1_polyA  U2    U2_polyA
ath-miR156a-5p  120   3         98    5
ath-miR166b     45    1         67    2
```

**Calculation Method:**

For ath-miR156a-5p in sample U1:

- Total reads = 120
- poly-A tail reads = 3 (unmodified)
- Suspected Nm modified = 120 - 3 = 117
- **Nm modification ratio = (120 - 3) / 120 = 97.5%**

**Python Calculation Script Example:**

```python
import pandas as pd

# Read merged results
df = pd.read_csv('results/mirna.merged.tsv', sep='\t')

# Calculate Nm modification ratio for each sample
for sample in ['U1', 'U2']:
    total_col = sample
    polya_col = f'{sample}_polyA'
    nm_ratio_col = f'{sample}_Nm_ratio'
    
    df[nm_ratio_col] = (df[total_col] - df[polya_col]) / df[total_col]
    df[nm_ratio_col] = df[nm_ratio_col].fillna(0)

# Save results
df.to_csv('results/mirna.merged.with_nm_ratio.tsv', sep='\t', index=False)

print(df[['Name', 'U1', 'U1_polyA', 'U1_Nm_ratio']].head())
```

**Output Example:**

```tsv
Name            U1    U1_polyA  U1_Nm_ratio
ath-miR156a-5p  120   3         0.975
ath-miR166b     45    1         0.978
ath-miR167a     89    2         0.978
```

### Complete Run Example

#### 1. Prepare Sample Sheet

```csv
sample_id,R1,R2
WT_rep1,/data/WT_rep1_1.fq.gz,/data/WT_rep1_2.fq.gz
WT_rep2,/data/WT_rep2_1.fq.gz,/data/WT_rep2_2.fq.gz
KO_rep1,/data/KO_rep1_1.fq.gz,/data/KO_rep1_2.fq.gz
KO_rep2,/data/KO_rep2_1.fq.gz,/data/KO_rep2_2.fq.gz
```

#### 2. Run Pipeline

```bash
cd /path/to/snipr-nm

nextflow run SNIPR-Nm.nf \
    --samplesheet samplesheet.csv \
    --mirna_fa ath_mirna_neo.fa \
    --outdir /home/username/srna_results \
    -resume
```

#### 3. Check Run Log

```log
N E X T F L O W  ~  version 23.04.0
Launching `SNIPR-Nm.nf` [sad_ardinghelli] DSL2 - revision: 8e3a7c3a8
------------------------------------------------------
executor >  local (6)
[aa/bc1234] process > CUTADAPT (WT_rep1)       [100%] 1 of 1, CPU: 8
[ab/cd5678] process > PEAR (WT_rep1)           [100%] 1 of 1, CPU: 4
[ac/ef9012] process > FASTQ2COUNT (WT_rep1)   [100%] 1 of 1
[ad/fg3456] process > MIRNA_COUNT (WT_rep1)   [100%] 1 of 1
[ae/gh7890] process > ANNOTATE_SAMPLE (WT_rep1) [100%] 1 of 1
[ef/ij1234] process > MERGE_SAMPLES           [100%] 1 of 1
[fg/kl5678] process > QC_SAMPLE (WT_rep1)     [100%] 1 of 1
[gh/mn9012] process > MULTIQC                 [100%] 1 of 1
==================================================
SNIPR-Nm run complete
Results : /home/username/srna_results
Merged  : /home/username/srna_results/mirna.merged.tsv
QC html : /home/username/srna_results/multiqc_report.html
==================================================
```

#### 4. Check Output Files

```bash
# Check output directory structure
ls -la /home/username/srna_results/

# Check merged results
head -n 10 /home/username/srna_results/mirna.merged.tsv

# Open QC report
open /home/username/srna_results/multiqc_report.html
```

---

## QC Metrics

`snipr_qc.py` calculates the following metrics on clean R1:

| Metric | Description | Reference Threshold |
|--------|-------------|---------------------|
| **miRNA 18–25nt proportion** | Proportion of miRNA-typical length reads | > 30% is good |
| **Adapter contamination rate** | Reads still containing 3' adapters | < 5% is good |
| **Low complexity reads** | High proportion of repetitive sequences (possibly rRNA/tRNA contamination) | < 20% |
| **Duplication rate** | Duplicate sequence proportion | Varies with library complexity |
| **Average quality** | Phred score | Q30+ is good |

> ⚠️ Thresholds are for reference only; actual judgment should be based on library construction method and experimental background.

---

## FAQ

### Q: Error "No such file or directory"
Check if paths in `samplesheet.csv` are **complete absolute paths** and files actually exist:
```bash
ls -la /complete/path/to/U1_1.fq.gz
```

### Q: QC data is all 0 or reports no data
Confirm `--outdir` uses **absolute path**:
```bash
--outdir /home/username/srna_results
```
Relative paths (`./results`) may cause Python scripts to fail to find JSON files in some cases.

### Q: Want to add more samples
Simply append rows in `samplesheet.csv`; no pipeline modification needed.

### Q: How to investigate when a step reports an error
```bash
nextflow run SNIPR-Nm.nf ... -resume
# Check .nextflow.log or .command.log in corresponding work directory
```

---

## References

Tang et al. "A high-throughput method SNIPR-Nm reveals small RNA methylation ratio." bioRxiv (2025). doi: 10.64898/2025.12.11.693694
