#!/usr/bin/env nextflow
/*
================================================================================
SNIPR-Nm  —  Nextflow Pipeline
Small RNA-seq analysis: Cutadapt → PEAR merge → miRNA counting + poly-A detection

Usage:
  # Local
  nextflow run SNIPR-Nm.nf \
      --samplesheet samplesheet.csv \
      --mirna_fa ath_mirna_neo.fa \
      --outdir ./results

  # LSF cluster
  nextflow run SNIPR-Nm.nf -c conf/lsf.config \
      --samplesheet samplesheet.csv \
      --mirna_fa ath_mirna_neo.fa \
      --outdir /path/to/results

  # Resume
  nextflow run SNIPR-Nm.nf -resume ...

Input:
  samplesheet.csv  — CSV with columns: sample_id, R1, R2
  <mirna_fa>      — FASTA of known + novel miRNA references
================================================================================
*/

nextflow.enable.dsl = 2

/*
────────────────────────────────────────────────────────────────────────────────
Parameters
────────────────────────────────────────────────────────────────────────────────
*/
params.samplesheet       = "samplesheet.csv"
params.mirna_fa          = "ath_mirna_neo.fa"
params.outdir            = "./results"
params.adaptor_r1        = "AGATCGGAAGAGCACA"
params.adaptor_r2        = "GATCGTCGGACTGTAGAA"
params.cutadapt_opts     = "-O 5 -m 17 -e 0.18 --discard-untrimmed"
params.pear_min_ovlp    = 17
params.cutadapt_threads = 8
params.pear_threads    = 4
params.bin_dir          = "${workflow.projectDir}/bin"

// 工具路径（推荐在 nextflow.config 或命令行中覆盖）
params.cutadapt_bin = "cutadapt"
params.pear_bin     = "pear"
params.tsvjoin_bin  = "tsv-join"
params.perl_bin     = "perl"
params.py_bin       = "python3"

/*
────────────────────────────────────────────────────────────────────────────────
Samplesheet channel
────────────────────────────────────────────────────────────────────────────────
*/
Channel
    .fromPath(params.samplesheet)
    .splitCsv(header: true)
    .map { row -> [
        row.sample_id,
        file(row.R1,  checkExists: true),
        file(row.R2,  checkExists: true)
    ] }
    .set { input_samples }

/*
────────────────────────────────────────────────────────────────────────────────
Process: CUTADAPT
────────────────────────────────────────────────────────────────────────────────
*/
process CUTADAPT {
    tag { sample_id }
    publishDir "${params.outdir}/${sample_id}", mode: 'copy',
        pattern: 'cutadapt.log'

    input:
        tuple val(sample_id), path(r1), path(r2)

    output:
        tuple val(sample_id),
              path("clean_R1.fq.gz"),
              path("clean_R2.fq.gz"),
              path("cutadapt.log"),   emit: cutadapt_out
        path "cutadapt.log",          emit: log

    script:
        """
        ${params.cutadapt_bin} \
            -a ${params.adaptor_r1} \
            -A ${params.adaptor_r2} \
            ${params.cutadapt_opts} \
            -j ${task.cpus} \
            -o clean_R1.fq.gz \
            -p clean_R2.fq.gz \
            $r1 $r2 \
            > cutadapt.log 2>&1
        """
}

/*
────────────────────────────────────────────────────────────────────────────────
Process: PEAR
────────────────────────────────────────────────────────────────────────────────
*/
process PEAR {
    tag { sample_id }

    input:
        tuple val(sample_id), path(clean_r1), path(clean_r2)

    output:
        tuple val(sample_id),
              path("assembled.assembled.fastq"),
              emit: assembled
        path "assembled.pear_stats.txt",
              emit: stats

    script:
        """
        ${params.pear_bin} -j ${task.cpus} \
              -f $clean_r1 \
              -r $clean_r2 \
              -n ${params.pear_min_ovlp} \
              -o assembled \
              -y 2G \
              > assembled.pear_stats.txt 2>&1
        """
}

/*
────────────────────────────────────────────────────────────────────────────────
Process: FASTQ2COUNT
────────────────────────────────────────────────────────────────────────────────
*/
process FASTQ2COUNT {
    tag { sample_id }
    publishDir "${params.outdir}/${sample_id}", mode: 'copy',
        pattern: "reads.count.tsv"

    input:
        tuple val(sample_id), path(assembled_fq)

    output:
        tuple val(sample_id),
              path("reads.count.tsv"),
              emit: counts

    script:
        """
        ${params.perl_bin} ${params.bin_dir}/fastq2count.pl \
            < $assembled_fq > reads.count.tsv
        """
}

/*
────────────────────────────────────────────────────────────────────────────────
Process: MIRNA_COUNT
────────────────────────────────────────────────────────────────────────────────
*/
process MIRNA_COUNT {
    tag { sample_id }

    input:
        tuple val(sample_id), path(counts_tsv)
        path mirna_fa

    output:
        tuple val(sample_id),
              path("mirna.count.tsv"),
              emit: mirna_counts

    script:
        """
        ${params.perl_bin} ${params.bin_dir}/mirna_count.pl \
            $counts_tsv $mirna_fa > mirna.count.tsv
        """
}

/*
────────────────────────────────────────────────────────────────────────────────
Process: ANNOTATE_SAMPLE
────────────────────────────────────────────────────────────────────────────────
*/
process ANNOTATE_SAMPLE {
    tag { sample_id }
    // Each sample gets a unique filename — prevents Nextflow "name collision"
    // when all sample outputs are later collected into MERGE_SAMPLES.
    publishDir "${params.outdir}/${sample_id}", mode: 'copy',
        pattern: "mirna.${sample_id}.annotated.tsv"

    input:
        tuple val(sample_id), path(mirna_tsv)

    output:
        tuple val(sample_id),
              path("mirna.${sample_id}.annotated.tsv"),
              emit: annotated

    script:
        """
        awk 'NR==1{print "Name\\t"sample"\\t"sample"_polyA"} {gsub(/^>/,""); print}' sample="${sample_id}" $mirna_tsv \
            > mirna.${sample_id}.annotated.tsv
        """
}

/*
────────────────────────────────────────────────────────────────────────────────
Process: MERGE_SAMPLES
────────────────────────────────────────────────────────────────────────────────
*/
process MERGE_SAMPLES {
    tag { "MERGE" }
    publishDir params.outdir, mode: 'copy', overwrite: true

    input:
        // Pass annotated TSV paths so Nextflow stages them into the work dir.
        // ANNOTATE_SAMPLE channel dependency ensures MERGE waits for all samples.
        path annotated_files
        path samplesheet

    output:
        path "mirna.merged.tsv", emit: merged_table

    script:
        // Stage all annotated files, ls to find them, merge in sample order.
        // Using shell $(ls ...) — works in Nextflow script: """ blocks.
        """
        TSVS=(\$(ls mirna.*.annotated.tsv 2>/dev/null | sort))
        n=\${#TSVS[@]}
        echo "Merging \$n samples"

        if [ \$n -eq 0 ]; then
            echo "ERROR: no sample files found" >&2; exit 1
        fi

        # Header + content from the first sample
        cp "\${TSVS[0]}" mirna.merged.tsv

        # Iteratively left-join each subsequent sample
        for ((i=1; i<\$n; i++)); do
            ${params.tsvjoin_bin} -H --data-fields "Name" --key-fields "Name" \\
                -f "\${TSVS[\$i]}" mirna.merged.tsv -a 2,3 > mirna.merged.tmp.tsv
            mv mirna.merged.tmp.tsv mirna.merged.tsv
        done

        echo "Done — merged \$n samples"
        """
}

/*
────────────────────────────────────────────────────────────────────────────────
Process: QC_SAMPLE
────────────────────────────────────────────────────────────────────────────────
*/
process QC_SAMPLE {
    tag { sample_id }
    publishDir "${params.outdir}/${sample_id}", mode: 'copy',
        pattern: "*.json"
    publishDir "${params.outdir}/${sample_id}", mode: 'copy',
        pattern: "*_qc_report.html"

    input:
        tuple val(sample_id), path(clean_r1), path(clean_r2), path(cutadapt_log)

    output:
        tuple val(sample_id),  path("${sample_id}_qc_metrics.json"),  emit: qc_json
        path "${sample_id}_qc_report.html",                          emit: qc_html
        path "${sample_id}_quality_by_pos.json",                    emit: qc_qual
        path "${sample_id}_length_dist.json",                       emit: qc_len

    script:
        """
        ${params.py_bin} ${params.bin_dir}/snipr_qc.py \
            --fq $clean_r1 \
            --sample $sample_id \
            --outdir . \
            --adapter ${params.adaptor_r1} \
            --min_len 10 \
            --max_len 50
        """
}

/*
────────────────────────────────────────────────────────────────────────────────
Process: MULTIQC
────────────────────────────────────────────────────────────────────────────────
*/
process MULTIQC {
    tag { "MULTIQC" }
    publishDir params.outdir, mode: 'copy', overwrite: true

    input:
        path qc_jsons
        path samplesheet
        val qc_outdir

    output:
        path "multiqc_report.html",  emit: report
        path "multiqc_data/",         emit: data_dir

    script:
        """
        ${params.py_bin} ${params.bin_dir}/aggregate_qc.py \
            --samplesheet $samplesheet \
            --qc_dir "${qc_outdir}" \
            --outdir .
        """
}

/*
────────────────────────────────────────────────────────────────────────────────
Workflow
────────────────────────────────────────────────────────────────────────────────
*/
workflow {
    // 1. Cutadapt
    CUTADAPT(input_samples)

    // 2. PEAR merge
    PEAR(CUTADAPT.out.cutadapt_out
        .map { sid, r1, r2, log -> [ sid, r1, r2 ] })

    // 3. FASTQ → count table
    FASTQ2COUNT(PEAR.out.assembled)

    // 4. miRNA counting + poly-A detection
    mirna_fa_file = file(params.mirna_fa, checkExists: true)
    MIRNA_COUNT(FASTQ2COUNT.out.counts, mirna_fa_file)

    // 5. Annotate per-sample TSV (unique filenames to avoid name collision)
    ANNOTATE_SAMPLE(MIRNA_COUNT.out.mirna_counts)

    // 6. Merge all samples
    //    Use broadcast channels: ANNOTATE_SAMPLE outputs are collected once
    //    (each sample has a unique filename, so no name collision).
    //    The collect() waits for all ANNOTATE_SAMPLE jobs to complete before
    //    MERGE_SAMPLES starts.
    MERGE_SAMPLES(
        ANNOTATE_SAMPLE.out.annotated
            .map { sid, tsv -> tsv }
            .collect(),
        file(params.samplesheet)
    )

    // 7. QC — [sid, clean_r1, clean_r2, cutadapt_log]
    // cutadapt_out already has [sid, r1, r2, log], just use it directly
    QC_SAMPLE(CUTADAPT.out.cutadapt_out)

    // 8. MULTIQC — waits for QC_SAMPLE to finish, then aggregates all QC reports
    // Resolve outdir to absolute path so Python can find JSONs from any work dir.
    def qc_outdir = params.outdir.startsWith('/') ? params.outdir : new File(params.outdir).absolutePath
    MULTIQC(
        QC_SAMPLE.out.qc_json
            .map { sid, f -> f }
            .collect(),
        file(params.samplesheet),
        qc_outdir
    )
}

/*
────────────────────────────────────────────────────────────────────────────────
Completion
────────────────────────────────────────────────────────────────────────────────
*/
workflow.onComplete {
    log.info ""
    log.info "=================================================="
    log.info "SNIPR-Nm run complete"
    log.info "Results : ${params.outdir}"
    log.info "Merged  : ${params.outdir}/mirna.merged.tsv"
    log.info "QC html : ${params.outdir}/multiqc_report.html"
    log.info "=================================================="
}
