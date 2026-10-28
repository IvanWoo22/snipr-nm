#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Quote paths and structured arguments. cutadapt_opts is an explicitly supplied shell option string.
def q(value) { "'" + value.toString().replace("'", "'\"'\"'") + "'" }

if (params.help) {
    log.info '''SNIPR-Nm: Cutadapt → PEAR → exact/polyA counts and staged QC
Usage: nextflow run SNIPR-Nm.nf --samplesheet samples.csv --mirna_fa reference.fa --outdir results
CSV: sample_id,R1,R2; optional R1_md5,R2_md5. Paths are relative to the CSV.
No species/reference is assumed. Every FASTA entry, including external controls, is retained.
Use -resume with the same work directory. See README_zh.md for resources and output definitions.'''
    exit 0
}
if (!params.mirna_fa) error 'An explicit --mirna_fa reference FASTA is required'
if (params.pear_min_ovlp < 1 || params.pear_min_length < 1 || params.qc_sample_reads < 1) {
    error 'Overlap, assembled length and QC sample size must be positive'
}

process PREPARE {
    publishDir params.outdir, mode: 'copy', pattern: 'reference_annotations.tsv'
    input:
    path sheet
    path reference
    path scripts, stageAs: 'tools/*'
    output:
    path 'validated_samples.tsv', emit: samples
    path 'reference_annotations.tsv', emit: annotations
    script:
    """
    ${q(params.py_bin)} tools/pipeline_utils.py prepare ${q(sheet)} ${q(reference)} validated_samples.tsv
    """
}

process INPUT_AUDIT {
    tag { sample_id }
    publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: 'input_manifest.json'
    input:
    tuple val(sample_id), path(r1, stageAs:'raw1/*'), path(r2, stageAs:'raw2/*'), val(md51), val(md52)
    path scripts, stageAs: 'tools/*'
    output:
    tuple val(sample_id), path(r1), path(r2), path('input_manifest.json')
    script:
    """
    ${q(params.py_bin)} tools/pipeline_utils.py audit input_manifest.json ${q(r1)} ${q(r2)} --md5 ${q(md51)} ${q(md52)}
    """
}

process CUTADAPT {
    tag { sample_id }
    publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: 'cutadapt.{log,json}'
    input:
    tuple val(sample_id), path(r1, stageAs:'raw1/*'), path(r2, stageAs:'raw2/*'), path(manifest)
    path scripts, stageAs: 'tools/*'
    output:
    tuple val(sample_id), path(r1), path(r2), path('clean_R1.fq.gz'), path('clean_R2.fq.gz'), path('cutadapt.json'), emit: clean
    path 'cutadapt.log'
    script:
    """
    ${q(params.py_bin)} tools/pipeline_utils.py verify ${q(manifest)} ${q(r1)} ${q(r2)}
    ${q(params.cutadapt_bin)} -a ${q(params.adaptor_r1)} -A ${q(params.adaptor_r2)} \
        ${params.cutadapt_opts} -j ${task.cpus} --json cutadapt.json \
        -o clean_R1.fq.gz -p clean_R2.fq.gz ${q(r1)} ${q(r2)} > cutadapt.log 2>&1
    ${q(params.py_bin)} tools/pipeline_utils.py verify ${q(manifest)} ${q(r1)} ${q(r2)}
    """
}

process PEAR {
    tag { sample_id }
    publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: 'pear.{log,json}'
    input:
    tuple val(sample_id), path(r1), path(r2), path(cutadapt_json)
    path scripts, stageAs: 'tools/*'
    output:
    tuple val(sample_id), path('assembled.assembled.fastq.gz'), path('pear.log'), emit: assembled
    path 'pear.json'
    path 'assembled.*.fastq.gz', emit: retained
    script:
    """
    ${q(params.py_bin)} tools/pipeline_utils.py pear --binary ${q(params.pear_bin)} \
        --r1 ${q(r1)} --r2 ${q(r2)} --cutadapt-json ${q(cutadapt_json)} --threads ${task.cpus} \
        --overlap ${params.pear_min_ovlp} --min-length ${params.pear_min_length} --memory ${q(params.pear_memory)}
    mv pear_intermediates.json pear.json
    """
}

process FASTQ2COUNT {
    tag { sample_id }
    publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: 'reads.count.tsv.gz'
    input:
    tuple val(sample_id), path(fastq)
    path scripts, stageAs: 'tools/*'
    output:
    tuple val(sample_id), path('reads.count.tsv.gz')
    script:
    """
    set -o pipefail
    ${q(params.perl_bin)} tools/fastq2count.pl ${q(fastq)} | gzip -1 -n > reads.count.tsv.gz
    """
}

process MIRNA_COUNT {
    tag { sample_id }
    publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: 'mirna.count.tsv'
    input:
    tuple val(sample_id), path(counts)
    path reference
    path scripts, stageAs: 'tools/*'
    output:
    tuple val(sample_id), path('mirna.count.tsv')
    script:
    """
    ${q(params.perl_bin)} tools/mirna_count.pl --header ${q(counts)} ${q(reference)} > mirna.count.tsv
    """
}

process ANNOTATE_SAMPLE {
    tag { sample_id }
    publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: 'mirna.*.annotated.tsv'
    input:
    tuple val(sample_id), path(counts)
    path reference
    path scripts, stageAs: 'tools/*'
    output:
    path "mirna.${sample_id}.annotated.tsv"
    script:
    """
    ${q(params.py_bin)} tools/pipeline_utils.py annotate ${q(sample_id)} ${q(counts)} ${q(reference)}
    """
}

process MERGE_SAMPLES {
    publishDir params.outdir, mode: 'copy', pattern: 'mirna.*.tsv'
    input:
    path tables
    path sheet
    path reference
    path scripts, stageAs: 'tools/*'
    output:
    path 'mirna.*.tsv'
    script:
    """
    ${q(params.py_bin)} tools/pipeline_utils.py merge ${q(sheet)} ${q(reference)} ${tables.collect{q(it)}.join(' ')}
    """
}

process QC_SAMPLE {
    tag { sample_id }
    publishDir "${params.outdir}/${sample_id}", mode: 'copy', pattern: "${sample_id}_*"
    input:
    tuple val(sample_id), path(raw1, stageAs:'raw1/*'), path(raw2, stageAs:'raw2/*'), path(clean1), path(clean2), path(cutadapt_json), path(merged), path(pear_log), path(counts)
    path scripts, stageAs: 'tools/*'
    output:
    path "${sample_id}_*"
    script:
    """
    ${q(params.py_bin)} tools/pipeline_utils.py qc --sample ${q(sample_id)} \
        --fastqs ${q(raw1)} ${q(raw2)} ${q(clean1)} ${q(clean2)} ${q(merged)} \
        --cutadapt-json ${q(cutadapt_json)} --pear-log ${q(pear_log)} --counts ${q(counts)} \
        --adapter-r1 ${q(params.adaptor_r1)} --adapter-r2 ${q(params.adaptor_r2)} --sample-reads ${params.qc_sample_reads}
    """
}

process AGGREGATE_QC {
    publishDir params.outdir, mode: 'copy'
    input:
    path reports
    path sheet
    path scripts, stageAs: 'tools/*'
    output:
    path 'snipr_qc_report.html'
    path 'snipr_qc_data'
    script:
    """
    ${q(params.py_bin)} tools/aggregate_qc.py --samplesheet ${q(sheet)} --outdir . --files ${reports.collect{q(it)}.join(' ')}
    """
}

workflow {
    sheet = file(params.samplesheet, checkExists:true)
    reference = file(params.mirna_fa, checkExists:true)
    // Stage source files as task inputs so -resume invalidates changed implementations.
    scripts = file("${projectDir}/bin/*.{py,pl}").sort { it.name }
    PREPARE(sheet, reference, scripts)
    samples = PREPARE.out.samples.splitCsv(header:true, sep:'\t').map { row ->
        tuple(row.sample_id, file(row.R1,checkExists:true), file(row.R2,checkExists:true), row.R1_md5 ?: '-', row.R2_md5 ?: '-')
    }
    INPUT_AUDIT(samples, scripts)
    CUTADAPT(INPUT_AUDIT.out, scripts)
    PEAR(CUTADAPT.out.clean.map{sid,raw1,raw2,r1,r2,json -> tuple(sid,r1,r2,json)}, scripts)
    FASTQ2COUNT(PEAR.out.assembled.map{sid,fq,log -> tuple(sid,fq)}, scripts)
    MIRNA_COUNT(FASTQ2COUNT.out, reference, scripts)
    ANNOTATE_SAMPLE(MIRNA_COUNT.out, reference, scripts)
    MERGE_SAMPLES(ANNOTATE_SAMPLE.out.collect().map { it.sort { f -> f.name } }, sheet, reference, scripts)
    QC_SAMPLE(CUTADAPT.out.clean.join(PEAR.out.assembled).join(FASTQ2COUNT.out), scripts)
    AGGREGATE_QC(QC_SAMPLE.out.flatten().collect().map { it.sort { f -> f.name } }, sheet, scripts)
}

workflow.onComplete {
    if (workflow.success) {
        log.info "SNIPR-Nm completed successfully: ${params.outdir}/mirna.merged.tsv and snipr_qc_report.html"
    } else {
        log.error 'SNIPR-Nm failed; published files may be partial. Inspect the failed task before resuming.'
    }
}
