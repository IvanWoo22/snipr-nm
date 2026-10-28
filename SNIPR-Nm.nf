#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

def q(value) { "'" + value.toString().replace("'", "'\"'\"'") + "'" }
def deps(names) { names.collect { file("${projectDir}/bin/${it}", checkExists:true) } }

if (params.help) {
    log.info '''SNIPR-Nm: reference-independent preprocessing, native metrics and exact/polyA counts
Raw: nextflow run SNIPR-Nm.nf --samplesheet samples.csv --mirna_fa reference.fa
Counts: nextflow run SNIPR-Nm.nf --counts_sheet checkpoints.csv --mirna_fa reference.fa
Modes: --intermediates checkpoint (default) | compact. Use -resume with the same work directory.
Counts CSV: sample_id,counts,manifest. See README_zh.md for profiles and recovery boundaries.'''
    exit 0
}
if (!params.mirna_fa) error 'An explicit --mirna_fa reference FASTA is required'
if (!(params.intermediates in ['checkpoint','compact'])) error 'intermediates must be checkpoint or compact'
if (params.pear_min_ovlp < 1 || params.pear_min_length < 1 || params.qc_sample_reads < 1 || params.count_max_unique < 1) error 'Bounds must be positive'

process PREPARE_SAMPLES {
    label 'light'
    input: path sheet; path scripts, stageAs:'tools/*'
    output: path 'validated_samples.tsv'
    script: """${q(params.py_bin)} tools/input_tools.py prepare ${q(sheet)}"""
}
process PREPARE_COUNTS {
    label 'light'
    input: path sheet; path scripts, stageAs:'tools/*'
    output: path 'validated_samples.tsv'
    script: """${q(params.py_bin)} tools/checkpoint_tools.py prepare ${q(sheet)}"""
}
process PREPARE_REFERENCE {
    label 'light'
    publishDir params.outdir, mode:'copy',overwrite:true, pattern:'reference_annotations.tsv'
    input: path fasta; path scripts, stageAs:'tools/*'
    output: path fasta, emit: reference; path 'reference_annotations.tsv'
    script: """${q(params.py_bin)} tools/export_tools.py reference ${q(fasta)}"""
}
process BUILD_SCANNER {
    label 'light'
    input: path source
    output: path 'fastq_metrics'
    script: """${q(params.cxx_bin)} -O3 -std=c++17 ${q(source)} -lz -o fastq_metrics"""
}
process INPUT_AUDIT {
    tag { sample_id }; label 'audit'
    publishDir "${params.outdir}/${sample_id}",mode:'copy',overwrite:true,pattern:'input_manifest.json'
    input:
    tuple val(sample_id),path(r1,stageAs:'raw1/*'),path(r2,stageAs:'raw2/*'),val(md51),val(md52)
    path scripts,stageAs:'tools/*'
    output: tuple val(sample_id),path(r1),path(r2),path('input_manifest.json')
    script: """${q(params.py_bin)} tools/input_tools.py audit input_manifest.json ${q(r1)} ${q(r2)} --md5 ${q(md51)} ${q(md52)}"""
}
process RAW_METRICS {
    tag { sample_id }; label 'qc_scan'
    input:
    tuple val(sample_id),path(r1,stageAs:'raw1/*'),path(r2,stageAs:'raw2/*'),path(manifest)
    path scanner; path scripts,stageAs:'tools/*'
    output: tuple val(sample_id),path('raw.metrics.json')
    script:
    """
    ${q(params.py_bin)} tools/stage_metrics.py pair --scanner ${q(scanner)} --prefix raw --r1 ${q(r1)} --r2 ${q(r2)} \
      --adapter-r1 ${q(params.adaptor_r1)} --adapter-r2 ${q(params.adaptor_r2)} --sample-reads ${q(params.qc_sample_reads)}
    """
}
process CUTADAPT {
    tag { sample_id }; label 'trim'
    publishDir "${params.outdir}/${sample_id}",mode:'copy',overwrite:true,pattern:'cutadapt.{log,json}'
    input:
    tuple val(sample_id),path(r1,stageAs:'raw1/*'),path(r2,stageAs:'raw2/*'),path(manifest)
    path scripts,stageAs:'tools/*'
    output:
    tuple val(sample_id),path('clean_R1.fq.gz'),path('clean_R2.fq.gz'),path('cutadapt.json'),emit:clean
    path 'cutadapt.log'
    script:
    """
    ${q(params.py_bin)} tools/input_tools.py verify ${q(manifest)} ${q(r1)} ${q(r2)}
    ${q(params.py_bin)} tools/input_tools.py trim --binary ${q(params.cutadapt_bin)} --r1 ${q(r1)} --r2 ${q(r2)} \
      --adapter-r1 ${q(params.adaptor_r1)} --adapter-r2 ${q(params.adaptor_r2)} --options=${q(params.cutadapt_opts)} --threads ${task.cpus}
    ${q(params.py_bin)} tools/input_tools.py verify ${q(manifest)} ${q(r1)} ${q(r2)}
    """
}
process CLEAN_METRICS {
    tag { sample_id }; label 'qc_scan'
    input:
    tuple val(sample_id),path(r1),path(r2),path(cutadapt_json)
    path scanner; path scripts,stageAs:'tools/*'
    output: tuple val(sample_id),path('clean.metrics.json')
    script:
    """
    ${q(params.py_bin)} tools/stage_metrics.py pair --scanner ${q(scanner)} --prefix clean --r1 ${q(r1)} --r2 ${q(r2)} \
      --adapter-r1 ${q(params.adaptor_r1)} --adapter-r2 ${q(params.adaptor_r2)} --sample-reads ${q(params.qc_sample_reads)}
    """
}
process PEAR {
    tag { sample_id }; label 'merge_reads'
    publishDir "${params.outdir}/${sample_id}",mode:'copy',overwrite:true,pattern:'pear.{log,json}'
    input:
    tuple val(sample_id),path(r1),path(r2),path(cutadapt_json)
    path scripts,stageAs:'tools/*'
    output:
    tuple val(sample_id),path('assembled.assembled.fastq.gz'),path('pear.log'),emit:assembled
    path 'pear.json'
    path 'assembled.*.fastq.gz',emit:retained
    script:
    """
    ${q(params.py_bin)} tools/pear_tools.py --binary ${q(params.pear_bin)} --r1 ${q(r1)} --r2 ${q(r2)} \
      --cutadapt-json ${q(cutadapt_json)} --threads ${task.cpus} --overlap ${q(params.pear_min_ovlp)} \
      --min-length ${q(params.pear_min_length)} --memory ${q(params.pear_memory)}
    mv pear_intermediates.json pear.json
    """
}
process MERGED_COUNT_METRICS {
    tag { sample_id }; label 'exact_count'
    input:
    tuple val(sample_id),path(fastq),path(pear_log)
    path scanner; path scripts,stageAs:'tools/*'
    output: tuple val(sample_id),path('reads.count.tsv.gz'),path('counting.json'),path('merged.metrics.json')
    script:
    """
    ${q(params.py_bin)} tools/stage_metrics.py merged --scanner ${q(scanner)} --fastq ${q(fastq)} \
      --adapter-r1 ${q(params.adaptor_r1)} --sample-reads ${q(params.qc_sample_reads)} --max-unique ${q(params.count_max_unique)}
    """
}
process ASSEMBLE_QC {
    tag { sample_id }; label 'light'
    publishDir "${params.outdir}/${sample_id}",mode:'copy',overwrite:true,pattern:"${sample_id}_*.json"
    input:
    tuple val(sample_id),path(raw),path(clean),path(counts),path(counting),path(merged),path(cutadapt_json),path(pear_log)
    path scripts,stageAs:'tools/*'
    output: tuple val(sample_id),path("${sample_id}_qc_metrics.json"),path("${sample_id}_length_dist.json"),path("${sample_id}_quality_by_pos.json")
    script:
    """
    ${q(params.py_bin)} tools/qc_assemble.py --sample ${q(sample_id)} --raw ${q(raw)} --clean ${q(clean)} \
      --merged ${q(merged)} --cutadapt ${q(cutadapt_json)} --pear ${q(pear_log)} --counting ${q(counting)}
    """
}
process SEAL_COUNTS {
    tag { sample_id }; label 'light'
    publishDir "${params.outdir}/${sample_id}",mode:'copy',overwrite:true,pattern:'{reads.count.tsv.gz,counts.manifest.json}'
    input:
    tuple val(sample_id),path(counts),path(counting),path(qc),path(lengths),path(qualities),path(manifest),path(cutadapt_json),path(pear_log)
    val parameters; path scanner; path scripts,stageAs:'tools/*'
    output: tuple val(sample_id),path(counts),path('counts.manifest.json')
    script:
    """
    ${q(params.py_bin)} tools/checkpoint_tools.py seal --sample ${q(sample_id)} --counts ${q(counts)} --counting ${q(counting)} \
      --qc ${q(qc)} --inputs ${q(manifest)} --parameters ${q(parameters)} --scanner ${q(scanner)}
    """
}
process COMPACT_SAMPLE {
    tag { sample_id }; label 'compact'
    publishDir "${params.outdir}/${sample_id}",mode:'copy',overwrite:true,pattern:'{reads.count.tsv.gz,counts.manifest.json,*_*.json,cutadapt.json,cutadapt.log,pear.log,intermediates.json}'
    input:
    tuple val(sample_id),path(r1,stageAs:'raw1/*'),path(r2,stageAs:'raw2/*'),path(manifest),path(raw_metrics)
    val parameters; path scanner; path scripts,stageAs:'tools/*'
    output:
    tuple val(sample_id),path('reads.count.tsv.gz'),path('counts.manifest.json'),emit:counts
    tuple val(sample_id),path("${sample_id}_qc_metrics.json"),path("${sample_id}_length_dist.json"),path("${sample_id}_quality_by_pos.json"),emit:qc
    path 'intermediates.json'
    path 'cutadapt.json'
    path 'cutadapt.log'
    path 'pear.log'
    script:
    """
    ${q(params.py_bin)} tools/compact_sample.py --sample ${q(sample_id)} --r1 ${q(r1)} --r2 ${q(r2)} \
      --manifest ${q(manifest)} --raw-metrics ${q(raw_metrics)} --parameters ${q(parameters)} --scanner ${q(scanner)} \
      --cutadapt ${q(params.cutadapt_bin)} --pear ${q(params.pear_bin)} --pear-memory ${q(params.pear_memory)} \
      --threads ${task.cpus} --sample-reads ${q(params.qc_sample_reads)} --max-unique ${q(params.count_max_unique)}
    """
}
process VERIFY_COUNTS {
    tag { sample_id }; label 'light'
    input: tuple val(sample_id),path(counts),path(manifest); path scripts,stageAs:'tools/*'
    output: tuple val(sample_id),path(counts),path(manifest)
    script: """${q(params.py_bin)} tools/checkpoint_tools.py verify ${q(sample_id)} ${q(counts)} ${q(manifest)}"""
}
process MIRNA_COUNT {
    tag { sample_id }; label 'light'
    publishDir "${params.outdir}/${sample_id}",mode:'copy',overwrite:true,pattern:'mirna.*.tsv'
    input: tuple val(sample_id),path(counts),path(manifest); path reference; path scripts,stageAs:'tools/*'
    output: path "mirna.${sample_id}.annotated.tsv", emit: annotated
    path 'mirna.count.tsv', emit: counts
    script:
    """
    ${q(params.perl_bin)} tools/mirna_count.pl --header ${q(counts)} ${q(reference)} > mirna.count.tsv
    ${q(params.py_bin)} tools/export_tools.py annotate ${q(sample_id)} mirna.count.tsv ${q(reference)}
    """
}
process MERGE_SAMPLES {
    label 'light'
    publishDir params.outdir,mode:'copy',overwrite:true,pattern:'mirna.*.tsv'
    input: path tables; path sheet; path reference; path scripts,stageAs:'tools/*'
    output: path 'mirna.*.tsv'
    script: """${q(params.py_bin)} tools/export_tools.py merge ${q(sheet)} ${q(reference)} ${tables.collect{q(it)}.join(' ')}"""
}
process RENDER_QC {
    tag { sample_id }; label 'light'
    publishDir "${params.outdir}/${sample_id}",mode:'copy',overwrite:true,pattern:"${sample_id}_qc_report.html"
    input: tuple val(sample_id),path(metrics),path(lengths),path(qualities); path scripts,stageAs:'tools/*'
    output: tuple val(sample_id),path("${sample_id}_qc_report.html")
    script: """${q(params.py_bin)} tools/render_sample_qc.py ${q(sample_id)} ${q(metrics)}"""
}
process AGGREGATE_QC {
    label 'light'
    publishDir params.outdir,mode:'copy',overwrite:true
    input: path reports; path sheet; path scripts,stageAs:'tools/*'
    output: path 'snipr_qc_report.html'; path 'snipr_qc_data'
    script: """${q(params.py_bin)} tools/aggregate_qc.py --samplesheet ${q(sheet)} --outdir . --files ${reports.collect{q(it)}.join(' ')}"""
}

workflow {
    reference=file(params.mirna_fa,checkExists:true)
    input_deps=deps(['input_tools.py','snipr_common.py'])
    export_deps=deps(['export_tools.py','snipr_common.py'])
    metrics_deps=deps(['stage_metrics.py','snipr_common.py'])
    checkpoint_deps=deps(['checkpoint_tools.py','snipr_common.py'])
    seal_deps=deps(['checkpoint_tools.py','snipr_common.py','input_tools.py','pear_tools.py','stage_metrics.py','qc_assemble.py'])
    PREPARE_REFERENCE(reference,export_deps)
    if(params.counts_sheet) {
        sheet=file(params.counts_sheet,checkExists:true)
        PREPARE_COUNTS(sheet,checkpoint_deps)
        ordered_sheet=PREPARE_COUNTS.out
        imported=ordered_sheet.splitCsv(header:true,sep:'\t').map {r -> tuple(r.sample_id,file(r.counts,checkExists:true),file(r.manifest,checkExists:true))}
        VERIFY_COUNTS(imported,checkpoint_deps)
        ready_counts=VERIFY_COUNTS.out
    } else {
        sheet=file(params.samplesheet,checkExists:true)
        PREPARE_SAMPLES(sheet,input_deps)
        ordered_sheet=PREPARE_SAMPLES.out
        samples=ordered_sheet.splitCsv(header:true,sep:'\t').map {r -> tuple(r.sample_id,file(r.R1,checkExists:true),file(r.R2,checkExists:true),r.R1_md5?:'-',r.R2_md5?:'-')}
        parameters=groovy.json.JsonOutput.toJson([adapter_r1:params.adaptor_r1,adapter_r2:params.adaptor_r2,cutadapt_opts:params.cutadapt_opts,pear_min_ovlp:params.pear_min_ovlp,pear_min_length:params.pear_min_length,qc_sample_reads:params.qc_sample_reads,count_max_unique:params.count_max_unique,intermediates:params.intermediates])
        BUILD_SCANNER(file("${projectDir}/src/fastq_metrics.cpp",checkExists:true))
        INPUT_AUDIT(samples,input_deps)
        RAW_METRICS(INPUT_AUDIT.out,BUILD_SCANNER.out,metrics_deps)
        if(params.intermediates=='compact') {
            compact_deps=deps(['compact_sample.py','input_tools.py','pear_tools.py','stage_metrics.py','qc_assemble.py','checkpoint_tools.py','snipr_common.py'])
            COMPACT_SAMPLE(INPUT_AUDIT.out.join(RAW_METRICS.out),parameters,BUILD_SCANNER.out,compact_deps)
            ready_counts=COMPACT_SAMPLE.out.counts
            numerical_qc=COMPACT_SAMPLE.out.qc
        } else {
            CUTADAPT(INPUT_AUDIT.out,input_deps)
            CLEAN_METRICS(CUTADAPT.out.clean,BUILD_SCANNER.out,metrics_deps)
            PEAR(CUTADAPT.out.clean,deps(['pear_tools.py','snipr_common.py']))
            MERGED_COUNT_METRICS(PEAR.out.assembled,BUILD_SCANNER.out,metrics_deps)
            cut_json=CUTADAPT.out.clean.map {sid,r1,r2,j -> tuple(sid,j)}
            pear_log=PEAR.out.assembled.map {sid,fq,log -> tuple(sid,log)}
            assembled_qc=RAW_METRICS.out.join(CLEAN_METRICS.out).join(MERGED_COUNT_METRICS.out).join(cut_json).join(pear_log)
            ASSEMBLE_QC(assembled_qc,deps(['qc_assemble.py','snipr_common.py']))
            numerical_qc=ASSEMBLE_QC.out
            seals=MERGED_COUNT_METRICS.out.map {sid,c,m,qc -> tuple(sid,c,m)}.join(numerical_qc).join(INPUT_AUDIT.out.map{sid,r1,r2,m -> tuple(sid,m)}).join(cut_json).join(pear_log)
            SEAL_COUNTS(seals,parameters,BUILD_SCANNER.out,seal_deps)
            ready_counts=SEAL_COUNTS.out
        }
        RENDER_QC(numerical_qc,deps(['render_sample_qc.py','report_html.py']))
        reports=numerical_qc.join(RENDER_QC.out).map {sid,m,l,q,h -> [m,l,q,h]}.flatten().collect().map{it.sort{f->f.name}}
        AGGREGATE_QC(reports,ordered_sheet,deps(['aggregate_qc.py','report_html.py']))
    }
    MIRNA_COUNT(ready_counts,PREPARE_REFERENCE.out.reference,deps(['mirna_count.pl','export_tools.py','snipr_common.py']))
    MERGE_SAMPLES(MIRNA_COUNT.out.annotated.collect().map{it.sort{f->f.name}},ordered_sheet,PREPARE_REFERENCE.out.reference,export_deps)
}
workflow.onComplete {
    if(workflow.success) log.info "SNIPR-Nm completed successfully: ${params.outdir}/mirna.merged.tsv"
    else log.error 'SNIPR-Nm failed; published files may be partial. Inspect failed tasks before resuming.'
}
