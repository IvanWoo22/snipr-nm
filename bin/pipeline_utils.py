#!/usr/bin/env python3
"""Validation, compression, count export and stage accounting for the Nextflow workflow."""
import argparse
import csv
import gzip
import hashlib
import json
import subprocess
from pathlib import Path
from snipr_common import references, samplesheet, fingerprint, write_json, write_html, open_text
from snipr_qc import compute_metrics


def prepare(sheet, fasta, out):
    rows = samplesheet(sheet, check_files=True)
    refs = references(fasta)
    with open(out, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=['sample_id','R1','R2','R1_md5','R2_md5'], delimiter='\t')
        writer.writeheader()
        for row in rows:
            for key in ('R1_md5','R2_md5'):
                value = row.get(key, '') or ''
                if value and (len(value) != 32 or any(c not in '0123456789abcdefABCDEF' for c in value)):
                    raise ValueError(f'{row["sample_id"]}: invalid {key}')
                row[key] = value.lower()
            writer.writerow({k:row.get(k,'') for k in writer.fieldnames})
    with open('reference_annotations.tsv','w',newline='') as fh:
        writer = csv.writer(fh,delimiter='\t'); writer.writerow(['Name','sequence'])
        writer.writerows(refs)


def audit(paths, expected, out):
    records = [fingerprint(p) for p in paths]
    for record, md5 in zip(records, expected):
        if md5 != '-' and record['md5'] != md5.lower():
            raise ValueError(f'MD5 mismatch: {record["path"]}')
    write_json(out, dict(schema_version=1, inputs=records))


def verify(manifest, paths):
    records = json.loads(Path(manifest).read_text())['inputs']
    if len(records) != len(paths):
        raise ValueError('Input manifest count mismatch')
    for old, path in zip(records, paths):
        now = fingerprint(path)
        if any(old[k] != now[k] for k in ('sha256','md5','bytes')):
            raise ValueError(f'Input changed since audit: {path}')


def compress_pear():
    """Only replace these four task-local outputs after a gzip round-trip checksum."""
    records = []
    for suffix in ('assembled','unassembled.forward','unassembled.reverse','discarded'):
        source = Path(f'assembled.{suffix}.fastq')
        if source.is_symlink() or not source.is_file():
            raise ValueError(f'Missing or nonlocal PEAR output: {source}')
        before = source.stat()
        temp = Path(str(source)+'.gz.partial'); dest = Path(str(source)+'.gz')
        original = hashlib.sha256()
        with source.open('rb') as inp, temp.open('wb') as target:
            with gzip.GzipFile(fileobj=target, mode='wb', compresslevel=1, mtime=0, filename='') as out:
                for block in iter(lambda: inp.read(4*1024*1024), b''):
                    original.update(block); out.write(block)
        restored = hashlib.sha256()
        with gzip.open(temp,'rb') as inp:
            for block in iter(lambda: inp.read(4*1024*1024), b''):
                restored.update(block)
        after = source.stat()
        if original.digest() != restored.digest() or (before.st_size,before.st_mtime_ns) != (after.st_size,after.st_mtime_ns):
            raise ValueError(f'Compression validation failed: {source}')
        temp.replace(dest)
        records.append(dict(file=str(dest), uncompressed_bytes=after.st_size,
                            compressed_bytes=dest.stat().st_size, uncompressed_sha256=original.hexdigest()))
        source.unlink()
    write_json('pear_intermediates.json', records)


def run_pear(args):
    clean_pairs = json.loads(Path(args.cutadapt_json).read_text())['read_counts']['output']
    if clean_pairs == 0:
        # PEAR 0.9.11 aborts on empty input. Preserve a valid all-zero sample.
        for suffix in ('assembled','unassembled.forward','unassembled.reverse','discarded'):
            Path(f'assembled.{suffix}.fastq').write_bytes(b'')
        Path('pear.log').write_text('PEAR not invoked: Cutadapt retained zero pairs.\n'
            'Assembled reads: 0\nDiscarded reads: 0\nNot assembled reads: 0\n')
    else:
        command = [args.binary,'-j',str(args.threads),'-f',args.r1,'-r',args.r2,
                   '-v',str(args.overlap),'-n',str(args.min_length),'-y',args.memory,'-o','assembled']
        with open('pear.log','w') as log:
            subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)
    compress_pear()


def load_counts(path, ref_names, expected_header):
    with open(path, newline='') as fh:
        rows = list(csv.reader(fh, delimiter='\t'))
    if not rows or rows[0] != expected_header:
        raise ValueError(f'Wrong count header: {path}')
    if [r[0] for r in rows[1:] if r] != ref_names:
        raise ValueError(f'Reference names/order mismatch: {path}')
    counts = []
    for row in rows[1:]:
        if len(row) != 3 or any(not v.isdigit() for v in row[1:]):
            raise ValueError(f'Invalid count row: {path}')
        counts.append([int(x) for x in row[1:]])
    return counts


def annotate(sample, counts, fasta):
    names = [name for name,_ in references(fasta)]
    values = load_counts(counts, names, ['Name','exact','polyA'])
    with open(f'mirna.{sample}.annotated.tsv','w',newline='') as fh:
        writer = csv.writer(fh,delimiter='\t'); writer.writerow(['Name',sample,sample+'_polyA'])
        writer.writerows([name,*v] for name,v in zip(names,values))


def merge(sheet, fasta, files):
    ids = [row['sample_id'] for row in samplesheet(sheet)]
    refs = references(fasta); names = [name for name,_ in refs]
    by_file = {Path(f).name: f for f in files}
    expected = {f'mirna.{sid}.annotated.tsv' for sid in ids}
    if len(by_file) != len(files) or set(by_file) != expected:
        raise ValueError('Missing, duplicate or unexpected sample count files')
    counts = {sid:load_counts(by_file[f'mirna.{sid}.annotated.tsv'],names,['Name',sid,sid+'_polyA']) for sid in ids}
    with open('mirna.merged.tsv','w',newline='') as fh:
        writer = csv.writer(fh,delimiter='\t')
        writer.writerow(['Name',*[col for sid in ids for col in (sid,sid+'_polyA')]])
        for index,name in enumerate(names):
            writer.writerow([name,*[v for sid in ids for v in counts[sid][index]]])
    for metric in ('exact','polyA','classified_total','no_tail_fraction'):
        with open(f'mirna.{metric}.tsv','w',newline='') as fh:
            writer = csv.writer(fh,delimiter='\t'); writer.writerow(['Name','sequence',*ids])
            for index,(name,sequence) in enumerate(refs):
                values = []
                for sid in ids:
                    e,p = counts[sid][index]; total = e+p
                    values.append({'exact':e,'polyA':p,'classified_total':total,
                                   'no_tail_fraction':e/total if total else 'NA'}[metric])
                writer.writerow([name,sequence,*values])


def qc(args):
    stages = {}
    qualities, lengths = {}, {}
    for label,path in zip(('raw_R1','raw_R2','clean_R1','clean_R2','merged'),args.fastqs):
        adapter = args.adapter_r2 if label.endswith('R2') else args.adapter_r1
        metric = compute_metrics(path, adapter, sample_reads=args.sample_reads)
        qualities[label] = metric.pop('quality_by_position')
        lengths[label] = metric.pop('length_distribution')
        stages[label] = metric
    cut = json.loads(Path(args.cutadapt_json).read_text())['read_counts']
    n = {k:v['total_reads'] for k,v in stages.items()}
    if n['raw_R1'] != n['raw_R2'] or n['raw_R1'] != cut['input']:
        raise ValueError('Raw paired read count / Cutadapt input mismatch')
    if n['clean_R1'] != n['clean_R2'] or n['clean_R1'] != cut['output']:
        raise ValueError('Clean paired read count / Cutadapt output mismatch')
    count_sum = 0
    with open_text(args.counts) as fh:
        for line in fh:
            seq,count = line.rstrip('\n').split('\t')
            if not count.isdigit(): raise ValueError('Invalid sequence count')
            count_sum += int(count)
    if count_sum != n['merged'] or n['merged'] > n['clean_R1']:
        raise ValueError('Merged read / sequence count conservation failed')
    # PEAR counts discarded pairs once; unmatched reads are one file per mate.
    import re
    log = Path(args.pear_log).read_text()
    pear_counts = {}
    for key in ('Assembled','Discarded','Not assembled'):
        match = re.search(r'^'+key+r' reads[ .]*:\s*([\d,]+)', log, re.M)
        if not match: raise ValueError(f'Missing PEAR accounting: {key}')
        pear_counts[key] = int(match[1].replace(',',''))
    if pear_counts['Assembled'] != n['merged'] or sum(pear_counts.values()) != n['clean_R1']:
        raise ValueError('PEAR pair accounting failed')
    result = dict(schema_version=2, sample_id=args.sample, stages=stages, pear_pairs=pear_counts,
                  sequence_count_sum=count_sum, accounting_pass=True)
    write_json(f'{args.sample}_qc_metrics.json',result)
    write_json(f'{args.sample}_quality_by_pos.json',qualities)
    write_json(f'{args.sample}_length_dist.json',lengths)
    headers = ['Stage','Reads','Bases','Mean quality','Q30 %','18-25 nt %','Sampled reads','Sample redundancy %']
    rows = [[label,*[m[k] for k in ('total_reads','total_bases','mean_quality','q30_pct',
             'mirna_18_25nt_pct','sampled_reads','sequence_redundancy_pct')]] for label,m in stages.items()]
    write_html(f'{args.sample}_qc_report.html',f'SNIPR QC: {args.sample}',headers,rows,
        'Full-file counts and global quality; detailed position quality, entropy, motif and redundancy use the first N reads. '
        'Raw/clean counts are reads per mate; merged counts are merged pairs. Sequence redundancy is not PCR duplication.')


def main():
    p = argparse.ArgumentParser(description=__doc__); sub = p.add_subparsers(dest='command',required=True)
    s=sub.add_parser('prepare'); s.add_argument('sheet');s.add_argument('fasta');s.add_argument('out')
    s=sub.add_parser('audit');s.add_argument('out');s.add_argument('paths',nargs=2);s.add_argument('--md5',nargs=2,default=['-','-'])
    s=sub.add_parser('verify');s.add_argument('manifest');s.add_argument('paths',nargs=2)
    sub.add_parser('compress-pear')
    s=sub.add_parser('pear');s.add_argument('--binary',required=True);s.add_argument('--r1',required=True);s.add_argument('--r2',required=True)
    s.add_argument('--cutadapt-json',required=True);s.add_argument('--threads',type=int,required=True)
    s.add_argument('--overlap',type=int,required=True);s.add_argument('--min-length',type=int,required=True);s.add_argument('--memory',required=True)
    s=sub.add_parser('annotate');s.add_argument('sample');s.add_argument('counts');s.add_argument('fasta')
    s=sub.add_parser('merge');s.add_argument('sheet');s.add_argument('fasta');s.add_argument('files',nargs='+')
    s=sub.add_parser('qc');s.add_argument('--sample',required=True);s.add_argument('--fastqs',nargs=5,required=True)
    s.add_argument('--cutadapt-json',required=True);s.add_argument('--pear-log',required=True);s.add_argument('--counts',required=True)
    s.add_argument('--adapter-r1',required=True);s.add_argument('--adapter-r2',required=True)
    s.add_argument('--sample-reads',type=int,default=100000)
    a=p.parse_args()
    if a.command=='prepare': prepare(a.sheet,a.fasta,a.out)
    elif a.command=='audit': audit(a.paths,a.md5,a.out)
    elif a.command=='verify': verify(a.manifest,a.paths)
    elif a.command=='compress-pear': compress_pear()
    elif a.command=='pear': run_pear(a)
    elif a.command=='annotate': annotate(a.sample,a.counts,a.fasta)
    elif a.command=='merge': merge(a.sheet,a.fasta,a.files)
    elif a.command=='qc': qc(a)


if __name__=='__main__': main()
