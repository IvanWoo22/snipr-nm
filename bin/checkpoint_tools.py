"""Versioned, reference-independent count checkpoints with strict import validation."""
import argparse,csv,gzip,hashlib,json,re
from pathlib import Path
from snipr_common import write_json

def validate_counts(path):
    digest=hashlib.sha256();total=unique=0;previous=None
    with gzip.open(path,'rb') as f:
        for line in f:
            digest.update(line)
            if not re.fullmatch(rb'[ACGTN]+\t[1-9][0-9]*\n',line):raise ValueError('Invalid checkpoint count row')
            seq,n=line.rstrip(b'\n').split(b'\t')
            if previous is not None and seq<=previous:raise ValueError('Checkpoint counts must be unique and sorted')
            previous=seq;total+=int(n);unique+=1
    with open(path,'rb') as f:compressed=hashlib.file_digest(f,'sha256').hexdigest()
    return dict(total_reads=total,unique_sequences=unique,content_sha256=digest.hexdigest(),compressed_sha256=compressed)

def seal(sample,counts,counting,qc,inputs,parameters,scanner):
    numbers=validate_counts(counts);meta=json.loads(Path(counting).read_text());metrics=json.loads(Path(qc).read_text())
    if any(meta.get(k)!=v for k,v in numbers.items()):raise ValueError('Count checkpoint metadata mismatch')
    if metrics.get('sample_id')!=sample or not metrics.get('accounting_pass') or metrics['sequence_count_sum']!=numbers['total_reads']:
        raise ValueError('Cannot seal checkpoint before QC passes')
    sources={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(Path(__file__).parent.glob('*.py'))}
    sources['fastq_metrics']=hashlib.sha256(Path(scanner).read_bytes()).hexdigest()
    manifest=dict(schema='snipr-counts-v1',sample_id=sample,counts=numbers,qc=metrics,
        inputs=json.loads(Path(inputs).read_text())['inputs'],parameters=parameters,
        implementation_sha256=sources,cutadapt=json.loads(Path('cutadapt.json').read_text()),pear_log=Path('pear.log').read_text(),
        counting_backend=meta,limitations='No read IDs, pairing or quality strings; not a replacement for raw FASTQ')
    write_json('counts.manifest.json',manifest)
    return manifest

def verify(sample,counts,manifest):
    m=json.loads(Path(manifest).read_text())
    if m.get('schema')!='snipr-counts-v1' or m.get('sample_id')!=sample:raise ValueError('Count checkpoint schema/sample mismatch')
    if m.get('qc',{}).get('accounting_pass') is not True or m['qc'].get('sample_id')!=sample:raise ValueError('Checkpoint QC not passed')
    n=validate_counts(counts)
    if n!=m['counts'] or n['total_reads']!=m['qc']['sequence_count_sum']:raise ValueError('Count checkpoint content/digest mismatch')
    if not m.get('inputs') or not m.get('parameters') or not m.get('implementation_sha256'):raise ValueError('Missing checkpoint provenance')
    return m

def prepare(sheet):
    with open(sheet,newline='',encoding='utf-8-sig') as f:
        reader=csv.DictReader(f)
        if not reader.fieldnames or not {'sample_id','counts','manifest'}<=set(reader.fieldnames):
            raise ValueError('Counts sheet needs sample_id,counts,manifest')
        if len(reader.fieldnames)!=len(set(reader.fieldnames)):
            raise ValueError('Duplicate counts samplesheet column')
        rows=list(reader)
    if not rows:raise ValueError('Empty counts samplesheet')
    seen=set();used=set()
    for r in rows:
        if None in r or any(r.get(k) is None for k in ('sample_id','counts','manifest')):raise ValueError('Counts sheet needs sample_id,counts,manifest')
        sid=r['sample_id']
        if not re.fullmatch('[A-Za-z0-9][A-Za-z0-9_.-]*',sid) or sid in seen or sid=='Name' or sid.endswith('_polyA'):raise ValueError('Invalid/duplicate sample ID')
        seen.add(sid)
        for key in ('counts','manifest'):
            if not r[key]:raise ValueError('Empty checkpoint path')
            p=Path(r[key]).expanduser();p=(Path(sheet).resolve().parent/p).resolve() if not p.is_absolute() else p.resolve()
            if not p.is_file() or p in used:raise ValueError('Missing or reused checkpoint file')
            used.add(p);r[key]=str(p)
    with open('validated_samples.tsv','w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=['sample_id','counts','manifest'],delimiter='\t');w.writeheader();w.writerows({k:r[k] for k in w.fieldnames} for r in rows)

def main():
    p=argparse.ArgumentParser();sub=p.add_subparsers(dest='action',required=True)
    s=sub.add_parser('prepare');s.add_argument('sheet')
    s=sub.add_parser('verify');s.add_argument('sample');s.add_argument('counts');s.add_argument('manifest')
    s=sub.add_parser('seal')
    for key in ('sample','counts','counting','qc','inputs','parameters','scanner'):s.add_argument('--'+key,required=True)
    a=p.parse_args()
    if a.action=='prepare':prepare(a.sheet)
    elif a.action=='verify':verify(a.sample,a.counts,a.manifest)
    else:seal(a.sample,a.counts,a.counting,a.qc,a.inputs,json.loads(a.parameters),a.scanner)
if __name__=='__main__':main()
