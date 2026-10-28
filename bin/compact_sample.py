"""Own all sample temporary FASTQs; retire only after counts/QC/checkpoint validation."""
import argparse,json
from pathlib import Path
from types import SimpleNamespace
from input_tools import verify,run_cutadapt
from pear_tools import run_pear
from stage_metrics import pair,merged
from qc_assemble import assemble
from checkpoint_tools import seal,verify as verify_checkpoint
from snipr_common import write_json

def run(a):
    params=json.loads(a.parameters)
    verify(a.manifest,[a.r1,a.r2])
    run_cutadapt(a.cutadapt,a.r1,a.r2,params['adapter_r1'],params['adapter_r2'],params['cutadapt_opts'],a.threads)
    verify(a.manifest,[a.r1,a.r2])
    pair(a.scanner,'clean_R1.fq.gz','clean_R2.fq.gz',params['adapter_r1'],params['adapter_r2'],a.sample_reads,'clean')
    run_pear(SimpleNamespace(binary=a.pear,r1='clean_R1.fq.gz',r2='clean_R2.fq.gz',cutadapt_json='cutadapt.json',
        threads=a.threads,overlap=params['pear_min_ovlp'],min_length=params['pear_min_length'],memory=a.pear_memory,compress=False))
    merged(a.scanner,'assembled.assembled.fastq',params['adapter_r1'],a.sample_reads,a.max_unique)
    assemble(a.sample,a.raw_metrics,'clean.metrics.json','merged.metrics.json','cutadapt.json','pear.log','counting.json')
    seal(a.sample,'reads.count.tsv.gz','counting.json',a.sample+'_qc_metrics.json',a.manifest,params,a.scanner)
    verify_checkpoint(a.sample,'reads.count.tsv.gz','counts.manifest.json')
    # Restrict deletion to regular files produced inside THIS task, never staged inputs.
    names=['clean_R1.fq.gz','clean_R2.fq.gz']+[f'assembled.{s}.fastq' for s in ('assembled','unassembled.forward','unassembled.reverse','discarded')]
    inventory=[]
    for name in names:
        p=Path(name)
        if p.is_symlink() or not p.is_file():raise ValueError('Refusing to retire missing/nonlocal temporary file: '+name)
        inventory.append(dict(path=name,bytes=p.stat().st_size,action='pending'))
    write_json('intermediates.json',inventory)
    for record in inventory:
        Path(record['path']).unlink();record['action']='removed_after_checkpoint_verification'
    write_json('intermediates.json',inventory)

def main():
    p=argparse.ArgumentParser()
    for key in ('sample','r1','r2','manifest','raw-metrics','parameters','scanner','cutadapt','pear','pear-memory'):p.add_argument('--'+key,required=True)
    for key in ('threads','sample-reads','max-unique'):p.add_argument('--'+key,type=int,required=True)
    a=p.parse_args();run(a)
if __name__=='__main__':main()
