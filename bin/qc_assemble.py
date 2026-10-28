"""Numerical QC assembly and conservation checks; consumes small files only."""
import argparse,json,re
from pathlib import Path
from snipr_common import write_json

def assemble(sample,raw,clean,merged,cutadapt,pear,counting):
    stages={**json.loads(Path(raw).read_text()),**json.loads(Path(clean).read_text()),'merged':json.loads(Path(merged).read_text())}
    if set(stages)!={'raw_R1','raw_R2','clean_R1','clean_R2','merged'}:raise ValueError('Missing QC stages')
    lengths={};qualities={}
    for label,metric in stages.items():
        lengths[label]=metric.pop('length_distribution');qualities[label]=metric.pop('quality_by_position')
        if sum(lengths[label].values())!=metric['total_reads'] or sum(int(k)*v for k,v in lengths[label].items())!=metric['total_bases']:
            raise ValueError('Invalid length/read/base totals')
    n={k:v['total_reads'] for k,v in stages.items()};cut=json.loads(Path(cutadapt).read_text())['read_counts']
    if n['raw_R1']!=n['raw_R2'] or n['raw_R1']!=cut['input']:raise ValueError('Raw paired read count / Cutadapt input mismatch')
    if n['clean_R1']!=n['clean_R2'] or n['clean_R1']!=cut['output']:raise ValueError('Clean paired read count / Cutadapt output mismatch')
    counts=json.loads(Path(counting).read_text())
    if counts['total_reads']!=n['merged'] or n['merged']>n['clean_R1']:raise ValueError('Merged read / sequence count conservation failed')
    pear_counts={};log=Path(pear).read_text()
    for key in ('Assembled','Discarded','Not assembled'):
        m=re.search(r'^'+key+r' reads[ .]*:\s*([\d,]+)',log,re.M)
        if not m:raise ValueError('Missing PEAR accounting: '+key)
        pear_counts[key]=int(m[1].replace(',',''))
    if pear_counts['Assembled']!=n['merged'] or sum(pear_counts.values())!=n['clean_R1']:raise ValueError('PEAR pair accounting failed')
    result=dict(schema_version=2,sample_id=sample,stages=stages,pear_pairs=pear_counts,sequence_count_sum=counts['total_reads'],accounting_pass=True)
    write_json(sample+'_qc_metrics.json',result);write_json(sample+'_length_dist.json',lengths);write_json(sample+'_quality_by_pos.json',qualities)
    return result

def main():
    p=argparse.ArgumentParser()
    for name in ('sample','raw','clean','merged','cutadapt','pear','counting'):p.add_argument('--'+name,required=True)
    a=p.parse_args();assemble(a.sample,a.raw,a.clean,a.merged,a.cutadapt,a.pear,a.counting)
if __name__=='__main__':main()
