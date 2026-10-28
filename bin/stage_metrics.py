"""Native stage scanner orchestration; no report rendering or reference dependencies."""
import argparse,gzip,hashlib,json,shutil,subprocess
from pathlib import Path
from snipr_common import write_json

def scan(scanner, fastq, adapter, sample_reads, output, counts=False, max_unique=1000000):
    command=[str(Path(scanner).resolve()),'--input',str(fastq),'--metrics',str(output),
        '--adapter',adapter,'--sample-reads',str(sample_reads)]
    if counts:command+=['--counts','reads.count.tsv','--max-unique',str(max_unique),'--count-meta','counting.json']
    subprocess.run(command,check=True)
    return json.loads(Path(output).read_text())

def pair(scanner,r1,r2,adapter1,adapter2,sample_reads,prefix):
    metrics={}
    for mate,path,adapter in [('R1',r1,adapter1),('R2',r2,adapter2)]:
        tmp=Path(f'.{prefix}_{mate}.metrics.json')
        metrics[prefix+'_'+mate]=scan(scanner,path,adapter,sample_reads,tmp);tmp.unlink()
    if metrics[prefix+'_R1']['total_reads']!=metrics[prefix+'_R2']['total_reads']:
        raise ValueError('Paired read count mismatch')
    write_json(prefix+'.metrics.json',metrics)
    return metrics

def merged(scanner,fastq,adapter,sample_reads,max_unique):
    metric=scan(scanner,fastq,adapter,sample_reads,'merged.metrics.json',True,max_unique)
    source=Path('reads.count.tsv');temp=Path('reads.count.tsv.gz.partial')
    with source.open('rb') as src,temp.open('wb') as target:
        with gzip.GzipFile(fileobj=target,mode='wb',compresslevel=1,mtime=0,filename='') as out:shutil.copyfileobj(src,out,1024*1024)
    # Verify the small persistent checkpoint before removing its plain temporary copy.
    with source.open('rb') as f:expected=hashlib.file_digest(f,'sha256').hexdigest()
    with gzip.open(temp,'rb') as f:observed=hashlib.file_digest(f,'sha256').hexdigest()
    if expected!=observed:raise ValueError('Count compression verification failed')
    temp.replace('reads.count.tsv.gz');source.unlink()
    meta=json.loads(Path('counting.json').read_text());meta['content_sha256']=observed
    meta['compressed_sha256']=hashlib.sha256(Path('reads.count.tsv.gz').read_bytes()).hexdigest()
    if meta['total_reads']!=metric['total_reads']:raise ValueError('Merged count/metric mismatch')
    write_json('counting.json',meta)
    return metric

def main():
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['pair','merged']);p.add_argument('--scanner',required=True)
    p.add_argument('--r1');p.add_argument('--r2');p.add_argument('--fastq');p.add_argument('--prefix',choices=['raw','clean'])
    p.add_argument('--adapter-r1',required=True);p.add_argument('--adapter-r2',default='')
    p.add_argument('--sample-reads',type=int,required=True);p.add_argument('--max-unique',type=int,default=1000000)
    a=p.parse_args()
    if a.mode=='pair':pair(a.scanner,a.r1,a.r2,a.adapter_r1,a.adapter_r2,a.sample_reads,a.prefix)
    else:merged(a.scanner,a.fastq,a.adapter_r1,a.sample_reads,a.max_unique)
if __name__=='__main__':main()
