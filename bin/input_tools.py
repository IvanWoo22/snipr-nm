"""Input identity checks; independent of QC, references and rendering."""
import argparse,csv,json,shlex,subprocess
from pathlib import Path
from snipr_common import samplesheet,fingerprint,write_json
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


def run_cutadapt(binary, r1, r2, adapter_r1, adapter_r2, options, threads):
    """Pass options as argv in both workflow modes; never evaluate shell syntax."""
    if threads < 1:
        raise ValueError('Cutadapt threads must be positive')
    command = [binary, '-a', adapter_r1, '-A', adapter_r2, *shlex.split(options),
               '-j', str(threads), '--json', 'cutadapt.json',
               '-o', 'clean_R1.fq.gz', '-p', 'clean_R2.fq.gz', str(r1), str(r2)]
    with open('cutadapt.log', 'w') as log:
        subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)


def main():
 p=argparse.ArgumentParser();s=p.add_subparsers(dest='action',required=True)
 a=s.add_parser('prepare');a.add_argument('sheet')
 a=s.add_parser('audit');a.add_argument('manifest');a.add_argument('paths',nargs=2);a.add_argument('--md5',nargs=2,default=['-','-'])
 a=s.add_parser('verify');a.add_argument('manifest');a.add_argument('paths',nargs=2)
 a=s.add_parser('trim')
 for key in ('binary','r1','r2','adapter-r1','adapter-r2','options'):a.add_argument('--'+key,required=True)
 a.add_argument('--threads',type=int,required=True)
 a=p.parse_args()
 if a.action=='prepare':
  rows=samplesheet(a.sheet,check_files=True)
  for row in rows:
   for key in ('R1_md5','R2_md5'):
    v=row.get(key,'') or ''
    if v and (len(v)!=32 or any(c not in '0123456789abcdefABCDEF' for c in v)):raise ValueError('Invalid MD5')
    row[key]=v.lower()
  with open('validated_samples.tsv','w',newline='') as f:
   w=csv.DictWriter(f,fieldnames=['sample_id','R1','R2','R1_md5','R2_md5'],delimiter='\t');w.writeheader();w.writerows({k:r.get(k,'') for k in w.fieldnames} for r in rows)
 elif a.action=='audit':audit(a.paths,a.md5,a.manifest)
 elif a.action=='verify':verify(a.manifest,a.paths)
 else:run_cutadapt(a.binary,a.r1,a.r2,a.adapter_r1,a.adapter_r2,a.options,a.threads)
if __name__=='__main__':main()
