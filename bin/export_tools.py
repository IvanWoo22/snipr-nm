"""Reference and matrix exports, independent of preprocessing and QC."""
import argparse,csv
from pathlib import Path
from snipr_common import references,samplesheet
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
    ids = [row['sample_id'] for row in csv.DictReader(open(sheet, newline=''), delimiter='\t')]
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


def main():
 p=argparse.ArgumentParser();s=p.add_subparsers(dest='action',required=True)
 a=s.add_parser('reference');a.add_argument('fasta')
 a=s.add_parser('annotate');a.add_argument('sample');a.add_argument('counts');a.add_argument('fasta')
 a=s.add_parser('merge');a.add_argument('sheet');a.add_argument('fasta');a.add_argument('files',nargs='+')
 a=p.parse_args()
 if a.action=='reference':
  with open('reference_annotations.tsv','w',newline='') as f:
   w=csv.writer(f,delimiter='\t');w.writerow(['Name','sequence']);w.writerows(references(a.fasta))
 elif a.action=='annotate':annotate(a.sample,a.counts,a.fasta)
 else:merge(a.sheet,a.fasta,a.files)
if __name__=='__main__':main()
