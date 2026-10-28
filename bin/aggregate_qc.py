#!/usr/bin/env python3
"""Aggregate staged SNIPR QC outputs. These are custom reports, not MultiQC."""
import argparse
import csv
import json
from pathlib import Path
import shutil
from report_html import write_html


def aggregate(sheet, files, outdir):
    with open(sheet, newline='') as fh:
        ids = [row['sample_id'] for row in csv.DictReader(fh, delimiter='\t' if str(sheet).endswith('.tsv') else ',')]
    if not ids or len(ids) != len(set(ids)): raise ValueError('Empty/duplicate sample order')
    indexed = {Path(f).name:Path(f) for f in files}
    suffixes = ('qc_metrics.json','quality_by_pos.json','length_dist.json','qc_report.html')
    expected = {f'{sid}_{suffix}' for sid in ids for suffix in suffixes}
    if len(indexed) != len(files) or set(indexed) != expected:
        raise ValueError(f'Missing, duplicate or unexpected QC files: missing={sorted(expected-set(indexed))}, extra={sorted(set(indexed)-expected)}')
    rows = []
    for sid in ids:
        metric = json.loads(indexed[f'{sid}_qc_metrics.json'].read_text())
        if metric.get('sample_id') != sid or metric.get('schema_version') != 2 or not metric.get('accounting_pass'):
            raise ValueError(f'Invalid QC schema/sample/accounting: {sid}')
        stages = metric['stages']
        lengths = json.loads(indexed[f'{sid}_length_dist.json'].read_text())
        qualities = json.loads(indexed[f'{sid}_quality_by_pos.json'].read_text())
        if set(stages) != {'raw_R1','raw_R2','clean_R1','clean_R2','merged'} or set(lengths) != set(stages) or set(qualities) != set(stages):
            raise ValueError(f'Missing QC stages: {sid}')
        for stage, stats in stages.items():
            if sum(lengths[stage].values()) != stats['total_reads']:
                raise ValueError(f'Length histogram total mismatch: {sid}/{stage}')
        raw,clean,merged = [stages[k]['total_reads'] for k in ('raw_R1','clean_R1','merged')]
        rows.append([sid,raw,clean,merged,100*clean/raw if raw else 'NA',
                     100*merged/clean if clean else 'NA',stages['merged']['q30_pct']])
    out = Path(outdir); data = out/'snipr_qc_data'; data.mkdir(parents=True,exist_ok=True)
    headers = ['Sample','Raw pairs','Clean pairs','Merged pairs','Trim retention %','Merge retention %','Merged Q30 %']
    with (data/'sample_summary.tsv').open('w',newline='') as fh:
        writer = csv.writer(fh,delimiter='\t');writer.writerow(headers);writer.writerows(rows)
    for filename,path in indexed.items():
        shutil.copyfile(path,data/filename)
    write_html(out/'snipr_qc_report.html','SNIPR-Nm batch QC',headers,rows,
               'All expected samples and read conservation checks passed. Detailed per-sample reports and JSONs are in snipr_qc_data/. '
               'This is the built-in SNIPR report, not a MultiQC report.')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--samplesheet',required=True)
    p.add_argument('--files',nargs='+',required=True,help='Staged metrics, length, position-quality and HTML files')
    p.add_argument('--outdir',required=True)
    a=p.parse_args();aggregate(a.samplesheet,a.files,a.outdir)


if __name__=='__main__':main()
