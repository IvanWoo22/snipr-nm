"""Render precomputed sample QC. Style changes do not invalidate numerical tasks."""
import argparse,json
from pathlib import Path
from report_html import write_html

def render(sample,metrics):
    result=json.loads(Path(metrics).read_text())
    if result['sample_id']!=sample or not result['accounting_pass']:raise ValueError('Invalid QC identity/accounting')
    columns=('total_reads','total_bases','mean_quality','q30_pct','mirna_18_25nt_pct','sampled_reads','sequence_redundancy_pct')
    rows=[[label,*[m[k] for k in columns]] for label,m in result['stages'].items()]
    write_html(sample+'_qc_report.html',f'SNIPR QC: {sample}',
        ['Stage','Reads','Bases','Mean quality','Q30 %','18-25 nt %','Sampled reads','Sample redundancy %'],rows,
        'Full-file counts and global quality; detailed position quality, entropy, motif and redundancy use the first N reads. '
        'Raw/clean counts are reads per mate; merged counts are merged pairs. Sequence redundancy is not PCR duplication.')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('sample');p.add_argument('metrics');a=p.parse_args();render(a.sample,a.metrics)
