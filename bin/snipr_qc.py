#!/usr/bin/env python3
"""Streaming QC: full-file length/quality counts; explicitly bounded sequence sampling."""
import argparse
from collections import Counter, defaultdict
import math
from pathlib import Path
from snipr_common import fastq_records, write_json, write_html


def compute_metrics(fq_path, adapter_seq, min_len=10, max_len=50, sample_reads=100000):
    if sample_reads < 1 or min_len < 1 or max_len < min_len:
        raise ValueError('Invalid QC bounds')
    lengths, qualities, sequences = Counter(), Counter(), Counter()
    positions = defaultdict(Counter)
    bases = Counter()
    low = adapter = total = sampled = 0
    for _, seq, qual in fastq_records(fq_path):
        total += 1
        lengths[len(seq)] += 1
        qualities.update(qual)
        bases.update(seq)
        if total > sample_reads:
            continue
        sampled += 1
        sequences[seq] += 1
        adapter += bool(adapter_seq and adapter_seq.upper() in seq)
        for pos, score in enumerate(qual):
            positions[pos][ord(score)-33] += 1
        dinuc = Counter(seq[i:i+2] for i in range(len(seq)-1)
                        if 'N' not in seq[i:i+2])
        n = sum(dinuc.values())
        entropy = -sum((c/n)*math.log2(c/n) for c in dinuc.values()) / 4 if n else 0
        low += entropy < 0.5
    nbases = sum(bases.values())
    pct = lambda n, d: 100*n/d if d else None
    median = None
    if total:
        targets = [(total-1)//2, total//2]
        cumulative, selected = 0, []
        for length, count in sorted(lengths.items()):
            selected.extend(length for target in targets if cumulative <= target < cumulative+count)
            cumulative += count
        median = sum(selected)/2
    metrics = dict(schema_version=2, total_reads=total, total_bases=nbases,
        reads_in_range=sum(c for l,c in lengths.items() if min_len <= l <= max_len),
        range_nt=[min_len,max_len], median_length=median,
        mean_quality=sum((ord(q)-33)*c for q,c in qualities.items())/nbases if nbases else None,
        q30_pct=pct(sum(c for q,c in qualities.items() if ord(q)-33 >= 30), nbases),
        gc_pct=pct(bases['G']+bases['C'], nbases),
        mirna_18_25nt_pct=pct(sum(c for l,c in lengths.items() if 18 <= l <= 25), total),
        sampled_reads=sampled, sampling='first N reads; not a random sample',
        sequence_redundancy_pct=pct(sampled-len(sequences), sampled),
        exact_adapter_motif_pct=pct(adapter, sampled),
        adapter_motif=adapter_seq, low_dinucleotide_entropy_pct=pct(low, sampled),
        low_entropy_definition='Shannon entropy of ACGT dinucleotides / 4 < 0.5; N-containing pairs excluded',
        length_distribution=dict(sorted(lengths.items())),
        quality_by_position={str(p+1): {'reads':sum(h.values()),
            'mean':sum(q*c for q,c in h.items())/sum(h.values()), 'histogram':dict(sorted(h.items()))}
            for p,h in sorted(positions.items())})
    return metrics


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--fq', required=True)
    p.add_argument('--sample', required=True)
    p.add_argument('--outdir', required=True)
    p.add_argument('--adapter', default='AGATCGGAAGAGCACA')
    p.add_argument('--min_len', type=int, default=10)
    p.add_argument('--max_len', type=int, default=50)
    p.add_argument('--sample-reads', type=int, default=100000)
    args = p.parse_args()
    out = Path(args.outdir); out.mkdir(parents=True, exist_ok=True)
    metrics = compute_metrics(args.fq, args.adapter, args.min_len, args.max_len, args.sample_reads)
    metrics['sample_id'] = args.sample
    write_json(out/f'{args.sample}_quality_by_pos.json', metrics.pop('quality_by_position'))
    write_json(out/f'{args.sample}_length_dist.json', metrics.pop('length_distribution'))
    write_json(out/f'{args.sample}_qc_metrics.json', metrics)
    write_html(out/f'{args.sample}_qc_report.html', f'SNIPR QC: {args.sample}', ['Metric','Value'],
               metrics.items(), 'Counts, length and global quality use all reads. Position quality, sequence redundancy, '
               'adapter motif and entropy use the first N reads. Empty denominators are null. '
               'Sequence redundancy is not PCR duplication. Exact adapter motifs do not estimate partial adapter contamination.')


if __name__ == '__main__':
    main()
