#!/usr/bin/env python3
"""
SNIPR-QC  —  Small RNA-seq per-sample quality control

Generates:
  {sample}.qc_metrics.json   — MultiQC-compatible JSON
  {sample}.qc_report.html    — detailed per-sample HTML report

Usage:
  python snipr_qc.py --fq <fastq.gz> --sample <name> --outdir <dir> \
      --adapter "AGATCGGAAGAGCACA"

Modules used: Bio.SeqIO (biopython), gzip
Install: pip install biopython
"""

import argparse
import json
import os
import sys
import gzip
import re
from collections import defaultdict
from io import StringIO

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
except ImportError:
    sys.stderr.write("ERROR: biopython not installed. Run: pip install biopython\n")
    sys.exit(1)


# ── Argument parsing ─────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(description="SNIPR-QC: Small RNA-seq QC")
    p.add_argument("--fq",        required=True,  help="Input FASTQ (.gz or plain)")
    p.add_argument("--sample",    required=True,  help="Sample name")
    p.add_argument("--outdir",    required=True,  help="Output directory")
    p.add_argument("--adapter",   default="AGATCGGAAGAGCACA",
                                            help="Expected adapter sequence (default: Illumina)")
    p.add_argument("--min_len",   type=int, default=10,
                                            help="Minimum read length to include (default: 10)")
    p.add_argument("--max_len",   type=int, default=50,
                                            help="Maximum read length (default: 50)")
    p.add_argument("--threads",   type=int, default=4,
                                            help="Parallel threads (default: 4)")
    return p.parse_args()


# ── FASTQ iteration ─────────────────────────────────────────────────────────

def iter_fastq(path):
    """Yield (name, seq, qual) tuples from a FASTQ file (gz or plain)."""
    opener = gzip.open if path.endswith(".gz") else open
    mode   = "rt"    if path.endswith(".gz") else "r"
    with opener(path, mode) as fh:
        rec = None
        buf = []
        for i, line in enumerate(fh):
            buf.append(line.rstrip("\n"))
            if len(buf) == 4:
                yield buf[0], buf[1], buf[3]
                buf = []


# ── Core metrics computation ─────────────────────────────────────────────────

def compute_metrics(fq_path, adapter_seq, min_len=10, max_len=50):
    """
    Single-pass FASTQ analysis.
    Returns a dict of computed metrics.
    """

    total_reads      = 0
    qual_by_pos      = defaultdict(list)   # pos -> [phred scores]
    gc_by_pos        = defaultdict(list)   # pos -> [0/1 GC indicator]
    len_counts       = defaultdict(int)    # length -> count
    seqs             = {}                  # seq -> count (for dup detection)
    adapter_hits      = 0
    n_bases_by_pos   = defaultdict(int)    # pos -> N count
    complexity_buckets = {"low": 0, "mid": 0, "high": 0}
    at_ratio_by_pos  = defaultdict(list)   # pos -> [AT fraction]

    adapter_3prime_pattern = re.compile(
        re.escape(adapter_seq) + r"+$", re.IGNORECASE
    )

    qual_offset = 33   # Illumina quality scores

    for name, seq, qual in iter_fastq(fq_path):
        total_reads += 1
        L = len(seq)

        # ── Length distribution ──────────────────────────────────────────────
        if min_len <= L <= max_len:
            len_counts[L] += 1
        elif L < min_len:
            len_counts[0] = len_counts.get(0, 0) + 1  # under-length bucket

        # ── Per-base quality and GC (first 50 bases) ───────────────────────
        clipped_qual = qual[:50] if len(qual) >= 50 else qual
        clipped_seq   = seq[:50]  if len(seq)  >= 50 else seq
        for i, (b, qc) in enumerate(zip(clipped_seq, clipped_qual)):
            phred = ord(qc) - qual_offset
            qual_by_pos[i].append(phred)
            gc_by_pos[i].append(1 if b in "GCgc" else 0)
            n_bases_by_pos[i] += (1 if b in "Nn" else 0)

        # ── Adapter contamination (3' end scan) ─────────────────────────────
        # strip trailing adapter-like runs
        trimmed_seq = adapter_3prime_pattern.sub("", seq)
        if len(trimmed_seq) < len(seq):
            adapter_hits += 1

        # ── Sequence deduplication (for dup rate) ─────────────────────────
        # Only count unique sequences up to a cap to save memory
        if total_reads <= 500_000:
            seqs[seq] = seqs.get(seq, 0) + 1

        # ── Complexity (DUST-like score approximated via dinucleotide entropy)
        if total_reads <= 500_000:
            if L < 4:
                complexity_buckets["low"] += 1
                continue
            dinucs = [seq[i:i+2] for i in range(L - 1)]
            unique_dn = len(set(dinucs)) / max(len(dinucs), 1)
            if unique_dn < 0.5:
                complexity_buckets["low"] += 1
            elif unique_dn < 0.7:
                complexity_buckets["mid"] += 1
            else:
                complexity_buckets["high"] += 1

    # ── Derived statistics ─────────────────────────────────────────────────────

    # Length stats
    total_in_range = sum(v for k, v in len_counts.items() if k > 0)
    median_len = None
    if total_in_range > 0:
        cum = 0
        half = total_in_range / 2
        for l in sorted(len_counts):
            if l <= 0:
                continue
            cum += len_counts[l]
            if cum >= half:
                median_len = l
                break

    # Quality stats per position
    qual_summary = {}
    for pos in sorted(qual_by_pos):
        scores = qual_by_pos[pos]
        qual_summary[pos + 1] = {
            "mean":   round(sum(scores) / len(scores), 2),
            "median": round(sorted(scores)[len(scores) // 2], 2),
            "q25":    round(sorted(scores)[len(scores) // 4], 2),
            "q75":    round(sorted(scores)[3 * len(scores) // 4], 2),
            "count":  len(scores),
        }

    # GC content per position
    gc_summary = {}
    for pos in sorted(gc_by_pos):
        vals = gc_by_pos[pos]
        gc_summary[pos + 1] = round(sum(vals) / len(vals), 3)

    # Mean quality per read
    mean_qual_per_read = []
    for name, seq, qual in iter_fastq(fq_path):
        scores = [ord(c) - qual_offset for c in qual]
        mean_qual_per_read.append(sum(scores) / len(scores) if scores else 0)
    mean_qual_overall = round(sum(mean_qual_per_read) / len(mean_qual_per_read), 2) \
        if mean_qual_per_read else 0

    # Adapter contamination rate
    adapter_rate = round(adapter_hits / max(total_reads, 1) * 100, 2)

    # Duplicate rate
    dup_reads  = sum(1 for v in seqs.values() if v > 1)
    dup_rate   = round(dup_reads / max(len(seqs), 1) * 100, 2) if seqs else 0

    # miRNA-typical length range (18–25 nt) fraction
    mirna_18_25 = sum(len_counts.get(l, 0) for l in range(18, 26))
    mirna_fraction = round(mirna_18_25 / max(total_in_range, 1) * 100, 2)

    # Overrepresented sequences (top 20)
    top_overrep = sorted(seqs.items(), key=lambda x: -x[1])[:20]
    overrep_seqs = [
        {"sequence": seq, "count": cnt, "fraction": round(cnt / max(total_reads, 1) * 100, 4)}
        for seq, cnt in top_overrep
    ]

    metrics = {
        "sample_id":       None,   # filled by caller
        "total_reads":     total_reads,
        "reads_in_range":  total_in_range,
        "median_length":   median_len,
        "mean_quality":    mean_qual_overall,
        "adapter_contamination_pct": adapter_rate,
        "duplicate_rate_pct":        dup_rate,
        "mirna_18_25nt_pct":         mirna_fraction,
        "low_complexity_pct":  round(complexity_buckets["low"] / max(total_reads, 1) * 100, 2),
        "complexity_buckets": complexity_buckets,
        "length_distribution": {str(k): v for k, v in sorted(len_counts.items())},
        "quality_by_position":  qual_summary,
        "gc_by_position":       gc_summary,
        "overrepresented_seqs":  overrep_seqs,
    }

    return metrics


# ── JSON output (MultiQC-compatible) ─────────────────────────────────────────

def write_multiqc_json(metrics, sample_id, outdir):
    """Write {sample}_qc_metrics.json for MultiQC aggregation."""

    data = {
        "sample_id": sample_id,
        "total_reads":         metrics["total_reads"],
        "reads_in_range":      metrics["reads_in_range"],
        "median_length":       metrics["median_length"],
        "mean_quality":        metrics["mean_quality"],
        "adapter_contamination_pct": metrics["adapter_contamination_pct"],
        "duplicate_rate_pct":  metrics["duplicate_rate_pct"],
        "mirna_18_25nt_pct":   metrics["mirna_18_25nt_pct"],
        "low_complexity_pct":  metrics["low_complexity_pct"],
    }

    path = os.path.join(outdir, f"{sample_id}_qc_metrics.json")
    with open(path, "w") as fh:
        json.dump(data, fh, indent=2)

    # Also write a per-position quality table that MultiQC can pick up
    qpath = os.path.join(outdir, f"{sample_id}_quality_by_pos.json")
    with open(qpath, "w") as fh:
        json.dump(metrics["quality_by_position"], fh, indent=2)

    # Length distribution
    lpath = os.path.join(outdir, f"{sample_id}_length_dist.json")
    with open(lpath, "w") as fh:
        json.dump(metrics["length_distribution"], fh, indent=2)

    return path


# ── HTML report ───────────────────────────────────────────────────────────────

def write_html_report(metrics, sample_id, outdir):
    """Generate a standalone HTML QC report."""

    # Length distribution chart data (JSON inline)
    len_data = [
        {"length": l, "count": c}
        for l, c in sorted(metrics["length_distribution"].items(), key=lambda x: int(x[0]))
        if int(l) > 0
    ]
    len_json = json.dumps(len_data)

    # Quality boxplot data
    qdata = []
    for pos, stats in sorted(metrics["quality_by_position"].items(), key=lambda x: int(x[0])):
        qdata.append({"pos": pos, **stats})
    qjson = json.dumps(qdata)

    # GC content data
    gc_data = [
        {"pos": p, "gc": v}
        for p, v in sorted(metrics["gc_by_position"].items(), key=lambda x: int(x[0]))
    ]
    gc_json = json.dumps(gc_data)

    total    = metrics["total_reads"]
    in_range = metrics["reads_in_range"]
    mirna_pct = metrics["mirna_18_25nt_pct"]

    html = f"""\
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>SNIPR-QC Report — {sample_id}</title>
<style>
  body  {{ font-family: -apple-system, Arial, sans-serif; margin: 2rem; background: #f8f9fa; }}
  .card {{ background: white; border-radius: 8px; padding: 1.5rem; margin-bottom: 1.5rem;
           box-shadow: 0 1px 4px rgba(0,0,0,.1); max-width: 900px; margin-inline: auto; }}
  h1, h2 {{ color: #2c3e50; }}
  h1     {{ border-bottom: 3px solid #3498db; padding-bottom: .5rem; }}
  .stats {{ display: grid; grid-template-columns: repeat(auto-fill, minmax(200px, 1fr)); gap: 1rem; }}
  .stat  {{ background: #f0f4f8; border-radius: 6px; padding: 1rem; text-align: center; }}
  .stat .val  {{ font-size: 1.8rem; font-weight: 700; color: #2c3e50; }}
  .stat .lab  {{ font-size: .8rem; color: #7f8c8d; margin-top: .25rem; }}
  .warn  {{ color: #e67e22; font-weight: 600; }}
  .good  {{ color: #27ae60; }}
  .bad   {{ color: #c0392b; }}
  canvas {{ max-width: 100%; }}
  table  {{ width: 100%; border-collapse: collapse; margin-top: 1rem; }}
  th, td {{ padding: .4rem .6rem; border-bottom: 1px solid #ecf0f1; text-align: left; }}
  th     {{ background: #2c3e50; color: white; }}
  tr:hover {{ background: #f5f7fa; }}
  .overrep {{ font-family: monospace; font-size: .85rem; }}
  .footer {{ text-align: center; color: #95a5a6; font-size: .8rem; margin-top: 2rem; }}
</style>
</head>
<body>
<div class="card">
<h1>SNIPR-QC Report — {sample_id}</h1>
</div>

<div class="card">
<h2>Summary Statistics</h2>
<div class="stats">
  <div class="stat">
    <div class="val">{total:,}</div>
    <div class="lab">Total reads</div>
  </div>
  <div class="stat">
    <div class="val">{in_range:,}</div>
    <div class="lab">Reads in length range ({metrics.get('median_length','?')} bp median)</div>
  </div>
  <div class="stat">
    <div class="val class="{'good' if float(mirna_pct) > 30 else 'warn'}">{mirna_pct}%</div>
    <div class="lab">miRNA-typical (18–25 nt)</div>
  </div>
  <div class="stat">
    <div class="val {'warn' if float(metrics['adapter_contamination_pct']) > 5 else 'good'}">{metrics['adapter_contamination_pct']}%</div>
    <div class="lab">Adapter contamination</div>
  </div>
  <div class="stat">
    <div class="val {'bad' if float(metrics['duplicate_rate_pct']) > 50 else 'good'}">{metrics['duplicate_rate_pct']}%</div>
    <div class="lab">Duplicate reads</div>
  </div>
  <div class="stat">
    <div class="val {'bad' if float(metrics['low_complexity_pct']) > 20 else 'good'}">{metrics['low_complexity_pct']}%</div>
    <div class="lab">Low-complexity reads</div>
  </div>
</div>
</div>

<div class="card">
<h2>Read Length Distribution</h2>
<canvas id="lenChart" height="200"></canvas>
<p style="font-size:.85rem;color:#7f8c8d;">
  Shaded area: typical miRNA range (18–25 nt). High abundance in this window suggests
  high-quality small RNA enrichment.
</p>
</div>

<div class="card">
<h2>Per-Base Quality (Phred score)</h2>
<canvas id="qualChart" height="200"></canvas>
</div>

<div class="card">
<h2>Per-Base GC Content</h2>
<canvas id="gcChart" height="150"></canvas>
</div>

<div class="card">
<h2>Overrepresented Sequences (Top 10)</h2>
<table>
<thead><tr><th>Sequence</th><th>Count</th><th>Fraction (%)</th></tr></thead>
<tbody>
{chr(10).join(
    f'<tr><td class="overrep">{s["sequence"]}</td>'
    f'<td>{s["count"]:,}</td>'
    f'<td>{s["fraction"]:.4f}%</td></tr>'
    for s in metrics["overrepresented_seqs"][:10]
)}
</tbody>
</table>
</div>

<div class="footer">
  Generated by SNIPR-QC — Small RNA-seq Quality Control
</div>

<script src="https://cdn.jsdelivr.net/npm/chart.js@4.4.0/dist/chart.umd.min.js"></script>
<script>
  // ── Length distribution ──────────────────────────────────────────────────
  const lenData = {len_json};
  new Chart(document.getElementById('lenChart'), {{
    type: 'bar',
    data: {{
      labels: lenData.map(d => d.length),
      datasets: [{{
        label: 'Read count',
        data: lenData.map(d => d.count),
        backgroundColor: lenData.map(d =>
          (d.length >= 18 && d.length <= 25) ? 'rgba(52,152,219,0.8)' : 'rgba(149,165,166,0.5)'
        ),
        borderWidth: 1,
      }}]
    }},
    options: {{
      responsive: true,
      plugins: {{ legend: {{ display: false }} }},
      scales: {{
        x: {{ title: {{ display: true, text: 'Read length (nt)' }}, ticks: {{ maxTicksLimit: 30 }} }},
        y: {{ title: {{ display: true, text: 'Read count' }}, beginAtZero: true }}
      }}
    }}
  }});

  // ── Per-base quality ──────────────────────────────────────────────────────
  const qData = {qjson};
  new Chart(document.getElementById('qualChart'), {{
    type: 'line',
    data: {{
      labels: qData.map(d => d.pos),
      datasets: [
        {{
          label: 'Mean Phred',
          data: qData.map(d => d.mean),
          borderColor: 'rgba(52,152,219,1)',
          backgroundColor: 'rgba(52,152,219,0.1)',
          fill: true,
          tension: 0.3,
        }},
        {{
          label: 'Q25 / Q75',
          data: qData.map(d => d.q75),
          borderColor: 'rgba(39,174,96,0.5)',
          backgroundColor: 'transparent',
          borderDash: [5,5],
          pointRadius: 0,
        }},
        {{
          label: 'Q25',
          data: qData.map(d => d.q25),
          borderColor: 'rgba(231,76,60,0.5)',
          backgroundColor: 'transparent',
          borderDash: [5,5],
          pointRadius: 0,
        }},
      ]
    }},
    options: {{
      responsive: true,
      plugins: {{ legend: {{ display: true, position: 'bottom' }} }},
      scales: {{
        x: {{ title: {{ display: true, text: 'Position (bp)' }}}},
        y: {{ title: {{ display: true, text: 'Phred score' }}, min: 0, max: 40 }}
      }}
    }}
  }});

  // ── GC content ─────────────────────────────────────────────────────────────
  const gcData = {gc_json};
  new Chart(document.getElementById('gcChart'), {{
    type: 'line',
    data: {{
      labels: gcData.map(d => d.pos),
      datasets: [{{
        label: 'GC fraction',
        data: gcData.map(d => d.gc),
        borderColor: 'rgba(155,89,182,1)',
        backgroundColor: 'rgba(155,89,182,0.1)',
        fill: true,
        tension: 0.3,
      }}]
    }},
    options: {{
      responsive: true,
      plugins: {{ legend: {{ display: false }} }},
      scales: {{
        x: {{ title: {{ display: true, text: 'Position (bp)' }}}},
        y: {{ title: {{ display: true, text: 'GC fraction' }}, min: 0, max: 1 }}
      }}
    }}
  }});
</script>
</body>
</html>"""

    path = os.path.join(outdir, f"{sample_id}_qc_report.html")
    with open(path, "w") as fh:
        fh.write(html)
    return path


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    metrics = compute_metrics(
        args.fq, args.adapter,
        min_len=args.min_len, max_len=args.max_len
    )
    metrics["sample_id"] = args.sample

    json_path = write_multiqc_json(metrics, args.sample, args.outdir)
    html_path = write_html_report(metrics, args.sample, args.outdir)

    print(f"  QC metrics JSON : {json_path}")
    print(f"  QC HTML report  : {html_path}")
    print(f"  Total reads     : {metrics['total_reads']:,}")
    print(f"  miRNA range %   : {metrics['mirna_18_25nt_pct']}%")
    print(f"  Adapter cont. % : {metrics['adapter_contamination_pct']}%")
    print(f"  Low complexity %: {metrics['low_complexity_pct']}%")


if __name__ == "__main__":
    main()
