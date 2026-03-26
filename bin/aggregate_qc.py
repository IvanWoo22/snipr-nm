#!/usr/bin/env python3
"""
aggregate_qc.py  —  Aggregate per-sample snipr_qc JSON files into a MultiQC-compatible report.

Usage:
  python3 aggregate_qc.py \
      --samplesheet samplesheet.csv \
      --qc_dir /path/to/results/{sample}/ \
      --outdir /path/to/results/

Outputs:
  {outdir}/multiqc_report.html
  {outdir}/multiqc_data/snipr_qc_all.txt
  {outdir}/multiqc_data/snipr_qual_{sample}.json   (copied)
  {outdir}/multiqc_data/snipr_len_{sample}.json    (copied)
"""

import argparse
import json
import os
import glob
import shutil
import sys


def parse_args():
    p = argparse.ArgumentParser(description="Aggregate SNIPR-QC JSON files into MultiQC report")
    p.add_argument("--samplesheet",  required=True, help="CSV: sample_id,R1,R2")
    p.add_argument("--qc_dir",      required=True, help="Results directory containing per-sample subdirs")
    p.add_argument("--outdir",       required=True, help="Output directory")
    return p.parse_args()


def load_samplesheet(path):
    """Return ordered list of sample IDs from samplesheet."""
    samples = []
    with open(path) as fh:
        header = fh.readline()
        for line in fh:
            line = line.strip()
            if not line:
                continue
            sample_id = line.split(",")[0]
            samples.append(sample_id)
    return samples


def build_general_stats_table(sample_ids, qc_dir):
    """Build the MultiQC general-stats tab-delimited table."""
    rows = []
    headers = [
        "Sample",
        "Total reads",
        "Reads in range",
        "Median length (bp)",
        "Mean quality",
        "Adapter contamination (%)",
        "Duplicate rate (%)",
        "miRNA 18-25nt (%)",
        "Low complexity (%)",
    ]
    rows.append("\t".join(headers))

    for sid in sample_ids:
        json_path = os.path.join(qc_dir, sid, f"{sid}_qc_metrics.json")
        if not os.path.exists(json_path):
            # Try older naming
            json_path = os.path.join(qc_dir, sid, "qc_metrics.json")

        if os.path.exists(json_path):
            d = json.load(open(json_path))
        else:
            d = {}

        row = [
            sid,
            str(d.get("total_reads", 0)),
            str(d.get("reads_in_range", 0)),
            str(d.get("median_length", "NA")),
            str(d.get("mean_quality", 0)),
            str(d.get("adapter_contamination_pct", 0)),
            str(d.get("duplicate_rate_pct", 0)),
            str(d.get("mirna_18_25nt_pct", 0)),
            str(d.get("low_complexity_pct", 0)),
        ]
        rows.append("\t".join(row))

    return "\n".join(rows) + "\n"


def copy_custom_content(sample_ids, qc_dir, out_data_dir):
    """Copy per-position quality and length distribution JSONs for MultiQC custom content."""
    for sid in sample_ids:
        sample_qc_dir = os.path.join(qc_dir, sid)

        for prefix in ("snipr_qual_", "snipr_len_"):
            src = os.path.join(sample_qc_dir, f"{sid}_quality_by_pos.json")
            if not os.path.exists(src):
                src = os.path.join(sample_qc_dir, "quality_by_pos.json")
            dst = os.path.join(out_data_dir, f"{prefix}{sid}.json")
            if os.path.exists(src):
                shutil.copy2(src, dst)

        # Copy per-sample HTML report
        html_src = os.path.join(sample_qc_dir, f"{sid}_qc_report.html")
        if os.path.exists(html_src):
            shutil.copy2(html_src, out_data_dir)


def build_combined_html(sample_ids, qc_dir, outdir):
    """Build a combined HTML report with per-sample metrics table."""

    # Read per-sample JSONs for table
    rows = []
    for sid in sample_ids:
        json_path = os.path.join(qc_dir, sid, f"{sid}_qc_metrics.json")
        if not os.path.exists(json_path):
            json_path = os.path.join(qc_dir, sid, "qc_metrics.json")
        if os.path.exists(json_path):
            d = json.load(open(json_path))
            mirna_pct = float(d.get("mirna_18_25nt_pct", 0))
            color = "#27ae60" if mirna_pct > 30 else "#e67e22"
            rows.append(f"""
        <tr>
          <td>{sid}</td>
          <td>{d.get('total_reads', 0):,}</td>
          <td>{d.get('reads_in_range', 0):,}</td>
          <td>{d.get('median_length', 'NA')}</td>
          <td>{d.get('mean_quality', 0):.1f}</td>
          <td>{d.get('adapter_contamination_pct', 0):.2f}%</td>
          <td>{d.get('duplicate_rate_pct', 0):.2f}%</td>
          <td style="color:{color};font-weight:600">{d.get('mirna_18_25nt_pct', 0):.2f}%</td>
          <td>{d.get('low_complexity_pct', 0):.2f}%</td>
          <td><a href="multiqc_data/{sid}_qc_report.html">report</a></td>
        </tr>""")

    rows_html = "\n".join(rows)

    html = f"""\
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>SNIPR-Nm QC — Aggregated Report</title>
<style>
  body  {{ font-family: -apple-system, Arial, sans-serif; margin: 2rem; background: #f8f9fa; }}
  .card {{ background: white; border-radius: 8px; padding: 1.5rem; margin-bottom: 1.5rem;
           box-shadow: 0 1px 4px rgba(0,0,0,.1); max-width: 1100px; margin-inline: auto; }}
  h1, h2 {{ color: #2c3e50; }}
  h1     {{ border-bottom: 3px solid #3498db; padding-bottom: .5rem; }}
  table  {{ width: 100%; border-collapse: collapse; margin-top: 1rem; }}
  th, td {{ padding: .5rem .6rem; border-bottom: 1px solid #ecf0f1; text-align: left; }}
  th     {{ background: #2c3e50; color: white; position: sticky; top: 0; }}
  tr:hover {{ background: #f5f7fa; }}
  .footer {{ text-align: center; color: #95a5a6; font-size: .8rem; margin-top: 2rem; }}
  a      {{ color: #3498db; }}
</style>
</head>
<body>
<div class="card">
  <h1>SNIPR-Nm QC — Aggregated Report</h1>
</div>

<div class="card">
  <h2>Per-Sample Summary</h2>
  <table>
    <thead>
      <tr>
        <th>Sample</th>
        <th>Total reads</th>
        <th>Reads in range</th>
        <th>Median len (bp)</th>
        <th>Mean quality</th>
        <th>Adapter (%)</th>
        <th>Duplicate (%)</th>
        <th>miRNA 18–25nt (%)</th>
        <th>Low complexity (%)</th>
        <th>Detail</th>
      </tr>
    </thead>
    <tbody>
{rows_html}
    </tbody>
  </table>
</div>

<div class="card">
  <h2>Key thresholds</h2>
  <ul>
    <li><strong>miRNA 18–25nt:</strong> &gt;30% indicates good small RNA enrichment</li>
    <li><strong>Adapter contamination:</strong> &lt;5% is desirable</li>
    <li><strong>Low complexity:</strong> &lt;20% is acceptable</li>
    <li><strong>Mean quality:</strong> Q30+ (Phred &ge; 30) is good</li>
  </ul>
</div>

<div class="footer">
  Generated by SNIPR-Nm aggregate_qc.py
</div>
</body>
</html>"""

    path = os.path.join(outdir, "multiqc_report.html")
    with open(path, "w") as fh:
        fh.write(html)
    return path


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    sample_ids = load_samplesheet(args.samplesheet)
    out_data = os.path.join(args.outdir, "multiqc_data")
    os.makedirs(out_data, exist_ok=True)

    # ── General stats table ────────────────────────────────────────────────
    table = build_general_stats_table(sample_ids, args.qc_dir)
    table_path = os.path.join(out_data, "snipr_qc_all.txt")
    with open(table_path, "w") as fh:
        fh.write(table)
    print(f"  General stats  : {table_path}")

    # ── Custom content JSONs ────────────────────────────────────────────────
    copy_custom_content(sample_ids, args.qc_dir, out_data)
    print(f"  Custom content copied to {out_data}")

    # ── Combined HTML report ────────────────────────────────────────────────
    report_path = build_combined_html(sample_ids, args.qc_dir, args.outdir)
    print(f"  Combined HTML  : {report_path}")
    print(f"\nAggregated {len(sample_ids)} samples.")


if __name__ == "__main__":
    main()
