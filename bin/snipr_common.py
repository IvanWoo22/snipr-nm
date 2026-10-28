"""Shared, dependency-free input validation for SNIPR-Nm command-line tools."""
import csv
import gzip
import hashlib
import html
import json
import re
from pathlib import Path


def open_text(path):
    return gzip.open(path, 'rt', encoding='ascii') if str(path).endswith('.gz') else open(path, encoding='ascii')


def fastq_records(path):
    with open_text(path) as fh:
        record = 0
        while True:
            name = fh.readline()
            if not name:
                return
            record += 1
            seq, plus, qual = (fh.readline().rstrip('\r\n') for _ in range(3))
            if (not name.startswith('@') or not name[1:].strip() or not plus.startswith('+')
                    or not re.fullmatch('[ACGTNacgtn]+', seq) or len(seq) != len(qual)
                    or not re.fullmatch('[!-~]+', qual)):
                raise ValueError(f'{path}: malformed/truncated FASTQ record {record}')
            yield name.rstrip('\r\n'), seq.upper(), qual


def references(path):
    records, names = [], set()
    name, sequence = None, []

    def append():
        if name is None:
            return
        seq = ''.join(sequence).upper().replace('U', 'T')
        if not name or '\t' in name or name in names:
            raise ValueError(f'{path}: empty, duplicate or tab-containing FASTA name: {name!r}')
        if len(seq) < 4 or not re.fullmatch('[ACGTN]+', seq):
            raise ValueError(f'{path}: invalid/short reference: {name}')
        records.append((name, seq))
        names.add(name)

    with open(path, encoding='utf-8') as fh:
        for line in fh:
            line = line.rstrip('\r\n')
            if not line.strip():
                continue
            if line.startswith('>'):
                append()
                name, sequence = line[1:], []
            elif name is None:
                raise ValueError(f'{path}: sequence before FASTA header')
            else:
                sequence.append(''.join(line.split()))
    append()
    if not records:
        raise ValueError(f'{path}: empty reference')
    return records


def samplesheet(path, check_files=False):
    rows, ids, paths = [], set(), set()
    with open(path, newline='', encoding='utf-8-sig') as fh:
        reader = csv.DictReader(fh)
        if not reader.fieldnames or not {'sample_id', 'R1', 'R2'} <= set(reader.fieldnames):
            raise ValueError('Samplesheet needs sample_id,R1,R2 columns')
        if len(reader.fieldnames) != len(set(reader.fieldnames)):
            raise ValueError('Duplicate samplesheet column')
        for row in reader:
            if None in row or any(row.get(k) is None for k in ('sample_id', 'R1', 'R2')):
                raise ValueError('Malformed samplesheet row')
            sid = row['sample_id']
            if not re.fullmatch('[A-Za-z0-9][A-Za-z0-9_.-]*', sid) or sid in ids:
                raise ValueError(f'Invalid/duplicate sample_id: {sid!r}')
            # Prevent collisions between sample and sample_polyA matrix columns.
            if sid.endswith('_polyA') or sid == 'Name':
                raise ValueError(f'Reserved sample_id: {sid}')
            ids.add(sid)
            for key in ('R1', 'R2'):
                if not row[key].strip():
                    raise ValueError(f'{sid}: empty {key}')
                p = Path(row[key]).expanduser()
                if not p.is_absolute():
                    p = Path(path).resolve().parent / p
                p = p.resolve()
                if p in paths:
                    raise ValueError(f'FASTQ reused in samplesheet: {p}')
                paths.add(p)
                if check_files and (not p.is_file() or p.stat().st_size == 0):
                    raise ValueError(f'Missing/empty FASTQ: {p}')
                row[key] = str(p)
            rows.append(row)
    if not rows:
        raise ValueError('Empty samplesheet')
    return rows


def fingerprint(path):
    path = Path(path)
    before = path.stat()
    sha, md5 = hashlib.sha256(), hashlib.md5()
    with path.open('rb') as fh:
        for block in iter(lambda: fh.read(4 * 1024 * 1024), b''):
            sha.update(block)
            md5.update(block)
    after = path.stat()
    if (before.st_size, before.st_mtime_ns, before.st_ctime_ns) != (after.st_size, after.st_mtime_ns, after.st_ctime_ns):
        raise ValueError(f'Input changed during checksum: {path}')
    return dict(path=str(path.resolve()), bytes=after.st_size, mtime_ns=after.st_mtime_ns,
                sha256=sha.hexdigest(), md5=md5.hexdigest())


def write_json(path, data):
    path = Path(path)
    tmp = path.with_name(path.name + '.partial')
    tmp.write_text(json.dumps(data, indent=2, allow_nan=False) + '\n')
    tmp.replace(path)


def write_html(path, title, headers, rows, note=''):
    esc = lambda value: html.escape(str(value))
    table = '<tr>' + ''.join(f'<th>{esc(h)}</th>' for h in headers) + '</tr>'
    table += ''.join('<tr>' + ''.join(f'<td>{esc(v)}</td>' for v in row) + '</tr>' for row in rows)
    Path(path).write_text('<!doctype html><html lang="en"><meta charset="utf-8">'
        f'<title>{esc(title)}</title><style>body{{font:15px system-ui;margin:2rem}}'
        'table{border-collapse:collapse}td,th{padding:.5rem;border:1px solid #ddd;text-align:left}'
        f'</style><h1>{esc(title)}</h1><p>{esc(note)}</p><table>{table}</table></html>')
