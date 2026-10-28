"""HTML presentation only; deliberately excluded from numerical task dependencies."""
import html
from pathlib import Path

def write_html(path, title, headers, rows, note=''):
    esc = lambda value: html.escape(str(value))
    table = '<tr>' + ''.join(f'<th>{esc(h)}</th>' for h in headers) + '</tr>'
    table += ''.join('<tr>' + ''.join(f'<td>{esc(v)}</td>' for v in row) + '</tr>' for row in rows)
    Path(path).write_text('<!doctype html><html lang="en"><meta charset="utf-8">'
        f'<title>{esc(title)}</title><style>body{{font:15px system-ui;margin:2rem}}'
        'table{border-collapse:collapse}td,th{padding:.5rem;border:1px solid #ddd;text-align:left}'
        f'</style><h1>{esc(title)}</h1><p>{esc(note)}</p><table>{table}</table></html>')
