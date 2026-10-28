"""PEAR execution and validated compression, independent of QC code."""
import argparse,gzip,hashlib,json,subprocess
from pathlib import Path
from snipr_common import write_json
def compress_pear():
    """Only replace these four task-local outputs after a gzip round-trip checksum."""
    records = []
    for suffix in ('assembled','unassembled.forward','unassembled.reverse','discarded'):
        source = Path(f'assembled.{suffix}.fastq')
        if source.is_symlink() or not source.is_file():
            raise ValueError(f'Missing or nonlocal PEAR output: {source}')
        before = source.stat()
        temp = Path(str(source)+'.gz.partial'); dest = Path(str(source)+'.gz')
        original = hashlib.sha256()
        with source.open('rb') as inp, temp.open('wb') as target:
            with gzip.GzipFile(fileobj=target, mode='wb', compresslevel=1, mtime=0, filename='') as out:
                for block in iter(lambda: inp.read(4*1024*1024), b''):
                    original.update(block); out.write(block)
        restored = hashlib.sha256()
        with gzip.open(temp,'rb') as inp:
            for block in iter(lambda: inp.read(4*1024*1024), b''):
                restored.update(block)
        after = source.stat()
        if original.digest() != restored.digest() or (before.st_size,before.st_mtime_ns) != (after.st_size,after.st_mtime_ns):
            raise ValueError(f'Compression validation failed: {source}')
        temp.replace(dest)
        records.append(dict(file=str(dest), uncompressed_bytes=after.st_size,
                            compressed_bytes=dest.stat().st_size, uncompressed_sha256=original.hexdigest()))
        source.unlink()
    write_json('pear_intermediates.json', records)

def run_pear(args):
    clean_pairs = json.loads(Path(args.cutadapt_json).read_text())['read_counts']['output']
    if clean_pairs == 0:
        # PEAR 0.9.11 aborts on empty input. Preserve a valid all-zero sample.
        for suffix in ('assembled','unassembled.forward','unassembled.reverse','discarded'):
            Path(f'assembled.{suffix}.fastq').write_bytes(b'')
        Path('pear.log').write_text('PEAR not invoked: Cutadapt retained zero pairs.\n'
            'Assembled reads: 0\nDiscarded reads: 0\nNot assembled reads: 0\n')
    else:
        command = [args.binary,'-j',str(args.threads),'-f',args.r1,'-r',args.r2,
                   '-v',str(args.overlap),'-n',str(args.min_length),'-y',args.memory,'-o','assembled']
        with open('pear.log','w') as log:
            subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)
    if getattr(args, "compress", True): compress_pear()


def main():
 p=argparse.ArgumentParser()
 for name in ('binary','r1','r2','cutadapt-json','memory'):p.add_argument('--'+name,required=True)
 for name in ('threads','overlap','min-length'):p.add_argument('--'+name,type=int,required=True)
 a=p.parse_args();run_pear(a)
if __name__=='__main__':main()
