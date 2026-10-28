#!/usr/bin/env perl
use strict;
use warnings FATAL => 'all';
use IO::Uncompress::Gunzip qw($GunzipError);

# One merged FASTQ, or stdin. Output sequence/count TSV, sorted by sequence.
die "Usage: $0 [merged.fastq[.gz]]\n" if @ARGV > 1;
my $path = shift @ARGV;
my $fh;
if (defined $path && $path =~ /\.gz$/) {
    $fh = IO::Uncompress::Gunzip->new($path, MultiStream => 1, Strict => 1)
        or die "$path: $GunzipError\n";
} elsif (defined $path) {
    open($fh, '<', $path) or die "$path: $!\n";
} else { $fh = *STDIN; }
my (%count, $record);
while (my $name = <$fh>) {
    ++$record;
    my ($seq, $plus, $qual) = (scalar <$fh>, scalar <$fh>, scalar <$fh>);
    die "Truncated FASTQ at record $record\n" unless defined $qual;
    s/\r?\n$// for ($name, $seq, $plus, $qual);
    die "Invalid FASTQ record $record\n" unless $name =~ /^@./ && $plus =~ /^\+/
        && $seq =~ /^[ACGTN]+$/ && length($qual) == length($seq) && $qual =~ /^[!-~]+$/;
    $count{$seq}++;
}
die "FASTQ gzip read failed: $GunzipError\n" if defined $path && $path =~ /\.gz$/ && $fh->error();
close($fh) or die "Cannot close FASTQ: $!\n";
for my $seq (sort keys %count) { print "$seq\t$count{$seq}\n"; }
