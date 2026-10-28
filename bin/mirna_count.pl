#!/usr/bin/env perl
use strict;
use warnings FATAL => 'all';
use IO::Uncompress::Gunzip qw($GunzipError);

# Preserve the historical CLI and headerless output by default.
# Memory scales with the reference, not every prefix of every query sequence.
my $header = @ARGV && $ARGV[0] eq '--header' ? shift @ARGV : '';
die "Usage: $0 [--header] <counts.tsv[.gz]> <reference.fa>\n" unless @ARGV == 2;
my ($counts_file, $mirna_fa) = @ARGV;
my (@refs, %names, %by_length);
open(my $fa, '<', $mirna_fa) or die "$mirna_fa: $!\n";
my ($name, $seq);
sub save_reference {
    return unless defined $name;
    die "Empty or invalid sequence for $name\n" unless $seq =~ /^[ACGTN]+$/;
    # The old prefix algorithm was defined for references >=4 bases.
    die "Reference shorter than 4 bases: $name\n" if length($seq) < 4;
    die "Duplicate reference name: $name\n" if $names{$name}++;
    push @refs, [$name, $seq];
    $by_length{length($seq)}{$seq} //= [0, 0];
}
while (<$fa>) {
    s/\r?\n$//;
    next if /^\s*$/;
    if (/^>(.*)$/) {
        save_reference();
        ($name, $seq) = ($1, '');
        die "Empty name or tab in FASTA header\n" if $name eq '' || $name =~ /\t/;
    } else {
        die "Sequence before FASTA header\n" unless defined $name;
        s/\s+//g;
        tr/acgtunACGTUN/ACGTTNACGTTN/;
        $seq .= $_;
    }
}
save_reference();
close($fa) or die "$mirna_fa: $!\n";
die "Empty reference FASTA\n" unless @refs;
my @lengths = sort {$a <=> $b} keys %by_length;
my $fh;
if ($counts_file =~ /\.gz$/) {
    $fh = IO::Uncompress::Gunzip->new($counts_file, MultiStream => 1, Strict => 1)
        or die "$counts_file: $GunzipError\n";
} else {
    open($fh, '<', $counts_file) or die "$counts_file: $!\n";
}
my $line = 0;
while (<$fh>) {
    ++$line;
    s/\r?\n$//;
    my ($read, $count) = split /\t/, $_, -1;
    die "Invalid count row $line\n" unless /\A[ACGTN]+\t[0-9]+\z/;
    for my $len (@lengths) {
        last if $len > length($read);
        my $entry = $by_length{$len}{substr($read, 0, $len)};
        next unless $entry;
        if (length($read) == $len) {
            $entry->[0] += $count;
        } else {
            my $tail = substr($read, $len);
            my $n = length($tail);
            my $a = ($tail =~ tr/A/A/);
            $entry->[1] += $count if $tail eq 'AA' || $tail eq 'AAA' || ($n > 3 && 4*$a >= 3*$n);
        }
    }
}
die "Count gzip read failed: $GunzipError\n" if $counts_file =~ /\.gz$/ && $fh->error();
close($fh) or die "Cannot close count input: $!\n";
print "Name\texact\tpolyA\n" if $header;
for my $ref (@refs) {
    my ($id, $sequence) = @$ref;
    print join("\t", $id, @{$by_length{length($sequence)}{$sequence}}), "\n";
}
