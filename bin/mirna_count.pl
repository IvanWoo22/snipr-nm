#!/usr/bin/env perl
use strict;
use warnings FATAL => 'all';

# Optimized miRNA counting with poly-A tail detection
#
# Algorithm:
#   1. Load all query sequences into a hash keyed by the full sequence
#   2. Group all reads by their reference prefix (the part matching the known miRNA)
#      Only keep groups that have at least one read sharing a prefix with some reference
#   3. For each reference in ath_mirna_neo.fa:
#        - Exact match: read sequence == reference sequence
#        - Poly-A tail:  read starts with reference + tail where:
#            tail length > 3  AND  A content >= 75%
#            OR tail is exactly AA or AAA
#   Complexity: O(n + m) — one pass to build the hash, one pass per reference
#               vs original O(n × m) where m = number of references
#
# Usage: perl mirna_count.pl <counts.tsv> <ath_mirna_neo.fa>
#   counts.tsv    — sequence\tcount lines from fastq2count.pl
#   ath_mirna_neo.fa — miRNA reference sequences (one line per sequence, name line starts with >)
#
# Output: tab-separated "mirna_name\texact_match_count\tpolyA_variant_count"

my $polyA_min_len = 3;
my $polyA_min_A_frac = 0.75;

die "Usage: $0 <counts.tsv> <mirna_fa>\n" unless @ARGV == 2;

my ($counts_file, $mirna_fa) = @ARGV;  # order: counts first, FASTA second (matching original)

# ── Step 1: Load all query sequences into %read_count ──────────────────────
my %read_count;
open( my $fh, '<', $counts_file );
while (<$fh>) {
    chomp;
    my ( $seq, $cnt ) = split /\t/, $_, 2;
    $read_count{$seq} //= 0;
    $read_count{$seq} += $cnt;
}
close($fh);

# ── Step 2: Pre-build a prefix-index of reads ─────────────────────────────
# For each read, compute all prefixes from length $polyA_min_len+1 up to full length
# and insert into %prefix_map{prefix}[] so we can look up all reads sharing
# the same prefix in O(1) instead of scanning all reads for each reference.
my %prefix_map;

for my $seq ( keys %read_count ) {
    my $len = length($seq);
    # Only index prefixes long enough to potentially match a reference
    # (minimum reference length is $polyA_min_len + 1 = 4)
    for my $pre_len ( $polyA_min_len + 1 .. $len ) {
        # Force a plain-string copy so the key is always a clean scalar
        (my $pre = substr( $seq, 0, $pre_len )) =~ /(.*)/;
        push @{ $prefix_map{$pre} }, $seq;
    }
}

# ── Step 3: Process each reference miRNA ──────────────────────────────────
open( my $fa_fh, '<', $mirna_fa );
while ( my $name = <$fa_fh> ) {
    chomp($name);
    $name =~ s/^>//;    # strip leading >
    chomp( my $ref = <$fa_fh> );
    my $ref_len = length($ref);

    my $exact    = 0;
    my $polyA    = 0;
    my %seen_seq;       # avoid double-counting if same seq appears via multiple prefixes

    # Collect candidate reads: any read whose prefix matches this reference
    # We look at all prefixes of the reference to find indexed reads.
    my %candidates;
    for my $pre_len ( $polyA_min_len + 1 .. $ref_len ) {
        my $pre = substr( $ref, 0, $pre_len );
        if ( exists $prefix_map{$pre} ) {
            for my $read_seq ( @{ $prefix_map{$pre} } ) {
                $candidates{$read_seq} = 1;
            }
        }
    }

    # Also always include the exact match
    $candidates{$ref} = 1 if exists $read_count{$ref};

    for my $seq ( keys %candidates ) {
        my $cnt = $read_count{$seq} // next;

        if ( $seq eq $ref ) {
            $exact += $cnt;
            $seen_seq{$seq} = 1;
        }
        elsif ( $seq =~ /^\Q$ref\E(.*)/s ) {
            my $tail = $1;
            next if length($tail) == 0;    # same as exact — already counted
            next if length($tail) <= $polyA_min_len
                       && $tail ne 'AA'
                       && $tail ne 'AAA';
            my $a_frac = eval {
                my $a_count =()= $tail =~ /A/g;
                $a_count / length($tail);
            };
            next unless defined $a_frac && $a_frac >= $polyA_min_A_frac;

            $polyA += $cnt unless $seen_seq{$seq}++;
        }
    }

    print "$name\t$exact\t$polyA\n";
}
close($fa_fh);
