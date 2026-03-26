#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

# Read paired FASTQ files (R1 and R2) where PEAR has merged them
# into a single assembled file. Merged reads have both reads overlapping
# in a single record, so we only need to parse the assembled file.
#
# Input:  interleaved or single FASTQ from PEAR merged reads
# Output: tab-separated "sequence\tcount" (one line per unique seq)

my %count;

while ( my $name  = <> ) {
    chomp( my $seq  = <> );
    <>;                      # quality header line (+)
    chomp( my $qual = <> ); # quality string — discarded

    $count{$seq}++;
}

for my $seq ( sort keys %count ) {
    print "$seq\t$count{$seq}\n";
}
