#!/usr/bin/env perl

# Select FASTQ records by QNAME.

# FASTQ records provided to STDIN are printed to STDOUT if and only if
# they occur in the given list of QNAMEs.

# Command line args:
#   [1] QNAME file with a QNAME to include on each line.

use warnings;
use strict;

################################################################
# Collect args

my $nArgs = 1;
unless (scalar(@ARGV) == $nArgs) {
    die "USAGE: $0 qname.txt < input.fastq > output.fastq";
}

my($qnameFile) = @ARGV;

################################################################
# Read list of QNAMEs into a hash

open(QNAME, $qnameFile) or die("ERROR: Can't open QNAME file '$qnameFile': $!");
my %qnames = ();
while (my $line = <QNAME>) {
    chomp($line);
    $qnames{$line}++;
}
close QNAME;

################################################################
# Read records from STDIN, printing those with QNAMEs on the list to STDOUT

while (my $line1 = <STDIN>) {
    my $record = $line1;
    $record .= <STDIN>;
    $record .= <STDIN>;
    $record .= <STDIN>;
    $line1 =~ /^@(\S+)\s/ or die("ERROR: Unexpected line '$line1'");
    my $qname = $1;
    $qname =~ s/\/[12]$//;
    if (exists $qnames{$qname}) {
        print $record;
    }
}

################################################################

exit 0;
