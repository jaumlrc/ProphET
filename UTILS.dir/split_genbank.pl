#!/bin/env perl
 
use strict;

my $usage = "\n\nsplit_genbank.pl <gb file> [<suffix (default=.gbk)>]\n\n";
die $usage if scalar( @ARGV ) != 1 && scalar( @ARGV ) != 2;

open ARQ, $ARGV[0];

my $suffix = ".gbk";
$suffix = $ARGV[1] if $#ARGV == 1;

my @todoArq = <ARQ>;
my $todoString = join '', @todoArq;

my @record = split 'LOCUS', $todoString;
shift @record;

foreach my $currRecord ( @record ){
  my ($fileName) = ( $currRecord =~ /\s+(\w+)/ );
  open OUT_ARQ, ">" . $fileName . $suffix;
  print OUT_ARQ "LOCUS" . $currRecord;
  close OUT_ARQ;
}

close ARQ;
