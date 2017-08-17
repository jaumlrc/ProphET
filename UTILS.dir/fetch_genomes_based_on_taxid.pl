#!/bin/env perl

use strict;
use LWP::Simple;
use XML::Simple;
#use Data::Dumper;

my ( $name, $outname, $url, $xml, $out, $count, $query_key, $webenv, $ids );
my @genomeId;
my @genomeId_nuccore;
my $base  = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
#my $limit = 'srcdb+refseq[prop]+AND+gene+in+chromosome[prop])';
my $limit = 'srcdb+refseq[prop]';

my $s = $ARGV[0];
my $debug = 0;

undef @genomeId;
$query_key = $webenv = '';
$s =~ s/ /+/g;

# ESearch
$url = $base . "esearch.fcgi?db=genome&term=txid" . $s . "[Organism:exp]&usehistory=y";
print "URL esearch 1: " . $url . "\n" if $debug;
$xml = get($url);

if ( $xml =~ /<Count>(\d+)<\/Count>/ ){
	$count = $1;
}

print "Number of genomes: $count\n" if $debug;

if ( $count > 20 ) {
	$url =
	    $base
	  . "esearch.fcgi?db=genome&term==txid"
	  . $s
	  . "[Organism:exp]&retmax=$count&usehistory=y";
	print "URL esearch 2: " . $url . "\n" if $debug;
	$xml = get($url);
}
while ( $xml =~ /<Id>(\d+?)<\/Id>/gs ) {
	my $curr_genome_id = $1;
	print "$curr_genome_id\n" if $debug;
	push( @genomeId, $curr_genome_id );
}
$ids = join( ',', @genomeId );

# ELink
$url = $base . "elink.fcgi?dbfrom=genome&db=nuccore&id=$ids&term=$limit&usehistory=y";
print "URL elink: $url\n" if $debug;
$xml       = get($url);

# create object
my $xmlIn = new XML::Simple( ForceArray => 1 );

# read XML file
my $xmlContent = $xmlIn->XMLin($xml);

my $id_nuccore = '';
foreach my $value ( @{$xmlContent->{'LinkSet'}->[0]->{'LinkSetDb'}->[0]->{'Link'}} ){
  $id_nuccore .= $value->{'Id'}->[0] . ",";
}

$id_nuccore =~ s/,$//;

# EFetch
$url = $base . "efetch.fcgi?db=nuccore&id=$id_nuccore&rettype=gb&retmode=text";
print "URL efetch: $url\n" if $debug;
$out = get($url);
	
open( OUT, ">$s.gb" );
print OUT $out;
close OUT;
sleep(5);
