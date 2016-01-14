#!/bin/env perl

use strict;
use LWP::Simple;

my ( $name, $outname, $url, $xml, $out, $count, $query_key, $webenv, $ids );
my @genomeId;
my $base  = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
#my $limit = 'srcdb+refseq[prop]+AND+gene+in+chromosome[prop])';
my $limit = 'srcdb+refseq[prop]';

my $s = $ARGV[0];

undef @genomeId;
$query_key = $webenv = '';
$s =~ s/ /+/g;

# ESearch
$url = $base . "esearch.fcgi?db=genome&term=txid" . $s . "[Organism:exp]&usehistory=y";
print "URL1: " . $url . "\n";
$xml = get($url);
$count = $1 if ( $xml =~ /<Count>(\d+)<\/Count>/ );
if ( $count > 20 ) {
	$url =
	    $base
	  . "esearch.fcgi?db=genome&term==txid"
	  . $s
	  . "[Organism:exp]&retmax=$count&usehistory=y";
	print "URL2: " . $url . "\n";
	$xml = get($url);
}
while ( $xml =~ /<Id>(\d+?)<\/Id>/gs ) {
	push( @genomeId, $1 );
}
$ids = join( ',', @genomeId );

# ELink
$url = $base . "elink.fcgi?dbfrom=genome&db=nuccore&cmd=neighbor_history&id=$ids&term=$limit&usehistory=y";
$xml       = get($url);
print "URL elink: $url\n";
$query_key = $1 if ( $xml =~ /<QueryKey>(\d+)<\/QueryKey>/ );
$webenv    = $1 if ( $xml =~ /<WebEnv>(\S+)<\/WebEnv>/ );

# EFetch
	$url = $base . "efetch.fcgi?db=nuccore&query_key=$query_key&WebEnv=$webenv&rettype=gb&retmode=text";
	print "URL efetch: $url\n";
	$out = get($url);
	
open( OUT, ">$s.gb" );
print OUT $out;
close OUT;
sleep(5);
