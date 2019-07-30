#!/usr/bin/env perl

use strict;
use LWP::Simple;
use XML::Simple;

#use Data::Dumper;

my ( $name, $outname, $url, $xml, $out, $count, $query_key, $webenv, $ids );
my @genomeId;
my @genomeId_nuccore;
my $base = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/';

#my $limit = 'srcdb+refseq[prop]+AND+gene+in+chromosome[prop])';
my $limit = 'srcdb+refseq[prop]';

my $tax_id              = $ARGV[0];
my $debug               = 0;
my $DOWNLOAD_INCREMENTS = 500;
my $delay = 1;

open( LOG, ">$tax_id.ncbi_utils.log" )
  or die "ERROR: Unable to write on file $tax_id.ncbi_utils.log\n";

undef @genomeId;
$query_key = $webenv = '';
$tax_id =~ s/ /+/g;

# ESearch
$url =
    $base
  . "esearch.fcgi?db=genome&term=txid"
  . $tax_id
  . "[Organism:exp]&usehistory=y";
print "URL esearch 1: " . $url . "\n" if $debug;
sleep($delay);
$xml = get($url);

if ( $xml =~ /<Count>(\d+)<\/Count>/ ) {
	$count = $1;
}

print "Number of records in Genome database: $count\n";

if ( $count > 20 ) {
	$url =
	    $base
	  . "esearch.fcgi?db=genome&term==txid"
	  . $tax_id
	  . "[Organism:exp]&retmax=$count&usehistory=y";
	print "URL esearch 2: " . $url . "\n" if $debug;
	sleep($delay);
	$xml = get($url);
}
while ( $xml =~ /<Id>(\d+?)<\/Id>/gs ) {
	my $curr_genome_id = $1;
	print "$curr_genome_id\n" if $debug;
	push( @genomeId, $curr_genome_id );
}

#-----------------------------------------------------------------------
# Converting genome Ids to nuccore

my $num_genomeId = scalar(@genomeId);
for ( my $ind = 0 ; $ind < $num_genomeId ; $ind += $DOWNLOAD_INCREMENTS ) {
	my $last = $ind + ( $DOWNLOAD_INCREMENTS - 1 );
	$last = $num_genomeId - 1 if $last >= $num_genomeId;

	print "Converting genomeids "
	  . ( $ind + 1 ) . " to "
	  . ( $last + 1 ) . "...\n";

	my $ids = join( ',', @genomeId[ $ind .. $last ] );

	# ELink
	$url = $base
	  . "elink.fcgi?dbfrom=genome&db=nuccore&id=$ids&term=$limit&usehistory=y";
	print "URL elink: $url\n" if $debug;
	sleep($delay);
	$xml = get($url);

	# create object
	my $xmlIn = new XML::Simple( ForceArray => 1 );

	# read XML file
	my $xmlContent = $xmlIn->XMLin($xml);

	foreach my $value (
		@{ $xmlContent->{'LinkSet'}->[0]->{'LinkSetDb'}->[0]->{'Link'} } )
	{
		my $curr_id = $value->{'Id'}->[0];
		push( @genomeId_nuccore, $curr_id );
		print $curr_id . "\n" if $debug;
	}

}
getc() if $debug;


#-----------------------------------------------------------------------

# Downloading genomes
$num_genomeId = scalar(@genomeId_nuccore);


print "Number of genomes under TaxID $tax_id: " . $num_genomeId . "\n";
print LOG "Number of genomes under TaxID $tax_id: " . $num_genomeId . "\n";
if ( -e "$tax_id.gb" ) {
	`rm $tax_id.gb`;
}

getc() if $debug;
sleep($delay);

for ( my $ind = 0 ; $ind < $num_genomeId ; $ind += $DOWNLOAD_INCREMENTS ) {
	my $last = $ind + ( $DOWNLOAD_INCREMENTS - 1 );
	$last = $num_genomeId - 1 if $last >= $num_genomeId;

	print "Downloading genomes "
	  . ( $ind + 1 ) . " to "
	  . ( $last + 1 ) . "...\n";

	my $ids = join( ',', @genomeId_nuccore[ $ind .. $last ] );

	# EFetch
	$url = $base . "efetch.fcgi?db=nuccore&id=$ids&rettype=gb&retmode=text";
	print "URL efetch: $url\n" if $debug;
	sleep($delay);
	$out = get($url);

	open( OUT, ">>$tax_id.gb" );
	print OUT $out;
	close OUT;
	sleep(10);
}

close(LOG);

getc() if $debug;


