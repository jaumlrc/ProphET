#!/usr/bin/env perl
use strict;

use FindBin;
use lib "$FindBin::Bin/GFFLib";


use GFFFile;
use GFFUtils;
use Bio::SeqIO;
use Bio::Seq;
use Bio::Tools::CodonTable;

my $usage =
    "\nusage: gff2gene_protein_seq.pl  <GFF file> <fasta in> "
  . "<either a codon table number or a tab delimited file* indicating contig/chrom and codon tables to be used> "
  . "<trasncript seq. file> <CDS seq. file> <peptide seq. file> [<extend transcript by X bp on both ends, usefull for making fake, fixed length UTRs>] [1=add gene name to header] \n\n";

$usage .= "Available codon tables:\n";

my $tables = Bio::Tools::CodonTable->tables;
while ( my ( $id, $name ) = each %{$tables} ) {
	$usage .= "$id = $name\n";
}
$usage .= "\n";

$usage .= "*File format:\n";
$usage .= "<chrom>\\t<codon table number>\\n\n";
$usage .= "\n";

die $usage if ( scalar(@ARGV) != 6 && scalar(@ARGV) != 7 && scalar(@ARGV) != 8 );

my $gff_filename               = $ARGV[0];
my $genome_filename            = $ARGV[1];
my $codon_table_file_or_number = $ARGV[2];
my $out_trans                  = $ARGV[3];
my $out_cds                    = $ARGV[4];
my $out_pep                    = $ARGV[5];
my $transcript_buffer          = $ARGV[6];
my $add_gene_name_to_header    = $ARGV[7];


my $single_codon_table = -1;
my %seq_codon_table;

if ( defined $transcript_buffer ) {
	die "\n\nERROR: the length of transcript extension should be >= 1\n\n"
	  if ( $transcript_buffer <= 0 );
}


# Testing if codon table parameter is a number
if( $codon_table_file_or_number =~ /^\d+$/ ){
	print "\n\nINFO: Using a single codon table to all chromosomes: $codon_table_file_or_number\n\n";
	$single_codon_table = $codon_table_file_or_number;

}else{
	print "\n\nINFO: Codon table FILE provided by user: $codon_table_file_or_number.\nUsing multiple codon tables\n";
	  
  	open CODON, "$codon_table_file_or_number" or die "\n\nERROR: Not able to open file: $codon_table_file_or_number!!!\n\n";
	my $cont_line = 1;
	while (<CODON>) {
		my $line = $_;
		chomp $line;

		my @cols = split "\t", $line;

		die "\n\nERROR: Error on line $cont_line: $line\n"
		  . "Every line should have the following format:\n"
		  . " <chrom>\\t<codon table>\n\n"
		  if ( scalar(@cols) != 2 or $cols[0] eq "" or $cols[1] eq "" );

		$seq_codon_table{ $cols[0] } = $cols[1];

		print STDERR
		  "Using codon table \'$cols[1]\' for chromosome \'$cols[0]\'\n";

		$cont_line++;
	}
	close(CODON);

	print STDERR "\n";

	
}


my $gffFile = GFFFile::new($gff_filename);

print "Reading GFF file...\n";
$gffFile->read();

# Read FASTA file
print "Reading FASTA file...\n";
my %fastaSeq;
my @fastaSeqOrder;
my $seq_in = Bio::SeqIO->new(
	'-file'   => $genome_filename,
	'-format' => "fasta"
);

while ( my $inseq = $seq_in->next_seq ) {
	$fastaSeq{ $inseq->id() } = $inseq->seq();
	push( @fastaSeqOrder, $inseq->id() );
}

$seq_in->close();

my $gffGenes = $gffFile->get_genes_hash();

# Ordering genes based on template name and start coord
my @gffGenesArray = values %{$gffGenes};
GFFUtils::sort_gene_arrays( \@gffGenesArray, 0 );

print "Number of genes in the GFF file: " . scalar(@gffGenesArray) . "\n";

open OUT_TRAN, ">$out_trans" or die "Unable to open file $out_trans to write\n";
open OUT_CDS,  ">$out_cds"   or die "Unable to open file $out_cds to write\n";
open OUT_PEP,  ">$out_pep"   or die "Unable to open file $out_pep to write\n";

my $cont_trans = 0;
my $cont_cds   = 0;

for my $currGene (@gffGenesArray) {
	my $gene_id = $currGene->get_id();
	my $chrom   = $currGene->get_chrom();
	my $strand  = $currGene->get_strand();

	die "Unable to find $chrom in the FASTA file $genome_filename"
	  if not defined $fastaSeq{$chrom};

	my $fastaSeqStr = $fastaSeq{$chrom};

	my $gffTranscripts = $currGene->get_transcripts_hash();
	for my $currTranscript ( values %{$gffTranscripts} ) {
		my $transcript_id = $currTranscript->get_id();
		$transcript_id = 
		$currTranscript->get_id() . " " . $currTranscript->get_name() if ( $add_gene_name_to_header == 1 );
		
		my $transcript_seq = "";
		for my $currExon ( @{ $currTranscript->get_exon_array() } ) {
			$transcript_seq .= substr(
				$fastaSeqStr,
				$currExon->get_start() - 1,
				$currExon->get_end() - $currExon->get_start() + 1
			);
		}

		if ( defined $transcript_buffer && $transcript_buffer >= 1 ) {
			my $curr_start = $currTranscript->get_start();
			my $curr_end   = $currTranscript->get_end();

			# Setting effective start = start - buffer.
			# Taking care of genes near beginning of sequence
			my $effective_start = $curr_start - $transcript_buffer;
			$effective_start = 1 if $effective_start < 1;

			# Setting effective end = end - buffer.
			# Taking care of genes near beginning of sequence
			my $effective_end = $curr_end + $transcript_buffer;
			$effective_end = length($fastaSeqStr)
			  if $effective_end > length($fastaSeqStr);

			# Adding buffer sequence to original transcript sequence
			$transcript_seq =
			    substr( $fastaSeqStr, $effective_start, $transcript_buffer - 1 )
			  . $transcript_seq
			  . substr( $fastaSeqStr, $curr_end, $effective_end - $curr_end );

		}

		my $cds_seq = "";
		for my $currCDS ( @{ $currTranscript->get_CDS_array() } ) {
			$cds_seq .= substr(
				$fastaSeqStr,
				$currCDS->get_start() - 1,
				$currCDS->get_end() - $currCDS->get_start() + 1
			);
		}

		if ( $strand eq "-" ) {
			$transcript_seq = reverse uc($transcript_seq);
			$transcript_seq =~ tr/AGCT/TCGA/;

			$cds_seq = reverse uc($cds_seq);
			$cds_seq =~ tr/AGCT/TCGA/;
		}

		# Tranlating CDS
		my $cdsSeqObj = Bio::Seq->new(
			-seq        => $cds_seq,
			-display_id => "temp",
			-alphabet   => "dna"
		);

		my $codon_table;

		if ( $single_codon_table != -1 ) {
			$codon_table = $single_codon_table;
		}
		else {
			$codon_table = $seq_codon_table{$chrom};
			die
			  "Error. Using file to define the codon table of each chromosome. "
			  . "But codon table for chromosome \'$chrom\' not found\n"
			  if not defined $codon_table;
		}

		my $pepObj  = $cdsSeqObj->translate( -codontable_id => $codon_table );
		my $pep_seq = $pepObj->seq();

		# Writing to files
		print OUT_TRAN ">$transcript_id\n$transcript_seq\n";
		$cont_trans++;

		if ( scalar @{ $currTranscript->get_CDS_array() } != 0 ) {
			print OUT_CDS ">$transcript_id\n$cds_seq\n";
			print OUT_PEP ">$transcript_id\n$pep_seq\n";
			$cont_cds++;
		}

	}
}
close(OUT_TRAN);
close(OUT_CDS);
close(OUT_PEP);

print "Number of transcripts sequences: $cont_trans\n";
print "Number of protein/CDS sequences: $cont_cds\n";
