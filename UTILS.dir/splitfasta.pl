#!/usr/bin/env perl

use strict;

use Bio::Seq;
use Bio::SeqIO;
use IO::Handle;

use Pod::Usage;
use Getopt::Long;
use FindBin;

=head1 NAME

Splits the fasta sequences in term number of sequences or number of files.

=head1 SYNOPSIS

usage: splitfasta.pl --fasta_in <file> [ --num_seq_per_file <number of sequences> |  --num_files <number of files> ] [--out_dir] --pre <prefix>


=head1 OPTIONS

B<--fasta_in> - FASTA input

B<--num_seq> - Number of sequences in each file B<(Optional)>

B<--num_files> - Number of file to be divided the original file B<(Optional)>

B<--pre> - prefix of the output files

B<--out_dir> - output directory B<(Optional)>

B<--help> - print this message B<(Optional)>

=head1 DESCRIPTION

$num_seq_per_file_per_file
=head1 CONTACT
 Gustavo C. Cerqueira (2015)
 cerca11@gmail.com
 gustavo@broadinstitute.org
=cut


my $fasta_in;
my $num_seq_per_file;
my $num_files;
my $pre;
my $help;
my $output_dir;


GetOptions(
'fasta_in=s'  => \$fasta_in,
'num_seq_per_file=s'   => \$num_seq_per_file,
'num_files=s' => \$num_files,
'out_dir=s' => \$output_dir,
'pre=s'       => \$pre,
'help!'       => \$help
);

if ( defined($help) ) {
	pod2usage( -verbose => 2, -exitval => 0 );
}


if ( not defined($fasta_in) ) {
	pod2usage(
	-message => "Error: Parameter --fasta_in is required !!!!\n\n",
	-verbose => 1,
	-exitval => 1,
	-output  => \*STDERR
	);
}

if ( not defined($pre) ) {
	pod2usage(
	-message => "Error: Parameter --pre is required !!!!\n\n",
	-verbose => 1,
	-exitval => 1,
	-output  => \*STDERR
	);
}

if ( defined($num_seq_per_file) && defined($num_files)) {
  pod2usage(
	-message => "Error: Either use the parameter --num_seq or --num_files!!!!\n\n",
	-verbose => 1,
	-exitval => 1,
	-output  => \*STDERR
	);
}

if ( not defined($output_dir) ) {
	$output_dir = '.';
}


if ( not defined($num_seq_per_file) && not defined($num_files)) {
  $num_seq_per_file = 1;
}


# Counting the number of sequences in the FASTA file
my $num_seq_in_orig_file = 0;
my $inSeqIO     = Bio::SeqIO->new(-file => $fasta_in, '-format' => 'Fasta');
while ( my $inSeq = $inSeqIO->next_seq() ){
	$num_seq_in_orig_file++ if( defined( $inSeq->seq() ) );
}
$inSeqIO->close();

# Calculate the number of sequences per each file
# given the number of files as result of split and the number of sequences
# in the original file
if( defined ( $num_files ) ){
    $num_seq_per_file  = int( $num_seq_in_orig_file / $num_files );
    $num_seq_per_file = 1 if $num_seq_per_file == 0;
}


# Generating files
my $count_seqs = 1;
my $count_files = 1;
my $seqs_written = 0; # Number of seqs. already read from the original final and written in one of the segmented files

my $out_file = "$output_dir/$pre.$count_files.fas";

mkdir($output_dir) unless(-d $output_dir);

open OUT, ">$out_file";

$inSeqIO     = Bio::SeqIO->new(-file => $fasta_in, '-format' => 'Fasta');
while ( my $inSeq = $inSeqIO->next_seq() ){

  if( !defined( $inSeq->seq() ) ){
    print STDERR  $inSeq->id() . " have an empty string as sequence. Sequence discarded\n";
  } else {
    print OUT ">" . $inSeq->id() . "\n" . $inSeq->seq() . "\n";
    $count_seqs++;
    $seqs_written++;
  }
  
  
  # if max number of seqs. per file was reached
  if( $count_seqs >= $num_seq_per_file ){


	# ... if the number of files was defined, meaning that the number of sequences is the directive of this script
	# close the current file and open another one
  	if( not defined($num_files) ){
  		close(OUT);
  		
  		last if ( $seqs_written == $num_seq_in_orig_file );
  		
	  	$count_files++;  	
		my $out_file = "$output_dir/$pre.$count_files.fas";
		open OUT, ">$out_file";
	  	$count_seqs = 1;
  	}else{

	  # ... if the number of seq per file was reached but this is not the last file
	  # close the current file and open another
	  if( $count_files != $num_files ){
  		close(OUT);

  		last if ( $seqs_written == $num_seq_in_orig_file );
  		
	  	$count_files++;  	
		my $out_file = "$output_dir/$pre.$count_files.fas";
		open OUT, ">$out_file";
	  	$count_seqs = 1;
  	
  	  }
  	}	
  }
}
close(OUT);
$inSeqIO->close();

print "Number of files created: $count_files\n";
