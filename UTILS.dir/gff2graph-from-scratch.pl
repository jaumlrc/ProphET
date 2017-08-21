#!/bin/env perl
use strict;

use FindBin;
use lib "$FindBin::Bin/GFFLib";

use GD;
use Bio::SeqIO;
use Bio::Graphics;
use Bio::SeqFeature::Generic;
use Getopt::Std;
use GFFFile;

my $MAX_E_VALUE    = 1000;
my $DEFAULT_WIDTH  = 5000; 

my $DEFAULT_COLOR  = "blue";
my $LEFT_BORDER    = 0;


my $gff_file     		= $ARGV[0];
my $fasta_file   		= $ARGV[1];
my $outputFile  			= $ARGV[2];
my $Blast_matches 		= $ARGV[3];			
my $Phage_coordinates 	= $ARGV[4];		

### Gerar a figura inicial e a Linha baseada no Fasta:

my $width = $DEFAULT_WIDTH;		## Sets the width of figure
my $overallMaxLength = 0;		## Stores the length of the largest scaffold. This will be used to define the ratio bp/pixel
my %fasta_seq;					## Stores the length of each sequence. Not being used currently, but it will be used in the next versions of ProphET
my @fasta_seqnames;

######################
# Opening FASTA file and Processing

my $inSeqIO = Bio::SeqIO->new( -file => $fasta_file, '-format' => 'Fasta' );
while ( my $inSeq = $inSeqIO->next_seq() ) {
	if ( !defined( $inSeq->seq() ) || length( $inSeq->seq() ) == 0 ) {
		print STDERR $inSeq->id()
		  . " have an empty string as sequence. Sequence discarded\n";
	}
	else {
		my $seq_length = $inSeq->length();
		my $seq_name   = $inSeq->id();

		push( @fasta_seqnames, $seq_name );	##  Array with the names of the scaffolds
		$fasta_seq{$seq_name}{len} = $seq_length;

		$overallMaxLength = $seq_length if ( $overallMaxLength < $seq_length );
	}
}
$inSeqIO->close();

print STDERR "Largest sequence in the FASTA files has $overallMaxLength\n";


######################
# Generate the base panel in which all the widgets will be rendered
my $panel = Bio::Graphics::Panel->new(  
	-length      => $overallMaxLength, 	## Sets the lenght in bp to the length of the largest scaffold
	-key_style	 => 'between',
	-image_class => 'GD::SVG',
	-width       => $width,				## sets the length in pixels
	-pad_left    => 40,					
	-pad_right   => 200, 
	-pad_top     => 40,
	-pad_bottom  => 40,
	);

######################
# Adds the X-axis arrow indicating the coordinates on the scaffold
my $full_length = Bio::SeqFeature::Generic->new(-start=>1,-end=> $overallMaxLength);
  $panel->add_track($full_length,
                    -glyph   => 'arrow',
                    -tick    => 2,
                    -fgcolor => 'black',
                    -double  => 4,
                    -key     => "Genome_sequence",
                    -height  => 15,
                    );

######################
# Adds the track that will contain all the CDSs in the scaffold

 my $track = $panel->add_track(-glyph => 'graded_segments',
                                -label  => 0,
                                -bgcolor => 'blue',
                                -key => "GFF CDSs:",
                                -fgcolor => 'black',
                                -bump     => +1,
                                -height   => 8,

  								);
  open (INPUT1,"$gff_file");

  while (<INPUT1>) { # reads GFF file
    chomp;
    next if /^\#/;  # ignore comments
    my($name,$start,$end) = split /\t+/;
    my $feature = Bio::SeqFeature::Generic->new(-display_name=> $name,
                                                -start        => $start,
                                                -end          => $end,
    );                                           
    $track->add_feature($feature);
  }
  close (INPUT1);





######################
# Add track that will contain the BLAST results against ProphET database
 my $track = $panel->add_track(-glyph => 'graded_segments',
                                -label  => 0,
                                -bgcolor => 'green',
                                -fgcolor => 'black',
                                -key => "Blast Matches:",
                                -fgcolor => 'black',
                                -font2color => 'red',                                
                                -bump     => +1,
                                -height   => 8,

  								);
  open (INPUT2,"$Blast_matches");

  while (<INPUT2>) { # reads blast results file
    chomp;
    next if /^\#/;  # ignore comments
    my($name,$start,$end) = split /\t+/;
    my $feature = Bio::SeqFeature::Generic->new(-display_name=> $name,
                                                -start        => $start,
                                                -end          => $end,
    );                                           
    $track->add_feature($feature);
  }
  close (INPUT2);

######################
# Add track that will contain the boundaries of the prophage predictions
my $track = $panel->add_track(-glyph => 'graded_segments',
                                -label  => 0,
                                -bgcolor => 'cyan',
                                -key => "Predicted Phages:",
                                -fgcolor => 'black',
                                -bump     => +10,
                                -height   => 8,

  								);
  open (INPUT3,"$Phage_coordinates");

  while (<INPUT3>) { # reads prophage coordinates
    chomp;
    next if /^\#/;  # ignore comments
    my($name,$start,$end) = split /\t+/;
    my $feature = Bio::SeqFeature::Generic->new(-display_name=> $name,
                                                -start        => $start,
                                                -end          => $end,
    );                                           
    $track->add_feature($feature);
  }
  close (INPUT3);

open OUTPUT, ">" . $outputFile; 
print OUTPUT $panel->svg;
close OUTPUT;



