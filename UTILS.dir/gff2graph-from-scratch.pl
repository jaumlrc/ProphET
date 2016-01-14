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
#my $DEFAULT_WIDTH  = 1500;
my $DEFAULT_WIDTH  = 5000; #retirei 2 zeros (era 500000) para diminuir a figura

my $DEFAULT_COLOR  = "blue";
my $LEFT_BORDER    = 0;


my $gff_file     = $ARGV[0];
my $fasta_file   = $ARGV[1];
my $outputFile   = $ARGV[2];
my $Blast_matches = $ARGV[3];			
my $Phage_coordnates = $ARGV[4];		

### Gerar a figura inicial e a Linha baseada no Fasta:

my $width = $DEFAULT_WIDTH;		## Joao - Tamanho da figura;
my $overallMaxLength = 0;		## Joao - COntador inicial do maior genoma do arquivo - para traçar a linha

my %fasta_seq;					## Joao - Hash inicializado apenas uma vez. Nao mais acessado.
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

		push( @fasta_seqnames, $seq_name );	## Joao - hash contendo o nome e tamanho das sequencias fasta
		$fasta_seq{$seq_name}{len} = $seq_length;

		$overallMaxLength = $seq_length if ( $overallMaxLength < $seq_length );
	}
}
$inSeqIO->close();

print STDERR "Largest sequence in the FASTA files has $overallMaxLength\n";


######################
# Gera o painel onde serão plotadas as figuras
my $panel = Bio::Graphics::Panel->new(  ## Joao - gera o painel onde serão geradas as figuras
	-length      => $overallMaxLength, 	## Joao - Coloca a linha referente ao maior tamanho da sequencia
	-key_style	 => 'between',
	-image_class => 'GD::SVG',
	-width       => $width,				## Joao - Usa o tamanho da imagem previamente setado
	-pad_left    => 40,					## Joao - Pads sao os deslocamentos das bordas
	-pad_right   => 200, 
	-pad_top     => 40,
	-pad_bottom  => 40,
	);
######################
# Gera a seta com base no fasta
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
# Gera os resultados do GFF

 my $track = $panel->add_track(-glyph => 'graded_segments',
                                -label  => 0,
                                -bgcolor => 'blue',
                                -key => "GFF CDSs:",
                                -fgcolor => 'black',
                                -bump     => +1,
                                -height   => 8,

  								);
  open (INPUT1,"<$ARGV[0]");

  while (<INPUT1>) { # read blast file
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
# Gera os resultados do BLAST
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
  open (INPUT2,"<$ARGV[3]");

  while (<INPUT2>) { # read blast file
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
# Gera os resultados do Phage prediction
my $track = $panel->add_track(-glyph => 'graded_segments',
                                -label  => 0,
                                -bgcolor => 'cyan',
                                -key => "Predicted Phages:",
                                -fgcolor => 'black',
                                -bump     => +10,
                                -height   => 8,

  								);
  open (INPUT3,"<$ARGV[4]");

  while (<INPUT3>) { # read blast file
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

open OUTPUT, ">" . $outputFile;   ### Joao - Aqui termina de gerar as figuras. O que esta abaixo eh subrotina
print OUTPUT $panel->svg;
close OUTPUT;



