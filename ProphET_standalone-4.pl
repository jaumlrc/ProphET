#!/bin/env perl


use strict;
use Pod::Usage;
use Getopt::Long;
use FindBin;

=head1 NAME

ProphET is a user friendly algorithm to identify prophages within prokaryote genomes.

=head1 SYNOPSIS

usage: ProphET_standalone.pl --fasta <file> --gff <file> --out <file> [--help]       


=head1 OPTIONS

B<--fasta> - Bacterial genome Fasta file

B<--gff> - Bacterial GFF file

B<--out> - output

B<--help> - print this message B<(Optional)>

=head1 DESCRIPTION

B<Important! The Fasta and GFF file MUST have the exact same ID:>
B<Ex.:>

Fasta:

>NC_005362.1

TTGTTTGATCTAGATAAATTTTGGCAATTTTTTAATGCTGAGATGAAAAAAAGCTACAGCACGGTTGCCT

ATAATGCTTGGTTTAAAAATACTAAACCAATTTCCTTTAATAAAAAGACAAAAGAAATGATAATCGCTGT

GFF:

NC_005362.1     .       gene    1       1365    .       +       .       ID=LJ_RS00005;Name=LJ_RS00005;

NC_005362.1     .       mRNA    1       1365    .       +       .       ID=LJ_RS00005.t01;Parent=LJ_RS00005;

NC_005362.1     .       exon    1       1365    .       +       .       ID=LJ_RS00005.t01-E1;Parent=LJ_RS00005.t01;

NC_005362.1     .       CDS     1       1365    .       +       0       ID=LJ_RS00005.p01;Parent=LJ_RS00005.t01;

 
=head1 CONTACT
 
 Joao Luis R. Cunha (2015)
 jaumlrc@gmail.com
 jaumlrc@broadinstitute.org
 
 Gustavo C. Cerqueira (2015)
 cerca11@gmail.com
 gustavo@broadinstitute.org
=cut





my $fasta;   #Fasta file
my $gff;     #GFF file
my $name;    #Output name
my $help;


GetOptions(     'fasta=s' => \$fasta,
                        'gff=s'                => \$gff,
                        'out=s'                                 => \$name,
                        'help!'                                 => \$help);

if( defined($help) ){
   pod2usage(-verbose => 2 ,-exitval => 0);
}

if( not defined($fasta) ){
   pod2usage( -message => "Error: Parameter --fasta is required !!!!\n\n", -verbose => 1,-exitval => 1, -output => \*STDERR);
}

if( not defined($gff) ){
   pod2usage(-message => "Error: Parameter --gff is required !!!!\n\n", -verbose => 1 ,-exitval => 1, -output => \*STDERR);
}

if( not defined($name) ){
   pod2usage(-message=>"Error: Parameter --out is required !!!!\n\n", -verbose => 1 ,-exitval => 1, -output => \*STDERR);
}


my $debug = 0;

my $usage ="ProphET_standalone.pl <fasta file> <gff file> <output_name> ";

my $images_location = "Phages-$name.intermediate.logs";
mkdir("$images_location", 0755);


##########
#Checking the input files
#Checking if there is only one fasta sequence
#Checking all the gffs IDs matches the fasta sequence

my $fasta_count = `grep -c '>' $fasta`;
if ($fasta_count > 1) {die "The input fasta file has more than one genome sequence. Check if you are not submitting a file with plasmids as well as the genome file\n"};

my $gff_ids_count = `awk '{print \$1}' $gff | sort -u | wc | awk '{print \$1}'`;
my $gff_ids_ids = `awk '{print \$1}' $gff | sort -u`;
if ($gff_ids_count > 1) {die "The input gff file has more than one genome id:\n$gff_ids_ids Check if you are not submitting a file with plasmids as well as the genome file\n"};

my $fasta_only_id = `grep '>' $fasta | sed s'/>//'`;
my $gff_only_id = `awk '{print \$1}' $gff | sort -u`;
if ($fasta_only_id ne $gff_only_id) {die "The Fasta ID does not match the GFF id:\n Fasta_ID:$fasta_only_id GFF_ID:$gff_only_id"};




##########
#Module 1 
#-Input: Genome fasta and GFF file
#-Obtaining the bacterial predicted proteome, 
#-Blast against our phage database,
#-formating the output to be used as input to the next steps,
#-Setting paths to third party programs
#-Output: Phage matches coordinates

my $UTILS_dir = "$FindBin::Bin/UTILS.dir";
chomp $UTILS_dir;

my $root_dir = "$FindBin::Bin/root.dir";
chomp $root_dir;

my $prophet_database = "$FindBin::Bin/PhrophET_phage_proteins_database.dir";
chomp $prophet_database;

my $emboss_extractseq_path = `grep 'Emboss_extractseq_path' $root_dir/Third_party_programs_paths.log | awk '{print \$2}'`;
my $blastall_path = `grep 'Blastall_path' $root_dir/Third_party_programs_paths.log | awk '{print \$2}'`;
my $bedtools_path = `grep 'Bedtools_path' $root_dir/Third_party_programs_paths.log | awk '{print \$2}'`;

chomp $emboss_extractseq_path;
chomp $blastall_path;
chomp $bedtools_path;

`$UTILS_dir/./gff2gene_protein_seq.pl $gff $fasta 11 $name.trans $name.cds $name.prot`;  #Script to generate the proteins fasta based on the genome fasta and gff
`$blastall_path -p blastp -d $prophet_database/Phage_proteins_without_ABC-t.db -i $name.prot -e 1e-5 -m8 -a8 -o $name.blast`; #Blast predicted protein BLAST agains our phage db
`cat $name.blast | awk '{OFS="\t";print \$0, \$1 }' \$1 | sort -k 1,1 -k 12,12rn | uniq -f 12 | awk '{OFS="\t";print \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8, \$9, \$10, \$11, \$12 }' > $name.blastt-best-matches`;  #Extracting only BLAST best matches
`cat $name.blastt-best-matches | awk '{print \$1}' > $name.blastt-best-matches-ids`; #Extracting the gene ID
`grep -wf $name.blastt-best-matches-ids $gff | grep 'CDS' | awk '{print \$1"\t"\$4"\t"\$5}' | sort -nk 3,3 > $name.matches-coords`; #Extracting the gene ID coordinates
`$UTILS_dir/./union.pl --in $name.matches-coords --seg_name 1 --seg_start 2 --seg_end 3 > $name.matches-coords-union`; #Union comand - Thera are some genes with overlap, so we use union to do not count some positions twice
`grep 'tRNA' $gff | awk '{print \$1"\t"\$4"\t"\$5}' > tRNA_file_coordnates`; #Extraction tRNA coordinates from the GFF
my $tRNAfile = "tRNA_file_coordnates" ;

##########
# Module 2 - Sliding window and phage predictions:
#-Input: Coordinates with phage matches 
#-Sliding window and intersect with phage matches - Windows of 10,000 with 1,000 increment
#-Estimates the Density of nucleotides with phage matches in a given window
#-Output: Phage matches for each window


my $blast = "$name.matches-coords-union"; #File with the initial and termina coordenates of the genes with matches with our Phage db
open(BLAST, "$name.matches-coords-union");
my @genes = <BLAST>;
close(BLAST);

open(RESULTS, ">$blast-window-output"); #File with the number of bases with phage matches in a given window
open(LOG, ">$blast-log-matches");       #Log-file with the coordinates of matches in a given window


my @last_pos = split(/[\t]/, $genes[$#genes]); #As we use a 10000b window, we stop the analysis 10000b after the last match with a phage gene
my $id_genoma = $last_pos[0];
my $final_pos = $last_pos[2];
my $final_real = $final_pos+10000;      
chomp $id_genoma;

for (my $i=10000; $i<= $final_real; $i=$i+1000) {       #Sliding window to intersect with the phage matches coordinates
	open(WINDOW, ">$blast-window");   # Window overwrites at each iteration
	my $coord_ini = $i-9999;        
	chomp $coord_ini;
	chomp $i;
	print WINDOW "$id_genoma\t$coord_ini\t$i\n";
	close(WINDOW);
	
	`$bedtools_path intersect -a $blast -b $blast-window > $blast-intersect-entre-phage-janela`; # Makes the intersect between the current window and our phage-matches coordinates

	open (INTERSECTED, "$blast-intersect-entre-phage-janela"); #Open the intersect file
	my @merged = <INTERSECTED>;
	close(INTERSECTED);

	print LOG "#Window\t$coord_ini\t$i\n"; #Stores the coordinates with phage matches in a given window in the log file
	
	my $gene_content = 0;
	for (my $j=0; $j<=$#merged; $j++) {           #Runs trough all the genes with phage matches in the current window 
		my @splitted = split(/[\t]/, $merged[$j]);
		my $initial = $splitted[1];
		my $final = $splitted[2];
		my $size = ($final+1)-$initial;   #Stores the size of the gene with phage matches (could be more than one per window)
		$gene_content = $gene_content+$size;  #Stores the number whole extent of genes with phage matches in a give window
                chomp $initial;
                chomp $final;
                print LOG "$id_genoma\t$initial\t$final\n"; #Log with the coordinates of the genes with phage matches in each window
        }

        print LOG "--------------------------------\n"; #Inter-window log separator
        print RESULTS "$id_genoma\t$coord_ini\t$i\t$gene_content\n"; #Tem output file - containing the number of phage-related nucleotides in each window

}

close(LOG);
close(RESULTS);

##########
# Module 3 - Cluster consecutive windows with phage matches
#-Input: Phage matches for each window
#-Group consecutive windows with phage matches given the phage content in a given window is higher than half of the window (5000b), generatin raw-clusters
#-Then, group raw-clusters with overlaps
#-Output - Raw phage coordinates prediction

open(INPUT, "$blast-window-output");
my @matches = <INPUT>;
close(INPUT);

open(TEMP, ">$blast-merged");

for (my $k=0; $k<=$#matches; $k++) {   #Loop to cluster consectutive windows with "phage content" higher than 5000pb, half of the window
    my @splitted = split(/[\t]/, $matches[$k]);
	my $coverage= $splitted[3];
	my $ini_coordenates= $splitted[1];
	my $fin_coordenates= $splitted[2];
	my $id = $splitted[0];

        if ($coverage > 5000){                          
			while ($coverage > 5000 && $k <$#matches) {     #while coverage>5000, group consecutive windows
	        my @splitted2 = split(/[\t]/, $matches[$k]);
        	$coverage= $splitted2[3];
        	$k ++;
        }

        my @splitted3 = split(/[\t]/, $matches[$k-2]); #Get the initial, final coordinates and size of the clustered windows
        my $fin_coordenates3= $splitted3[2];
        my $size = ($fin_coordenates3+1) - $ini_coordenates; #### Nao tenho usado o size para nada...
        chomp $ini_coordenates;
        chomp $fin_coordenates3;
        chomp $id;
        chomp $size;
 
 	print TEMP "$id\t$ini_coordenates\t$fin_coordenates3\n";
        }
}
close(TEMP);

my @unnited_gustavo =`$UTILS_dir/./union.pl --in $blast-merged --seg_name 1 --seg_start 2 --seg_end 3 > $blast-merged-1 `; #Merge regions raw-clusters with overlap in sequences

##########
# Module 4 - Trimming raw phages with less than 8 genes with phage matches
#-Input: Raw phage coordinates prediction
#-Check in the raw phage has at least 8 genes with phage matches
#-Output: raw phages that have at least 8 genes with phage matches

open(MERGED, "$blast-merged-1") or die "couldnt open merged1\n";
my @merged1 = <MERGED>;
close(MERGED);

open(BLAST, "$blast") or die "couldnt open file_united";
my @blast=<BLAST>;
close(BLAST);

open(TEMP2, ">$blast-merged-united-2");


 for (my $l=0; $l<= $#merged1; $l++) {
	my @P_positions= split(/[\t]/, $merged1[$l]);
        my $P_start = $P_positions[1];
        my $P_end = $P_positions[2];
        my $contador = 0;

       for (my $m=0; $m <= $#blast; $m++) {
       	    my @B_positions = split(/[\t]/, $blast[$m]);
            my $B_start = $B_positions[1];
            my $B_end = $B_positions[2];
            
	    if ($B_start > $P_start && $B_end < $P_end)  {
		$contador++;
            }                                     
       	    if ($contador == 8) {
            print TEMP2 "$merged1[$l]";
            last;
            }
        }      
}
close(TEMP2);

##########
# Module 5 - Trimming raw phages borders by tRNAs or last gene with phage match
#-Input: Raw phage coordinates prediction with at least 8 genes
#-Trimns the phage border to the last gene with phage matches and then searches 3kb upstream and downstream for tRNA genes, extending accordingly
#-Output: Polished final phage prediction


open(PHAGES, "$blast-merged-united-2");
my @phages = <PHAGES>;
close(PHAGES);

open(MATCHES, "$name.matches-coords-union");
my @matches = <MATCHES>;
close(MATCHES);

open(TRNA, "$tRNAfile");
my @tRNA = <TRNA>;
close(TRNA);

open(OUTLOG, ">$name-trna.log");

open(PHAGESFINAL, ">$name-phages-final-prediction");

for (my $n=0; $n<=$#phages; $n++) {			#Searches for the last gene with phage match in the beggining and end of the raw phage
	my @phage_coords = split(/[\t]/, $phages[$n]);
	my $phage_id = $phage_coords[0];
	my $phage_initial = $phage_coords[1];
	my $phage_final = $phage_coords[2];

	my @inside_initial =();
	my @inside_final =();

	for (my $p=0; $p<=$#matches; $p++) {
	my @matches_coords = split(/[\t]/, $matches[$p]);
	my $matches_id = $matches_coords[0];
	my $matches_initial = $matches_coords[1];
	my $matches_final = $matches_coords[2];

		if ($matches_final > $phage_initial && $matches_final < $phage_final) {
			push(@inside_initial, $matches_coords[1]);
		}	
	
		if ($matches_initial > $phage_initial && $matches_initial < $phage_final) {
			push(@inside_final, $matches_coords[2]);
		}	
	}

	$phage_initial = $inside_initial[0];
	my $number = $#inside_final;
	$phage_final = $inside_final[$number];
	my $phage_beginnin = $phage_initial; # Para poder alterar com o tRNA
	my $phage_ending = $phage_final; # Para poder alterar com o tRNA

	my $inital_minus = $phage_initial - 3000;
	my $initial_plus = $phage_initial + 3000;
	my $final_minus = $phage_final - 3000;
	my $final_plus = $phage_final + 3000;

	my $distancia = 3000; #distancia inicial para permitir o tRNA
	my $distancia2 = 3000;
	my $diff = 3000;  #diferenca para mudar com o tRNA
	my $diff2 = 3000;

	for (my $o=0; $o<=$#tRNA; $o++) { 	# Checks for tRNA genes near phage borders
		my @tRNA_coords = split(/[\t]/, $tRNA[$o]);
		my $trna_id = $tRNA_coords[0];
		my $trna_initial = $tRNA_coords[1];
		my $trna_final = $tRNA_coords[2];

		if ($trna_initial > $inital_minus && $trna_final < $initial_plus ) {
			if ($trna_initial < $phage_initial && $trna_final < $phage_initial) {
				$diff= $phage_initial - $trna_final;

				if ($distancia > $diff ) {
					$distancia = $diff;

					$phage_beginnin = $trna_initial;
                }
			}
			elsif ($trna_initial > $phage_initial && $trna_final > $phage_initial) {
				$diff = $trna_initial - $phage_initial;
				if ($distancia > $diff) {
					$distancia = $diff;
					$phage_beginnin = $trna_initial;
				}
			}	

			elsif ($trna_initial < $phage_initial && $trna_final > $phage_initial) {
				$diff =0;
				$phage_beginnin = $trna_initial;
			}

			print OUTLOG $tRNA[$o];
		}

		elsif ($trna_initial > $final_minus && $trna_final < $final_plus ) {
                        if ($trna_initial < $phage_final && $trna_final < $phage_final) {
				$diff2 = $phage_final - $trna_final;
				if ($distancia2 > $diff2) {
					$distancia2 = $diff2;
					$phage_ending = $trna_final;
				}
			}
			
			elsif ($trna_initial > $phage_final && $trna_final > $phage_final ) {
				$diff2 = $trna_initial - $phage_final;
				if ($distancia2 > $diff2) {
					$distancia2 = $diff2;
					$phage_ending = $trna_final;
				}
			}
			
			elsif ($trna_initial < $phage_final && $trna_final > $phage_final) {
				$diff2=0;
				$phage_ending = $trna_final;
			}
		
			print OUTLOG $tRNA[$o];
		}
		
	}



$phage_initial = $phage_beginnin;
$phage_final = $phage_ending;

        chomp $phage_id;
	chomp $phage_beginnin;
	chomp $phage_ending;

	print PHAGESFINAL "$phage_id\t$phage_beginnin\t$phage_ending\n";

}

##########
# Module 6 - Generate the image file and organizes the outputs file
#-Input: Final phage prediction, genome fasta, gff file
#-Generates the image output file
#-Output: Image output file



`cat $gff | grep 'CDS' | awk '{print \$1"\t"\$4"\t"\$5}' > $gff-ultraformated`;
`$UTILS_dir/./gff2graph-from-scratch.pl $gff-ultraformated $fasta $name-phages.svg $name.matches-coords-union $name-phages-final-prediction`;

#open(FINALPROGRAM, ">$name-run-sucessfuly-completed.log");
#print FINALPROGRAM "The Run was sucessfull";
#close(FINALPROGRAM);

open(PHAGES, "$name-phages-final-prediction");
my @fagos = <PHAGES>;
close(PHAGES);

mkdir("$name-Phage_predictions.dir", 0755);

for (my $z=0; $z<= $#fagos; $z++) {
	my @separated = split(/[\t]/, $fagos[$z]);
	my $namef = $separated[0];
	my $initf = $separated[1];
	my $finaf = $separated[2];

	my $number = $z+1;
	
	chomp $namef;
	chomp $initf;
	chomp $finaf;
	chomp $number;
	chomp $fasta;

	`$emboss_extractseq_path -sequence $fasta -regions \"$initf-$finaf\" -osdbname2 $namef-$number $namef-$number` ;
	`mv $namef-$number $name-Phage_predictions.dir`;
}

`cp $name-phages.svg $images_location`;
`cp $name-phages-final-prediction $images_location`;


`mv $name-trna.log $images_location`; #ok
`mv $name.matches-coords-union-log-matches $images_location`; #ok
`mv $name.matches-coords-union-merged-united-2 $images_location`; #ok
`mv $name.matches-coords-union-window-output $images_location`; #ok
`mv $name.blast $images_location`; 
`mv $name.blastt-best-matches $images_location`;
`mv $name.cds $images_location`;
`mv $name.prot $images_location`;
`mv $name.trans $images_location`;
`mv tRNA_file_coordnates $images_location`;

`rm $name.matches-coords-union-window`;
`rm $name.matches-coords-union-merged-1`;
`rm $name.matches-coords`;
`rm $name.matches-coords-union`;
`rm $name.blastt-best-matches-ids`;
`rm $name.matches-coords-union-merged`;
`rm $name.matches-coords-union-intersect-entre-phage-janela`;
`rm $gff-ultraformated`;
#`rm error.log`;

`mv $images_location $name-Phage_predictions.dir`;
`mv $name-phages.svg $name-Phage_predictions.dir`;
`mv $name-phages-final-prediction $name-Phage_predictions.dir`;
