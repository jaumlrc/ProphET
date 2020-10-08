#!/usr/bin/env perl

use strict;
use Pod::Usage;

use Getopt::Long;
use FindBin;

use lib "$FindBin::RealBin/UTILS.dir/GFFLib"; 

use GFFFile;


=head1 NAME

ProphET is a user friendly algorithm to identify prophages in bacterial genomes.

=head1 SYNOPSIS

usage: ProphET_standalone.pl --fasta_in <file> --gff_in <file> --outdir <string> [--evalue_cutoff <number>] [--window_size <number>]> [--grid] [--help]


=head1 OPTIONS

B<--fasta_in> - Bacterial genome FASTA file

B<--gff_in> - Bacterial GFF file

B<--outdir> - output directory

B<--grid> - Use UGER for BLAST jobs (Currently only works in the Broad Institute UGER grid system) B<(Optional)>

B<--evalue_cutoff> - Sets the BLAST evalue cutoff. Please do not use scientific notation (default=0.00001) B<(Optional)>

B<--window_size> - Sets the size of the sliding windown used to calculate gene density (minimum value = 1000, default = 10000) B<(Optional)>

B<--help> - print this and some additional info. about FASTA and GFF input format B<(Optional)>

=head1 DESCRIPTION

B<Important! The FASTA and GFF file MUST have the exact scaffold/chrom IDs.
The GFF file should have the format described by Sequence Ontology consortium:
https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md >

B<Ex.:>

 FASTA:

 >NC_005362.1
 TTGTTTGATCTAGATAAATTTTGGCAATTTTTTAATGCTGAGATGAAAAAAAGCTACAGCACGGTTGCCT
 ATAATGCTTGGTTTAAAAATACTAAACCAATTTCCTTTAATAAAAAGACAAAAGAAATGATAATCGCTGT

GFF:

 NC_005362.1     .       gene    1       1365    .       +       .       ID=LJ_RS00005;Name=LJ_RS00005;
 NC_005362.1     .       mRNA    1       1365    .       +       .       ID=LJ_RS00005.t01;Parent=LJ_RS00005;
 NC_005362.1     .       exon    1       1365    .       +       .       ID=LJ_RS00005.t01-E1;Parent=LJ_RS00005.t01;
 NC_005362.1     .       CDS     1       1365    .       +       0       ID=LJ_RS00005.p01;Parent=LJ_RS00005.t01;

=head1 CONTACT

 Joao Luis R. Cunha (2017)
 jaumlrc@gmail.com

 Gustavo C. Cerqueira (2017)
 cerca11@gmail.com
=cut

my $fasta_in;    #Fasta file
my $gff_in;      #GFF file
my $gff_trna;    #GFF file with tRNA coordinates
my $outdir;      #Output name

my $evalue_cutoff;
my $window_size;

my $help;
my $grid;

GetOptions(
	'fasta_in=s'       => \$fasta_in,
	'gff_in=s'         => \$gff_in,
	'gff_trna=s'       => \$gff_trna,
	'outdir=s'         => \$outdir,
	'grid'             => \$grid,
	'evalue_cutoff=s'   => \$evalue_cutoff,
	'window_size=s'     => \$window_size,
	'help!'            => \$help
);

if ( defined($evalue_cutoff) ) {	
    if ( $evalue_cutoff !~ /^[0-9.]+$/ ) {
	    pod2usage(
	        -message => "Error: Parameter --evalue_cutoff should be a number !!!!\n\n",
	        -verbose => 1,
	        -exitval => 1,
	        -output  => \*STDERR
	    );
    }
}else{
	$evalue_cutoff = 0.00001;
}


if ( defined($window_size) ) {
    if ( $window_size !~ /^[0-9]+$/ || $window_size < 1000) {
        pod2usage(
            -message => "Error: Parameter --window_size should be an integer value higher or equal to 1000  !!!!\n\n",
            -verbose => 1,
            -exitval => 1,
            -output  => \*STDERR
        );
    }
}else{
    $window_size = 10000;
}


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

if ( not defined($gff_in) ) {
	pod2usage(
		-message => "Error: Parameter --gff_in is required !!!!\n\n",
		-verbose => 1,
		-exitval => 1,
		-output  => \*STDERR
	);
}

if ( not defined($outdir) ) {
	pod2usage(
		-message => "Error: Parameter --outdir is required !!!!\n\n",
		-verbose => 1,
		-exitval => 1,
		-output  => \*STDERR
	);
}

# Checking if users has provided
my $separate_gff_trna = 0;
if ( defined($gff_trna) ) {
	$separate_gff_trna = 1;
}

my $debug                  = 0;
my $blast_grid_output_directory = "$outdir/blast_grid";

mkdir( "$outdir", 0755 );

# Number of jobs issued when using the grid
my $number_of_jobs = 20;


# Setting path of some important directories
my $UTILS_DIR = "$FindBin::Bin/UTILS.dir";
my $CONFIG_DIR = "$FindBin::Bin/config.dir";
my $OBA_DIR  = "/cil/shed/apps/internal/OBA/GridSubmissions";


# Retrieving from configuration file path to auxiliar applications
my $EMBOSS_EXTRACTSEQ_PATH =
`grep 'Emboss_extractseq_path' $CONFIG_DIR/Third_party_programs_paths.log | awk '{printf \$2}'`;

my $BLAST_PATH =
`grep 'Blastall_path' $CONFIG_DIR/Third_party_programs_paths.log | awk '{printf \$2}'`;

my $BEDTOOLS_PATH =
`grep 'Bedtools_path' $CONFIG_DIR/Third_party_programs_paths.log | awk '{printf \$2}'`;


# Database with prophage proteins
my $PROPHET_DB_DIR = "$FindBin::Bin/PhrophET_phage_proteins_database.dir";

my $usage = "ProphET_standalone.pl <fasta file> <gff file> <output_name> ";


#########
#Processing the input files and separating in one fasta per GFF

# Get the scaffold IDs from the gff


my $gff_handler = GFFFile::new($gff_in);
$gff_handler->read();
my @scaffold_ids = $gff_handler->get_chrom_names();

#################
# Print scaffolds in the GFF

print "\nProcessing the following scaffolds/chromosomes:\n";
map { chomp $_ } @scaffold_ids;
map { print STDERR "$_\n" } @scaffold_ids;

`perl $UTILS_DIR/fasta2line $fasta_in > $outdir/fasta.line`;

print STDERR "\n";

# Array storing BLAST cmds and output files
my @cmds;
my %output_files;    # One set per each scaffold

# Iterate through each scaffold/chromosome
foreach my $scaff_chrom (@scaffold_ids) {

	print STDERR "Processing scaffold/chromosome: $scaff_chrom ...\n";

	my $intermediate_files_dir = "$outdir/$scaff_chrom";
	mkdir( "$intermediate_files_dir", 0755 );

	my $curr_fasta = "$intermediate_files_dir/$scaff_chrom.fasta";
	my $curr_gff   = "$intermediate_files_dir/$scaff_chrom.gff";
	

	# Generate a GFF and FASTA per each scaffold/chromosome
	`grep -v "^#" $gff_in | awk 'id == \$1 {print \$0}' id=$scaff_chrom $gff_in > $curr_gff`;
	`awk 'id == \$2 {print ">"\$2"\\n"\$1}' id=$scaff_chrom  $outdir/fasta.line > $curr_fasta`;

	##########
	#Checking the input files
	#Checking all the gffs IDs matches the fasta sequence

	my $num_seqs = `grep -c '>' $curr_fasta | awk '{print \$1}'`;
	die "ERROR: The file $curr_fasta has either more than one sequence or no sequence." if ( $num_seqs != 1 );
		

	my $gff_ids_count =
	  `awk '{print \$1}' $curr_gff | sort -u | wc | awk '{print \$1}'`;
	my $gff_ids_ids = `awk '{print \$1}' $curr_gff | sort -u`;
	if ( $gff_ids_count > 1 ) {
		die
			"ERROR: The input gff file has more than one genome id:\n" . 
			"$gff_ids_ids Check if you are not submitting a file with plasmids as well as the genome file\n";
	}


	my $fasta_only_id = `grep '>' $curr_fasta | sed s'/>//'`;
	my $gff_only_id   = `awk '{print \$1}' $curr_gff | sort -u`;
	if ( $fasta_only_id ne $gff_only_id ) {
		die
			"ERROR: The FASTA headers do not match the GFF sequence id:\n FASTA:$fasta_only_id GFF_SEQ_ID:$gff_only_id";
	}
	
	
	#Script to generate the proteins fasta based on the genome fasta and gff
	print STDERR "Generating file containing protein and gene sequence...\n";
`$UTILS_DIR/./gff2gene_protein_seq.pl $curr_gff $curr_fasta 11 $intermediate_files_dir/$scaff_chrom.trans $intermediate_files_dir/$scaff_chrom.cds $intermediate_files_dir/$scaff_chrom.prot`;


	# BLAST predicted protein against our phage db
	if( $grid ){
		split_blast(
			$scaff_chrom,
"$BLAST_PATH -p blastp -d $PROPHET_DB_DIR/Phage_proteins_without_ABC-t.db -e $evalue_cutoff -m8 -a8",
			"$intermediate_files_dir/$scaff_chrom.prot",
			\$number_of_jobs,
			\@cmds,
			\@{ $output_files{$scaff_chrom} },
			$outdir,
			$blast_grid_output_directory
		);

	# BLAST locally
	}else{
		print STDERR
		  "BLASting protein sequences against phage proteins db (e-value <= $evalue_cutoff) ...\n";
		`$BLAST_PATH -p blastp -d $PROPHET_DB_DIR/Phage_proteins_without_ABC-t.db -i $intermediate_files_dir/$scaff_chrom.prot -e $evalue_cutoff -m8 -a8 -o $intermediate_files_dir/$scaff_chrom.blast`;	
	}
}

`rm $outdir/fasta.line`;

if( $grid ){
	
	# Issue BLAST jobs in the grid and wait...
	execute_blast_on_grid( $outdir );
	
	# Merge BLAST results
	foreach my $scaff_chrom (@scaffold_ids){
		my $intermediate_files_dir = "$outdir/$scaff_chrom";
		my $curr_blast = "$intermediate_files_dir/$scaff_chrom.blast";

		if ($debug) {
			foreach my $file ( @{ $output_files{$scaff_chrom} } ) {
				print "Files to merge: $file\n" if $debug;
			}
		}
		merge_blast( \@{ $output_files{$scaff_chrom} }, $curr_blast );
	}
}	

foreach my $scaff_chrom (@scaffold_ids) {

	my $intermediate_files_dir = "$outdir/$scaff_chrom";
	
	my $curr_fasta = "$intermediate_files_dir/$scaff_chrom.fasta";
	my $curr_gff   = "$intermediate_files_dir/$scaff_chrom.gff";
	my $curr_blast = "$intermediate_files_dir/$scaff_chrom.blast";
	my $curr_blast_best_matches = "$intermediate_files_dir/$scaff_chrom.blastt-best-matches";
	my $curr_blast_best_matches_ids = "$intermediate_files_dir/$scaff_chrom.blastt-best-matches-ids";
	my $curr_blast_matches_coords = "$intermediate_files_dir/$scaff_chrom.matches-coords";
	my $curr_blast_union = "$intermediate_files_dir/$scaff_chrom.matches-coords-union";
	my $tRNAfile = "$intermediate_files_dir/$scaff_chrom.tRNA_file_coordinates";
	
		
	# Skip scaffold if BLAST result is empty
	next if ( -z "$curr_blast" );

	#Extracting only BLAST best matches
	print STDERR "Parsing results...\n";
`cat $curr_blast | awk '{OFS="\t";print \$0, \$1 }' \$1 | sort -k 1,1 -k 12,12rn | uniq -f 12 | awk '{OFS="\t";print \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8, \$9, \$10, \$11, \$12 }' > $curr_blast_best_matches`;

	#Extracting the gene ID
`cat $curr_blast_best_matches | awk '{print \$1}' > $curr_blast_best_matches_ids`;

	#Extracting the gene ID coordinates
`grep -wf $curr_blast_best_matches_ids $curr_gff | grep 'CDS' | awk '{print \$1"\t"\$4"\t"\$5}' | sort -nk 3,3 > $curr_blast_matches_coords`;

	#Union comand - There are some genes with overlap, so we use union so we do not count some positions twice
	print STDERR "Collapsing results...\n";
`$UTILS_DIR/./union.pl --in $curr_blast_matches_coords --seg_name 1 --seg_start 2 --seg_end 3 > $curr_blast_union`;

	#Extraction tRNA coordinates from the GFF
	print STDERR "Extracting tRNA records...\n";
	my $gff_input_selector;

	if ($separate_gff_trna) {
		$gff_input_selector = $gff_trna;
	}
	else {
		$gff_input_selector = $curr_gff;
	}
	`awk '\$3=="tRNA"{print \$0}' $gff_input_selector | awk '{print \$1"\t"\$4"\t"\$5}' > $tRNAfile`;

	

	my $blast_window_output = compute_density_sliding_windows( $scaff_chrom, $curr_blast, $curr_blast_union, $intermediate_files_dir, $window_size );

	my $blast_merged_final = group_consecutive_windows( $scaff_chrom, $blast_window_output, $intermediate_files_dir, $window_size );
	
	my $blast_merged_united = trim_phages( $scaff_chrom, $curr_blast_union, $blast_merged_final, $intermediate_files_dir );
	
	my ($ref_arr_phages_found , $phages_curr_scaffold ) = 
		trim_phages_by_trna_or_coding_gene( $scaff_chrom, $curr_blast_union, $blast_merged_united, $tRNAfile, $intermediate_files_dir, $outdir );
	
	render_phage( $scaff_chrom, $ref_arr_phages_found, $curr_fasta, $curr_gff, $phages_curr_scaffold, $curr_blast_union, $outdir );
	
	slice_out_phage_seq( $scaff_chrom, $ref_arr_phages_found, $curr_fasta, $outdir ); 
	
}

# Copy phage_db stats to results dir
`cp $PROPHET_DB_DIR/phage_db.summary.stats $outdir/phage_db.summary.stats`;


exit(0);

#########
# Sliding window and phage predictions:
#-Input: Coordinates with phage matches
#-Sliding window and intersect with phage matches - Windows of 10,000 with 1,000 increment
#-Estimates the Density of nucleotides with phage matches in a given window
#-Output: Phage matches for each window
	
sub compute_density_sliding_windows{
	
	my ( $scaff_chrom, $blast, $blast_union, $intermediate_files_dir, $window_size  ) = @_;

    print STDERR "Computing density of prophage genes for each sliding window (window size = $window_size bp)...\n";
	
	my $blast_window_output = "$intermediate_files_dir/$scaff_chrom.blast.window.output";
	my $blast_log_matches = "$intermediate_files_dir/$scaff_chrom.blast.log.matches";
	my $blast_window = "$intermediate_files_dir/$scaff_chrom.blast.window";
	my $blast_intersect_entre_phage_janela = "$intermediate_files_dir/$scaff_chrom.entre.phage.janela";
	
	
	open( BLAST, "$blast_union" )
	  or die "Unable to open file $blast_union!\n";
	my @genes = <BLAST>;
	close(BLAST);

	open( RESULTS, ">$blast_window_output" )
	  or die "Unable to write on file $blast_window_output!\n";
	;    #File with the number of bases with phage matches in a given window
	
	open( LOG, ">$blast_log_matches" )
	  or die "Unable to write on file $blast_log_matches!\n";
	;    #Log-file with the coordinates of matches in a given window

	my @last_pos =
	  split( /[\t]/, $genes[$#genes] )
	  ; #As we use a $window_size window, we stop the analysis $window_size after the last match with a phage gene
	my $id_genoma  = $last_pos[0];
	my $final_pos  = $last_pos[2];
	my $final_real = $final_pos + $window_size;
	chomp $id_genoma;

	for ( my $i = $window_size ; $i <= $final_real ; $i = $i + 1000 )
	{    #Sliding window to intersect with the phage matches coordinates
		open( WINDOW, ">$blast_window" )
		  or die "Unable to write on file $blast_window!\n"
		  ;    # Window overwrites at each iteration
		my $coord_ini = $i - ( $window_size - 1 );
		chomp $coord_ini;
		chomp $i;
		print WINDOW "$id_genoma\t$coord_ini\t$i\n";
		close(WINDOW);

		`$BEDTOOLS_PATH intersect -a $blast_union -b $blast_window > $blast_intersect_entre_phage_janela`
		  ; # Makes the intersect between the current window and our phage-matches coordinates

		open( INTERSECTED, "$blast_intersect_entre_phage_janela" )
		  or die "Unable to open file $blast_intersect_entre_phage_janela!\n";
		;    #Open the intersect file
		my @merged = <INTERSECTED>;
		close(INTERSECTED);

		print LOG "#Window\t$coord_ini\t$i\n"
		  ; #Stores the coordinates with phage matches in a given window in the log file

		my $gene_content = 0;
		for ( my $j = 0 ; $j <= $#merged ; $j++ )
		{    #Runs trough all the genes with phage matches in the current window
			my @splitted = split( /[\t]/, $merged[$j] );
			my $initial  = $splitted[1];
			my $final    = $splitted[2];
			my $size =
			  ( $final + 1 ) -
			  $initial
			  ; #Stores the size of the gene with phage matches (could be more than one per window)
			$gene_content =
			  $gene_content + $size
			  ; #Stores the number whole extent of genes with phage matches in a give window
			chomp $initial;
			chomp $final;
			print LOG "$id_genoma\t$initial\t$final\n"
			  ; #Log with the coordinates of the genes with phage matches in each window
		}

		print LOG
		  "--------------------------------\n";    #Inter-window log separator
		print RESULTS "$id_genoma\t$coord_ini\t$i\t$gene_content\n"
		  ; #Tem output file - containing the number of phage-related nucleotides in each window

	}

	close(LOG);
	close(RESULTS);	
	return $blast_window_output;
}



##########
# Group consecutive windows with phage matches
#-Input: Phage matches for each window
#-Group consecutive windows with phage matches given the phage content in a given window is higher than half of the window (5000b), generatin raw-clusters
#-Then, group raw-clusters with overlaps
#-Output - Raw phage coordinates prediction

sub group_consecutive_windows {
	
	print STDERR "Grouping consecuting windows containing the putative prophage...\n";
	
	
	my ( $scaff_chrom, $blast_window_output, $intermediate_files_dir, $window_size) = @_;
	
	my $half_window = int( $window_size / 2 );

	my $blast_merged = "$intermediate_files_dir/$scaff_chrom.blast.merged";
	my $blast_merged_final = "$intermediate_files_dir/$scaff_chrom.blast.merged.final";

	open( INPUT, "$blast_window_output" )
	  or die "Unable to open file $blast_window_output!\n";
	my @matches = <INPUT>;
	close(INPUT);

	open( TEMP, ">$blast_merged" )
	  or die "Unable to write on file $blast_merged!\n";

	for ( my $k = 0 ; $k <= $#matches ; $k++ )
	{ #Loop to cluster consectutive windows with "phage content" higher than $half_window, half of the window
		my @splitted        = split( /[\t]/, $matches[$k] );
		my $coverage        = $splitted[3];
		my $ini_coordinates = $splitted[1];
		my $end_coordinates = $splitted[2];
		my $id              = $splitted[0];

		if ( $coverage > $half_window ) {
			while ( $coverage > $half_window && $k < $#matches )
			{    #while coverage>$half_window, group consecutive windows
				my @splitted2 = split( /[\t]/, $matches[$k] );
				$coverage = $splitted2[3];
				$k++;
			}

			my @splitted3 =
			  split( /[\t]/, $matches[ $k - 2 ] )
			  ; #Get the initial, final coordinates and size of the clustered windows
			my $end_coordinates3 = $splitted3[2];
			my $size =
			  ( $end_coordinates3 + 1 ) -
			  $ini_coordinates;    #### Nao tenho usado o size para nada...
			chomp $ini_coordinates;
			chomp $end_coordinates3;
			chomp $id;
			chomp $size;

			print TEMP "$id\t$ini_coordinates\t$end_coordinates3\n";
		}
	}
	close(TEMP);
	
	#Merge regions raw-clusters with overlap in sequences
	`$UTILS_DIR/union.pl --in $blast_merged --seg_name 1 --seg_start 2 --seg_end 3 > $blast_merged_final`;

	return $blast_merged_final;
}

##########
# Trimming raw phages with less than 8 genes with phage matches
#-Input: Raw phage coordinates prediction
#-Check in the raw phage has at least 8 genes with phage matches
#-Output: raw phages that have at least 8 genes with phage matches

sub trim_phages {
	
	print STDERR "Trimming prophage...\n";
	
	my ($scaff_chrom, $blast_union, $blast_merged_final, $intermediate_files_dir ) = @_;
	
	my $blast_merged_united = "$intermediate_files_dir/$scaff_chrom.blast.merged.united";

	
	open( MERGED, "$blast_merged_final" ) or die "couldnt open $blast_merged_final\n";
	my @merged = <MERGED>;
	close(MERGED);

	open( BLAST, "$blast_union" ) or die "couldnt open $blast_union";
	my @blast = <BLAST>;
	close(BLAST);

	open( TEMP2, ">$blast_merged_united" )
	  or die "Unable to write on file $blast_merged_united!\n";

	for ( my $l = 0 ; $l <= $#merged ; $l++ ) {
		my @P_positions = split( /[\t]/, $merged[$l] );
		my $P_start     = $P_positions[1];
		my $P_end       = $P_positions[2];
		my $contador    = 0;

		for ( my $m = 0 ; $m <= $#blast ; $m++ ) {
			my @B_positions = split( /[\t]/, $blast[$m] );
			my $B_start     = $B_positions[1];
			my $B_end       = $B_positions[2];

			if ( $B_start > $P_start && $B_end < $P_end ) {
				$contador++;
			}
			if ( $contador == 8 ) {
				print TEMP2 "$merged[$l]";
				last;
			}
		}
	}
	close(TEMP2);
	
	return $blast_merged_united;
}


##########
# Trimming raw phages borders by tRNAs or last gene with a phage match
#-Input: Phage coordinates
#-Trims the phage border to the last gene with and then searches 3kb upstream and downstream for tRNA genes and extends the coordinates accordingly
#-Output: Polished final phage prediction

sub trim_phages_by_trna_or_coding_gene {
	
	print STDERR "Trimming prophage based on tRNA or last gene in the last and first window ...\n";
	

	my ($scaff_chrom, $blast_union, $blast_merged_united, $tRNAfile, $intermediate_files_dir, $outdir ) = @_;
	
	my $trna_log = "$intermediate_files_dir/$scaff_chrom.trna.log";
	my $phages_coord = "$outdir/phages_coords";
	my $phages_curr_scaffold = "$intermediate_files_dir/$scaff_chrom.phages_coords";


	open( PHAGES, "$blast_merged_united" )
	  or die "Unable to open file $blast_merged_united!\n";
	my @phages = <PHAGES>;
	close(PHAGES);
	map( chomp, @phages); 

	open( MATCHES, "$blast_union" )
	  or die "Unable to open file $blast_union!\n";
	my @matches = <MATCHES>;
	close(MATCHES);
	map( chomp, @matches);

	open( TRNA, "$tRNAfile" ) or die "Unable to open file $tRNAfile!\n";
	my @tRNA = <TRNA>;
	close(TRNA);
	map( chomp, @tRNA);

	open( OUTLOG, ">$trna_log" )
	  or die "Unable to write on file $trna_log!\n";

	open( PHAGESFINAL, ">>$phages_coord" )
	  or die "Unable to write on file $phages_coord!\n";

	open( PHAGES_CURR_SCAFFOLD, ">$phages_curr_scaffold" )
	  or die "Unable to write on file $phages_curr_scaffold!\n";

	my @phages_found;

	for ( my $n = 0 ; $n <= $#phages ; $n++ )
	{ #Searches for the last gene with phage match in the beggining and end of the raw phage
		my @phage_coords  = split( /[\t]/, $phages[$n] );
		
		my $phage_scaff_chrom = $phage_coords[0];
		my $phage_id      = $n + 1;
		my $phage_initial = $phage_coords[1];
		my $phage_final   = $phage_coords[2];

		my @inside_initial = ();
		my @inside_final   = ();

		for ( my $p = 0 ; $p <= $#matches ; $p++ ) {
			my @matches_coords  = split( /[\t]/, $matches[$p] );
			my $matches_id      = $matches_coords[0];
			my $matches_initial = $matches_coords[1];
			my $matches_final   = $matches_coords[2];

			if (   $matches_final > $phage_initial
				&& $matches_final < $phage_final )
			{
				push( @inside_initial, $matches_coords[1] );
			}

			if (   $matches_initial > $phage_initial
				&& $matches_initial < $phage_final )
			{
				push( @inside_final, $matches_coords[2] );
			}
		}

		$phage_initial = $inside_initial[0];
		my $number = $#inside_final;
		$phage_final = $inside_final[$number];
		my $phage_begining = $phage_initial;    # Para poder alterar com o tRNA
		my $phage_ending   = $phage_final;      # Para poder alterar com o tRNA

		my $inital_minus = $phage_initial - 3000;
		my $initial_plus = $phage_initial + 3000;
		my $final_minus  = $phage_final - 3000;
		my $final_plus   = $phage_final + 3000;

		my $distancia  = 3000;    #distancia inicial para permitir o tRNA
		my $distancia2 = 3000;
		my $diff       = 3000;    #diferenca para mudar com o tRNA
		my $diff2      = 3000;

		for ( my $o = 0 ; $o <= $#tRNA ; $o++ )
		{                         # Checks for tRNA genes near phage borders
			my @tRNA_coords  = split( /[\t]/, $tRNA[$o] );
			my $trna_id      = $tRNA_coords[0];
			my $trna_initial = $tRNA_coords[1];
			my $trna_final   = $tRNA_coords[2];

			if (   $trna_initial > $inital_minus
				&& $trna_final < $initial_plus )
			{
				if (   $trna_initial < $phage_initial
					&& $trna_final < $phage_initial )
				{
					$diff = $phage_initial - $trna_final;

					if ( $distancia > $diff ) {
						$distancia = $diff;

						$phage_begining = $trna_initial;
					}
				}
				elsif ($trna_initial > $phage_initial
					&& $trna_final > $phage_initial )
				{
					$diff = $trna_initial - $phage_initial;
					if ( $distancia > $diff ) {
						$distancia      = $diff;
						$phage_begining = $trna_initial;
					}
				}

				elsif ($trna_initial < $phage_initial
					&& $trna_final > $phage_initial )
				{
					$diff           = 0;
					$phage_begining = $trna_initial;
				}

				print OUTLOG $tRNA[$o];
			}

			elsif ($trna_initial > $final_minus
				&& $trna_final < $final_plus )
			{
				if (   $trna_initial < $phage_final
					&& $trna_final < $phage_final )
				{
					$diff2 = $phage_final - $trna_final;
					if ( $distancia2 > $diff2 ) {
						$distancia2   = $diff2;
						$phage_ending = $trna_final;
					}
				}

				elsif ($trna_initial > $phage_final
					&& $trna_final > $phage_final )
				{
					$diff2 = $trna_initial - $phage_final;
					if ( $distancia2 > $diff2 ) {
						$distancia2   = $diff2;
						$phage_ending = $trna_final;
					}
				}

				elsif ($trna_initial < $phage_final
					&& $trna_final > $phage_final )
				{
					$diff2        = 0;
					$phage_ending = $trna_final;
				}

				print OUTLOG $tRNA[$o];
			}

		}

		$phage_initial = $phage_begining;
		$phage_final   = $phage_ending;

		#chomp $phage_id;
		#chomp $phage_begining;
		#chomp $phage_ending;

		push @phages_found,
		  {
		  	scaff_chrom => $phage_scaff_chrom,
			id    => $phage_id,
			start => $phage_begining,
			end   => $phage_ending
		  };
		print PHAGESFINAL "$phage_scaff_chrom\t$phage_id\t$phage_begining\t$phage_ending\n";
		print PHAGES_CURR_SCAFFOLD
		  "$phage_id\t$phage_begining\t$phage_ending\n";

	}

	close(PHAGESFINAL);
	close(PHAGES_CURR_SCAFFOLD);
	
	return (\@phages_found, $phages_curr_scaffold);
}

#########
# Generate the image file 
#-Input: Final phage prediction, genome fasta, gff file
#-Generates the image output file
#-Output: Image output file

sub render_phage {
	
	print STDERR "Rendering graph depicting prophage genes ...\n";
	
	
	my ( $scaff_chrom, $ref_arr_phages_found, $fasta, $gff, $phages_coords, $blast_union, $intermediate_files_dir ) = @_; 
	
	my $gff_ultraformated = "$intermediate_files_dir/$scaff_chrom.gff_ultraformated";
	my $svg = "$intermediate_files_dir/$scaff_chrom.svg";


	if ( scalar(@{$ref_arr_phages_found}) != 0 ) {

		`cat $gff | grep 'CDS' | awk '{print \$1"\t"\$4"\t"\$5}' > $gff_ultraformated`;
		`$UTILS_DIR/./gff2graph-from-scratch.pl $gff_ultraformated $fasta $svg $blast_union $phages_coords`;
		`rm -r $gff_ultraformated`;

		#open(FINALPROGRAM, ">$name-run-sucessfuly-completed.log");
		#print FINALPROGRAM "The Run was sucessfull";
		#close(FINALPROGRAM);

	}
}

sub slice_out_phage_seq{

	my ( $scaff_chrom, $ref_arr_phages_found, $fasta, $outdir ) = @_; 
	
	if ( scalar(@{$ref_arr_phages_found}) != 0 ) {
		foreach my $curr_phage (@{$ref_arr_phages_found}) {
			my $scaff_chrom = $curr_phage->{scaff_chrom};
			my $id = $curr_phage->{id};
			my $start = $curr_phage->{start};
			my $end = $curr_phage->{end};

			my $cmd = "$EMBOSS_EXTRACTSEQ_PATH -sequence $fasta -regions \"$start-$end\" -osdbname2 phage_$id:$start-$end $outdir/$scaff_chrom.phage_$id.fas";
			#print STDERR "$cmd\n";			 
			`$cmd`;
		}
	}
}


#############
# Executing BLAST cmds if using the grid
sub execute_blast_on_grid {
	my ( $outdir ) = @_;

	print STDERR "BLASting protein sequences against phage proteins db...\n";

	open BLAST_CMDS, ">$outdir/blast.cmds"
	  or die "Unable to write on file $outdir/blast.cmds!\n";

	foreach my $cmd (@cmds) {
		print "Writing on file BLAST command: $cmd" if $debug;
		print BLAST_CMDS $cmd;
	}

	close(BLAST_CMDS);

	my $cmd = "$OBA_DIR/run_cmds_on_grid.py $outdir/blast.cmds";
	print "$cmd\n" if $debug;

	my $output = `$cmd`;
	print $output if $debug;

	if ( $output !~ "No failed commands" ) {
		die
"\nERROR: Unable to complete succesfully all BLAST commands!! Please check for errors in the UGER log files\n\n";
	}
	else {
		print "No failed commands.\n";
	}
}


sub split_blast {
	my (
		$batch_name,       $blast_cmd_base, $input_file,
		$ref_num_of_parts, $refArrCmds,     $refArrOutFiles,
		$outdir, $blastDirectory
	) = @_;

	print "Spliting FASTA file $input_file...\n";

	my $cmd_split_fasta =
"$UTILS_DIR/splitfasta.pl --fasta_in $input_file --num_files $$ref_num_of_parts --pre $batch_name --out_dir $blastDirectory";

	print $cmd_split_fasta . "\n" if $debug;

	my $output = `$cmd_split_fasta`;

	my @cols = split " ", $output;
	my $num_files_created = pop @cols;

	print STDERR "Number of files created: $num_files_created\n" if $debug;

	# Adjusting the value of created files
	$$ref_num_of_parts = $num_files_created;

	print STDERR $output . "\n" if $debug;

	for (
		my $file_index = 1 ;
		$file_index <= $$ref_num_of_parts ;
		$file_index++
	  )
	{
		my $outputFile =
		  $blastDirectory . "/" . "$batch_name.$file_index.blast";
		my $inputFile = $blastDirectory . "/" . "$batch_name.$file_index.fas";
		my $cmd       = $blast_cmd_base . " -i $inputFile -o $outputFile\n";
		push( @$refArrCmds,     $cmd );
		push( @$refArrOutFiles, $outputFile );
		print "Generating BLAST cmd: $cmd" if $debug;
	}
}

sub merge_blast {
	print "Merging BLAST files...\n";

	my ( $refArrOutFiles, $output_file ) = @_;

	my $all_files_to_combine = join " ", @{$refArrOutFiles};
	my $cmd = "cat $all_files_to_combine > $output_file";
	print $cmd . "\n" if $debug;
	`$cmd`;
}
