#!/usr/bin/env perl
use strict;

use Pod::Usage;
use Getopt::Long;
use File::Path;
use File::Copy;

=head1 NAME

INSTALL.pl

=head1 SYNOPSIS

INSTALL.pl
	[  --update_db_only --reparse ] |

=head1 OPTIONS

B<--update_db_only> - Only update the database of known prophages B<(Optional)>.

B<--reparse> - Only reformat the database of known prophages B<(Optional)>.

B<--help> - prints the usage information. B<(Optional)>

=head1 DESCRIPTION

=head1 CONTACT
 Gustavo C. Cerqueira (2018)
 cerca11@gmail.com
 gcerqueira@pgdx.com
=cut

my ($update_db_only, $phage_families_file, $reparse);
my $help;
my $temp = "ProphET_install_temp.dir";
my $database_dir = "ProphET_phage_proteins_database.dir";
my $config_dir = "../config.dir";
my $default_phage_families_file = "Prophages_names_sem_Claviviridae_Guttaviridae-TxID";


GetOptions(
    'update_db_only'	=> \$update_db_only,
    'phage_families_file=s' => \$phage_families_file,
    'reparse' => \$reparse,
    'help'			=> \$help
);

if( defined($help) ){
   pod2usage(-verbose => 1 ,-exitval => 2);
}


goto DOWNLOADING_DB if defined( $update_db_only );

goto PARSING_DB if defined( $reparse );

#-----------------------------------------

print "Looking for required programs in the enviroment PATH...\n";
my $config_file = "./config.dir/Third_party_programs_paths.log";

open(LOGS, ">$config_file" ) or die "Unable to write on config file: $config_file\n";

my $emboss_extractseq = `which extractseq`;
die "\nERROR: Unable to find \"extractseq\", EMBOSS suite\n" if( $emboss_extractseq eq '' );
chomp $emboss_extractseq;
print "\tFound EMBOSS extractseq: $emboss_extractseq\n";

my $blastp = `which blastp`;
die "\nERROR: Unable to find \"blastp\", BLAST suite\n\n" if( $blastp eq '' );
chomp $blastp;
print "\tFound blastp: $blastp\n";

my $formatdb = `which makeblastdb`;
die "\nERROR: Unable to find \"makeblastdb\", BLAST suite\n\n" if( $formatdb eq '' );
chomp $formatdb;
print "\tFound makeblastdb: $formatdb\n";

my $bedtools = `which bedtools`;
die "\nERROR: Unable to find \"bedtools\"\n\n" if( $bedtools eq '' );
chomp $bedtools;
print "\tFound bedtools: $bedtools\n";

#-----------------------------------------
print "Saving program paths in $config_file ...\n";
print LOGS "Emboss_extractseq_path\t$emboss_extractseq\n";
print LOGS "Blastp_path\t$blastp\n";
print LOGS "Formatdb_path\t$formatdb\n";
print LOGS "Bedtools_path\t$bedtools\n";
close(LOGS);


#-----------------------------------------
print "Looking for required Perl libraries...\n";

my $output = system("perl -e 'use Bio::Perl;'");
die "\nERROR: Unable to find Perl module Bio::Perl\n\n" if( $output );

$output = system("perl -e 'use LWP::Simple;'");
die "\nERROR: Unable to find Perl module LWP::Simple\n\n" if( $output );

$output = system("perl -e 'use XML::Simple;'");
die "\nERROR: Unable to find Perl module XML::Simple\n\n" if( $output );

$output = system("perl -e 'use GD;'");
die "\nERROR: Unable to find Perl module GD\n\n" if( $output );


#-----------------------------------------
print "Downloading GFFLib ...\n";
#$output = system("svn --force export https://github.com/gustavo11/GFFLib/trunk UTILS.dir/GFFLib");
$output = system("rm -rf UTILS.dir/GFFLib");
$output = system("git clone https://github.com/gustavo11/GFFLib.git UTILS.dir/GFFLib");
die "ERROR: Unable to download GFFLib from github\n\n" if( $output );



DOWNLOADING_DB:

#-----------------------------------------
print "Creating database directory...\n";

if( -e $database_dir ){
	my $datestring = localtime();
	$datestring =~ s/ /_/g;
	my $src = $database_dir;
	my $dst = "$database_dir.$datestring.bak";

    move( $src, $dst  )
        || die("ERROR: Unable to move directory $src to $dst!");
}

mkdir($database_dir, 0755) or
die "ERROR: Unable to create directory $database_dir\n";


#-----------------------------------------
print "Creating database temp directory ...\n";

if( -e $temp ){
	rmtree( $temp )
		|| die("ERROR: Unable to remove directory $temp!");
}

mkdir($temp, 0755) or
die "ERROR: Unable to create directory $temp.\n";

chdir "$temp" or
	die "ERROR: Unable to enter directory $temp\n";

#-----------------------------------------

print "Downloading Phage sequences ...\n";

my $eff_phage_families_file;
if ( defined( $phage_families_file) ){
	$eff_phage_families_file = $config_dir . "/" . $phage_families_file;
}else{
	$eff_phage_families_file = $config_dir . "/" . $default_phage_families_file;
}

if ( ! ( -e $eff_phage_families_file ) ){
	die "ERROR: Phage families file $eff_phage_families_file does not exist!";
}

$output = system("../UTILS.dir/extrair_ncbi_prophage_families.pl $eff_phage_families_file");
die "ERROR: Unable to execute extrair_ncbi_prophage_families.pl\n\n" if( $output );

`perl ../UTILS.dir/obtain_prot_with_annot_seq.pl $eff_phage_families_file > Phage_proteins_pre_raw.db`;

PARSING_DB:

if (defined( $reparse )){
    chdir "$temp" or
	die "ERROR: Unable to enter directory $temp\n";
}
#-----------------------------------------
print "Formating sequences ...\n";
`sed s'/[*]//g' Phage_proteins_pre_raw.db > Phage_proteins_pre_raw_without_stop.db `; # Remove asterisks representing STOP codons
`perl ../UTILS.dir/script_remover_vazios.pl Phage_proteins_pre_raw_without_stop.db > Phage_proteins_raw_pregrep.db`;
print "checking for misformated sequences ...\n";
# remove problematic salmonella items by search for lines not starting with > and not starting with a amino acid code
# Note that this clears out any equences starting with a "-"
# I dont think that is a problem, but I may be thinking of anything
`cat Phage_proteins_raw_pregrep.db | grep -B 1 "^[^>ABCDEFGHIJKLMNOPQRSTUVWXYZ]" > bad_lines`;
# use the bad lines file (-f) literally (-F)  to get inverse (-v)
# if the file is non, empty, that is
`if [ -s bad_lines ] ;then grep -v -F -f bad_lines Phage_proteins_raw_pregrep.db > Phage_proteins_raw.db ; else cp Phage_proteins_raw_pregrep.db Phage_proteins_raw.db ; fi`;


print "Formatting blast database ...\n";

`makeblastdb -dbtype prot -in  Phage_proteins_raw.db`;
if ($? == -1) {
    die "ERROR: Unable to execute formatdb!\n";
}

`../UTILS.dir/fasta2line Phage_proteins_raw.db > Phage_proteins_raw.line`;


#-----------------------------------------
print "Removing ABC-Transporters ...\n";

# Retrieve ABC transporter from database based on their annotation
`grep -wf ../config.dir/ABC_transporters_to_grep.txt Phage_proteins_raw.line | sort -u | awk '{print ">"\$2"\\n"\$1}' > ABC_transporters_seqs.fasta`;

# system("export BLASTPATH=\$\(dirname `which blastp`)");
# print "\$BLASTPATH\n";
# BLAST those ABC transporters against the rest of the database
`blastp -db Phage_proteins_raw.db -query ABC_transporters_seqs.fasta -evalue 1e-5 -outfmt 6 -out ABC_trans_BLAST_matches`;
if ($? == -1) {
    die "ERROR: Unable to execute blastp!\n";
}

# Retrieve the ID of matches against ABC transporters
`cat ABC_trans_BLAST_matches | awk '{print \$2}' | sort -u > IDs_Matches_com_ABC_transporters`;

# Remove those matches from the database
if ( -z 	"IDs_Matches_com_ABC_transporters" && -z "ABC_transporters_seqs.fasta"){
	`awk '{print ">"\$2"\\n"\$1}' Phage_proteins_raw.line > Phage_proteins_without_ABC-t.db`
}elsif( -z 	"IDs_Matches_com_ABC_transporters" ){
	die "ERROR: Incomplete or corrupted BLAST search for ABC transporters!!\n";
}else{
	`grep -vf IDs_Matches_com_ABC_transporters Phage_proteins_raw.line | awk '{print ">"\$2"\\n"\$1}' > Phage_proteins_without_ABC-t.db`;
}

#-----------------------------------------
print "Finalizing phage database...\n";
`cp Phage_proteins_without_ABC-t.db ../$database_dir/Phage_proteins_without_ABC-t_pre.db`;
`cp phage_db.summary.stats ../$database_dir`;
chdir "../$database_dir";
print "doublechecking for empty/misformatted sequences\n";
`cat Phage_proteins_without_ABC-t_pre.db | grep -B 1 "^[^>ABCDEFGHIJKLMNOPQRSTUVWXYZ]" > bad_lines`;
`if [ -s bad_lines ] ;then grep -v -F -f bad_lines Phage_proteins_without_ABC-t_pre.db > Phage_proteins_without_ABC-t.db ; else cp Phage_proteins_without_ABC-t_pre.db Phage_proteins_without_ABC-t.db ; fi`;
# `mv tmp.txt Phage_proteins_raw.db`;
print "Formatting final phage database...\n";
`makeblastdb -dbtype prot -in  Phage_proteins_without_ABC-t.db`;


#-----------------------------------------
chdir "../";
#`rm -rf ProphET_install_temp.dir`;
print "Installation completed!\n";


exit(0);
