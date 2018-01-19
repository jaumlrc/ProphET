#!/usr/bin/env perl
use strict;

use Pod::Usage;
use Getopt::Long;


=head1 NAME

INSTALL.pl

=head1 SYNOPSIS

INSTALL.pl 
	[  --update_db_only ] | 

=head1 OPTIONS

B<--update_db_only> - Only update the database of known prophages B<(Optional)>.

B<--help> - prints the usage information. B<(Optional)>

=head1 DESCRIPTION

=head1 CONTACT
 Gustavo C. Cerqueira (2018)
 cerca11@gmail.com
 gcerqueira@pgdx.com
=cut

my ($update_db_only);
my $help;

GetOptions(	'update_db_only'	=> \$update_db_only,
			'help'			=> \$help );
		
if( defined($help) ){
   pod2usage(-verbose => 1 ,-exitval => 2);
} 
		
goto DOWNLOADING_DB if defined( $update_db_only );

#-----------------------------------------

print "Looking for required programs in the enviroment PATH...\n";
my $config_file = "./config.dir/Third_party_programs_paths.log";

open(LOGS, ">$config_file" ) or die "Unable to write on config file: $config_file\n";

my $emboss_extractseq = `which extractseq`; 
die "\nERROR: Unable to find \"extractseq\", EMBOSS suite\n" if( $emboss_extractseq eq '' );
chomp $emboss_extractseq;
print "\tFound EMBOSS extractseq: $emboss_extractseq\n";

my $blastall = `which blastall`;
die "\nERROR: Unable to find \"blastall\", BLAST suite\n\n" if( $blastall eq '' );
chomp $blastall;
print "\tFound blastall: $blastall\n";

my $formatdb = `which formatdb`;
die "\nERROR: Unable to find \"formatdb\", BLAST suite\n\n" if( $formatdb eq '' );
chomp $formatdb;
print "\tFound blastall: $formatdb\n";

my $bedtools = `which bedtools`;
die "\nERROR: Unable to find \"bedtools\"\n\n" if( $bedtools eq '' );
chomp $bedtools;
print "\tFound bedtools: $bedtools\n";

#-----------------------------------------
print "Saving program paths in $config_file ...\n";
print LOGS "Emboss_extractseq_path\t$emboss_extractseq\n";
print LOGS "Blastall_path\t$blastall\n";
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
$output = system("svn --force export https://github.com/gustavo11/GFFLib/trunk UTILS.dir/GFFLib");
die "ERROR: Unable to download GFFLib from github\n\n" if( $output );


#-----------------------------------------
print "Creating database directory...\n";

my $database_dir = "PhrophET_phage_proteins_database.dir";

if( !( -e $database_dir ) ){
	mkdir($database_dir, 0755) or 
	die "ERROR: Unable to create directory $database_dir\n";
}


#-----------------------------------------

DOWNLOADING_DB:

print "Downloading Phage sequences ...\n";
my $temp = "ProphET_install_temp.dir";

if( !( -e $temp ) ){
	mkdir($temp, 0755) or 
	die "ERROR: Unable to create directory $temp. If $temp already exists, please remove it.\n";
}
	
chdir "$temp" or 
	die "ERROR: Unable to enter directory $temp\n";
		
$output = system("../UTILS.dir/extrair_ncbi_prophage_families.pl ../config.dir/Prophages_names_sem_Claviviridae_Guttaviridae-TxID");
	
die "ERROR: Unable to execute extrair_ncbi_prophage_families.pl\n\n" if( $output );

`perl ../UTILS.dir/obtain_prot_with_annot_seq.pl ../config.dir/Prophages_names_sem_Claviviridae_Guttaviridae-TxID > Phage_proteins_pre_raw.db`;

#-----------------------------------------
print "Formating sequences ...\n";
`sed s'/[*]//g' Phage_proteins_pre_raw.db > Phage_proteins_pre_raw_without_stop.db `; # Remove asterisks representing STOP codons
`perl ../UTILS.dir/script_remover_vazios.pl Phage_proteins_pre_raw_without_stop.db > Phage_proteins_raw.db`;
`../UTILS.dir/fasta2line Phage_proteins_raw.db > Phage_proteins_raw.line`;


#-----------------------------------------
print "Removing ABC-Transporters ...\n";
`grep -wf ../config.dir/ABC_transporters_to_grep.txt Phage_proteins_raw.line | sort -u | awk '{print ">"\$2"\\n"\$1}' > ABC_transporters_seqs.fasta`;
`$formatdb -p T -i Phage_proteins_raw.db`;
`blastall -p blastp -d Phage_proteins_raw.db -i ABC_transporters_seqs.fasta -e 1e-5 -m8 -o ABC_trans_BLAST_matches`;
`cat ABC_trans_BLAST_matches | awk '{print \$2}' | sort -u > IDs_Matches_com_ABC_transporters`;
`grep -vf IDs_Matches_com_ABC_transporters Phage_proteins_raw.line | awk '{print ">"\$2"\\n"\$1}' > Phage_proteins_without_ABC-t.db`;

#-----------------------------------------


`cp Phage_proteins_without_ABC-t.db ../$database_dir`;
chdir "../$database_dir";
`formatdb -p T -i Phage_proteins_without_ABC-t.db`;

#-----------------------------------------
chdir "../";
`rm -rf ProphET_install_temp.dir`;
print "Installation completed!\n";

