#!/bin/env perl
use strict;


open(LOGS, "./config.dir/Third_party_programs_paths.log");

my $emboss_extractseq = `which extractseq`; 
die "\nERROR: Unable to find \"extractseq\", EMBOSS suite\n\n" if( $emboss_extractseq eq '' );

my $blastall = `which blastall`;
die "\nERROR: Unable to find \"blastall\", BLAST suite\n\n" if( $blastall eq '' );

my $bedtools = `which bedtools`;
die "\nERROR: Unable to find \"bedtools\"\n\n" if( $bedtools eq '' );

chomp $emboss_extractseq;
chomp $blastall;
chomp $bedtools;

print LOGS "Emboss_extractseq_path\t$emboss_extractseq\n";
print LOGS "Blastall_path\t$blastall\n";
print LOGS "Bedtools_path\t$bedtools\n";
close(LOGS);


`svn --force export https://github.com/gustavo11/GFFLib/trunk UTILS.dir/GFFLib`;


print "Downloading Phage sequences\n--------------------------------------------------\n";

mkdir("ProphET_instal_temp.dir", 0755);
chdir "ProphET_instal_temp.dir";
`perl ../UTILS.dir/extrair_ncbi_prophage_families.pl ../config.dir/Prophages_names_sem_Claviviridae_Guttaviridae-TxID`;
`perl ../UTILS.dir/obtain_prot_with_annot_seq.pl ../config.dir/Prophages_names_sem_Claviviridae_Guttaviridae-TxID >Phage_proteins_pre_raw.db`;
print "Formating sequences\n--------------------------------------------------\n";
`sed s'/[*]//g' Phage_proteins_pre_raw.db > Phage_proteins_pre_raw_without_stop.db `; # Remove asterisks representing STOP codons
`perl ../UTILS.dir/script_remover_vazios.pl Phage_proteins_pre_raw_without_stop.db > Phage_proteins_raw.db`;
`../UTILS.dir/./fasta2line Phage_proteins_raw.db > Phage_proteins_raw.line`;

print "Removing ABC-Transporters\n--------------------------------------------------\n";
`grep -wf ../config.dir/ABC_transporters_to_grep.txt Phage_proteins_raw.line | sort -u | awk '{print ">"\$2"\\n"\$1}' > ABC_transporters_seqs.fasta`;
`formatdb -p T -i Phage_proteins_raw.db`;
`blastall -p blastp -d Phage_proteins_raw.db -i ABC_transporters_seqs.fasta -e 1e-5 -m8 -o ABC_trans_BLAST_matches`;
`cat ABC_trans_BLAST_matches | awk '{print \$2}' | sort -u > IDs_Matches_com_ABC_transporters`;
`grep -vf IDs_Matches_com_ABC_transporters Phage_proteins_raw.line | awk '{print ">"\$2"\\n"\$1}' > Phage_proteins_without_ABC-t.db`;

mkdir("../PhrophET_phage_proteins_database.dir", 0755);
`cp Phage_proteins_without_ABC-t.db ../PhrophET_phage_proteins_database.dir`;
chdir "../PhrophET_phage_proteins_database.dir";
`formatdb -p T -i Phage_proteins_without_ABC-t.db`;

chdir "../";
`rm -rf ProphET_instal_temp.dir`;
print "Database generation complete\n--------------------------------------------------\n";

