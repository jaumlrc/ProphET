#!user/bin/perl
use strict;
#Obtem e junta todos os arquivos do proteoma de virus gerados pelo scrip do gustavo

my $input = $ARGV[0];     #Tabela contendo o nome das familias
open (INPUT, "$input");
my @N_seq = <INPUT>;
close(INPUT);
#open (RESULTS, ">Phage_proteins_pre_raw.db");


my @proteinseq = ();
for (my $j=0; $j<=$#N_seq; $j++) {
        my @temp1 = split (/\t/, $N_seq[$j]);
        my $name = $temp1[0];
	chdir "$name.dir";
	my @temp_seq = `cat *.fasta.prot`;
	push (@proteinseq, @temp_seq);
	chdir "../";
}

print @proteinseq;
#chomp @proteinseq;
#print RESULTS @proteinseq;
#close(RESULTS);
