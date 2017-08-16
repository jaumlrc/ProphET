#!user/bin/perl
use strict;
#Utiliza como entrada um arquivo contento o nome da família do virus e o Taxonomy ID do NCBI separados por tab, e utiliza o script do Gustavo para obter todas as informacoes do genoma.

my $input = $ARGV[0];
open (INPUT, "$input");
my @N_seq = <INPUT>;
close(INPUT);

for (my $i=0; $i<=$#N_seq; $i++) {
	my @temp1 = split (/\t/, $N_seq[$i]); 
	my $name = $temp1[0];
	my $Tx_id = $temp1[1];
	mkdir ("$name.dir", 0755);
	chdir "$name.dir";
	`../../UTILS.dir/./retrieve_proteins.sh $name $Tx_id report.txt`;
	chdir "../";
}
