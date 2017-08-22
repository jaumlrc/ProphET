#!/usr/bin/env perl

use strict;
# Retrieves and concatenate the protein sequences retrieved from Genbank (NCBI)
#
# This script seems to be unecessary after the latest modifications to
# retrieve_proteins.sh. I explain...
#
# Earlier retrieve_proteins.sh produced the proteome file of each family.
# But now, for the sake of improving speed, one single file with the proteome
# of all families is being generated 

my $input = $ARGV[0];     # File containing the name of bacteriophage families
open (INPUT, "$input");
my @N_seq = <INPUT>;
close(INPUT);


my @proteinseq = ();
for (my $j=0; $j<=$#N_seq; $j++) {
        my @temp1 = split (/\t/, $N_seq[$j]);
        my $name = $temp1[0];
	chdir "$name.dir";
	my @temp_seq = `cat all.prot.fas`;
	push (@proteinseq, @temp_seq);
	chdir "../";
}

print @proteinseq;
