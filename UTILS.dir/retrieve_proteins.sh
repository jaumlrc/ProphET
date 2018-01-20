#!/bin/bash

family=$1
txid=$2


###################################
## Retrieving sequence from Genbank
echo "Retrieving representatives of virus family $family, TaxID $txid ..." 

../../UTILS.dir/fetch_genomes_based_on_taxid.pl $txid


##########################################################
# Extract mat_peptide (segments of a polyprotein) and CDS

echo "Extracting segments of polyproteins and coding sequences..."
extractfeat $txid.gb -sense 1 -type 'mat_peptide|CDS' -describe 'product|protein_id|locus_tag' -outseq $txid.sense -auto
perl -pe 's/(>[\W\w]+?) /\1 + /g' $txid.sense > $txid.fasta  
 
# Extract anti-sense
extractfeat $txid.gb -sense -1 -type 'mat_peptide|CDS' -describe 'product|protein_id|locus_tag' -outseq $txid.antisense -auto
perl -pe 's/(>[\W\w]+?) /\1 - /g' $txid.antisense >> $txid.fasta


#################################################
# Extract and report all features in the files
echo "Extracting all other protein coding features ..."
extractfeat $txid.gb -describe 'product|protein_id|locus_tag' -outseq $txid.all_features -auto
grep ">" $txid.all_features | sed 's/>/ /'  | awk '{print $2,$3}' > all_features
awk '{print $2}' all_features | sort | uniq > all_feature_types


#################################################
# Retrieving the product of each CDS or mat_peptide
echo "Retrieving the featured product for each CDS or mat_peptide ..."

if [ -f feature_product.txt ]
then 
    rm feature_product.txt
fi

grep ">" $txid.fasta | perl -lane '$line = $_; ($feature) = ( $line =~ /^>([\w\W]+?)\s/ ); ($product) = ( $line =~ /product=\"(.+?)\"/ ); print "$feature\t$product"' >> feature_product.txt


##############################################
# Translate genes
echo "Translating genes ..."
transeq $txid.fasta $txid.prot -auto

# Remove "*" STOP codons from each sequence
echo "Removing * representing STOP codons ..."
../../UTILS.dir/fasta2line $txid.prot | perl -lane '$seq = $F[0]; $seq =~ s/\*$//; print $seq . "\t" . join " ", @F[1..$#F]' | ../../UTILS.dir/line2fasta > all.prot.fas

# Number of CDS and mat_peptides
num_CDS_mat=`grep -c ">" all.prot.fas`
echo "Number of coding features under TaxID $txid : $num_CDS_mat"
echo "Number of coding features under TaxID $txid : $num_CDS_mat" >> $txid.prot.log

