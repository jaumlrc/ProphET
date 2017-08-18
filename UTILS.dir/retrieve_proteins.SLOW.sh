#!/bin/bash

family=$1
txid=$2


###################################
## Retrieving sequence from Genbank
echo "Retrieving representatives of virus family $family, TaxID $txid ..." 

../../UTILS.dir/fetch_genomes_based_on_taxid.pl $txid

# Split genbank file
../../UTILS.dir/split_genbank.pl $txid.gb

# Number of genbank files
num_gb_files=`ls -l *.gbk |  wc | awk '{print $1}'`
echo "Number of sequences retrieved from genbank: $num_gb_files" 

## List organisms and definition associated to each genbank file
if [ -f gb_description.txt ]
then 
    rm gb_description.txt
fi


echo "Generating a summary of each genome..."
for gb_seq in `ls -l *.gbk | awk '{print $9}'`
 do org=`awk '$1 == "ORGANISM" {print $0}' $gb_seq | sed 's/ORGANISM//' | sed 's/^ \+//' | sed 's/ /_/g'`
 def=`awk '$1 == "DEFINITION" {print $0}' $gb_seq | sed 's/DEFINITION//' | sed 's/^ \+//'`
 segment=`perl -lane 'print $seg if ( ($seg) = ( $_ =~ /\/segment=\"([\w\W]+?)\"/ ) )' $gb_seq`
 chromosome=`perl -lane 'print $seg if ( ($seg) = ( $_ =~ /\/chromosome=\"([\w\W]+?)\"/ ) )' $gb_seq` 
 
 if [[ ("$segment" == "") && ("$chromosome" != "")  ]]; then segment=$chromosome; fi
 if [[ ("$chromosome" == "") && ("$segment" != "")  ]]; then segment=$segment; fi
 if [[ ("$segment" == "") && ("$chromosome" == "")  ]]; then segment=ND; fi
 
 segment_print=`echo $segment | sed 's/[sS]egment //' | sed 's/ [sS]egment//'`
 gb=`echo $gb_seq | sed 's/\.gbk//'`
 
 echo -e "$gb\t$org\t$segment_print\t$def" >> gb_description.txt  
done

# Number of species
num_species=`awk '{print $2}' gb_description.txt | sort | uniq | wc | awk '{print $1}'`
echo "Number of species: $num_species" 


##########################################################
# Extract mat_peptide (segments of a polyprotein) and CDS

echo "Extracting segments of polyproteins and coding sequences..."
for i in `ls -l *.gbk | awk '{print $9}'`
 do filename_wo_gb=`echo $i | sed 's/\.gbk//g'`
 echo "Working on $i ..."
 # Extract sense
 extractfeat $i -sense 1 -type 'mat_peptide|CDS' -describe 'product|protein_id|locus_tag' -outseq $filename_wo_gb.sense -auto
 
 perl -pe 's/(>[\W\w]+?) /\1 + /g' $filename_wo_gb.sense > $filename_wo_gb.fasta  
 
 # Extract anti-sense
 extractfeat $i -sense -1 -type 'mat_peptide|CDS' -describe 'product|protein_id|locus_tag' -outseq $filename_wo_gb.antisense -auto
 
 perl -pe 's/(>[\W\w]+?) /\1 - /g' $filename_wo_gb.antisense >> $filename_wo_gb.fasta

done

#################################################
# Extract and report all features in the files
echo "Extracting all other protein coding features ..."
for i in `ls -l *.gbk | awk '{print $9}'`
 do filename_wo_gb=`echo $i | sed 's/\.gbk//g'`
 echo "Working on $i ..."
 # Extract sense
 extractfeat $i -describe 'product|protein_id|locus_tag' -outseq $filename_wo_gb.all_features -auto
done
grep ">" *.all_features | sed 's/>/ /'  | awk '{print $2,$3}' > all_features
awk '{print $2}' all_features | sort | uniq > all_feature_types


#################################################
# Retrieving the product of each CDS or mat_peptide
echo "Retrieving the featured product for each CDS or mat_peptide ..."

if [ -f feature_product.txt ]
then 
    rm feature_product.txt
fi

for i in *.fasta
 do grep ">" $i | perl -lane '$line = $_; ($feature) = ( $line =~ /^>([\w\W]+?)\s/ ); ($product) = ( $line =~ /product=\"(.+?)\"/ ); print "$feature\t$product"' >> feature_product.txt
done


##############################################
# Translate genes
echo "Translating genes ..."
for i in *.fasta; do transeq $i $i.prot -auto; done

# Combine all proteins
echo "Gathering all protein sequences in a single file ..."
cat *.prot > all.prot.fas

# Remove "*" STOP codons from each sequence
echo "Removing * representing STOP codons ..."
../../UTILS.dir/fasta2line all.prot.fas | perl -lane '$seq = $F[0]; $seq =~ s/\*$//; print $seq . "\t" . $F[1]' | ../../UTILS.dir/line2fasta > tmp
mv tmp all.prot.fas


# Number of CDS and mat_peptides
num_CDS_mat=`grep -c ">" all.prot.fas`
echo "Number of features: $num_CDS_mat"

