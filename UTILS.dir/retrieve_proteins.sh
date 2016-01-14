#!/bin/bash

family=$1
txid=$2
report_file=$3


###################################
## Retrieving sequence from Genbank
echo "Retrieving representatives of virus family: $txid" 

`../../UTILS.dir/./fetch_genomes_based_on_taxid.pl $txid` #adicionei o ./ e retirei ~gustavo/devel/robots/robot_genbank/scripts/

# Split genbank file
# for i in *.gb; do split_genbank.pl $i; done
`../../UTILS.dir/./split_genbank.pl $txid.gb` #adicionei o ./

# Number of genbank files
num_gb_files=`ls -l *.gb | grep -v original.gb | wc | awk '{print $1}'`
echo "    Number of sequences retrieved from genbank: " $num_gb_files >> $report_file

# List organisms and definition associated to each genbank file
rm gb_description.txt
for gb_seq in `ls -l *.gb | grep -v original.gb | awk '{print $9}'`
 do org=`awk '$1 == "ORGANISM" {print $0}' $gb_seq | sed 's/ORGANISM//' | sed 's/^ \+//' | sed 's/ /_/g'`
 def=`awk '$1 == "DEFINITION" {print $0}' $gb_seq | sed 's/DEFINITION//' | sed 's/^ \+//'`
 segment=`perl -lane 'print $seg if ( ($seg) = ( $_ =~ /\/segment=\"([\w\W]+?)\"/ ) )' $gb_seq`
 chromosome=`perl -lane 'print $seg if ( ($seg) = ( $_ =~ /\/chromosome=\"([\w\W]+?)\"/ ) )' $gb_seq` 
 
 if [[ ("$segment" == "") && ("$chromosome" != "")  ]]; then segment=$chromosome; fi
 if [[ ("$chromosome" == "") && ("$segment" != "")  ]]; then segment=$segment; fi
 if [[ ("$segment" == "") && ("$chromosome" == "")  ]]; then segment=ND; fi
 
 segment_print=`echo $segment | sed 's/[sS]egment //' | sed 's/ [sS]egment//'`
 gb=`echo $gb_seq | sed 's/\.gb//'`
 
 echo -e "$gb\t$org\t$segment_print\t$def" >> gb_description.txt  
done

# Number of species
num_species=`awk '{print $2}' gb_description.txt | sort | uniq | wc | awk '{print $1}'`
echo "    Number of species: " $num_species >> $report_file


##########################################################
# Extract mat_peptide (segments of a polyprotein) and CDS

for i in `ls -l *.gb | grep -v original.gb | awk '{print $9}'`
 do filename_wo_gb=`echo $i | sed 's/\.gb//g'`
 echo "Working on $i ..."
 # Extract sense
 extractfeat $i -sense 1 -type 'mat_peptide|CDS' -describe 'product|protein_id|locus_tag' -outseq $filename_wo_gb.sense  
 
 perl -pe 's/(>[\W\w]+?) /\1 + /g' $filename_wo_gb.sense > $filename_wo_gb.fasta  
 
 # Extract anti-sense
 extractfeat $i -sense -1 -type 'mat_peptide|CDS' -describe 'product|protein_id|locus_tag' -outseq $filename_wo_gb.antisense
 
 perl -pe 's/(>[\W\w]+?) /\1 - /g' $filename_wo_gb.antisense >> $filename_wo_gb.fasta

done

#################################################
# Extract and report all features in the files
for i in `ls -l *.gb | grep -v original.gb | awk '{print $9}'`
 do filename_wo_gb=`echo $i | sed 's/\.gb//g'`
 echo "Working on $i ..."
 # Extract sense
 extractfeat $i -describe 'product|protein_id|locus_tag' -outseq $filename_wo_gb.all_features
done
grep ">" *.all_features | sed 's/>/ /'  | awk '{print $2,$3}' > all_features
awk '{print $2}' all_features | sort | uniq > all_feature_types


#################################################
# Retrieving the product of each CDS or mat_peptide
rm feature_product.txt
for i in *.fasta
 do grep ">" $i | perl -lane '$line = $_; ($feature) = ( $line =~ /^>([\w\W]+?)\s/ ); ($product) = ( $line =~ /product=\"(.+?)\"/ ); print "$feature\t$product"' >> feature_product.txt
done


##############################################
# Translate genes
for i in *.fasta; do transeq $i $i.prot; done

# Combine all proteins
cat *.prot > all.prot.fas

# Remove "*" STOP codons from each sequence
`../../UTILS.dir/./fasta2line all.prot.fas | perl -lane '$seq = $F[0]; $seq =~ s/\*$//; print $seq . "\t" . $F[1]' | ../../UTILS.dir/./line2fasta > tmp` ## Adicionei o ./
mv tmp all.prot.fas


# Number of CDS and mat_peptides
num_CDS_mat=`grep -c ">" all.prot.fas`
echo "    Number of features: $num_CDS_mat" >> $report_file
