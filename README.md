
<a href="http://www.biorxiv.org/content/early/2017/08/16/176750"><h2>ProphET, Prophage Estimation Tool: a standalone prophage sequence prediction tool with self-updating reference database.</h2></a>

João L. Reis-Cunha<sup>1,2</sup>, Daniella C. Bartholomeu<sup>2</sup>, Ashlee M. Earl<sup>1</sup>,  Bruce W. Birren<sup>1</sup>, Gustavo C. Cerqueira<sup>1</sup>


<a href="http://www.biorxiv.org/content/early/2017/08/16/176750">Manuscript draft in BioRxiv</a>

------

<sup>1</sup> Broad Institute of Harvard and MIT, Cambridge, Massachusetts, United States

<sup>2</sup> Instituto de Ciências Biológicas, Universidade Federal de Minas Gerais, Brazil

------
<h4>Contact</h4>

<a href="mailto:jaumlrc@gmail.com">jaumlrc@gmail.com</a>

gustavo@broadinstitute.org

------
<h4>Required libraries and programs:</h4>

Broad users don't need to install any of the of programs and libraries listed below. If you are **Broadie** please follow the instructions on [README_BROAD_USERS.md](README_BROAD_USERS.md)  before installing and running ProphET.

* EMBOSS suite

* BEDTools suite

* BLAST

* Perl module Bio::Perl

* Perl module Bio::Graphics

* Perl module LWP::Simple

* Perl module XML::Simple

* Perl module GD

* Perl moduel GD::SVG



------
<h4>ProphET installation:</h4>

To either install ProphET or to update ProphET bacteriophage database please execute the following command from ProphET's home directory:
```
$ ./INSTALL.pl
```

This will search for required libraries, set the paths of required programs and download from Genbank (NCBI) all genomes associated to 16 families of bacteriophages
(listed in [config.dir/Prophages_names_sem_Claviviridae_Guttaviridae-TxID](config.dir/Prophages_names_sem_Claviviridae_Guttaviridae-TxID) ).


Some warnings will be issued during the setup of ProphET DB. See some examples below: 
```
Warning: bad /anticodon value '(pos:complement(13054..13056),aa:Met,seq:cat)'
Warning: NC_022920: Bad value '(pos:complement(13054..13056),aa:Met,seq:cat)' for tag '/anticodon'
```
Those warnings refer to unexpected format for coordinates of tRNA features and they won't affect the execution.


**If the script fails and reports missing Perl modules/libraries, please folow the instrucions on file  [README_INSTALLING_PERL_MODULES.md](README_INSTALLING_PERL_MODULES.md) on how to install those.**
 

------

<h4>Testing installation:</h4>

From ProphET's home directory execute the following command:
```
$ ./ProphET_standalone.pl --fasta test.fasta --gff_in test.gff --outdir test
```
The execution should take ~ 5 minutes.

Three putative prophages should be reported and its coordinates indicated in the file *test/phages_coords*:
```
NC_005362.1     1       327710  378140
NC_005362.1     2       502194  519268
NC_005362.1     3       1292553 1330556
```

The nucleotide sequence of each prophage can be found in:
```
test/NC_005362.1.phage_1.fas
test/NC_005362.1.phage_2.fas
test/NC_005362.1.phage_3.fas
```

A simple diagram depicting all coding genes in the bacterial genome, coding genes with significant matches to phage genes and the location of predicted prophages can be found in:
```
test/NC_005362.1.svg
```


------

<h4>Usage:</h4>

```
NAME
       ProphET is a user friendly algorithm to identify prophages within prokaryote genomes.

SYNOPSIS
       usage: ProphET_standalone.pl --fasta_in <file> --gff_in <file> 
       --outdir <string> [--grid] [--gff_trna <file> ] [--help]

OPTIONS
       --fasta_in - Bacterial genome Fasta file

       --gff_in - Bacterial GFF file

       --gff_trna - Optional parameter, in case the tRNAs are reported in a separate GFF 
       please provide it here <(Optional)>

       --outdir - output directory

       --grid - Use UGER for BLAST jobs 
       (Currently only works in the Broad Institute UGER grid system)

       --help - print this message (Optional)
```

