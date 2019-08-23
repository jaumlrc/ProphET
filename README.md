
<a href="http://www.biorxiv.org/content/early/2017/08/16/176750"><h2>ProphET, Prophage Estimation Tool: a standalone prophage sequence prediction tool with self-updating reference database.</h2></a>

João L. Reis-Cunha<sup>1,2</sup>, Daniella C. Bartholomeu<sup>2</sup>, Ashlee M. Earl<sup>1</sup>,  Bruce W. Birren<sup>1</sup>, Gustavo C. Cerqueira<sup>1</sup>


<a href="http://www.biorxiv.org/content/early/2017/08/16/176750">Manuscript draft in BioRxiv</a>

------

<sup>1</sup> Broad Institute of Harvard and MIT, Cambridge, Massachusetts, United States (2017)

<sup>2</sup> Instituto de Ciências Biológicas, Universidade Federal de Minas Gerais, Brazil

------
<h4>Contact</h4>

<a href="mailto:jaumlrc@gmail.com">jaumlrc@gmail.com</a>

gustavo@broadinstitute.org

------

ProphET is an open source software developed to be used in the Linux platform. Users are free to bundle executables and modify the script. However, we do not guarantee the efficiency and precision of predictions if modifications were performed in the script and dependencies.

------


<h4>Required libraries and programs:</h4>

Broad users don't need to install any of the of programs and libraries listed below. If you are **Broadie** please follow the instructions on [README_BROAD_USERS.md](README_BROAD_USERS.md)  before installing and running ProphET.

* EMBOSS suite

* BEDTools suite

* BLAST

* Perl module Bio::Perl

* Perl module SVG

* Perl module GD

* Perl moduel GD::SVG

* Perl module Bio::Graphics

* Perl module LWP::Simple

* Perl module XML::Simple

* Perl module Mozilla::CA

* Perl module LWP::Protocol::https

------

<h4>ProphET Dependencies:</h4>

ProphET requires that BLAST (legacy) EMBOSS and BEDTools are installed and added to the enviroment variable PATH.

**BLAST:**

BLAST legacy can be downloaded from the following link:
ftp://ftp.ncbi.nlm.nih.gov/blast/executables/legacy.NOTSUPPORTED/2.2.26/

Download BLAST legacy using wget for Linux:
```
$ wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/legacy.NOTSUPPORTED/2.2.26/blast-2.2.26-x64-linux.tar.gz
```

Unpack files using tar:
```
$ tar -xvzf blast-2.2.26-x64-linux.tar.gz
```

Change directory to where blastall and formatdb executables are: 
```
$ cd blast-2.2.26/bin
```

Print the full path to current folder using pwd:
```
$ pwd -P
```
example output from pwd: 
```
/data/usr/BLAST/blast-2.2.26/bin
```

Add BLAST binary folder to the $PATH environment variable:
```
PATH=<path_to_blastall/bin>:$PATH
```
example:
```
PATH=/data/usr/BLAST/blast-2.2.26/bin:$PATH
```

Test if blastall and formatdb commands are on the path:
```
$ blastall --> Will print the blastall help information
$ formatdb --> will print "[formatdb] ERROR: No database name was specified
```

**EMBOSS**

EMBOSS suite installation instructions can be obtained from this link: 
http://emboss.sourceforge.net/download/

**BEDTools**

BEDTools suite installation instructions can be obtained from this link:
https://bedtools.readthedocs.io/en/latest/content/installation.html

**Perl modules/libraries:**

If the script fails and reports missing Perl modules/libraries, please follow the instructions on file  [README_INSTALLING_PERL_MODULES.md](README_INSTALLING_PERL_MODULES.md) on how to install those.


**Adding third party programs to the $PATH enviroment**

EMBOSS and BEDTools folders also have to be added to the $PATH environment prior to run ProphET INSTALL.pl installation. This can be done using the following command:
```
PATH=<path/to/EMBOSS/>:$PATH
PATH=<path/to/BEDTools>$PATH
```
Example:
```
PATH=/home/bin/EMBOSS-6.3.1/emboss/:$PATH
PATH=/home/bin/bedtools:$PATH
```
 
 
**ProphET was tested with the following versions of third-party dependencies:**

-blast-legacy-2.2.26
ftp://ftp.ncbi.nlm.nih.gov/blast/executables/legacy.NOTSUPPORTED/2.2.26/blast-2.2.26-x64-linux.tar.gz
-blast-legacy-2.2.9 
ftp://ftp.ncbi.nlm.nih.gov/blast/executables/legacy.NOTSUPPORTED/2.2.9/blast-2.2.9-amd64-linux.tar.gz
-bedtools v2.23.0 
https://github.com/arq5x/bedtools2/releases/download/v2.23.0/bedtools-2.23.0.tar.gz
-bedtools-v2.28.0 
https://github.com/arq5x/bedtools2/releases/download/v2.28.0/bedtools-2.28.0.tar.gz
-EMBOSS-6.3.1 
ftp://emboss.open-bio.org/pub/EMBOSS/old/6.3.1/EMBOSS-6.3.1.tar.gz
-EMBOSS-6.6.0 
ftp://emboss.open-bio.org/pub/EMBOSS/EMBOSS-6.6.0.tar.gz

 
------

<h4>ProphET installation:</h4>

To install ProphET and download bacteriophage database please execute the following command from ProphET's home directory:
```
$ ./INSTALL.pl
```

This will search for required libraries, set the paths of required programs and download from Genbank (NCBI) all genomes associated to 16 families of bacteriophages
(listed in [config.dir/Prophages_names_sem_Claviviridae_Guttaviridae-TxID](config.dir/Prophages_names_sem_Claviviridae_Guttaviridae-TxID) ).

IMPORTANT: Please ensure that the third-party programs blastall, EMBOSS and BEDTools were added to the environment variable PATH, or the installation will crash. See “ProphET Dependencies” section for more instructions.


Some warnings will be issued during the setup of ProphET DB. See some examples below: 
```
Warning: bad /anticodon value '(pos:complement(13054..13056),aa:Met,seq:cat)'
Warning: NC_022920: Bad value '(pos:complement(13054..13056),aa:Met,seq:cat)' for tag '/anticodon'
```
Those warnings refer to unexpected format for coordinates of tRNA features and they won't affect the execution.

------

<h4>Testing installation:</h4>

From ProphET's home directory execute either the following command (GFF file containing both coding genes and tRNAs):
```
$ ./ProphET_standalone.pl --fasta test.fasta --gff_in test.gff --outdir test
```

The execution should take ~ 5 minutes.

Two putative prophages should be reported and their coordinates indicated in the file *test/phages_coords*:
```
FORMAT:
<scaffold>  <#prophage> <genomic.start.coord> <genomic.end.coord>

CONTENT:
NC_005362.1     1       327710  378140
NC_005362.1     2       1292553 1330556
```

Small differences between the coordinates reported above and the coordinates obtained by your first test run of ProphET are expected.
This is due to changes in the database of known prophage proteins, which is updated on each installation of ProphET. 


The nucleotide sequence of each prophage can be found in:
```
test/NC_005362.1.phage_1.fas
test/NC_005362.1.phage_2.fas
```

The program also renders a simple diagram depicting all coding genes in the bacterial genome, coding genes with significant matches to phage genes and the location of predicted prophages:
```
test/NC_005362.1.svg
```

------

<h4>Updating phage database:</h4>

To update the bacteriophage database with the latest sequences deposited at Genbank, please execute the following command from ProphET's home directory:

```
$ ./INSTALL.pl --update_db_only
```

The current database will be backed up as 
```
PhrophET_phage_proteins_database.dir.<current date and time>.bak
```

and an updated database will be saved at:
```
PhrophET_phage_proteins_database.dir
```

All instances of the prophage DB (current and backups) include a file  reporting the download date and stats.: `phage_db.summary.stats`. This file is copied to the results directory of every ProphET execution to enable auditing and reproducibility of results.


------
<h4>Before running ProphET in your favorite bacterial genome</h4>

* Check if the GFF file that will be provided to ProphET has the format specified by [The Sequence Ontology Consortium](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md)

* If your GFF does not meet those specifications, a converter is provided as part of GFFLib (package installed during Prophet setup): 
```
GFFLib/gff_rewrite --input <GFF input> -output <GFF output> --add_missing_features
```

* The GFF converter will not work for all cases. If you happen to encounter one of those, please issue a ticket reporting that.

* Check if all sequences IDs in the FASTA file (header of each sequence) matches perfectly the source field in the GFF file (first column of the GFF) and vice-versa.

------

<h4>Usage:</h4>

```
    ProphET_standalone.pl --fasta_in <file> --gff_in <file> --outdir
    <string> [--grid] [--gff_trna <file> ] [--help]

Options:
    --fasta_in - Bacterial genome Fasta file

    --gff_in - Bacterial GFF file

    --gff_trna - Optional parameter, in case the tRNAs are reported in a
    separate GFF please provide it here <(Optional)>

    --outdir - output directory

    --grid - Use UGER for BLAST jobs (Currently only works in the Broad
    Institute UGER grid system) (Optional)

    --help - print this and some additional info. about FASTA and GFF input
    format (Optional)
```

