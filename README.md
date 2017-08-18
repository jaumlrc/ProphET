
<h2>ProphET, Prophage Estimation Tool: a standalone prophage sequence prediction tool with self-updating reference database.</h2>

João L. Reis-Cunha<sup>1,2</sup>, Daniella C. Bartholomeu<sup>2</sup>, Ashlee M. Earl<sup>1</sup>,  Bruce W. Birren<sup>1</sup>, Gustavo C. Cerqueira<sup>1</sup>

------

<sup>1</sup> Broad Institute of Harvard and MIT, Cambridge, Massachusetts, United States

<sup>2</sup> Instituto de Ciências Biológicas, Universidade Federal de Minas Gerais, Brazil

------
<h4>Contact</h4>
jaumlrc@gmail.com

gustavo@broadinstitute.org

------
<h4>Required libraries and programs:</h4>

* EMBOSS suite

* BEDTools suite

* BLAST

* Perl module Bio::Perl

* Perl module LWP::Simple

* Perl module XML::Simple

* Perl module GD



------
<h4>Installing CPAN:</h4>

Perl modules can be installed using CPAN. Please first certify that CPAN is installed and configured by issuing the command below:

``$ perl -MCPAN -e shell``


* If the command above returns the prompt ```cpan[1]>``` or similar prompt then CPAN is already configured. So quit the cpan shell by typing: ```cpan[1]> quit```

* If the command returns a text saying that 'CPAN requires configuration...' follow the steps for automatic configuration. Select the default option in every question. Quit CPAN after the configuration is done by typing: ```cpan[1]> quit```

* If the command returns: ```Can't locate CPAN.pm in @INC (@INC contains:... ``` then you will need Administrative privileges to install CPAN either using apt-get ```sudo  apt-get install build-essential``` or yum ```sudo yum install perl-CPAN```

<h4>Installing Perl modules:</h4>

Now you are ready to install the required Perl modules. Issue the following commands:

```
$ perl -MCPAN -e 'install Bio::Perl'
$ perl -MCPAN -e 'install LWP::Simple'
$ perl -MCPAN -e 'install XML::Simple'
$ perl -MCPAN -e 'install GD'
```

------
<h4>ProphET installation:</h4>

To either install ProphET or to update ProphET bacteriophage database please execute the following command from ProphET's home directory:
```
$ ./INSTALL.pl
```

This will set the paths of required programs and download from Genbank (NCBI) all genomes associated to 16 families of bacteriophages
(listed in *config.dir/Prophages_names_sem_Claviviridae_Guttaviridae-TxID*).
 

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
ProphET_standalone.pl --fasta_in <file> --gff_in <file> --outdir
      <string> [--grid] [--gff_trna <file> ] [--help]

  Options:
    --fasta_in - Bacterial genome FASTA file

    --gff_in - Bacterial GFF file

    --gff_trna - Optional parameter, in case the tRNA are reported in a
    separate GFF please provide it here <(Optional)>

    --outdir - output directory

    --grid - Use UGER for BLAST jobs (Currently only works in the Broad Institute UGER grid system)

    --help - print this message (Optional)
```

