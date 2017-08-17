
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


------
<h4>Instalation:</h4>

To either install ProphET or to update ProphET bacteriophage database please execute the following command from ProphET home directory:
```
$ ./INSTALL.pl
```

This will set the paths of required programs and download from Genbank (NCBI) all genomes associated to 16 families of bacteriophages
(listed in *config.dir/Prophages_names_sem_Claviviridae_Guttaviridae-TxID*).
 

------
<h4>Execution:</h4>

```
Usage:
    usage: ProphET_standalone.pl --fasta_in <file> --gff_in <file> --outdir
      <string> [--grid] [--gff_trna <file> ] [--help]

  Options:s
    --fasta - Bacterial genome FASTA file

    --gff_in - Bacterial GFF file

    --gff_trna - Optional parameter, in case the tRNA are reported in a
    separate GFF please provide it here <(Optional)>

    --outdir - output directory

    --grid - Use UGER for BLAST jobs

    --help - print this message (Optional)
```
