
** ProphET, Prophage Estimation Tool: a standalone prophage sequence prediction tool with self-updating reference database. **

João L. Reis-Cunha1,2*, Bruce W. Birren1, Daniella C. Bartholomeu2, Gustavo C. Cerqueira1*


1 Broad Institute of Harvard and MIT, Cambridge, Massachusetts, United States
2 Instituto de Ciências Biológicas, Universidade Federal de Minas Gerais, Brazil



  Usage:
    usage: ProphET_standalone.pl --fasta_in <file> --gff_in <file> --outdir
      <string> [--grid] [--gff_trna <file> ] [--help]

  Options:
    --fasta - Bacterial genome Fasta file

    --gff_in - Bacterial GFF file

    --gff_trna - Optional parameter, in case the tRNA are reported in a
    separate GFF please provide it here <(Optional)>

    --outdir - output directory

    --grid - Use UGER for BLAST jobs

    --help - print this message (Optional)
