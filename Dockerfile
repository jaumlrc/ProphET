FROM        ubuntu:14.04
MAINTAINER  Nick Waters nickp60@gmail.com
RUN export LC_ALL=en_US.UTF-8
RUN export LANG=en_US.UTF-8
RUN locale-gen en_US.UTF-8
RUN apt-get update
RUN apt-get install \
	build-essential \
	ncbi-blast+ \
	emboss \
	git \
	expat \
	libexpat1-dev \
	libxml-parser-perl \
	libgd2-xpm-dev \
	cpanminus \
        perl-doc \
	bedtools --yes


# RUN curl ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.7.1/ncbi-blast-2.7.1+-x64-linux.tar.gz -o ncbi-blast-2.7.1+-x64-linux.tar.gz ;  tar zxvpf ncbi-blast-2.7.1+-x64-linux.tar.gz ; mv ncbi*/bin/* /usr/bin/
RUN git clone https://github.com/nIckp60/ProphET.git
RUN cpanm Module::CPANfile
WORKDIR ProphET
RUN cpanm Data::Stag
RUN cpanm IO::String
RUN cpanm Bio::Perl
RUN cpanm SVG
RUN cpanm GD
RUN cpanm GD::SVG
RUN cpanm Bio::Graphics
RUN cpanm LWP::Simple
RUN cpanm XML::SAX::Expat
RUN cpanm XML::Simple
RUN cpanm Mozilla::CA
RUN ./INSTALL.pl
RUN git pull https://github.com/nIckp60/ProphET.git
# COPY ./ProphET_standalone.pl ./ProphET_standalone.pl
RUN ./ProphET_standalone.pl --help
RUN ./ProphET_standalone.pl --fasta_in test.fasta --gff_in test.gff --outdir tmp --cores 2
RUN which perl
ENTRYPOINT [ "/usr/bin/perl", "/ProphET/ProphET_standalone.pl" ]
