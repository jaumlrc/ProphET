FROM        ubuntu:14.04
MAINTAINER  Nick Waters nickp60@gmail.com
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

RUN cpanm Carton

# RUN curl ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.7.1/ncbi-blast-2.7.1+-x64-linux.tar.gz -o ncbi-blast-2.7.1+-x64-linux.tar.gz ;  tar zxvpf ncbi-blast-2.7.1+-x64-linux.tar.gz ; mv ncbi*/bin/* /usr/bin/
RUN git clone https://github.com/nickp60/ProphET.git
RUN cpanm Module::CPANfile
WORKDIR ProphET
RUN carton install
RUN carton exec ./INSTALL.pl
RUN git pull
# COPY ./ProphET_standalone.pl ./ProphET_standalone.pl
RUN carton exec ./ProphET_standalone.pl --help
RUN carton exec ./ProphET_standalone.pl --fasta_in test.fasta --gff_in test.gff --outdir tmp --cores 2
# ENTRYPOINT [ "carton", "exec", "/ProphET/ProphET_standalone.pl" ]
