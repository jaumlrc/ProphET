FROM        ubuntu:14.04
MAINTAINER  Nick Waters nickp60@gmail.com
RUN apt-get update
RUN apt-get install \
	build-essential \
	# ncbi-blast+ \ # this doesnt include the legacy_blast.pl :(
	emboss \
	git \
	expat \
	libexpat1-dev \
	libxml-parser-perl \
	libgd2-xpm-dev \
	cpanminus \
	bedtools --yes

RUN cpanm Carton

RUN curl ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.7.1/ncbi-blast-2.7.1+-x64-linux.tar.gz -o ncbi-blast-2.7.1+-x64-linux.tar.gz ;  tar zxvpf ncbi-blast-2.7.1+-x64-linux.tar.gz ; mv ncbi*/bin/* /usr/bin/
RUN cachebuster=c9s53b3he git clone https://github.com/nickp60/ProphET.git
RUN cpanm Module::CPANfile
RUN cd ProphET && carton install
RUN cd ProphET &&  carton exec ./INSTALL.pl

WORKDIR ProphET
RUN carton exec ./ProphET_standalone.pl --help
ENTRYPOINT [ "carton", "exec", "ProphET_standalone.pl" ]
