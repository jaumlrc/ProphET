FROM        ubuntu:14.04
# FROM        perl:latest
MAINTAINER  Nick Waters nickp60@gmail.com
RUN cachebuster=c953b3f apt-get update
RUN apt-get install \
	build-essential \
	# ncbi-blast+ \
	emboss \
	git \
	expat \
	libexpat1-dev \
	libxml-parser-perl \
	libgd2-xpm-dev \
	cpanminus \
	bedtools --yes
# RUN curl -L http://cpanmin.us | perl - App::cpanminus
RUN cpanm Carton
# RUN cpanm --installdeps .

RUN cachebuster=c953b32 git clone https://github.com/nickp60/ProphET.git
# we gotta manually install the GFFlib repo, as it doesnt have a Makefile.pl
# RUN cachebuster=b953b36 git clone https://github.com/gustavo11/GFFLib.git UTILS.dir/GFFLib
RUN cachebuster=c953b3d curl ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.7.1/ncbi-blast-2.7.1+-x64-linux.tar.gz -o ncbi-blast-2.7.1+-x64-linux.tar.gz ;  tar zxvpf ncbi-blast-2.7.1+-x64-linux.tar.gz ; mv ncbi*/bin/* /usr/bin/
RUN cpanm Module::CPANfile
RUN which blastn
RUN which legacy_blast.pl
RUN cachebuster=b953b3g cd ProphET && carton install
	#--deployment
RUN cachebuster=b953b36s cd ProphET &&  ./INSTALL.pl

WORKDIR ProphET
CMD carton exec ./ProphET_standalone.pl
