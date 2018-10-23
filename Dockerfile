FROM        ubuntu:14.04
# FROM        perl:latest
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
	cpanminus --yes
# RUN curl -L http://cpanmin.us | perl - App::cpanminus
RUN cpanm Carton
# RUN cpanm --installdeps .

RUN cachebuster=c953b35 git clone https://github.com/nickp60/ProphET.git
# we gotta manually install the GFFlib repo, as it doesnt have a Makefile.pl
RUN cachebuster=b953b36 git clone https://github.com/gustavo11/GFFLib.git UTILS.dir/GFFLib
RUN cpanm Module::CPANfile
RUN cd ProphET && carton install
	#--deployment
RUN  cd ProphET &&  ./INSTALL.pl

WORKDIR ProphET
CMD carton exec ./ProphET_standalone.pl
