FROM        perl:latest
MAINTAINER  Nick Waters nickp60@gmail.com

RUN apt-get install curl
RUN curl -L http://cpanmin.us | perl - App::cpanminus
RUN cpanm Carton
# RUN cpanm --installdeps .

RUN cachebuster=b953b35 git clone https://github.com/nickp60/ProphET.git
# we gotta manually install the GFFlib repo, as it doesnt have a Makefile.pl
RUN cachebuster=b953b36 git clone https://github.com/gustavo11/GFFLib.git UTILS.dir/GFFLib
RUN cd ProphET && carton install --deployment
RUN ./INSTALL.pl

WORKDIR ProphET
CMD carton exec ./ProphET_standalone.pl
