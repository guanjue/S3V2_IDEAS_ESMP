FROM ubuntu:latest

COPY . /app
WORKDIR /app

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends build-essential r-base r-cran-randomforest python3.6 python3-pip python3-setuptools python3-dev git time bash wget make bedtools vim

RUN pip3 install numpy
RUN pip3 install scipy
RUN pip3 install multiprocess

RUN Rscript -e "install.packages('data.table')"
RUN Rscript -e "install.packages('doParallel')"
RUN Rscript -e "install.packages('foreach')"

RUN git clone https://github.com/guanjue/S3V2_IDEAS_ESMP.git

RUN cd /app && mkdir /app/gsl && cd /app/gsl && wget ftp://ftp.gnu.org/gnu/gsl/gsl-2.2.1.tar.gz && tar -zxvf gsl-2.2.1.tar.gz && cd /app/gsl/gsl-2.2.1 && mkdir /app/gsl/gsl && ./configure --prefix=/app/gsl/gsl && make && make check && make install 
