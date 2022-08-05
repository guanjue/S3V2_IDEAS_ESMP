### start a AWS EC2 instance with Debian system (at least 4CUP and 32GB mem)

### 



### 
sudo apt-get update
### install 
sudo apt-get install -y --no-install-recommends build-essential r-base r-cran-randomforest python3.6 python3-pip python3-setuptools python3-dev git time bash wget make bedtools vim

### install GSL
mkdir /home/admin/gsl
cd /home/admin/gsl/ && mkdir /home/admin/gsl/gsl && cd /home/admin/gsl/gsl && wget ftp://ftp.gnu.org/gnu/gsl/gsl-2.2.1.tar.gz && tar -zxvf gsl-2.2.1.tar.gz && cd /home/admin/gsl/gsl/gsl-2.2.1 && mkdir /home/admin/gsl/gsl/gsl && ./configure --prefix=/home/admin/gsl/gsl/gsl && make && make check && make install 
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/admin/gsl/gsl/gsl/lib/

### install python lib
pip3 install numpy
pip3 install scipy
pip3 install multiprocess
### install R lib
R
install.packages('data.table')
install.packages('doParallel')
install.packages('foreach')

### install S3V2-IDEAS
mkdir /home/admin/app
mkdir /home/admin/app/git
cd /home/admin/app/git
git clone https://github.com/guanjue/S3V2_IDEAS_ESMP.git

### 
