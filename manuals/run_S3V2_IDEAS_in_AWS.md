# Running S3V2-IDEAS in operating system set in AWS

### Outline
**[(1) Install gsl python and R and used library](#Install-gsl-python-and-R-and-used-library.)**<br>
#####
**[(2) Running S3V2-IDEAS in AWS](#Running-S3V2-IDEAS-in-AWS)**<br>


####
## S3V2-IDEAS require to use the AWS Debian system with the default settings to run
## Install gsl python and R and used library
```
### start a AWS EC2 instance with Debian system (at least 4CUP and 32GB mem)
### login AWS instance
ssh -i ~/.ssh/your_key_pem_file.pem admin@your_instance_open_IP_address

### apt-get install packages 
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
```

## 
## Running S3V2-IDEAS in AWS
```
### git clone S3V2-IDEAS
mkdir /home/admin/app
mkdir /home/admin/app/git
cd /home/admin/app/git
git clone https://github.com/guanjue/S3V2_IDEAS_ESMP.git

### test run
cd /home/admin/app/git/S3V2_IDEAS_ESMP/test_data
time bash run_S3V2_IDEAS_ESMP_AWS.sh
```


###### Then, users need to modify the "metadata.forEScall_AWS.txt" and "run_S3V2_IDEAS_ESMP_AWS.sh" files to run the pipeline on their own datasets.

###### The instruction for setting up the "metadata.forEScall_AWS.txt" and "run_S3V2_IDEAS_ESMP_AWS.sh" files can be found in the following page:
**[Inputs-ParameterSettings-RunningSteps-Outputs for S3V2_IDEAS_ESMP](https://github.com/guanjue/S3V2_IDEAS_ESMP/blob/master/manuals/inoutput_for_S3V2_IDEAS_pipeline.md)**<br>




