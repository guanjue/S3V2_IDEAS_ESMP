## Prerequisites and S3V2_IDEAS_ESMP installation

### 
### Operating System are as follows:
```
The pipeline needs to be run on a linux-64 operating system with at least 8GB memory to run. 
For users without access to the linux-64 system, people can use the Docker to set up a Operating System environment to run the pipeline.
We have tested this pipeline in different linux-64 operating systems. 
The pipeline can be run on, but not limited to, the following two version of the linux-64 systems:

(1)
LSB Version:	:base-4.0-amd64:base-4.0-noarch:core-4.0-amd64:core-4.0-noarch:graphics-4.0-amd64:graphics-4.0-noarch:printing-4.0-amd64:printing-4.0-noarch
Distributor ID:	RedHatEnterpriseServer
Description:	Red Hat Enterprise Linux Server release 6.10 (Santiago)
Release:	6.10
Codename:	Santiago

(2)
LSB Version:	:core-4.1-amd64:core-4.1-noarch
Distributor ID:	CentOS
Description:	CentOS Linux release 7.8.2003 (Core)
Release:	7.8.2003
Codename:	Core

```

### Dependencies are as follows:
```
(1) python3 (https://www.python.org/download/releases/3.0/)
(2) python dependencies: 'numpy', 'scipy', 'multiprocess'
(3) R (https://www.r-project.org/)
(4) R dependencies: 'data.table', 'doParallel', 'foreach'
The pipeline should work in the latest version of R. We have also tested other versions of R and the pipeline should also work in the following versions: 3.5.2; 3.6.2; 4.0.0
(5) gawk
(6) bedtools: The instruction about bedtools can be found in the following link(https://bedtools.readthedocs.io/en/latest/content/installation.html)
(7) GSL 2.2.1: The instruction about GSL can be found in the following link (https://www.gnu.org/software/gsl/manual/gsl-ref.html)
After installing the GSL system, users need to add the ~/gsl/lib into the LD_LIBRARY_PATH
```

### Installing Dependencies and Setting Up Environment. [(Running the S3V2_IDEAS pipeline in Docker)](https://github.com/guanjue/S3V2_IDEAS_ESMP/blob/master/run_S3V2_IDEAS_in_Docker.md)
```
### Use 'apt-get' to install python3, R, bedtools and some basic softwares in Linux-64 system.
apt-get update && apt-get install -y --no-install-recommends build-essential r-base r-cran-randomforest python3.6 python3-pip python3-setuptools python3-dev git time bash wget make bedtools vim

### python3 packages
pip3 install numpy --user
pip3 install scipy --user
pip3 install multiprocess --user

### R packages
Rscript -e "install.packages('data.table')"
Rscript -e "install.packages('doParallel')"
Rscript -e "install.packages('foreach')"

### GSL 2.2.1: 
# The Admin permission is required to install the package. 
# However, it is often not available for users using the public linux servers. 
# Please contact the Administrator of server to see how to use this package.
cd /where_you_install_gsl_2_2_1/ && mkdir /where_you_install_gsl_2_2_1/gsl && cd /where_you_install_gsl_2_2_1/gsl && wget ftp://ftp.gnu.org/gnu/gsl/gsl-2.2.1.tar.gz && tar -zxvf gsl-2.2.1.tar.gz && cd /where_you_install_gsl_2_2_1/gsl/gsl-2.2.1 && mkdir /where_you_install_gsl_2_2_1/gsl/gsl && ./configure --prefix=/where_you_install_gsl_2_2_1/gsl/gsl && make && make check && make install 

### Clone the S3V2_IDEAS_ESMP from GitHub.
cd  /The_folder_you_want_to_install_S3V2_IDEAS_ESMP_package/
git clone https://github.com/guanjue/S3V2_IDEAS_ESMP.git

### Add the /where_you_install_gsl_2_2_1/gsl/gsl/lib into the LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/where_you_install_gsl_2_2_1/gsl/gsl/lib/

# The S3V2_IDEAS_ESMP pipeline should be able to run.
```

