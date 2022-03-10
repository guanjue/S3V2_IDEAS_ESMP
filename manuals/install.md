## Prerequisites and S3V2_IDEAS_ESMP installation

### 
### Operating System are as follows:
```
The pipeline needs to be run on a linux-64 operating system with at least 8GB memory to run. 
For users without access to the linux-64 system, people can use the Docker to set up a Operating System environment to run the pipeline.
We have tested this pipeline in different linux-64 operating systems. 
The pipeline can be run on, but not limited to, the following two version of the linux-64 systems:

Noted: This pipeline was designed to incorporate multiple epigenomic dataset across multiple cell types. 
This type of analysis is typically run on the high-performance servers. 
So, running the pipeline in HPC server is suggested.

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
##### The pipeline will generated several whole genome matrix which require large memory.
##### For example, a system has 7 epigenetic features in 21 cell types in human, we used 30GB memories and 4 threads to run the S3V2-IDEAS pipeline. It will take around 25 hours to be done.
##### Other systems and the computing time that we have tried have been shown in the following table:

<img src="https://github.com/guanjue/S3V2_IDEAS_ESMP/blob/master/figures/computing_resource_timing.png" width="600"/>


### Dependencies are as follows:
```
(1) python3 (https://www.python.org/download/releases/3.0/)
(2) python dependencies: 'numpy_1.13.3', 'scipy_0.19.1', 'multiprocess_0.70.11.1' & later versions
(3) R (https://www.r-project.org/): The pipeline should work in the latest version of R. We have also tested other versions of R and the pipeline should also work in the following versions: 3.5.2; 3.6.2; 4.0.0
(4) R dependencies: 'doParallel_1.0.16', 'foreach_1.4.8', 'data.table_1.12.0' & later versions
(5) gawk
(6) bedtools (v2.27.1 and above): The instruction about bedtools can be found in the following link(https://bedtools.readthedocs.io/en/latest/content/installation.html)
(7) GSL 2.2.1: The instruction about GSL can be found in the following link (https://www.gnu.org/software/gsl/manual/gsl-ref.html)
After installing the GSL system, users need to add the ~/gsl/lib into the LD_LIBRARY_PATH
```

### Conda Installing Dependencies and Setting Up Environment. 
```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --add channels mvdbeek
conda update conda

### conda create s3v2_ideas env
conda create -n s3v2_ideas gsl=2.2.1 python=3 r=3.6 bedtools ucsc_tools numpy scipy multiprocess r-doParallel r-foreach r-data.table

### conda activate s3v2_ideas env
conda activate s3v2_ideas

```


### Manually Installing Dependencies and Setting Up Environment. 
```
### Use 'apt-get' to install python3, R 3.5.2, bedtools (v2.27.1) and some basic softwares in Linux-64 system.
apt-get update && apt-get install -y --no-install-recommends build-essential r-base r-cran-randomforest python3.6 python3-pip python3-setuptools python3-dev git time bash wget make bedtools vim

### python3 packages
pip3 install numpy --user
pip3 install scipy --user
pip3 install multiprocess --user
### Version of the python modules works for our operating system:
### numpy_1.13.3 scipy_0.19.1     multiprocess_0.70.11.1
### Other versions especially the latest version of the python modules should works as well


### R packages
Rscript -e "install.packages('data.table')"
Rscript -e "install.packages('doParallel')"
Rscript -e "install.packages('foreach')"
### Version of the R libraries works for our operating system:
### doParallel_1.0.16 foreach_1.4.8     data.table_1.12.0
### Other versions especially the latest version of the R libraries should works as well

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

## If the installations have problems cannot be solve, you can try to run the package in Docker with preset environment:
### [(Running the S3V2_IDEAS pipeline in Docker)](https://github.com/guanjue/S3V2_IDEAS_ESMP/blob/master/manuals/run_S3V2_IDEAS_in_Docker.md)

