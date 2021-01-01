## Prerequisites and S3V2_IDEAS_ESMP installation

### 
### Operating System are as follows:
```
The pipeline needs to be run on a linux-64 operating system. For users without access to the linux-64 system, people can use the Docker to set up a Operating System environment to run the pipeline.

We have tested this pipeline in different linux-64 operating systems. The pipeline can be run on, but not limited to, the following two version of the linux-64 systems:

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
(4) R dependencies: 'data.table', 'doParallel', 'foreach', 'LSD', 'pheatmap'
The pipeline should work in the latest version of R. We have also tested other versions of R and the pipeline should also work in the following versions: 3.5.2; 3.6.2; 4.0.0
(5) gawk
(6) bedtools (https://bedtools.readthedocs.io/en/latest/content/installation.html)

```

### Installing S3V2_IDEAS_ESMP package
#### Clone the github repository 
```
cd /where_user_clone_the_S3norm_GitHub/
git clone https://github.com/guanjue/S3V2_IDEAS_ESMP.git
```
#### Install dependency: 
##### If some of dependencies were not installed, user can use the following command to install them
```
###### For python dependencies, they can be installed by the following scripts
pip install --upgrade pip --user
pip install --upgrade numpy --user
pip install --upgrade scipy --user
pip install multiprocess --user

R packages
install.packages('data.table')
install.packages('doParallel')
install.packages('foreach')
install.packages('LSD')
install.packages('pheatmap')


bedtools
######(https://bedtools.readthedocs.io/en/latest/content/installation.html)

IDEAS requires GSL 2.2.1 and python.
The instruction about GSL can be found in the following link:

https://www.gnu.org/software/gsl/manual/gsl-ref.html
Add the ~/gsl/lib into the LD_LIBRARY_PATH


IDEAS also requires UCSC utilities
IDEAS already include the required utilities in the package. But if user is using different system, please replace the UCSC utilities by the version of user's system.
The follow link includes the UCSC utilities for other systems.

http://hgdownload.soe.ucsc.edu/admin/exe/


