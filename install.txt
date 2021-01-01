## Prerequisites and S3V2_IDEAS_ESMP installation

### S3V2_IDEAS_ESMP dependencies are as follows:
```
python3 (https://www.python.org/download/releases/3.0/)
python dependencies: 'numpy', 'scipy', 'multiprocess'
R (https://www.r-project.org/)
R dependencies: 'data.table', 'doParallel', 'foreach', 'LSD', 'pheatmap'
gawk
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


