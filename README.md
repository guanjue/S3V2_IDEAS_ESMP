# S3V2_IDEAS_ESMP

## In the package, it can first normalize average read counts signal by S3norm ver2 and then use IDEAS to either to do genome segmentations or to call master peaks list across multiple datasets. This package will use bigWig file as input files. 


<img src="https://github.com/guanjue/S3V2_IDEAS_ESMP/blob/master/figures/overall_package.png" width="800"/>

##### Figure 1. The overview of S3V2-IDEAS package. There are two main steps in the S3V2-IDEAS package: (A) the data normalization and denoising step by S3V2 and (B and C) the data integration step by IDEAS. The data integration step has two modes. (B) In the epi-genetic state mode, multiple epigenetic features can be integrated into the epigenetic states model (D). (C) In the signal intensity state mode, the signal of one epigenetic feature can be clustered into different signal intensity states (E). (F) An additional master peak list can be extracted from the signal intensity state tracks in multiple cell types. 

#####################################################################################

## Table of Contents
**[(1) Prerequisites and S3V2_IDEAS_ESMP installation](#Prerequisites-and-S3V2_IDEAS_ESMP-installation)**<br>
#####
**[(2) Inputs for S3V2_IDEAS_ESMP](#Inputs-for-S3V2_IDEAS_ESMP)**<br>
#####
**[(3) How to run S3V2_IDEAS_ESMP package](#How-to-run-S3V2_IDEAS_ESMP-package)**<br>
#####
**[(4) Outputs of S3V2_IDEAS_ESMP](#Outputs-of-S3V2_IDEAS_ESMP)**<br>
#####
**[(5) Contacts and References](#Contacts-and-References)**<br>
#####

#####################################################################################

python packages
pip install numpy --user
pip install scipy --user
pip install multiprocess --user

R packages
install.packages('data.table')


## Prerequisites and S3V2_IDEAS_ESMP installation
### S3V2_IDEAS_ESMP dependencies are as follows:
#### python3 (https://www.python.org/download/releases/3.0/)
#### python dependencies: 'numpy', 'scipy', 'multiprocess'
#### R (https://www.r-project.org/)
#### R dependencies: 'data.table', 'doParallel', 'foreach', 'LSD', 'pheatmap'
#### gawk

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


```


#####################################################################################

## Inputs for S3V2_IDEAS_ESMP
#### (1) The input metadata list for S3norm
##### It only needs one input metadata file which tells the package where are input bigwig files. An example of the metadata is in the "metadata.for_master_peak_calls.txt" file with 4 columns (columns are separated by tab):
##### 1st column: cell type name (!!!The cell type name should not have "." in it!!!)
##### 2nd column: epigenetic feature
##### 3rd column: cell type id
##### 4th column: absolute path to the IP bigwig files
##### 5th column: absolute path to the CONTROL bigwig files (If there is no control signal track, this column can be leave as empty)
```
>>> head metadata.for_master_peak_calls.txt
CMP	H3K27ac	r1	/storage/home/gzx103/scratch/test_S3V2/S3V2_IDEAS_ESMP/test_data/input_bw_files/CMP.H3K27ac.r1.chr11.bw	/storage/home/gzx103/scratch/test_S3V2/S3V2_IDEAS_ESMP/test_data/input_bw_files/CMP.H3K27ac.r1.chr11.ctrl.bw
CMP	H3K4me1	r1	/storage/home/gzx103/scratch/test_S3V2/S3V2_IDEAS_ESMP/test_data/input_bw_files/CMP.H3K4me1.r1.chr11.bw	/storage/home/gzx103/scratch/test_S3V2/S3V2_IDEAS_ESMP/test_data/input_bw_files/CMP.H3K4me1.r1.chr11.ctrl.bw
CMP	H3K4me3	r1	/storage/home/gzx103/scratch/test_S3V2/S3V2_IDEAS_ESMP/test_data/input_bw_files/CMP.H3K4me3.r1.chr11.bw	/storage/home/gzx103/scratch/test_S3V2/S3V2_IDEAS_ESMP/test_data/input_bw_files/CMP.H3K4me3.r1.chr11.ctrl.bw
ERY_fl	H3K27ac	r1	/storage/home/gzx103/scratch/test_S3V2/S3V2_IDEAS_ESMP/test_data/input_bw_files/ERY_fl.H3K27ac.r1.chr11.bw	/storage/home/gzx103/scratch/test_S3V2/S3V2_IDEAS_ESMP/test_data/input_bw_files/ERY_fl.H3K27ac.r1.chr11.ctrl.bw
ERY_fl	H3K4me1	r1	/storage/home/gzx103/scratch/test_S3V2/S3V2_IDEAS_ESMP/test_data/input_bw_files/ERY_fl.H3K4me1.r1.chr11.bw	/storage/home/gzx103/scratch/test_S3V2/S3V2_IDEAS_ESMP/test_data/input_bw_files/ERY_fl.H3K4me1.r1.chr11.ctrl.bw
ERY_fl	H3K4me3	r1	/storage/home/gzx103/scratch/test_S3V2/S3V2_IDEAS_ESMP/test_data/input_bw_files/ERY_fl.H3K4me3.r1.chr11.bw	/storage/home/gzx103/scratch/test_S3V2/S3V2_IDEAS_ESMP/test_data/input_bw_files/ERY_fl.H3K4me3.r1.chr11.ctrl.bw
```
#####################################################################################

## How to run S3V2_IDEAS_ESMP package
#### We prepared a testing dataset in the "test_data/" folder.
#### User can use the 'test_data/run_S3V2_IDEAS_ESMP.sh' script to run S3V2_IDEAS_ESMP package on these datasets.
##### After perparing the input data, user just need to set the parameters in "run_S3V2_IDEAS_ESMP.sh" to run S3V2_IDEAS_ESMP.
##### Users just need to change abosolute path in the follow parameters to run the package on the test dataset: "script_dir=, output_dir=, metadata=, GENOMESIZES=, BLACK="
##### For the testing datasets, they include three histone marks (H3K4me3, H3K4me1, H3K27ac) in two cell types (CMP and ERY_fl) in the VISION project (http://usevision.org). The S3V2-IDEAS package can finish running on this set of data within 2.5 hours by using 50GB with 4 threads.
```
### required inputs
###### the absolute path to the bin folder
script_dir='/storage/home/gzx103/scratch/test_S3V2/S3V2_IDEAS_ESMP/bin/'
###### your output folder
output_dir='/storage/home/gzx103/scratch/test_S3V2/S3V2_IDEAS_ESMP/test_data/'
###### the absolute path to the your modified "metadata.for_master_peak_calls.txt" file
metadata='/storage/home/gzx103/scratch/test_S3V2/S3V2_IDEAS_ESMP/test_data/metadata.forEScall.txt'
###### The output name
id_name='test_S3V2_IDEAS_package'

###### genome
GENOME='mm10'
###### genome size (can be found in the "S3V2_IDEAS_ESMP/genomesize/" folder)
GENOMESIZES='/storage/home/gzx103/scratch/test_S3V2/S3V2_IDEAS_ESMP/genomesize/mm10.chrom.chr11.fortest.sizes'
###### blacklist (can be found in the "S3V2_IDEAS_ESMP/blacklist/" folder)
BLACK='/storage/home/gzx103/scratch/test_S3V2/S3V2_IDEAS_ESMP/blacklist/mm10-blacklist.v2.bed'

###### number of threads in system
###### When the number of threads is too large, the multi-threads in python may fail. So it is more stable to keep it below 4. 
threads=4
###### bin size of the signal resolution
bin_size=200
###### email address
email='your_email@xxx.edu'

###### other parameters 
get_sigtrack='T'
normalization='T'
get_bw='T'
run_ideas='T'
local_bg_bin=5
cap_sig=16
### User can use the "other_parafile" parameter to incorporate previous epigenetic state model
### We provided two epigenetic state models with 8/7 epigenetic features that can be found in the "prior_ES_models/" folder
other_parafile='F'
IDEAS_track_link='http://your_acess_link_that_can_be_used_for_track_hub_in_genome_browser/'

cd $output_dir
time python $script_dir/S3V2_IDEAS_package.py \
-u $get_sigtrack -v $normalization -y $get_bw -z $run_ideas \
-s $script_dir -o $output_dir -g $GENOME -c $GENOMESIZES -b $BLACK \
-i $metadata -d $id_name -e $email -t $threads -w $IDEAS_track_link -x $other_parafile \
-l $bin_size -n $local_bg_bin -a $cap_sig

```

##### The package will generated several whole genome matrix which require large memory such as:
##### For our analysis in 20 cell types in mouse, we usually ask 50GB memories and 4 threads to run the package.
```
#PBS -l nodes=1:ppn=4
#PBS -l walltime=50:00:00
#PBS -j oe
#PBS -l pmem=50gb
```

##### Then Run:
```
time bash run_S3V2_IDEAS_ESMP.sh
```

#####################################################################################

## Outputs of S3V2_IDEAS_ESMP
### All outputs will be saved in the user provided "$output_dir".
##### The epigenetic state genome segmentation (multiple epigenetic features) and master peak list (one epigenetic feature) will be saved in the following folder: "your_output_name_IDEAS_output/"
```
ls -ltrh test_S3V2_IDEAS_package_IDEAS_output/
```
##### The S3norm normalized average read counts will be save in the "test_S3V2_IDEAS_package_bws_RC" folder
##### The -log10(p-value) based on S3norm normalized average read counts will be save in the "test_S3V2_IDEAS_package_bws_NBP" folder
##### The signal composition of the epigenetic state will be the "test_S3V2_IDEAS_package.pdf"
##### The genome segmentation will be saved in "test_S3V2_IDEAS_package_IDEAS_output/Tracks/" folder
##### If there is one epigenetic feature, a master peak list will be saved as the "test_S3V2_IDEAS_package.cCRE.M.bed" file



#####################################################################################

## Contacts and References
#### Contacts: 
##### gzx103@psu.edu

#### S3V2-IDEAS
Manuscript in preparation
#### S3norm
Guanjue, Xiang, Cheryl A. Keller, ..., Yu Zhang, Ross C. Hardison. "S3norm: simultaneous normalization of sequencing depth and signal-to-noise ratio in epigenomic data." Nucleic Acids Research 48.8 (2020): e43-e43.
#### IDEAS
Yu Zhang, ..., Ross C. Hardison. (2016). "Jointly characterizing epigenetic dynamics across multiple human cell types." Nucleic acids research, 44(14), 6721-6731.
#### Vision paper
Guanjue, Xiang, Cheryl A. Keller, ..., Yu Zhang, Ross C. Hardison. "An integrative view of the regulatory and transcriptional landscapes in mouse hematopoiesis." Genome Research (2020).





