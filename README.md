# S3V2_IDEAS_ESMP

## In the pipeline, it can first normalize average read counts signal by S3norm ver2 and then use IDEAS to either to do genome segmentations or to call master peaks list across multiple datasets. This pipeline will use bigWig file as input files. 


<img src="https://github.com/guanjue/S3V2_IDEAS_ESMP/blob/master/figures/overall_pipeline.png" width="800"/>

##### Figure 1. The overview of S3V2-IDEAS pipeline. There are two main steps in the S3V2-IDEAS pipeline: (A) the data normalization and denoising step by S3V2 and (B and C) the data integration step by IDEAS. The data integration step has two modes. (B) In the epi-genetic state mode, multiple epigenetic features can be integrated into the epigenetic states model (D). (C) In the signal intensity state mode, the signal of one epigenetic feature can be clustered into different signal intensity states (E). (F) An additional master peak list can be extracted from the signal intensity state tracks in multiple cell types. 

#####################################################################################

## Table of Contents
**[(1) Prerequisites and S3V2_IDEAS_ESMP installation](#Prerequisites-and-S3V2_IDEAS_ESMP-installation)**<br>
#####
**[(2) Inputs for S3V2_IDEAS_ESMP](#Inputs-for-S3V2_IDEAS_ESMP)**<br>
#####
**[(3) How to run S3V2_IDEAS_ESMP pipeline](#How-to-run-S3V2_IDEAS_ESMP-pipeline)**<br>
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
#### python/2.7 (https://www.python.org/downloads/release/python-2716/)
#### python dependencies: 'numpy', 'scipy', 'multiprocess'
#### R (https://www.r-project.org/)
#### R dependencies: 'data.table', 'doParallel', 'foreach', 'LSD', 'pheatmap'
#### gawk

### Installing S3V2_IDEAS_ESMP pipeline
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

IDEAS requires GSL 2.2.1 and python 2.7.
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
##### It only needs one input metadata file which tells the pipeline where are input bigwig files. An example of the metadata is in the "metadata.for_master_peak_calls.txt" file with 4 columns (columns are separated by tab):
##### 1st column: cell type name (!!!The cell type name should not have "." in it!!!)
##### 2nd column: epigenetic feature
##### 3rd column: cell type id
##### 4th column: absolute path to the IP bigwig files
##### 5th column: absolute path to the CONTROL bigwig files (If there is no control signal track, this column can be leave as empty)
```
>>> head metadata.for_master_peak_calls.txt
Cell1	ATAC	01	/storage/home/gzx103/scratch/S3V2norm_compare/hg38_cCREs/bw/c10_chr16_read.ATAC.bw
Cell1	H3K4me3	01	/storage/home/gzx103/scratch/S3V2norm_compare/hg38_EP/bw/c11_chr16_read.H3K4me3.bw
Cell2	ATAC	02	/storage/home/gzx103/scratch/S3V2norm_compare/hg38_cCREs/bw/c11_chr16_read.ATAC.bw
Cell2	H3K4me3	02	/storage/home/gzx103/scratch/S3V2norm_compare/hg38_EP/bw/c11_chr16_read.H3K4me3.bw
```

#####################################################################################

## How to run S3V2_IDEAS_ESMP pipeline
#### Use 'run_S3V2_IDEAS_ESMP.sh' to run S3norm pipeline.
##### After perparing the input data, user just need to set the parameters in "run_S3V2_IDEAS_ESMP.sh" to run S3V2_IDEAS_ESMP.
#####
```
### required inputs
###### the absolute path to the bin folder
script_dir='/storage/home/gzx103/group/software/S3V2_IDEAS_ESMP/bin/'
###### your output folder
output_dir='/storage/home/gzx103/scratch/S3V2norm_compare/S3V2_cCRE_pipeline_test/'
###### the absolute path to the your modified "metadata.for_master_peak_calls.txt" file
metadata='/storage/home/gzx103/scratch/S3V2norm_compare/S3V2_cCRE_pipeline_test/input_files/metadata.for_master_peak_calls.txt'
###### The output name
id_name='test_S3V2_IDEAS_ESMP_pipeline'


###### genome
GENOME='hg38'
###### genome size (can be found in the "S3V2_IDEAS_ESMP/genomesize/" folder)
GENOMESIZES='/storage/home/gzx103/group/software/S3V2_IDEAS_ESMP/genomesize/hg38.chrom.chr16.fortest.sizes'
###### blacklist (can be found in the "S3V2_IDEAS_ESMP/blacklist/" folder)
BLACK='/storage/home/gzx103/group/software/S3V2_IDEAS_ESMP/blacklist/hg38-blacklist.v2.bed'

###### number of threads in system
threads=4
###### bin size of the signal resolution
bin_size=200
###### email address
email='your_email@xxx.edu'

###### User can use the "other_parafile" parameter to incorporate previous epigenetic state model
###### We provided two epigenetic state models with 8/7 epigenetic features that can be found in the "prior_ES_models/" folder
other_parafile='F'

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
ls -ltrh test_S3V2_IDEAS_ESMP_pipeline_IDEAS_output/
```
##### The S3norm normalized average read counts will be save in the "test_S3V2_IDEAS_cCRE_pipeline_bws_RC" folder
##### The -log10(p-value) based on S3norm normalized average read counts will be save in the "test_S3V2_IDEAS_cCRE_pipeline_bws_NBP" folder
##### The signal composition of the epigenetic state will be the "test_S3V2_IDEAS_cCRE_pipeline.pdf"
##### The genome segmentation will be saved in "test_S3V2_IDEAS_ESMP_pipeline_IDEAS_output/Tracks/" folder
##### If there is one epigenetic feature, a master peak list will be saved as the "test_S3V2_IDEAS_ESMP_pipeline.cCRE.M.bed" file



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





