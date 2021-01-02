# S3V2_IDEAS_ESMP

## In the package, it can first normalize average read counts signal by S3norm ver2 and then use IDEAS to either to do genome segmentations or to call master peaks list across multiple datasets. This package will use bigWig file as input files. 


<img src="https://github.com/guanjue/S3V2_IDEAS_ESMP/blob/master/figures/overall_pipeline.png" width="800"/>

##### Figure 1. The overview of S3V2-IDEAS package. There are two main steps in the S3V2-IDEAS package: (A) the data normalization and denoising step by S3V2 and (B and C) the data integration step by IDEAS. The data integration step has two modes. (B) In the epi-genetic state mode, multiple epigenetic features can be integrated into the epigenetic states model (D). (C) In the signal intensity state mode, the signal of one epigenetic feature can be clustered into different signal intensity states (E). (F) An additional master peak list can be extracted from the signal intensity state tracks in multiple cell types. 

#####################################################################################

## Table of Contents
**[(1) Prerequisites and S3V2_IDEAS_ESMP installation](https://github.com/guanjue/S3V2_IDEAS_ESMP/blob/master/install.md)**<br>
**[(2) Inputs for S3V2_IDEAS_ESMP](https://github.com/guanjue/S3V2_IDEAS_ESMP/blob/master/input_for_S3V2_IDEAS_pipeline.md)**<br>
**[(2) Outputs for S3V2_IDEAS_ESMP](https://github.com/guanjue/S3V2_IDEAS_ESMP/blob/master/output_for_S3V2_IDEAS_pipeline.md)**<br>


#####################################################################################

## Outputs of S3V2_IDEAS_ESMP
### All outputs will be saved in the user provided "$output_dir".
##### The epigenetic state genome segmentation (multiple epigenetic features) and master peak list (one epigenetic feature) will be saved in the following folder: "your_output_name_IDEAS_output/"
```
ls -ltrh test_S3V2_IDEAS_pipeline_IDEAS_output/
total 21M
-rw-rw---- 1 gzx103 gzx103_collab 2.3K Sep  8 04:00 test_S3V2_IDEAS_pipeline.profile0
-rw-rw---- 1 gzx103 gzx103_collab 2.0K Sep  8 04:00 test_S3V2_IDEAS_pipeline.para0
-rw-rw---- 1 gzx103 gzx103_collab  22K Sep  8 04:02 test_S3V2_IDEAS_pipeline.cluster
-rw-rw---- 1 gzx103 gzx103_collab  21M Sep  8 04:02 test_S3V2_IDEAS_pipeline.state
-rw-rw---- 1 gzx103 gzx103_collab 2.1K Sep  8 04:02 test_S3V2_IDEAS_pipeline.para
-rw-rw---- 1 gzx103 gzx103_collab 1.1K Sep  8 04:02 test_S3V2_IDEAS_pipeline.profile
-rw-rw---- 1 gzx103 gzx103_collab 4.9K Sep  8 04:02 test_S3V2_IDEAS_pipeline.pdf
drwxrws--- 2 gzx103 gzx103_collab 4.0K Sep  8 04:02 Tracks
-rw-rw---- 1 gzx103 gzx103_collab 1.8K Sep  8 04:02 log.txt

##### (1) The S3V2 normalized average read counts will be save in the "test_S3V2_IDEAS_package_bws_RC" folder
##### (2) The -log10(p-value) based on S3norm normalized average read counts will be save in the "test_S3V2_IDEAS_package_bws_NBP" folder
##### (3) The signal composition of the epigenetic state will be the "test_S3V2_IDEAS_package.pdf"
##### (4) The genome segmentation will be saved in "test_S3V2_IDEAS_package_IDEAS_output/Tracks/" folder. These bigBed files can be loaded into UCSC genome browser as track hub by the "hub_test_S3V2_IDEAS_pipeline.txt" file in the "Tracks/" folder.
##### (5) If there is one epigenetic feature, a master peak list will be saved as the "test_S3V2_IDEAS_package.cCRE.M.bed" file
```


#####################################################################################

## Contacts and References
#### Contacts: 
##### gzx103@psu.edu

#### S3V2-IDEAS
Guanjue Xiang, Belinda M. Giardine, Shaun Mahony, Yu Zhang, Ross C Hardison. "Snapshot: clustering and visualizing epigenetic history during cell differentiation." bioRxiv (2020): .
#### S3norm
Guanjue, Xiang, Cheryl A. Keller, ..., Yu Zhang, Ross C. Hardison. "S3norm: simultaneous normalization of sequencing depth and signal-to-noise ratio in epigenomic data." Nucleic Acids Research 48.8 (2020): e43-e43.
#### IDEAS
Yu Zhang, ..., Ross C. Hardison. (2016). "Jointly characterizing epigenetic dynamics across multiple human cell types." Nucleic acids research, 44(14), 6721-6731.
#### Vision paper
Guanjue, Xiang, Cheryl A. Keller, ..., Yu Zhang, Ross C. Hardison. "An integrative view of the regulatory and transcriptional landscapes in mouse hematopoiesis." Genome Research (2020).





