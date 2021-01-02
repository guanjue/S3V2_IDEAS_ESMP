# S3V2_IDEAS_ESMP

## In the package, it can first normalize average read counts signal by S3norm ver2 and then use IDEAS to either to do genome segmentations or to call master peaks list across multiple datasets. This package will use bigWig file as input files. 


<img src="https://github.com/guanjue/S3V2_IDEAS_ESMP/blob/master/figures/overall_pipeline.png" width="600"/>

##### Figure 1. The overview of S3V2-IDEAS package. There are two main steps in the S3V2-IDEAS package: (A) the data normalization and denoising step by S3V2 and (B and C) the data integration step by IDEAS. The data integration step has two modes. (B) In the epi-genetic state mode, multiple epigenetic features can be integrated into the epigenetic states model (D). (C) In the signal intensity state mode, the signal of one epigenetic feature can be clustered into different signal intensity states (E). (F) An additional master peak list can be extracted from the signal intensity state tracks in multiple cell types. 

#####################################################################################

## Table of Contents
## 
**[(1) Prerequisites and S3V2_IDEAS_ESMP Installation](https://github.com/guanjue/S3V2_IDEAS_ESMP/blob/master/install.md)**<br>
## 
**[(2) Inputs and Outputs for S3V2_IDEAS_ESMP](https://github.com/guanjue/S3V2_IDEAS_ESMP/blob/master/inoutput_for_S3V2_IDEAS_pipeline.md)**<br>
## 
**[(3) Running S3V2-IDEAS in operating system set in Docker](https://github.com/guanjue/S3V2_IDEAS_ESMP/blob/master/run_S3V2_IDEAS_in_Docker.md)**<br>

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





