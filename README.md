## S3V2_IDEAS_ESMP: a package for normalizing, denoising and integrating epigenomic datasets across different cell types

#### The S3V2-IDEAS package first normalizes average read counts by S3V2 and then uses IDEAS to either do genome segmentations or call master peaks list across multiple datasets. This package uses bigWig files as the inputs. 


<img src="https://github.com/guanjue/S3V2_IDEAS_ESMP/blob/master/figures/overall_pipeline.png" width="600"/>

##### Figure 1. The overview of S3V2-IDEAS package. There are two main steps in the S3V2-IDEAS package: (A) the data normalization and denoising step by S3V2 and (B and C) the data integration step by IDEAS. The data integration step has two modes. (B) In the epi-genetic state mode, multiple epigenetic features can be integrated into the epigenetic states model (D). (C) In the signal intensity state mode, the signal of one epigenetic feature can be clustered into different signal intensity states (E). (F) An additional master peak list can be extracted from the signal intensity state tracks in multiple cell types. 

#####################################################################################

## Table of Contents
## 
**[(1) Prerequisites and S3V2_IDEAS_ESMP Installation](https://github.com/guanjue/S3V2_IDEAS_ESMP/blob/master/manuals/install.md)**<br>
## 
**[(2) Inputs-ParameterSettings-RunningSteps-Outputs for S3V2_IDEAS_ESMP](https://github.com/guanjue/S3V2_IDEAS_ESMP/blob/master/manuals/inoutput_for_S3V2_IDEAS_pipeline.md)**<br>
## 
**[(3) Running S3V2-IDEAS in operating system set in Docker](https://github.com/guanjue/S3V2_IDEAS_ESMP/blob/master/manuals/run_S3V2_IDEAS_in_Docker.md)**<br>

#####################################################################################
 

## Contacts and References
#### Contacts: 
##### gzx103@psu.edu

#### S3V2-IDEAS
Guanjue Xiang, Belinda M Giardine, Shaun Mahony, Yu Zhang, Ross C Hardison, S3V2-IDEAS: a package for normalizing, denoising and integrating epigenomic datasets across different cell types, Bioinformatics, 2021;, btab148, https://doi.org/10.1093/bioinformatics/btab148
#### S3norm
Guanjue Xiang, Cheryl A Keller, Belinda Giardine, Lin An, Qunhua Li, Yu Zhang, Ross C Hardison, S3norm: simultaneous normalization of sequencing depth and signal-to-noise ratio in epigenomic data, Nucleic Acids Research, Volume 48, Issue 8, 07 May 2020, Page e43, https://doi.org/10.1093/nar/gkaa105
#### IDEAS genome segmentation
Yu Zhang, Lin An, Feng Yue, Ross C Hardison, Jointly characterizing epigenetic dynamics across multiple human cell types, Nucleic Acids Research, Volume 44, Issue 14, 19 August 2016, Pages 6721–6731, https://doi.org/10.1093/nar/gkw278
#### Vision project 
###### IDEAS analysis with 8 epigenetic features in 20 mouse hematopoietic cell types & Updates in IDEAS
Guanjue, Xiang, Cheryl A. Keller, ..., Yu Zhang, Ross C. Hardison. "An integrative view of the regulatory and transcriptional landscapes in mouse hematopoiesis." Genome research. 2020 Mar 1;30(3):472-84, https://doi.org/10.1101/gr.255760.119





