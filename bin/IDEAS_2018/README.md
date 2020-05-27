# IDEAS: Integrative and Discriminative Epigenome Annotation System

#### Advanced sequencing technologies have generated a plethora of data for many chromatin marks in multiple tissues and cell types, yet there is lack of a generalized tool for optimal utility of those data. A major challenge is to quantitatively model the epigenetic dynamics across both the genome and many cell types for understanding their impacts on differential gene regulation and disease. We introduce IDEAS, an integrative and discriminative epigenome annotation system, for jointly characterizing epigenetic landscapes in many cell types and detecting differential regulatory regions. A key distinction between our method and existing state-of-the-art algorithms is that IDEAS integrates epigenomes of many cell types simultaneously in a way that preserves the position-dependent and cell type-specific information at fine scales, thereby greatly improving segmentation accuracy and producing comparable annotations across cell types. 
#### (Zhang, Yu, Lin An, Feng Yue, and Ross C. Hardison. "Jointly characterizing epigenetic dynamics across multiple human cell types." Nucleic acids research 44, no. 14 (2016): 6721-6731.)


<img src="https://github.com/guanjue/IDEAS_2018/blob/master/example_figures/f1_IDEAS_mechanism.png" width="800"/>

##### Figure 1. Illustration of IDEAS model. The IDEAS method borrows locus specific information across cell types to improve accuracy, and simultaneously accounts for local cell type relationships for inferring cell type-specific activities. In particular, to infer epigenetic state at a given locus in a target cell type, IDEAS uses the currently inferred states in other cell types at the same locus, but only those cell types showing similar local epigenetic landscapes with the target cell type, as priors to improve inference. The local window (dashed box) is dynamically determined by Markov chains, and cell types are clustered within the local window for their relationships with the target cell. The entire process is iterative with all the unknowns (epigenetic states and local cell type clustering) updated until convergence. The final segmentation is then colored using an automatic coloring script for visualization in browser.



<img src="https://github.com/guanjue/IDEAS_2018/blob/master/example_figures/f2_roadmap_result.png" width="800"/>

##### Figure 2. Inferred chromatin states in 127 cell types. (A) Mean epigenetic signal in the IDEAS inferred states (red labeled) and the ChromHMM inferred states (black labeled in brackets). Color key for each state is shown under the heatmap. Percentage of each state in the genome is shown in parenthesis. IDEAS states that do not have a one-to-one mapping with ChromHMM’s states are marked by asterisk. (B) Reproducibility of segmentation by IDEAS between three independent runs using the original program (blue) and the proposed training pipeline (yellow). Each box shows the agreement of segmentation between two runs, measured by adjusted rand index between the inferred chromatin states within matched cell types. Adjusted rand index is a standardized statistics of similarity between two clustering results, which corrects for chance and accounts for different numbers of clusters. (C) Segmentation example by IDEAS and ChromHMM in 127 cell types at genes CIITA and CLEC16A. Blowups highlight some differences between the two maps. Color keys of chromatin states are defined in (A). From: Accurate and reproducible functional maps in 127 human cell types via 2D genome segmentation Nucleic Acids Res. 2017;45(17):9823-9836. doi:10.1093/nar/gkx659


## Install IDEAS
#### Clone the github repository 
```
git clone https://github.com/guanjue/IDEAS_2018.git
```
###### IDEAS requires GSL 2.2.1 and python 2.7.
###### The instruction about GSL can be found in the following link:
https://www.gnu.org/software/gsl/manual/gsl-ref.html
###### Add the ~/gsl/lib into the LD_LIBRARY_PATH

###### IDEAS also requires UCSC utilities
###### IDEAS already include the required utilities in the package. But if user is using different system, please replace the UCSC utilities by the version of user's system.
###### The follow link includes the UCSC utilities for other systems.
http://hgdownload.soe.ucsc.edu/admin/exe/


## Input data
##### The parameter file for IDEAS. 
```
run_IDEAS.parafile
>>> head -100 run_IDEAS.parafile 
id= test_IDEAS          #job id, also used as output file names
email= giardine@bx.psu.edu
thread= 32              #number of threads to be used for parallelization

prepmat= 0              #1: preprocess data, 0: for data already processed for ideas
build= mm10             #hg19, hg38, mm9, mm10, not used if bedfile is specified
prenorm= 0              #1: normalize data (assumed 100Million reads in total), 0: do not normalize
bed= mm10.noblack_list.bin      #user specified windows
sig= mean               #mean: mean signal per window, max: max signal per window

ideas= 1                #1: run ideas, 0: not run ideas
train= 50                #number of random starts, used to select states, 0: no training
trainsz= 500000
log2= 0                 #take log2(x+num), 0: do not take log2
cap= 16                 #maximum signal is capped at 16
norm= 0                 #1: standardize by mean and std, 0: no normalization
num_state= 0            #specify number of states for the model, 0: let program determine
num_start= 100          #specify number of states at the initialization stage
minerr= 0.5             #minimum standard deviation in each state, usually between (0,1]
#otherpara= /gpfs/group/yzz2/default/scratch/roadmap_analysis/impute/bin_12mark_1e-4.para0
smooth= 0               #make states more homogeneous along genome, 0: original ideas
burnin= 20              #number of burnins, include both sampling and maximization
sample= 5               #number of steps for maximization, 1 may be fine
#split= mm10.noblack_list.bin.inv    #specify an interval file, ideas will run on different intervals separately. The name of interval file is $bed'.inv'
impute= None            #specify which marks to be imputed; or All or None
maketrack= 1            #1: make custom tracks for browser visual, 0: no tracks
#statefiles= /storage/home/gzx103/scratch/gtex_encode/bams/entex_data_output_0_16lim_ideas_01/ideas_state_filelist.txt  #only needed if ideas was not run; separate file names by ","
#hubURL= "http://bx.psu.edu/~yuzhang/tmp/"      #URL where the custom tracks will be stored
#mycolor= 255,0,0;255,255,0;0,255,0;0,0,255;50,50,50    #rgb color for each mark, semicolon delimited
#statecolor= /gpfs/group/yzz2/default/scratch/roadmap_analysis/impute/statecolort.txt                   #rgb color of each state
#statename= statename.txt               #state names
#cellinfo= cellinfo.txt #cell type information, order of cell types will be the same in browser, 4 columns: cell type id as shown in state files, cell type short label to be shown in browser, cell type long label, cell type text color
```
##### Usually, user only needs to change the following parameters in the parameter file:
```
thread= 32				#number of threads to be used for parallelization
build= mm10				#hg19, hg38, mm9, mm10, not used if bedfile is specified
bed= mm10.noblack_list.bin		#user specified windows
split= mm10.noblack_list.bin.inv	#specify an interval file, ideas will run on different intervals separately
cap= 16					#maximum signal is capped at 16
impute= None                            #specify which marks to be imputed; or All or None; If user wants to keep the imputed signal, set it as 'All'
```


##### The bin file for IDEAS: each column is separated by whitespace
###### 1st column: chromosome; 
###### 2nd column: bin start coordinate; 
###### 3rd column: bin end coordinate; 
###### 4th column: bin id 
```
mm10.noblack_list.bin
>>> head mm10.noblack_list.bin
chr1 0 200 R1
chr1 200 400 R2
chr1 400 600 R3
chr1 600 800 R4
chr1 800 1000 R5
chr1 1000 1200 R6
chr1 1200 1400 R7
chr1 1400 1600 R8
chr1 1600 1800 R9
chr1 1800 2000 R10
......
```

##### The input file list: each column is separated by whitespace
###### 1st column: cell type name; 
###### 2nd column: mark name; 
###### 3rd column: input file and its absolution path
```
run_IDEAS.input
>>> head run_IDEAS.input 
ERY_ad atac /storage/home/gzx103/group/software/IDEAS/IDEAS_2018/test_data/run_IDEAS_input/ERY_ad.atac.1M.txt
MEP atac /storage/home/gzx103/group/software/IDEAS/IDEAS_2018/test_data/run_IDEAS_input/MEP.atac.1M.txt
ERY_ad h3k27ac /storage/home/gzx103/group/software/IDEAS/IDEAS_2018/test_data/run_IDEAS_input/ERY_ad.h3k27ac.1M.txt
MEP h3k27ac /storage/home/gzx103/group/software/IDEAS/IDEAS_2018/test_data/run_IDEAS_input/MEP.h3k27ac.1M.txt
......
```

##### The signal file: one signal column 
###### 1st column: signal of the mark in the cell
###### User can convert the bigwig file into the signal file by using the ucsc bigWigAverageOverBed utility. 
#####!!!After using the bigWigAverageOverBed, the row should be reordered so that the order of rows match with the mm10.noblack_list.bin!!!
###### The UCSC utility can be downloaded from (http://hgdownload.soe.ucsc.edu/admin/exe/)
```
cell_mark=NK_atac
bigWigAverageOverBed $cell_mark'.bw' whole_genome_bin.bed $cell_mark'.bw.tab'
sort -k1,1 $cell_mark'.bw.tab' | cut -f5 > $cell_mark'.bw.tab.sig'
``` 
```
>>> head /storage/home/gzx103/group/software/IDEAS/IDEAS_2018/test_data/run_IDEAS_input/ERY_ad.atac.1M.txt
0
0
0
0
0
......
```


## Run IDEAS
##### (1) copy the 'run_IDEAS.sh' & 'run_IDEAS.parafile' into the working directory
```
cp ~/group/software/IDEAS/IDEAS_2018/run_IDEAS.sh working_dir/
cp ~/group/software/IDEAS/IDEAS_2018/run_IDEAS.parafile working_dir/
cp ~/group/software/IDEAS/IDEAS_2018/run_IDEAS.input working_dir/
<<<<<<< HEAD
=======
# also make sure all the paths in the "run_IDEAS.input" are absolute path.
>>>>>>> 01138efd09dcd5db00916be265069eba00889ffe

```
##### (2) change the cell type name, mark name and file location in the 'run_IDEAS.input' file
##### The input file list: each column is separated by whitespace
###### 1st column: cell type name; 
###### 2nd column: mark name; 
###### 3rd column: input file and its absolution path
```
run_IDEAS.input
>>> head run_IDEAS.input
ERY_ad atac /storage/home/gzx103/group/software/IDEAS/IDEAS_2018/test_data/run_IDEAS_input/ERY_ad.atac.1M.txt
MEP atac /storage/home/gzx103/group/software/IDEAS/IDEAS_2018/test_data/run_IDEAS_input/MEP.atac.1M.txt
ERY_ad h3k27ac /storage/home/gzx103/group/software/IDEAS/IDEAS_2018/test_data/run_IDEAS_input/ERY_ad.h3k27ac.1M.txt
MEP h3k27ac /storage/home/gzx103/group/software/IDEAS/IDEAS_2018/test_data/run_IDEAS_input/MEP.h3k27ac.1M.txt
......
```

##### (3) change the following parameters in the 'run_IDEAS.sh' file:
###### script_dir='absolute path to the IDEAS_2018 dir'
###### output_dir='absolute path to the output directory'
###### binfile='name of bin file'
```
>>> head -100 run_IDEAS.sh 
###### run IDEAS
######
### cp script in the directory
IDEAS_job_name=run_IDEAS
script_dir=/storage/home/gzx103/group/software/IDEAS/IDEAS_2018/
output_dir=/storage/home/gzx103/group/software/IDEAS/IDEAS_2018/test_data/run_IDEAS_result/
<<<<<<< HEAD
binfile=mm10.noblack_list.bin
=======
binfile=mm10.noblack_list.bin	#(Absolute path is required if file isn't under the working direactory)
>>>>>>> 01138efd09dcd5db00916be265069eba00889ffe

### make output directory
mkdir -p $output_dir
### cp scripts to the working directory
cp -r $script_dir'bin' ./
cp -r $script_dir'data' ./
### get genome inv file
time python $script_dir'bin/bed2inv.py' -i $binfile -o $binfile'.inv'
### run IDEAS
time Rscript bin/runme.R run_IDEAS.input run_IDEAS.parafile $output_dir
### rm tmp files
rm $output_dir*tmp*
### get heatmap
time Rscript bin/get_heatmap.R $output_dir$IDEAS_job_name'.para0' FALSE ./bin/createGenomeTracks.R
```

##### (4) change the following parameters in the parameter file:
```
thread= 32				#number of threads to be used for parallelization
build= mm10				#hg19, hg38, mm9, mm10, not used if bedfile is specified
bed= mm10.noblack_list.bin		#user specified windows. (Absolute path is required if file isn't under the working direactory)
```

##### (5) use 'run_IDEAS.sh' script to run IDEAS
```
time bash run_IDEAS.sh
```



## Output results for test data
### All output files will be saved to the following directory:
```
output_dir=/storage/home/gzx103/group/software/IDEAS/IDEAS_2018/test_data/run_IDEAS_result/
```

## The heatmap for IDEAS epigenetic state
<img src="https://github.com/guanjue/IDEAS_2018/blob/master/example_figures/f3_vision_result.png" width="800"/>

##### Figure 3. The epigenetic state inferred by IDEAS. Each row represents one epigenetic state. Each column represents one epigenetic mark. The color density represent the average signal of all genome loci with the corresponding epigenetic state. The dark blue represents high average signal. The white represent the low average signal. 

## The genome browser track (bigbed format) for IDEAS epigenetic state will be saved in the subdirectory (named as Tracks/) in the output directory: 
```
track_dir=/storage/home/gzx103/group/software/IDEAS/IDEAS_2018/test_data/run_IDEAS_result/Tracks/
```

## References

##### Zhang, Yu, Lin An, Feng Yue, and Ross C. Hardison. "Jointly characterizing epigenetic dynamics across multiple human cell types." Nucleic acids research 44, no. 14 (2016): 6721-6731.
##### Zhang, Yu, and Ross C. Hardison. "Accurate and reproducible functional maps in 127 human cell types via 2D genome segmentation." Nucleic acids research 45, no. 17 (2017): 9823-9836.


