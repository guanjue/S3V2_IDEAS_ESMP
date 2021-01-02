
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
##### Users just need to change abosolute path in the follow parameters to run the package on the test dataset: "script_dir=, output_dir=, metadata=, GENOMESIZES=, BLACK=".
##### Users also need to change the abosolute path to the bigWig files in the "metadata.forEScall.txt" file.
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
#####################################################################################
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