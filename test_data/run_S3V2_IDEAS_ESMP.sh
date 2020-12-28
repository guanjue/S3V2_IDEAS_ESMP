### required inputs
###### the absolute path to the bin folder
script_dir='/app/S3V2_IDEAS_ESMP/bin/'
###### your output folder
output_dir='/app/S3V2_IDEAS_ESMP/test_data/outputs/'
###### the absolute path to the your modified "metadata.for_master_peak_calls.txt" file
metadata='/app/S3V2_IDEAS_ESMP/test_data/metadata.forEScall.txt'
###### The output name
id_name='test_S3V2_IDEAS_pipeline'

###### genome
GENOME='mm10'
###### genome size (can be found in the "S3V2_IDEAS_ESMP/genomesize/" folder)
GENOMESIZES='/app/S3V2_IDEAS_ESMP/genomesize/mm10.chrom.chr11.fortest.sizes'
###### blacklist (can be found in the "S3V2_IDEAS_ESMP/blacklist/" folder)
BLACK='/app/S3V2_IDEAS_ESMP/blacklist/mm10-blacklist.v2.bed'

###### number of threads in system
###### When the number of threads is too large, the multi-threads in python3 may fail. So it is more stable to keep it below 4. 
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

[[ -d $output_dir ]] || mkdir $output_dir
cd $output_dir
time python3 $script_dir/S3V2_IDEAS_pipeline.py \
-u $get_sigtrack -v $normalization -y $get_bw -z $run_ideas \
-s $script_dir -o $output_dir -g $GENOME -c $GENOMESIZES -b $BLACK \
-i $metadata -d $id_name -e $email -t $threads -w $IDEAS_track_link -x $other_parafile \
-l $bin_size -n $local_bg_bin -a $cap_sig

