### required inputs
script_dir='/storage/home/gzx103/group/software/S3V2_IDEAS_ESMP/bin/'
output_dir='/storage/home/gzx103/scratch/S3V2norm_compare/S3V2_cCRE_pipeline_test/'
metadata='/storage/home/gzx103/scratch/S3V2norm_compare/S3V2_cCRE_pipeline_test/input_files/metadata.for_master_peak_calls.txt'
id_name='test_S3V2_IDEAS_cCRE_pipeline'

GENOME='hg38'
GENOMESIZES='/storage/home/gzx103/group/software/S3V2_IDEAS_ESMP/genomesize/hg38.chrom.chr16.fortest.sizes'
BLACK='/storage/home/gzx103/group/software/S3V2_IDEAS_ESMP/blacklist/hg38-blacklist.v2.bed'

### other parameters 
get_sigtrack='T'
normalization='T'
get_bw='T'
run_ideas='T'
threads=4
bin_size=200
local_bg_bin=5
cap_sig=16
email='gzx103@psu.edu'
other_parafile='F'
IDEAS_track_link='http://bx.psu.edu/~gzx103/tmp/'

time python $script_dir/S3V2_IDEAS_pipeline.py \
-u $get_sigtrack -v $normalization -y $get_bw -z $run_ideas \
-s $script_dir -o $output_dir -g $GENOME -c $GENOMESIZES -b $BLACK \
-i $metadata -d $id_name -e $email -t $threads -w $IDEAS_track_link -x $other_parafile \
-l $bin_size -n $local_bg_bin -a $cap_sig

