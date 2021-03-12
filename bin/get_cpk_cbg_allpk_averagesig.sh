mk=$1
script_dir1=$2
meta_data=$3

cat $meta_data | awk -F '\t' -v OFS='\t' -v used_mk=$mk '{if ($2==used_mk) print $1"."$2"."$3".ip.idsort.bedgraph"}' > $mk'.file_list_tmp1'
cat $meta_data | awk -F '\t' -v OFS='\t' -v used_mk=$mk '{if ($2==used_mk) print $1"."$2"."$3".ctrl.idsort.bedgraph"}' > $mk'.file_list_tmp2'
paste $mk'.file_list_tmp1' $mk'.file_list_tmp2' > $mk'.file_list.txt'
### get cpk cbg
time Rscript $script_dir1/get_common_pk_p.R $mk'.file_list.txt' $mk'_commonpkfdr01_z' 0.1 z 0.95 1.0
### get allpk
cat $mk'_commonpkfdr01_z.cbg.txt' | awk '{if ($1!=0) print 0; else print 1}' > $mk'_commonpkfdr01_z.allpk.txt'
### get average signal
time Rscript $script_dir1/get_average_sig.R $mk'.file_list.txt' $mk'.average_sig.bedgraph'

