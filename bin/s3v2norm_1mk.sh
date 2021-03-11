#!/bin/bash

## print a message and exit
function die () {
        echo $1 >> $log_file
        exit 100
}

#input parameters
#python '+script_dir+'/s3v2norm.pipemultithreads_IDEAS.py'+' -n '+str(threads)+' -e '+str(local_bg_bin)+' -t '+mk[0]+'.file_list.S3V2.txt'+' -k '+mk[0]+'_commonpkfdr01_z.cpk.txt'+' -g '+mk[0]+'_commonpkfdr01_z.cbg.txt'+' -s '+script_dir+' -i '+'F'+' -l '+str(0.0001)
#possible input parameters and their defaults
reference_method='max1'
threads_num=4
exp_win='5'
fdr_thresh='0.1'
rank_lim='0.001'
upperlim='500'
lowerlim='0'
p_method='z'
common_pk_binary='0'
common_bg_binary='0'
allpk_binary='F'
file_list='0'
script_folder='0'
log_file='0'

echo "Reading parameters for s3v2norm.sh"
#read actual named input parameters
while getopts ":n:e:f:l:a:b:p:k:g:i:s:t:r:x:" opt; do
  case $opt in
    n) threads_num="$OPTARG"
    ;;
    e) exp_win="$OPTARG"
    ;;
    f) fdr_thresh="$OPTARG"
    ;;
    l) rank_lim="$OPTARG"
    ;;
    a) upperlim="$OPTARG"
    ;;
    b) lowerlim="$OPTARG"
    ;;
    p) p_method="$OPTARG"
    ;;
    k) common_pk_binary="$OPTARG"
    ;;
    g) common_bg_binary="$OPTARG"
    ;;
    i) allpk_binary="$OPTARG"
    ;;
    s) script_folder="$OPTARG"
    ;;
    t) file_list="$OPTARG"
    ;;
    r) reference_method="$OPTARG"
    ;;
    x) log_file="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done

#check required params stx
if [ "$script_folder" == '0' ] || [ "$file_list" == '0' ] || [ "$log_file" == '0' ]; then
   echo "Error missing parameter, -s -t -x are required" 
   exit 100;
fi

echo 'Get S3norm normalized read counts......' >> $log_file
# cpk_file cbg_file allpk_binary_file?
#output_tar_file celltype.mark.S3V2.bedgraph

while read -r line; do
   #line = input_tar_file input_ref_file input_ct input_mk
   read -a f <<< "$line"
   outfile="${f[2]}.${f[3]}.S3V2.bedgraph"
   command="Rscript ${script_folder}/s3v2norm_IDEAS_1mk.R ${f[0]} ${f[1]} ${outfile} ${fdr_thresh} ${exp_win} F F ${rank_lim} ${upperlim} ${lowerlim} ${p_method} ${common_pk_binary} ${common_bg_binary} ${allpk_binary}"
   echo $command >> ${log_file}
   $command &
   let i++
   if [ $i -ge $threads_num ]; then
      wait
      i=0
   fi
done < $file_list 

wait
echo "Finished s3v2norm.sh" >> $log_file

exit 0

