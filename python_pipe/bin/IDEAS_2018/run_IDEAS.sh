###### run IDEAS
######
### cp script in the dir
IDEAS_job_name=run_IDEAS
script_dir=/storage/home/gzx103/group/software/IDEAS/IDEAS_2018/
working_dir=/storage/home/gzx103/group/software/IDEAS/IDEAS_2018/test_data/
output_dir=/storage/home/gzx103/group/software/IDEAS/IDEAS_2018/test_data/run_IDEAS_result/
binfile=mm10_noblacklist_200bin.bin

cd $working_dir
### make output dir
if [ -d $output_dir ]; then rm -r $output_dir; mkdir $output_dir; else mkdir $output_dir; fi
### cp scripts to the analysis dir
if [ -d bin ]; then rm -r bin; cp -r $script_dir'bin' ./ ; else cp -r $script_dir'bin' ./ ; fi
if [ -d data ]; then rm -r data; cp -r $script_dir'data' ./ ; else cp -r $script_dir'data' ./ ; fi
### get genome inv file
time python ./bin/bed2inv.py -i $binfile -o $binfile'.inv'
### run IDEAS
time Rscript ./bin/runme.R $IDEAS_job_name'.input' $IDEAS_job_name'.parafile' $output_dir
### rm tmp files
rm $output_dir*tmp*
### get heatmap
time Rscript bin/get_heatmap.R $output_dir$IDEAS_job_name'.para' FALSE ./bin/createGenomeTracks.R
