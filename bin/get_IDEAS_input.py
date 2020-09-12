#module load python/2.7
import os
from subprocess import call

################################################################################################
### s3norm
def get_IDEAS_input(id_input, script_dir, working_dir, output_dir, binfile, email_input, thread_input, build_input, cap_input, hubURL_input, otherpara_input, uniq_mk_num):
	### write .sh file
	sh_file_name = working_dir+'/'+id_input+'.sh'
	sh_file=open(sh_file_name,'w')
	sh_file.write('IDEAS_job_name='+id_input+'\n')
	sh_file.write('script_dir='+script_dir+'\n')
	sh_file.write('working_dir='+working_dir+'\n')
	sh_file.write('output_dir='+output_dir+'\n')
	sh_file.write('binfile='+binfile+'\n')
	sh_file.write('\n')
	sh_file.write('cd $working_dir'+'\n')
	sh_file.write('if [ -d $output_dir ]; then rm -r $output_dir; mkdir $output_dir; else mkdir $output_dir; fi'+'\n')
	sh_file.write('if [ -d bin ]; then rm -r bin; cp -r $script_dir/bin ./ ; else cp -r $script_dir/bin ./ ; fi'+'\n')
	sh_file.write('if [ -d data ]; then rm -r data; cp -r $script_dir/data/ ./ ; else cp -r $script_dir/data ./ ; fi'+'\n')
	sh_file.write('time Rscript ./bin/runme.R $IDEAS_job_name\'.input\' $IDEAS_job_name\'.parafile\' $output_dir'+'\n')
	sh_file.close()

	### write .parafile file
	parafile_file_name = working_dir+'/'+id_input+'.parafile'
	parafile_file=open(parafile_file_name,'w')
	parafile_file.write('id= '+id_input+'\n')
	parafile_file.write('email= '+email_input+'\n')
	parafile_file.write('thread= '+thread_input+'\n')
	parafile_file.write('prepmat= 0'+'\n')
	parafile_file.write('build= '+build_input+'\n')
	parafile_file.write('prenorm= 0'+'\n')
	parafile_file.write('bed= '+binfile+'\n')
	parafile_file.write('sig= mean'+'\n')
	parafile_file.write('ideas= 1'+'\n')
	parafile_file.write('train= 50'+'\n')
	parafile_file.write('trainsz= 500000'+'\n')
	parafile_file.write('log2= 0'+'\n')
	parafile_file.write('cap= '+cap_input+'\n')
	parafile_file.write('norm= 0'+'\n')
	if (uniq_mk_num=='1'):
		parafile_file.write('num_state= 4'+'\n')
		parafile_file.write('num_start= 4'+'\n')
	else:
		parafile_file.write('num_state= 0'+'\n')
		parafile_file.write('num_start= 100'+'\n')
	parafile_file.write('minerr= 0.5'+'\n')
	parafile_file.write('smooth= 0'+'\n')
	parafile_file.write('burnin= 20'+'\n')
	parafile_file.write('sample= 5'+'\n')
	parafile_file.write('impute= None'+'\n')
	parafile_file.write('maketrack= 1'+'\n')
	parafile_file.write('hubURL= "'+hubURL_input+'"\n')
	if otherpara_input!='F':
		parafile_file.write('otherpara= '+otherpara_input+'\n')
	if (uniq_mk_num=='1'):
		parafile_file.write('otherpara= '+script_dir+'/prior/'+'bp24_v2_cCRE_max1.para'+'\n')
	parafile_file.close()

############################################################################
### time python bin/get_IDEAS_input.py -i run_IDEAS -s '/storage/home/gzx103/group/software/IDEAS/IDEAS_2018/' -w '/storage/home/gzx103/group/software/IDEAS/IDEAS_2018/test_data/' -o '/storage/home/gzx103/group/software/IDEAS/IDEAS_2018/test_data/run_IDEAS_result/' -b 'mm10_noblacklist_200bin.bin' -e 'gzx103' -t 8 -g hg38 -c 16 -u http://bx.psu.edu/~gzx103/tmp/ -p F

import getopt
import sys
def main(argv):
	### read user provided parameters
	opts, args = getopt.getopt(argv,"hi:s:w:o:b:e:t:g:c:u:p:n:")
	for opt,arg in opts:
		if opt=="-h":
			print('time python get_IDEAS_input.py -i id_input -s script_dir -w working_dir -o output_dir -b binfile -e email_input -t thread_input -g build_input -c cap_input -u hubURL_input -p otherpara_input')
			return()
		elif opt=="-i":
			id_input=str(arg.strip())
		elif opt=="-s":
			script_dir=str(arg.strip())
		elif opt=="-w":
			working_dir=str(arg.strip())
		elif opt=="-o":
			output_dir=str(arg.strip())
		elif opt=="-b":
			binfile=str(arg.strip())
		elif opt=="-e":
			email_input=str(arg.strip())
		elif opt=="-t":
			thread_input=str(arg.strip())
		elif opt=="-g":
			build_input=str(arg.strip())
		elif opt=="-c":
			cap_input=str(arg.strip())
		elif opt=="-u":
			hubURL_input=str(arg.strip())
		elif opt=="-p":
			otherpara_input=str(arg.strip())
		elif opt=="-n":
			uniq_mk_num=str(arg.strip())

	######### run get IDEAS inputs
	print('get IDEAS sh & parafile file.......')
	get_IDEAS_input(id_input, script_dir, working_dir, output_dir, binfile, email_input, thread_input, build_input, cap_input, hubURL_input, otherpara_input, uniq_mk_num)

if __name__=="__main__":
	main(sys.argv[1:])



