import os
import numpy as np
import subprocess
from subprocess import call
from multiprocessing import Pool

################################################################################################
### read 2d array
def read2d_array(filename,dtype_used):
	import numpy as np
	data=open(filename,'r')
	data0=[]
	for records in data:
		tmp = [x.strip() for x in records.split('\t')]
		data0.append(tmp)
	data0 = np.array(data0,dtype=dtype_used)
	data.close()
	return(data0)

def run_S3V2norm(script_folder, input_tar_file, reference_name, output_tar_file, fdr_thresh, exp_win, rank_lim, upperlim, lowerlim, pmethod, cpk_file, cbg_file):
	step1=call('python '+script_folder+'/s3norm.py '+input_tar_file+' '+reference_name+' '+output_tar_file+' '+fdr_thresh+' '+exp_win+' F F '+rank_lim+' '+upperlim+' '+lowerlim+' '+pmethod+' '+cpk_file+' '+cbg_file, shell=True)
	return(0)

def run_S3norm(script_folder, input_tar_file, reference_name, output_name, NTmethod, B_init, fdr_thresh, rank_lim, upperlim, lowerlim, p_method, common_pk_binary, common_bg_binary):
	step1=call('python '+script_folder+'/s3norm.py -r '+reference_name+' -t '+input_tar_file+' -o '+output_name+' -m '+NTmethod+' -i '+B_init+' -f '+fdr_thresh+' -l '+rank_lim+' -a '+upperlim+' -b '+lowerlim+' -p '+p_method+' -k '+common_pk_binary+' -g '+common_bg_binary, shell=True)

def run_QTnorm(script_folder, input_tar_file, reference_name, output_name, upperlim, lowerlim):
	step1=call('Rscript '+script_folder+'/QTnorm_rc.R '+reference_name+' '+input_tar_file+' '+output_name+' '+upperlim+' '+lowerlim, shell=True)

def run_TSnorm(script_folder, input_tar_file, reference_name, output_name, upperlim, lowerlim):
	step1=call('Rscript '+script_folder+'/TSnorm_rc.R '+reference_name+' '+input_tar_file+' '+output_name+' '+upperlim+' '+lowerlim, shell=True)

def run_MAnorm(script_folder, input_tar_file, reference_name, output_name, upperlim, lowerlim, common_pk_binary):
	step1=call('Rscript '+script_folder+'/MAnorm_rc.R '+reference_name+' '+input_tar_file+' '+output_name+' '+upperlim+' '+lowerlim+' '+common_pk_binary, shell=True)

def get_bw(script_folder, input_tar_file_tmp, genome_size, output_name):
	step1=call('sort -k1,1 -k2,2n '+input_tar_file_tmp+' > '+input_tar_file_tmp+'.sort.bedgraph', shell=True)
	#step2=call(script_folder+'/bedGraphToBigWig '+input_tar_file_tmp+'.sort.bedgraph'+' '+genome_size+' '+output_name, shell=True)
	#step3=call('rm '+input_tar_file_tmp+'.sort.bedgraph', shell=True)

def get_nbp(script_folder, sig_file_tmp, input_tar_file_tmp, output_file_tmp):
	step1=call('Rscript '+script_folder+'/negative_binomial_neglog10p.R '+sig_file_tmp+' '+input_tar_file_tmp+' '+output_file_tmp, shell=True)


def get_merged_macs2_pk(input_fdr_bedgraph, ouput_pk, thresh):
	step1=call('macs2 bdgpeakcall -i '+input_fdr_bedgraph+' -o '+ouput_pk+'.tmp.txt'+' -c '+thresh, shell=True)
	step2=call('tail -n+2 '+ouput_pk+'.tmp.txt'+' | cut -f1,2,3 > '+ouput_pk, shell=True)
	step3=call('rm '+ouput_pk+'.tmp.txt', shell=True)

def get_macs2_pk(ct, mk, thresh):
	###
	input_1 = ct+'.'+mk+'.meanrc.txt.bedgraph.nbp.bedgraph.fdr.bedgraph'
	output_1 = ct+'.'+mk+'.RAW.bedgraph.nbp.fdr.macs2.bed'
	step1=get_merged_macs2_pk(input_1, output_1, thresh)
	###
	input_1 = ct+'.'+mk+'.TS.bedgraph.nbp.bedgraph.fdr.bedgraph'
	output_2 = ct+'.'+mk+'.TS.bedgraph.nbp.fdr.macs2.bed'
	step2=get_merged_macs2_pk(input_1, output_2, thresh)
	###
	input_1 = ct+'.'+mk+'.MA.bedgraph.nbp.bedgraph.fdr.bedgraph'
	output_3 = ct+'.'+mk+'.MA.bedgraph.nbp.fdr.macs2.bed'
	step3=get_merged_macs2_pk(input_1, output_3, thresh)
	###
	input_1 = ct+'.'+mk+'.QT.bedgraph.nbp.bedgraph.fdr.bedgraph'
	output_4 = ct+'.'+mk+'.QT.bedgraph.nbp.fdr.macs2.bed'
	step4=get_merged_macs2_pk(input_1, output_4, thresh)
	###
	input_1 = ct+'.'+mk+'.S3.bedgraph.s3norm.bedgraph.nbp.bedgraph.fdr.bedgraph'
	output_5 = ct+'.'+mk+'.S3.bedgraph.nbp.fdr.macs2.bed'
	step5=get_merged_macs2_pk(input_1, output_5, thresh)
	###
	input_1 = ct+'.'+mk+'.S3V2.bedgraph.nbp.bedgraph.fdr.bedgraph'
	output_6 = ct+'.'+mk+'.S3V2.bedgraph.nbp.fdr.macs2.bed'
	step6=get_merged_macs2_pk(input_1, output_6, thresh)
	###
	step1=call('cat '+output_1+' '+output_2+' '+output_3+' '+output_4+' '+output_5+' '+output_6+' | sort -k1,1 -k2,2n > '+ct+'.'+mk+'.pooled.macs2pk.bed', shell=True)
	step2=call('bedtools merge -i '+ct+'.'+mk+'.pooled.macs2pk.bed'+' > '+ct+'.'+mk+'.merged.macs2pk.bed', shell=True)


def S3norm_pipeline(threads_num, script_folder, file_list, thresh):
	### step 2
	print('Get S3norm normalized read counts......')
	### get parameters
	sig_tar_file = read2d_array(file_list, str)[:,0]
	input_tar_file = read2d_array(file_list, str)[:,1]
	input_ct = read2d_array(file_list, str)[:,2]
	input_mk = read2d_array(file_list, str)[:,3]	
	pool1 = Pool(threads_num)
	len_input_tar_file = len(input_tar_file)
	pool1_para_thresh = np.repeat(thresh, len_input_tar_file)
	### run QTnorm
	pool1_paras = zip(input_ct, input_mk, pool1_para_thresh)
	pool1.starmap(get_macs2_pk, pool1_paras)
	pool1.close()
	pool1.join()



############################################################################
### time python ../src/s3norm.py -r reference.bedgraph -t target.bedgraph -o target -m non0mean -i 2.0 -f 0.05 -l 0.001 -a 100 -b 0 -p z -k 0 -g 0
### time python ../src/s3norm.py -r reference.bedgraph -t target.bedgraph -o target

### time python ../src/s3norm.py -r reference.bedgraph -t target.bedgraph -o target -m non0mean -i 2.0 -f 0.05 -l 0.001 -a 100 -b 0 -p z -k 0 -g 0
### time python ../src/S3norm_pipeline.py -s /Users/universe/Documents/2018_BG/S3norm/src/ -t file_list.txt

import getopt
import sys
def main(argv):
	### read user provided parameters
	opts, args = getopt.getopt(argv,"hn:s:i:t:")
	for opt,arg in opts:
		if opt=="-h":
			print('time python ../src/S3norm_pipeline.py -s script_folder -t input_file_list -r reference_method -m (Method for matching peaks and background: non0mean, non0median, mean, median) -i initial_B -f FDR_thresh -l rank_lim_p -a upperlimit -b lowerlimit -p (p-value_method: neglog10p, z) -k common_pk_binary (0 for nocommon_pk; common_pk_binary.txt) -g common_bg_binary (0 for nocommon_pk; common_bg_binary.txt)')
			return()	
		elif opt=="-n":
			threads_num=int(arg.strip())
		elif opt=="-s":
			script_folder=str(arg.strip())				
		elif opt=="-i":
			file_list=str(arg.strip())
		elif opt=="-t":
			thresh=str(arg.strip())

	############ Default parameters
	###### required parameters
	try:
		print('User provide script_folder: -s '+str(script_folder))
		print('User provide input_file_list: -i '+str(file_list))
		print('User provide input_file_list: -t '+str(thresh))
	except NameError:
		print('Missing required parameter(s): time python ../src/S3norm_pipeline.py -s /Users/universe/Documents/2018_BG/S3norm/src/ -t file_list.txt -g mm10.chrom.size')	
		return()	
	###
	###### optional parameters
	try:
		print('User provide threads_num: -n '+str(threads_num))
		if float(threads_num)+0!=float(threads_num):
			print('-n number of cpu')
			return()
	except NameError:
		print('Default threads_num: -n 4')
		threads_num = 4

	######### run s3norm
	print('start S3norm.......')
	S3norm_pipeline(threads_num, script_folder, file_list, thresh)

if __name__=="__main__":
	main(sys.argv[1:])

