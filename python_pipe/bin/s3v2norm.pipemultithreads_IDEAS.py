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

def run_S3V2norm(script_folder, input_tar_file, reference_name, output_tar_file, fdr_thresh, exp_win, rank_lim, upperlim, lowerlim, pmethod, cpk_file, cbg_file, allpk_binary_file):
	step2=call('Rscript '+script_folder+'/s3v2norm_IDEAS.R '+input_tar_file+' '+reference_name+' '+output_tar_file+' '+fdr_thresh+' '+exp_win+' F F '+rank_lim+' '+upperlim+' '+lowerlim+' '+pmethod+' '+cpk_file+' '+cbg_file+' '+allpk_binary_file, shell=True)
	return(0)

def check(a,b,c):
	print(a)
	print(b)
	print(c)
	return(0)

def S3norm_pipeline(threads_num, exp_win, fdr_thresh, rank_lim, upperlim, lowerlim, p_method, common_pk_binary, common_bg_binary, allpk_binary, script_folder, file_list, reference_method):
	### step 2
	print('Get S3norm normalized read counts......')
	### get parameters
	input_tar_file = read2d_array(file_list, str)[:,0]
	input_ref_file = read2d_array(file_list, str)[:,1]
	input_ct = read2d_array(file_list, str)[:,2]
	input_mk = read2d_array(file_list, str)[:,3]
	pmethod = 'z'
	pool1 = Pool(threads_num)
	len_input_tar_file = len(input_tar_file)
	pool1_para_script_dir = np.repeat(script_folder, len_input_tar_file)
	pool1_para_input_tar = input_tar_file
	pool1_para_input_ref = input_ref_file
	#pool1_para_output_tar = np.core.defchararray.add(input_ct, np.repeat('.',len_input_tar_file), input_mk, np.repeat('.bedgraph',len_input_tar_file))
	pool1_para_output_tar = map('.'.join, zip(input_ct, input_mk, np.repeat('S3V2.bedgraph',len_input_tar_file)))
	pool1_para_fdr_thresh = np.repeat(fdr_thresh, len_input_tar_file)
	pool1_para_exp_win = np.repeat(exp_win, len_input_tar_file)
	pool1_para_rank_lim = np.repeat(rank_lim, len_input_tar_file)
	pool1_para_upperlim = np.repeat(upperlim, len_input_tar_file)
	pool1_para_lowerlim = np.repeat(lowerlim, len_input_tar_file)
	pool1_para_pmethod = np.repeat(pmethod, len_input_tar_file)
	pool1_para_cpk = np.repeat(common_pk_binary, len_input_tar_file)
	pool1_para_cbg = np.repeat(common_bg_binary, len_input_tar_file)
	pool1_para_allpk = np.repeat(allpk_binary, len_input_tar_file)
	### run S3V2
	pool1_paras = zip(pool1_para_script_dir, pool1_para_input_tar, pool1_para_input_ref, pool1_para_output_tar, pool1_para_fdr_thresh, pool1_para_exp_win, pool1_para_rank_lim, pool1_para_upperlim, pool1_para_lowerlim, pool1_para_pmethod, pool1_para_cpk, pool1_para_cbg, pool1_para_allpk)
	pool1.starmap(run_S3V2norm, pool1_paras)
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
	opts, args = getopt.getopt(argv,"hn:e:f:l:a:b:p:k:g:i:s:t:r:")
	for opt,arg in opts:
		if opt=="-h":
			print('time python ../src/S3norm_pipeline.py -s script_folder -t input_file_list -r reference_method -m (Method for matching peaks and background: non0mean, non0median, mean, median) -i initial_B -f FDR_thresh -l rank_lim_p -a upperlimit -b lowerlimit -p (p-value_method: neglog10p, z) -k common_pk_binary (0 for nocommon_pk; common_pk_binary.txt) -g common_bg_binary (0 for nocommon_pk; common_bg_binary.txt)')
			return()	
		elif opt=="-n":
			threads_num=int(arg.strip())
		elif opt=="-e":
			exp_win=str(arg.strip())
		elif opt=="-f":
			fdr_thresh=str(arg.strip())
		elif opt=="-l":
			rank_lim=str(arg.strip())
		elif opt=="-a":
			upperlim=str(arg.strip())
		elif opt=="-b":
			lowerlim=str(arg.strip())
		elif opt=="-p":
			p_method=str(arg.strip())
		elif opt=="-k":
			common_pk_binary=str(arg.strip())
		elif opt=="-g":
			common_bg_binary=str(arg.strip())
		elif opt=="-i":
			allpk_binary=str(arg.strip())
		elif opt=="-s":
			script_folder=str(arg.strip())				
		elif opt=="-t":
			file_list=str(arg.strip())
		elif opt=="-r":
			reference_method=str(arg.strip())

	############ Default parameters
	###### required parameters
	try:
		print('User provide script_folder: -s '+str(script_folder))
		print('User provide input_file_list: -t '+str(file_list))
	except NameError:
		print('Missing required parameter(s): time python ../src/S3norm_pipeline.py -s /Users/universe/Documents/2018_BG/S3norm/src/ -t file_list.txt')	
		return()	
	###
	###### optional parameters
	try:
		print('User provide reference_method: -r '+str(reference_method))
		if reference_method!='mean' and reference_method!='median' and reference_method!='max1' and reference_method!='median1':
			print('-r (mean, median, max1, median1)')
			return()
	except NameError:
		print('Default reference_method: -m max1')
		reference_method = 'max1'
	###
	try:
		print('User provide threads_num: -n '+str(threads_num))
		if float(threads_num)+0!=float(threads_num):
			print('-n number of cpu')
			return()
	except NameError:
		print('Default threads_num: -n 4')
		threads_num = 4
	###
	try:
		print('User provide exp_win: -e '+str(exp_win))
	except NameError:
		print('Default exp_win: -e 5')
		exp_win = '5'
	###
	try:
		print('User provide fdr_thresh: -f '+str(fdr_thresh))
		if float(fdr_thresh)<0 or float(fdr_thresh)>1:
			print('-f (0.000001~1.0)')
			return()
	except NameError:
		print('Default fdr_thresh: -f 0.1')
		fdr_thresh = '0.1'
	###
	try:
		print('User provide rank_lim: -l '+str(rank_lim))
		if float(rank_lim)<0 or float(rank_lim)>1:
			print('-l (0.000001~1.0)')
			return()
	except NameError:
		print('Default rank_lim: -l 0.001')
		rank_lim = '0.001'
	###
	try:
		print('User provide upperlim: -a '+str(upperlim))
	except NameError:
		print('Default upperlim: -a 100000')
		upperlim = '100000'
	###
	try:
		print('User provide lowerlim: -b '+str(lowerlim))
	except NameError:
		print('Default lowerlim: -b 0')
		lowerlim = '0'
	###
	try:
		print('User provide p_method: -p '+str(p_method))
		if p_method!='z' and p_method!='neglog10p':
			print('-p (neglog10p, z)')
			return()
	except NameError:
		print('Default p_method: -p z')
		p_method = 'z'
	###
	try:
		print('User provide common_pk_binary: -k '+str(common_pk_binary))
	except NameError:
		print('Default common_pk_binary: -k 0')
		common_pk_binary = '0'
	###
	try:
		print('User provide common_bg_binary: -g '+str(common_bg_binary))
	except NameError:
		print('Default common_bg_binary: -g 0')
		common_bg_binary = '0'
	###

	try:
		print('User provide allpk_binary: -i '+str(allpk_binary))
	except NameError:
		print('Default allpk_binary: -i F')
		allpk_binary = 'F'
	
	######### run s3norm
	print('start S3norm.......')
	S3norm_pipeline(threads_num, exp_win, fdr_thresh, rank_lim, upperlim, lowerlim, p_method, common_pk_binary, common_bg_binary, allpk_binary, script_folder, file_list, reference_method)


if __name__=="__main__":
	main(sys.argv[1:])

