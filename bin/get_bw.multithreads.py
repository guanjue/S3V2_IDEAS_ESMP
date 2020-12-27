import os
import os.path
import numpy as np
import subprocess
from subprocess import call
from multiprocessing import Pool

### read_metadata
def read_metadata(filename,dtype_used):
	import numpy as np
	data=open(filename,'r')
	data0=[]
	for records in data:
		tmp = [x.strip() for x in records.split('\t')]
		if len(tmp)==4:
			tmp.append('')
		data0.append(tmp)
	data0 = np.array(data0,dtype=dtype_used)
	data.close()
	return(data0)

### get bw
def run_get_bw(script_folder, GENOMESIZES, bedgraph_file, OUTDIR):	
	a=call('sort -k1,1 -k2,2n '+bedgraph_file+' > '+bedgraph_file+'.pksort.bedgraph', shell=True)
	b=call(script_folder+'/bedGraphToBigWig '+bedgraph_file+'.pksort.bedgraph'+' '+GENOMESIZES+' '+bedgraph_file+'.bw', shell=True)
	c=call('rm '+bedgraph_file+'.pksort.bedgraph', shell=True)
	c=call('mv '+bedgraph_file+'.bw'+' '+OUTDIR, shell=True)
	return(0)

### get bw pipeline
def get_bw(script_folder, OUTDIR, threads_num, GENOMESIZES, metadata):
	### 1
	print('create the output directory if necessary......')
	if not os.path.exists(OUTDIR+'_NBP'):
		os.makedirs(OUTDIR+'_NBP')
	if not os.path.exists(OUTDIR+'_RC'):
		os.makedirs(OUTDIR+'_RC')
	### 2
	print('get bw......')
	print('get meta_info')
	meta_info = read_metadata(metadata, str)
	para_len = meta_info.shape[0]
	ct_list = meta_info[:,0]
	mk_list = meta_info[:,1]
	id_list = meta_info[:,2]
	uniq_mk_list = np.unique(mk_list)
	para_len_uniq_mk = len(uniq_mk_list)
	print('for S3V2 RC')
	script_folder_list = np.repeat(script_folder, para_len)
	GENOMESIZES_list = np.repeat(GENOMESIZES, para_len)
	bedgraph_file_list = map(''.join, zip(ct_list, np.repeat('_',para_len), id_list, np.repeat('.',para_len), mk_list, np.repeat('.S3V2.bedgraph',para_len)))
	OUTDIR_list = np.repeat(OUTDIR+'_RC', para_len)
	pool0 = Pool(threads_num)
	pool0_paras = zip(script_folder_list, GENOMESIZES_list, bedgraph_file_list, OUTDIR_list)
	pool0.starmap(run_get_bw, pool0_paras)
	pool0.close()
	pool0.join()
	print('for average signal')
	script_folder_list = np.repeat(script_folder, para_len_uniq_mk)
	GENOMESIZES_list = np.repeat(GENOMESIZES, para_len_uniq_mk)
	bedgraph_file_list = map(''.join, zip(uniq_mk_list, np.repeat('.average_sig.bedgraph.S3V2.ave.bedgraph',para_len_uniq_mk)))
	OUTDIR_list = np.repeat(OUTDIR+'_RC', para_len_uniq_mk)
	pool0 = Pool(threads_num)
	pool0_paras = zip(script_folder_list, GENOMESIZES_list, bedgraph_file_list, OUTDIR_list)
	pool0.starmap(run_get_bw, pool0_paras)
	pool0.close()
	pool0.join()
	print('for S3V2 -log10 p-value')
	script_folder_list = np.repeat(script_folder, para_len)
	GENOMESIZES_list = np.repeat(GENOMESIZES, para_len)
	bedgraph_file_list = map(''.join, zip(ct_list, np.repeat('_',para_len), id_list, np.repeat('.',para_len), mk_list, np.repeat('.S3V2.bedgraph.NBP.bedgraph',para_len)))
	OUTDIR_list = np.repeat(OUTDIR+'_NBP', para_len)
	pool0 = Pool(threads_num)
	pool0_paras = zip(script_folder_list, GENOMESIZES_list, bedgraph_file_list, OUTDIR_list)
	pool0.starmap(run_get_bw, pool0_paras)
	pool0.close()
	pool0.join()
	print('for average signal')
	script_folder_list = np.repeat(script_folder, para_len_uniq_mk)
	GENOMESIZES_list = np.repeat(GENOMESIZES, para_len_uniq_mk)
	bedgraph_file_list = map(''.join, zip(uniq_mk_list, np.repeat('.average_sig.bedgraph.S3V2.ave.bedgraph.NBP.bedgraph',para_len_uniq_mk)))
	OUTDIR_list = np.repeat(OUTDIR+'_NBP', para_len_uniq_mk)
	pool0 = Pool(threads_num)
	pool0_paras = zip(script_folder_list, GENOMESIZES_list, bedgraph_file_list, OUTDIR_list)
	pool0.starmap(run_get_bw, pool0_paras)
	pool0.close()
	pool0.join()
	print('get bw......Done')




############################################################################

#script_dir='/storage/home/gzx103/scratch/S3V2norm_compare/scripts/'
#OUTDIR='/storage/home/gzx103/scratch/S3V2norm_compare/hg38_IDEAS_state_SC/bws'
#threads_num=4
#GENOMESIZES='/storage/home/gzx103/scratch/S3V2norm_compare/hg38_IDEAS_state_SC/GRCh38_basic_chroms.sizes.tsv'
#metadata='/storage/home/gzx103/scratch/S3V2norm_compare/hg38_IDEAS_state_SC/metadata.for.test.txt'

#time python3 get_bw.multithreads.py -s $script_dir -o $OUTDIR -t $threads_num -g $GENOMESIZES -i $metadata

import getopt
import sys
def main(argv):
	### read user provided parameters
	opts, args = getopt.getopt(argv,"hs:o:t:g:i:")
	for opt,arg in opts:
		if opt=="-h":
			print('time python3 get_bw.multithreads.py -s $script_dir -o $OUTDIR -t $threads_num -g $GENOMESIZES -i $metadata')
			return()	
		elif opt=="-s":
			script_folder=str(arg.strip())
		elif opt=="-o":
			OUTDIR=str(arg.strip())
		elif opt=="-t":
			threads_num=int(arg.strip())
		elif opt=="-g":
			GENOMESIZES=str(arg.strip())
		elif opt=="-i":
			metadata=str(arg.strip())

	############ Default parameters
	###### required parameters
	try:
		print('User provide script_folder: -s '+str(script_folder))
		print('User provide script_folder: -o '+str(OUTDIR))
		print('User provide input_file_list: -t '+str(threads_num))
		print('User provide input_file_list: -g '+str(GENOMESIZES))
		print('User provide input_file_list: -i '+str(metadata))
	except NameError:
		print('Missing required parameter(s):')	
		return()

	######### run get_signal_track
	print('start get_signal_track.......')
	get_bw(script_folder, OUTDIR, threads_num, GENOMESIZES, metadata)




if __name__=="__main__":
	main(sys.argv[1:])

