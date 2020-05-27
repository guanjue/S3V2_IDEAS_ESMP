import os
import os.path
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

### get signal track
def run_getsignal_track(script_folder, ip_file, ctrl_file, bed_file_noid, bed_file_withid, ct, mk, id):
	### get ip signal
	a=call(script_folder+'/bigWigAverageOverBed '+ip_file+' '+bed_file_withid+' '+ct+'.'+mk+'.'+id+'.ip.tab', shell=True)
	### get ctrl signal
	if ctrl_file!='':
		b=call(script_folder+'/bigWigAverageOverBed '+ctrl_file+' '+bed_file_withid+' '+ct+'.'+mk+'.'+id+'.ctrl.tab', shell=True)
		b=call('cut -f5 '+ct+'.'+mk+'.'+id+'.ctrl.tab'+' > '+ct+'.'+mk+'.'+id+'.ctrl.tab.txt', shell=True)
		b=call('rm '+ct+'.'+mk+'.'+id+'.ctrl.tab', shell=True)
	else:
		print(ct+'.'+mk+'.'+id+': no control')
		b=call('cat '+bed_file_withid+' | awk \'{print 1}\' '+' > '+ct+'.'+mk+'.'+id+'.ctrl.tab.txt', shell=True)
	### get id sort signal txt
	c=call('paste '+ct+'.'+mk+'.'+id+'.ip.tab'+' '+ct+'.'+mk+'.'+id+'.ctrl.tab.txt'+' | sort -k1,1n | cut -f5,7 '+' > '+ct+'.'+mk+'.'+id+'.ip_ctrl.idsort.txt', shell=True)
	c=call('paste '+bed_file_noid+' '+ct+'.'+mk+'.'+id+'.ip_ctrl.idsort.txt'+' | cut -f1,2,3,4 > '+ct+'.'+mk+'.'+id+'.ip.idsort.bedgraph', shell=True)
	c=call('paste '+bed_file_noid+' '+ct+'.'+mk+'.'+id+'.ip_ctrl.idsort.txt'+' | cut -f1,2,3,5 > '+ct+'.'+mk+'.'+id+'.ctrl.idsort.bedgraph', shell=True)
	### rm tmp files
	d=call('rm '+ct+'.'+mk+'.'+id+'.ip_ctrl.idsort.txt'+' '+ct+'.'+mk+'.'+id+'.ip.tab'+' '+ct+'.'+mk+'.'+id+'.ctrl.tab.txt', shell=True)
	return(0)

### get signal track pipeline
def get_signal_track(script_folder, OUTDIR, threads_num, GENOMESIZES, bin_size, BLACK, metadata):
	### 1
	print('create the output directory if necessary......')
	if not os.path.exists(OUTDIR):
		os.makedirs(OUTDIR)
	os.chdir(OUTDIR)
	### 2
	print('create bed file......')
	a=call('bedtools makewindows -g '+GENOMESIZES+' -w '+str(bin_size)+' > '+OUTDIR+'/windows.bed', shell=True)
	b=call('bedtools subtract -a '+OUTDIR+'/windows.bed'+' -b '+BLACK+' > '+OUTDIR+'/windowsNoBlack.bed', shell=True)
	c=call('Rscript '+script_folder+'/add_id.R'+' '+OUTDIR+'/windowsNoBlack.bed'+' '+OUTDIR+'/windowsNoBlack.withid.bed', shell=True)
	d=call('cut -f1,2,3 '+OUTDIR+'/windowsNoBlack.withid.bed'+' > '+OUTDIR+'/windowsNoBlack.noid.bed', shell=True)
	print('create bed file......Done')
	### 3
	print('get signal track......')
	print('get meta_info')
	meta_info = read_metadata(metadata, str)
	para_len = meta_info.shape[0]
	ct_list = meta_info[:,0]
	mk_list = meta_info[:,1]
	id_list = meta_info[:,2]
	ip_list = meta_info[:,3]
	ctrl_list = meta_info[:,4]
	script_folder_list = np.repeat(script_folder, para_len)
	bed_file_noid_list = np.repeat(OUTDIR+'/windowsNoBlack.noid.bed', para_len)
	bed_file_withid_list = np.repeat(OUTDIR+'/windowsNoBlack.withid.bed', para_len)
	print('Pool')
	pool0 = Pool(threads_num)
	pool0_paras = zip(script_folder_list, ip_list, ctrl_list, bed_file_noid_list, bed_file_withid_list, ct_list, mk_list, id_list)
	pool0.starmap(run_getsignal_track, pool0_paras)
	pool0.close()
	pool0.join()
	print('get signal track......Done')




############################################################################

#script_dir='/storage/home/gzx103/scratch/S3V2norm_compare/scripts/'
#OUTDIR='/storage/home/gzx103/scratch/S3V2norm_compare/hg38_IDEAS_state_SC/'
#threads_num=4
#GENOMESIZES='/storage/home/gzx103/scratch/S3V2norm_compare/hg38_IDEAS_state_SC/GRCh38_basic_chroms.sizes.tsv'
#bin_size=200
#BLACK='/gpfs/group/rch8/legacy/group/genomes/hg38/blacklist.lifted.hg38.bed'
#metadata='/storage/home/gzx103/scratch/S3V2norm_compare/hg38_IDEAS_state_SC/metadata.for.test.txt'

#time python get_signal_track.py -s $script_dir -o $OUTDIR -t $threads_num -g $GENOMESIZES -l $bin_size -b $BLACK -i $metadata

import getopt
import sys
def main(argv):
	### read user provided parameters
	opts, args = getopt.getopt(argv,"hs:o:t:g:l:b:i:")
	for opt,arg in opts:
		if opt=="-h":
			print('time python get_signal_track.py -s $script_dir -o $OUTDIR -t $threads_num -g $GENOMESIZES -l $bin_size -b $BLACK -i $metadata')
			return()	
		elif opt=="-s":
			script_folder=str(arg.strip())
		elif opt=="-o":
			OUTDIR=str(arg.strip())
		elif opt=="-t":
			threads_num=int(arg.strip())
		elif opt=="-g":
			GENOMESIZES=str(arg.strip())
		elif opt=="-l":
			bin_size=str(arg.strip())
		elif opt=="-b":
			BLACK=str(arg.strip())
		elif opt=="-i":
			metadata=str(arg.strip())

	############ Default parameters
	###### required parameters
	try:
		print('User provide script_folder: -s '+str(script_folder))
		print('User provide script_folder: -o '+str(OUTDIR))
		print('User provide input_file_list: -t '+str(threads_num))
		print('User provide input_file_list: -g '+str(GENOMESIZES))
		print('User provide input_file_list: -l '+str(bin_size))
		print('User provide input_file_list: -b '+str(BLACK))
		print('User provide input_file_list: -i '+str(metadata))
	except NameError:
		print('Missing required parameter(s):')	
		return()

	######### run get_signal_track
	print('start get_signal_track.......')
	get_signal_track(script_folder, OUTDIR, threads_num, GENOMESIZES, bin_size, BLACK, metadata)




if __name__=="__main__":
	main(sys.argv[1:])

