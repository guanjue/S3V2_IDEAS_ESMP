import os
import numpy as np
import subprocess
from subprocess import call
import sys

def creat2log(info, log_file):
	log = open(log_file, 'w')
	log.write(info+'\n')
	log.close()

def add2log(info, log_file):
	log = open(log_file, 'a')
	log.write(info+'\n')
	log.close()

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

def read1d_array(filename,dtype_used):
	import numpy as np
	data=open(filename,'r')
	data0=[]
	for records in data:
		tmp = records.strip()
		data0.append(tmp)
	data0 = np.array(data0,dtype=dtype_used)
	data.close()
	return(data0)

def check(file_check):
	try:
		if not os.path.exists(file_check):
			sys.exit("something wrong with NBP value! *NBP.bedgraph file is not there")
		read1d_array(file_check, float)
		print('ok')
	except ValueError:
		sys.exit("something wrong with NBP value!")

def S3V2_IDEAS_pipeline(get_sigtrack, normalization, get_bw, run_ideas, script_dir, OUTDIR, GENOME, GENOMESIZES, BLACK, metadata, bin_size, local_bg_bin, id_name, email, threads, cap_sig, IDEAS_track_link, other_parafile):
	####################################
	### start!
	log_file=id_name+'.log.txt'
	creat2log('Start!', log_file)

	####################################
	### get signal track
	if get_sigtrack=='T':
		add2log('get signal track......', log_file)
		#a=call('python3 '+script_dir+'/get_signal_track.multithreads.py -s '+script_dir+' -o '+OUTDIR+' -t '+str(threads)+' -g '+GENOMESIZES+' -l '+str(bin_size)+' -b '+BLACK+' -i '+metadata, shell=True)
		#a=call(script_dir+'/prepBwFiles.sh '+OUTDIR+' '+GENOMESIZES+' '+BLACK+' '+metadata+' '+script_dir+' '+str(threads)+' '+log_file+' '+str(bin_size), shell=True)
		a=call(script_dir+'/prepBwFiles.sh '+OUTDIR+' '+GENOMESIZES+' '+BLACK+' '+metadata+' '+script_dir+' '+str(threads)+' '+log_file+' '+str(bin_size), shell=True)
		add2log('get signal track......Done', log_file)
	else:
		add2log('skip get signal track', log_file)

	####################################
	### normalization
	if normalization=='T':
		add2log('normalization......', log_file)
		# 1: get cpk cbg allpk average_sig
		add2log('get cpk cbg allpk average_sig......', log_file)
		a=call('cat '+metadata+' | cut -f2 | sort -u > '+id_name+'.mk_list.txt', shell=True)
		mk_list=id_name+'.mk_list.txt'
		mks = read2d_array(mk_list,str)
		uniq_mk_num = len(mks)
		for mk in mks:
			add2log(mk[0], log_file)
			a=call('bash '+script_dir+'/get_cpk_cbg_allpk_averagesig.sh '+mk[0]+' '+script_dir+' '+metadata, shell=True)
		add2log('get cpk cbg allpk average_sig......Done', log_file)
		# 2: S3norm average across marks
		add2log('S3norm average across marks......', log_file)
		for mk in mks:
			print(mk[0])
			add2log(mk[0], log_file)
			if uniq_mk_num==1:
				add2log('One mk mode...', log_file)
				#a=call('Rscript '+script_dir+'/get_max_median1.R'+' '+'max1'+' '+mk[0]+'.file_list_tmp1'+' '+mk[0]+'.average_sig.max1.bedgraph', shell=True)
				#a=call('python3 '+script_dir+'/s3norm_1mk.py'+' -r '+script_dir+'/prior/H3K27ac.average_sig.sample200k.seed2019.bedgraph'+' -t '+mk[0]+'.average_sig.bedgraph'+' -o '+mk[0]+'.average_sig.max1.bedgraph.S3.bedgraph'+' -c T', shell=True)
				a=call('cp '+mk[0]+'.average_sig.bedgraph'+' '+mk[0]+'.average_sig.bedgraph.S3.bedgraph', shell=True)
			else:
				add2log('multiple mks mode...', log_file)
				a=call('python3 '+script_dir+'/s3norm.py'+' -r '+script_dir+'/prior/H3K27ac.average_sig.sample200k.seed2019.bedgraph'+' -t '+mk[0]+'.average_sig.bedgraph'+' -o '+mk[0]+'.average_sig.bedgraph.S3.bedgraph'+' -c T', shell=True)
		add2log('S3norm average across marks......Done', log_file)
		# 3: S3V2 across samples
		add2log('S3V2 across samples......', log_file)
		for mk in mks:
			print(mk[0])
			add2log(mk[0], log_file) 
			if uniq_mk_num==1:
				add2log('S3V2norm each sample 1 mk mode......', log_file)
				### get input file list
				a=call('cat '+mk[0]+'.file_list_tmp1'+' | awk -F \'.\' -v OFS=\'\\t\' -v average='+mk[0]+'.average_sig.bedgraph.S3.bedgraph \'{print $0, average, $1"_"$3, $2}\' > '+mk[0]+'.file_list.S3V2.txt', shell=True)
				### S3V2 normalization
				#b=call(script_dir+'/s3v2norm.sh'+' -n '+str(threads)+' -e '+str(local_bg_bin)+' -t '+mk[0]+'.file_list.S3V2.txt'+' -k '+mk[0]+'_commonpkfdr01_z.cpk.txt'+' -g '+mk[0]+'_commonpkfdr01_z.cbg.txt'+' -s '+script_dir+' -i '+'F'+' -l '+str(0.0001)+' -x '+log_file, shell=True)
				b=call(script_dir+'/s3v2norm.sh'+' -n '+str(threads)+' -e '+str(local_bg_bin)+' -t '+mk[0]+'.file_list.S3V2.txt'+' -k '+mk[0]+'_commonpkfdr01_z.cpk.txt'+' -g '+mk[0]+'_commonpkfdr01_z.cbg.txt'+' -s '+script_dir+' -i '+'F'+' -l '+str(0.0001)+' -x '+log_file, shell=True)
			else:
				add2log('S3V2norm each sample multiple mks mode......', log_file)
				### get input file list
				a=call('cat '+mk[0]+'.file_list_tmp1'+' | awk -F \'.\' -v OFS=\'\\t\' -v average='+mk[0]+'.average_sig.bedgraph.S3.bedgraph \'{print $0, average, $1"_"$3, $2}\' > '+mk[0]+'.file_list.S3V2.txt', shell=True)
				### S3V2 normalization
				#b=call(script_dir+'/s3v2norm.sh'+' -n '+str(threads)+' -e '+str(local_bg_bin)+' -t '+mk[0]+'.file_list.S3V2.txt'+' -k '+mk[0]+'_commonpkfdr01_z.cpk.txt'+' -g '+mk[0]+'_commonpkfdr01_z.cbg.txt'+' -s '+script_dir+' -i '+mk[0]+'_commonpkfdr01_z.allpk.txt'+' -l '+str(0.0001)+' -x '+log_file, shell=True)
				b=call(script_dir+'/s3v2norm.sh'+' -n '+str(threads)+' -e '+str(local_bg_bin)+' -t '+mk[0]+'.file_list.S3V2.txt'+' -k '+mk[0]+'_commonpkfdr01_z.cpk.txt'+' -g '+mk[0]+'_commonpkfdr01_z.cbg.txt'+' -s '+script_dir+' -i '+mk[0]+'_commonpkfdr01_z.allpk.txt'+' -l '+str(0.0001)+' -x '+log_file, shell=True)
			### get S3V2 ave
			add2log('get S3V2 ave......', log_file)
			a=call('cat '+metadata+' | awk -F \'\\t\' -v OFS=\'\\t\' -v used_mk='+mk[0]+' \'{if ($2==used_mk) print $1"_"$3"."$2".S3V2.bedgraph"}\' > '+mk[0]+'.tmp1.list.txt', shell=True)
			b=call('cat '+metadata+' | awk -F \'\\t\' -v OFS=\'\\t\' -v used_mk='+mk[0]+' \'{if ($2==used_mk) print $1"."$2"."$3".ctrl.idsort.bedgraph.norm.bedgraph"}\' > '+mk[0]+'.tmp2.list.txt', shell=True)
			c=call('paste '+mk[0]+'.tmp1.list.txt'+' '+mk[0]+'.tmp2.list.txt'+' > '+mk[0]+'.getave_nbp.list.txt', shell=True)
			d=call('rm '+mk[0]+'.tmp1.list.txt'+' '+mk[0]+'.tmp2.list.txt', shell=True)
			e=call('Rscript '+script_dir+'/get_average_sig.R'+' '+mk[0]+'.getave_nbp.list.txt'+' '+mk[0]+'.average_sig.bedgraph.S3V2.ave.bedgraph', shell=True)
			add2log('get S3V2 ave......Done', log_file)
		add2log('S3V2 across samples......Done', log_file)
		# 4: normalize controls
		add2log('normalize controls......', log_file)
		a=call('cat '+metadata+' | awk -F \'\\t\' -v OFS=\'\\t\' \'{print $1"."$2"."$3".ctrl.idsort.bedgraph"}\' > all.ctrl.list.txt', shell=True)
		b=call('Rscript '+script_dir+'/non0scale.R all.ctrl.list.txt', shell=True)
		add2log('normalize controls......Done', log_file)
		# release some space
		a=call('rm *.ip.idsort.bedgraph', shell=True)
		b=call('rm *.ctrl.idsort.bedgraph', shell=True)
		# 5: Get NBP
		add2log('Get NBP......', log_file)
		for mk in mks:
			print(mk[0])
			add2log(mk[0], log_file) 
			if uniq_mk_num==1:
				add2log('get NBP 1 mk mode......', log_file)
				### get NBP
				a=call('Rscript '+script_dir+'/global_nbp_NB_cm.1mk.R '+mk[0]+'.getave_nbp.list.txt'+' '+mk[0]+'.average_sig.bedgraph.S3V2.ave.bedgraph'+' '+mk[0]+'_commonpkfdr01_z.cbg.txt', shell=True)
			else:
				add2log('get NBP multiple mks mode......', log_file)
				### get NBP
				a=call('Rscript '+script_dir+'/global_nbp_NB_cm.R '+mk[0]+'.getave_nbp.list.txt'+' '+mk[0]+'.average_sig.bedgraph.S3.bedgraph'+' '+mk[0]+'.average_sig.bedgraph.S3V2.ave.bedgraph'+' '+mk[0]+'_commonpkfdr01_z.cbg.txt', shell=True)
		add2log('Get NBP......Done', log_file)
		# 6: Get IDEAS input
		add2log('Get IDEAS input......', log_file)
		print('create the IDEAS input directory if necessary......')
		if not os.path.exists(id_name+'_IDEAS_input_NB/'):
			os.makedirs(id_name+'_IDEAS_input_NB/')
		# cut 4th column
		for mk in mks:
			all_file = read2d_array(mk[0]+'.getave_nbp.list.txt', str)
			a=call('while read -r IP CTRL; do cut -f4 $IP\'.NBP.bedgraph\' > $IP\'.NBP.txt\'; mv $IP\'.NBP.txt\' '+id_name+'_IDEAS_input_NB/; done < '+mk[0]+'.getave_nbp.list.txt', shell=True)
			file_check = id_name+'_IDEAS_input_NB/'+all_file[0,0]+'.NBP.txt'
			### check if .NBP.txt is there
			check(file_check)
			###
			if uniq_mk_num==1:
				b=call('cut -f4 '+mk[0]+'.average_sig.bedgraph.S3V2.ave.bedgraph.NBP.bedgraph'+' > '+mk[0]+'.average_sig.bedgraph.S3V2.bedgraph.NBP.txt', shell=True)
				b=call('mv '+mk[0]+'.average_sig.bedgraph.S3V2.bedgraph.NBP.txt '+id_name+'_IDEAS_input_NB/', shell=True)
		add2log('Get IDEAS input......Done', log_file)
	else:
		add2log('skip normalization......', log_file)

	####################################
	### get_bw, redo as bash
	if get_bw=='T':
		add2log('get_bw......', log_file)
		#a=call('python3 '+script_dir+'/get_bw.multithreads.py -s '+script_dir+' -o '+OUTDIR+id_name+'_bws'+' -t '+str(threads)+' -g '+GENOMESIZES+' -i '+metadata, shell=True)
		#a=call(script_dir+'/get_bws.sh '+script_dir+' '+OUTDIR+id_name+'_bws'+' '+str(threads)+' '+GENOMESIZES+' '+log_file, shell=True)
		a=call(script_dir+'/get_bws.sh '+script_dir+' '+OUTDIR+id_name+'_bws'+' '+str(threads)+' '+GENOMESIZES+' '+log_file, shell=True)
		add2log('get_bw......Done', log_file)
		### rm bedgraphs
		a=call('rm *.bedgraph.norm.bedgraph', shell=True)
		b=call('rm *.S3V2.bedgraph.NBP.bedgraph', shell=True)
		c=call('rm *.S3V2.ave.bedgraph.NBP.bedgraph', shell=True)
		d=call('rm *.bedgraph.S3.bedgraph', shell=True)
		e=call('rm *.S3V2.bedgraph', shell=True)
		f=call('rm *.S3V2.ave.bedgraph', shell=True)
		g=call('rm *.average_sig.bedgraph', shell=True)

	####################################
	### run_ideas
	if run_ideas=='T':
		add2log('run_ideas......', log_file)
		# get bin .txt
		a=call('cat '+metadata+' | cut -f2 | sort -u > '+id_name+'.mk_list.txt', shell=True)
		mk_list=id_name+'.mk_list.txt'
		mks = read2d_array(mk_list,str)
		uniq_mk_num = len(mks)
		print(uniq_mk_num)
		a=call('cat '+OUTDIR+'/windowsNoBlack.withid.bed'+' | awk -F \'\\t\' -v OFS=\' \' \'{print $1,$2,$3,$4}\' > '+OUTDIR+'/windowsNoBlack.withid.IDEASbins.txt', shell=True)
		# get sh & parafile
		b=call('python3 '+script_dir+'/get_IDEAS_input.py'+' -i '+id_name+' -s '+script_dir+'/IDEAS_2018/'+' -w '+OUTDIR+' -o '+OUTDIR+'/'+id_name+'_IDEAS_output/'+' -b '+OUTDIR+'/windowsNoBlack.withid.IDEASbins.txt'+' -e '+email+' -t '+str(threads)+' -g '+GENOME+' -c '+str(cap_sig)+' -u '+IDEAS_track_link+' -p '+other_parafile+' -n '+str(uniq_mk_num), shell=True)
		# get input file
		c=call('cat '+metadata+' | awk -F \'\\t\' -v OFS=\' \' -v S3norm_NBP_dir='+OUTDIR+'/'+id_name+'_IDEAS_input_NB/'+' \'{print $1"_"$3,$2,S3norm_NBP_dir$1"_"$3"."$2".S3V2.bedgraph.NBP.txt"}\' | sort -u > '+OUTDIR+'/'+id_name+'.input', shell=True)
		# add average track for input
		if uniq_mk_num==1:
			d=call('head -1 '+metadata+' | awk -F \'\\t\' -v OFS=\' \' -v S3norm_NBP_dir='+OUTDIR+'/'+id_name+'_IDEAS_input_NB/'+' \'{print "AVERAGE",$2,S3norm_NBP_dir$2".average_sig.bedgraph.S3V2.bedgraph.NBP.txt"}\' | sort -u > '+OUTDIR+'/'+id_name+'.input.tmp', shell=True)
			d=call('cat '+OUTDIR+'/'+id_name+'.input >> '+OUTDIR+'/'+id_name+'.input.tmp', shell=True)
			d=call('head -1 '+metadata+' | awk -F \'\\t\' -v OFS=\' \' -v S3norm_NBP_dir='+OUTDIR+'/'+id_name+'_IDEAS_input_NB/'+' \'{print "AVERAGE",$2"_1",S3norm_NBP_dir$2".average_sig.bedgraph.S3V2.bedgraph.NBP.txt"}\' | sort -u >> '+OUTDIR+'/'+id_name+'.input.tmp', shell=True)
			d=call('cat '+metadata+' | awk -F \'\\t\' -v OFS=\' \' -v S3norm_NBP_dir='+OUTDIR+'/'+id_name+'_IDEAS_input_NB/'+' \'{print $1"_"$3,$2"_1",S3norm_NBP_dir$1"_"$3"."$2".S3V2.bedgraph.NBP.txt"}\' | sort -u >> '+OUTDIR+'/'+id_name+'.input.tmp', shell=True)
			d=call('mv '+OUTDIR+'/'+id_name+'.input.tmp'+' '+OUTDIR+'/'+id_name+'.input', shell=True)
		### run IDEAS
		e=call('time bash '+id_name+'.sh', shell=True)
		### get cCRE list
		if uniq_mk_num==1:
			#redo as bash
			f=call('python3 '+script_dir+'/get_cCRE.pipeline.py -s '+script_dir+' -i '+id_name+' -w '+OUTDIR+'/'+id_name+'_IDEAS_output/'+' -l '+str(bin_size), shell=True)
	else:
		add2log('skip run_ideas......', log_file)

	####################################



############################################################################

import getopt
import sys
def main(argv):
	### read user provided parameters
	opts, args = getopt.getopt(argv,"hu:v:y:z:s:o:g:c:b:i:l:n:d:e:t:a:w:x:")
	for opt,arg in opts:
		if opt=="-h":
			print('time python3 $script_dir/S3V2_IDEAS_pipeline.py -u $get_sigtrack -v $normalization -y $get_bw -z $run_ideas -s $script_dir -o $OUTDIR -g $GENOME -c $GENOMESIZES -b $BLACK -i $metadata -d $id_name -e $email -t $threads -w $IDEAS_track_link -x $other_parafile -l $bin_size -n $local_bg_bin -a $cap_sig')
			return()	
		elif opt=="-u":
			get_sigtrack=str(arg.strip())
		elif opt=="-v":
			normalization=str(arg.strip())
		elif opt=="-y":
			get_bw=str(arg.strip())
		elif opt=="-z":
			run_ideas=str(arg.strip())
		elif opt=="-s":
			script_dir=str(arg.strip())
		elif opt=="-o":
			OUTDIR=str(arg.strip())
		elif opt=="-g":
			GENOME=str(arg.strip())
		elif opt=="-c":
			GENOMESIZES=str(arg.strip())
		elif opt=="-b":
			BLACK=str(arg.strip())
		elif opt=="-i":
			metadata=str(arg.strip())
		elif opt=="-l":
			bin_size=int(arg.strip())
		elif opt=="-n":
			local_bg_bin=int(arg.strip())				
		elif opt=="-d":
			id_name=str(arg.strip())
		elif opt=="-e":
			email=str(arg.strip())
		elif opt=="-t":
			threads=int(arg.strip())
		elif opt=="-a":
			cap_sig=float(arg.strip())
		elif opt=="-w":
			IDEAS_track_link=str(arg.strip())
		elif opt=="-x":
			other_parafile=str(arg.strip())

	############ Default parameters
	###### required parameters
	try:
		print('User provide script_dir: -s '+str(script_dir))
		print('User provide OUTDIR: -o '+str(OUTDIR))
		print('User provide GENOME: -g '+str(GENOME))
		print('User provide GENOMESIZES: -c '+str(GENOMESIZES))
		print('User provide BLACK: -b '+str(BLACK))
		print('User provide metadata: -i '+str(metadata))
		print('User provide id_name: -d '+str(id_name))
	except NameError:
		print('Missing required parameter(s): time python3 ../src/S3norm_pipeline.py -s /Users/universe/Documents/2018_BG/S3norm/src/ -t file_list.txt')	
		return()

	###
	###### optional parameters
	try:
		print('User provide get_sigtrack: -u '+str(get_sigtrack))
		if get_sigtrack=='':
			print('-u T or F')
			return()
	except NameError:
		print('Default get_sigtrack: -u T')
		get_sigtrack = 'T'
	try:
		print('User provide normalization: -v '+str(normalization))
		if normalization=='':
			print('-v T or F')
			return()
	except NameError:
		print('Default normalization: -v T')
		normalization = 'T'
	try:
		print('User provide get_bw: -y '+str(get_bw))
		if get_bw=='':
			print('-y T or F')
			return()
	except NameError:
		print('Default get_bw: -y T')
		get_bw = 'T'
	try:
		print('User provide run_ideas: -z '+str(run_ideas))
		if run_ideas=='':
			print('-z T or F')
			return()
	except NameError:
		print('Default run_ideas: -z T')
		run_ideas = 'T'


	###
	try:
		print('User provide IDEAS_track_link: -w '+str(IDEAS_track_link))
		if IDEAS_track_link=='':
			print('no IDEAS_track_link provided')
			return()
	except NameError:
		print('Default IDEAS_track_link: -w http://accessible_link_for_genome_browser_track/')
		IDEAS_track_link='http://accessible_link_for_genome_browser_track/'

	###
	try:
		print('User provide email: -e '+str(email))
		if email=='':
			print('no email provided')
			return()
	except NameError:
		print('Default email: -e *@*.*')
		email='*@*.*'

	###
	try:
		print('User provide other_parafile: -x '+str(other_parafile))
		if other_parafile=='':
			print('no other_parafile provided')
			return()
	except NameError:
		print('Default other_parafile: -x F')
		other_parafile='F'

	###
	try:
		print('User provide bin_size: -l '+str(bin_size))
		if float(bin_size)<=0:
			print('-l integer_greater_than_0')
			return()
	except NameError:
		print('Default bin_size: -l 200')
		bin_size = 200
	###
	try:
		print('User provide local_bg_bin: -n '+str(local_bg_bin))
		if float(local_bg_bin)<=0:
			print('-n integer_greater_than_0')
			return()
	except NameError:
		print('Default local_bg_bin: -n 5')
		local_bg_bin = 5
	###
	try:
		print('User provide threads: -t '+str(threads))
		if float(threads)<=0:
			print('-n integer_greater_than_0')
			return()
	except NameError:
		print('Default threads: -t 2')
		threads = 2
	###
	try:
		print('User provide cap_sig: -a '+str(cap_sig))
		if float(cap_sig)<=0:
			print('-n integer_greater_than_0')
			return()
	except NameError:
		print('Default cap_sig: -t 16')
		cap_sig = 16
	###

	######### run s3norm
	print('start S3V2_IDEAS_pipeline.......')
	S3V2_IDEAS_pipeline(get_sigtrack, normalization, get_bw, run_ideas, script_dir, OUTDIR, GENOME, GENOMESIZES, BLACK, metadata, bin_size, local_bg_bin, id_name, email, threads, cap_sig, IDEAS_track_link, other_parafile)



if __name__=="__main__":
	main(sys.argv[1:])

