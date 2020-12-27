import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import norm

################################################################################################
### read 2d array
def read2d_array(filename,dtype_used):
	import numpy as np
	data=open(filename,'r')
	data0=[]
	for records in data:
		tmp = [x.strip() for x in records.split(' ')]
		data0.append(tmp)
	data0 = np.array(data0,dtype=dtype_used)
	data.close()
	return data0

################################################################################################
### write 2d matrix
def write2d_array(array,output):
	r1=open(output,'w')
	for records in array:
		for i in range(0,len(records)-1):
			r1.write(str(records[i])+'\t')
		r1.write(str(records[len(records)-1])+'\n')
	r1.close()

################################################################################################
### bed2inv
def bed2inv(input_bed, output_inv):
	data0 = read2d_array(input_bed, str)
	data_chr_dict = {}
	data_chr_vec = []
	data_chr_inv = []

	### loop bed
	for bed_info in data0:
		chrom = bed_info[0]
		if not (chrom in data_chr_dict):
			data_chr_dict[chrom] = 1
			data_chr_vec.append(chrom)
		else:
			data_chr_dict[chrom] = data_chr_dict[chrom] + 1

	### get inv file
	start_num = 0
	for chrom in data_chr_vec:
		chrom_bin_num = data_chr_dict[chrom]
		end_num = start_num+chrom_bin_num
		inv_vec = [chrom, start_num, end_num]
		data_chr_inv.append(inv_vec)
		start_num = end_num

	data_chr_inv = np.array(data_chr_inv)
	write2d_array(data_chr_inv,output_inv)



############################################################################
#time python3 bed2inv.py -i input.bed -o output.inv

import getopt
import sys
def main(argv):
	try:
		opts, args = getopt.getopt(argv,"hi:o:")
	except getopt.GetoptError:
		print 'time python3 bed2inv.py -i input.bed -o output.inv'
		sys.exit(2)

	for opt,arg in opts:
		if opt=="-h":
			print 'time python3 bed2inv.py -i input.bed -o output.inv'
			sys.exit()
		elif opt=="-i":
			input_bed=str(arg.strip())
		elif opt=="-o":
			output_inv=str(arg.strip())

	bed2inv(input_bed, output_inv)

if __name__=="__main__":
	main(sys.argv[1:])


