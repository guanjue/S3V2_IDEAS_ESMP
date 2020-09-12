library('data.table')
args = commandArgs(trailingOnly=TRUE)

jobname = args[1]
#jobname = 'run_IDEAS_S3V2_mean'

### state matrix
d = data.frame( fread(paste(jobname, '.state', sep='')) )

### read para file
dpara = data.frame(fread(paste(jobname, '.para', sep='')))
state_sig = dpara[,2]/dpara[,1]
state_sig_ordered = c(0:(length(state_sig)-1))[order(state_sig)]

### get bin bed
bed_info = d[,2:4]
for (i in 5:(dim(d)[2]-1)){
	print(i)
	output_name = colnames(d)[i]
	### get bin bed
	used_col = d[,i]
	bin321 = bed_info[used_col!=state_sig_ordered[1],]
	bin32 = bed_info[(used_col!=state_sig_ordered[1])&(used_col!=state_sig_ordered[2]),]
	bin3 = bed_info[(used_col!=state_sig_ordered[1])&(used_col!=state_sig_ordered[2])&(used_col!=state_sig_ordered[3]),]
	#bin3 = bed_info[(used_col!=state_sig_ordered[1])&(used_col!=state_sig_ordered[2])&(used_col!=state_sig_ordered[3]),]
	### write output
	output_name3 = paste(jobname, '.', output_name, '.3.bed', sep='')
	write.table(bin3, output_name3, quote=F, col.names=F, row.names=F, sep='\t')
	output_name32 = paste(jobname, '.', output_name, '.32.bed', sep='')
	write.table(bin32, output_name32, quote=F, col.names=F, row.names=F, sep='\t')
	output_name321 = paste(jobname, '.', output_name, '.321.bed', sep='')
	write.table(bin321, output_name321, quote=F, col.names=F, row.names=F, sep='\t')
}

### get ct list
ct_list = colnames(d)[-c(1:4, dim(d)[2])]
ct_list = cbind(c('AVERAGE', ct_list[ct_list!='AVERAGE']))
write.table(ct_list, paste(jobname, '.ct.list.txt', sep=''), quote=F, col.names=F, row.names=F, sep='\t')


