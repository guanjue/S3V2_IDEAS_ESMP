args = commandArgs(trailingOnly=TRUE)

non0scale = function(sig, refmean, refsd){
	sig_non0 = sig[sig>0]
	sig_non0_mean = mean(sig_non0)
	sig_non0_sd = sd(sig_non0)
	B_tmp = refsd/sig_non0_sd
	A_tmp = refmean - sig_non0_mean/sig_non0_sd*refsd
	signorm = sig * B_tmp + A_tmp
	signorm[signorm<0] = 0
	return(signorm)
}

file_list_file = args[1]#'H3K36me3.ctrl.list.txt'
file_list = read.table(file_list_file, header=F)

### get d1 sig
file_tmp = toString(file_list[1,1])
d1 = read.table(file_tmp, header=F, sep='\t')
dbed = d1[,1:3]
#dmat = d1[,4]
bin_num = dim(dbed)[1]
bin_num_used = 100000


dmat_s = matrix(0, nrow=bin_num_used, ncol=dim(file_list)[1])

used_id = sample(bin_num, bin_num_used)
dmat_s[,1]=d1[used_id,4]

### get dmat
for (i in 2:dim(file_list)[1]){
	file_tmp = toString(file_list[i,1])
	print(file_tmp)
	dtmp = read.table(file_tmp, header=F, sep='\t')[used_id,4]
	#dmat = cbind(dmat, dtmp)
	dmat_s[,i]=dtmp
	rm(dtmp)
}

### get average non0 mean sd
dmat_s_non0 = dmat_s[dmat_s!=0]
average_non0mean = mean(dmat_s_non0)
average_non0sd = sd(dmat_s_non0)


print(average_non0mean)
print(average_non0sd)
### non0scale
for (i in 1:dim(file_list)[1]){
	file_tmp = toString(file_list[i,1])
	output_name_tmp = paste(toString(file_list[i,1]), '.norm.bedgraph', sep='')
	dmat_sigi = read.table(file_tmp, header=F, sep='\t')[,4]
	if (length(unique(dmat_sigi))!=1){
		dmat_sigi_norm = non0scale(dmat_sigi, average_non0mean, average_non0sd)
	}else{
		dmat_sigi_norm = dmat_sigi
	}
	fwrite(cbind(dbed, dmat_sigi_norm), output_name_tmp, row.names=F, col.names=F, quote=F, sep='\t')
	rm(dmat_sigi)
}


