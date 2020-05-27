library(data.table)
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
d1 = as.data.frame(fread(file_tmp))
dbed = d1[,1:3]
#dmat = d1[,4]

dmat = matrix(0, nrow=dim(dbed)[1], ncol=dim(file_list)[1])
dmat[,1]=d1[,4]

### get dmat
for (i in 2:dim(file_list)[1]){
	file_tmp = toString(file_list[i,1])
	print(file_tmp)
	dtmp = as.data.frame(fread(file_tmp))[,4]
	#dmat = cbind(dmat, dtmp)
	dmat[,i]=dtmp
}

### get average non0 mean sd
dmat_non0 = dmat[dmat!=0]
average_non0mean = mean(dmat_non0)
average_non0sd = sd(dmat_non0)


print(average_non0mean)
print(average_non0sd)
### non0scale
for (i in 1:dim(file_list)[1]){
	output_name_tmp = paste(toString(file_list[i,1]), '.norm.bedgraph', sep='')
	dmat_sigi = dmat[,i]
	if (length(unique(dmat_sigi))!=1){
		dmat_sigi_norm = non0scale(dmat_sigi, average_non0mean, average_non0sd)
	}else{
		dmat_sigi_norm = dmat_sigi
	}
	fwrite(cbind(dbed, dmat_sigi_norm), output_name_tmp, row.names=F, col.names=F, quote=F, sep='\t')
}


