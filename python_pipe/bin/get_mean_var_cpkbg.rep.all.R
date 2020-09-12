### r
library(data.table)

### functions
mse = function(x1,x2){
mse_x12 = mean((x1-x2)^2)
return(mse_x12)
}

r2 = function(x1,x2){
sstot = sum((x1-mean(x1))^2)
ssres = sum((x1-x2)^2)
r2_x12 = 1-ssres/sstot
return(r2_x12)
}

get_pk_bg_non0mean_non0var_mse_r2 = function(f1_file,f2_file, cpk_used, cbg_used){
	### read data
	f1 = as.data.frame(fread(f1_file))[,4]
	f2 = as.data.frame(fread(f2_file))[,4]
	### cpk
	f1cpk = f1[cpk_used]
	f2cpk = f2[cpk_used]
	f1cpknon0_log2 = log2(f1cpk[f1cpk!=0])
	f2cpknon0_log2 = log2(f2cpk[f2cpk!=0])
#	f1cpknon0 = (f1cpk[f1cpk>0])
#	f2cpknon0 = (f2cpk[f2cpk>0])
	r2cpk = max(r2(log2(f1cpk+1),log2(f2cpk+1)),r2(log2(f2cpk+1),log2(f1cpk+1)))
	### cbg
	f1cbg = f1[cbg_used]
	f2cbg = f2[cbg_used]
#	f1cbgnon0 = log2(f1cbg[f1cbg!=0])
#	f2cbgnon0 = log2(f2cbg[f2cbg!=0])
	f1cbgnon0 = (f1cbg[f1cbg>0])
	f2cbgnon0 = (f2cbg[f2cbg>0])
	r2cbg = max(r2(log2(f1cbg+1),log2(f2cbg+1)),r2(log2(f2cbg+1),log2(f1cbg+1)))
	### MSE
#	cpk_tmp = c(mean(f1cpknon0),mean(f2cpknon0), sd(f1cpknon0),sd(f2cpknon0), mse(log2(f1cpk+1),log2(f2cpk+1)))
	cpk_tmp = c(mean(f1cpknon0_log2),mean(f2cpknon0_log2), sd(f1cpknon0_log2),sd(f2cpknon0_log2), mse(log2(f1cpk+1),log2(f2cpk+1)))

#	cbg_tmp = c(mean(f1cbgnon0),mean(f2cbgnon0), sd(f1cbgnon0),sd(f2cbgnon0), mse(log2(f1cbg+1),log2(f2cbg+1)))
	cbg_tmp = c(log2(mean(f1cbgnon0)),log2(mean(f2cbgnon0)), log2(sd(f1cbgnon0)), log2(sd(f2cbgnon0)), mse(log2(f1cbg+1),log2(f2cbg+1)))
	non0mean_non0var_mse_tmp = c( cpk_tmp, cbg_tmp, mse(log2(f1+1),log2(f2+1)), r2cpk, r2cbg, max(r2(log2(f1+1),log2(f2+1)),r2(log2(f2+1),log2(f1+1))) )
	print(cpk_tmp)
	print(cbg_tmp)
	return(non0mean_non0var_mse_tmp)
}




mk_list = c('atac', 'ctcf', 'h3k4me3', 'h3k27ac')
#mk_list = c('ctcf', 'h3k4me3', 'h3k27ac')
mk_list = c('h3k27ac')

for (mk in mk_list){
RAW_mat = c()
TS_mat = c()
MA_mat = c()
QT_mat = c()
S3_mat = c()
S3V2_mat = c()
### file names
ref_file=paste("/storage/home/gzx103/scratch/S3V2norm_compare/merged/", mk, ".average_sig.bedgraph", sep='')
cpk_file=paste("/storage/home/gzx103/scratch/S3V2norm_compare/merged/", mk, "_commonpkfdr01_nb.cpk.txt", sep='')
cbg_file=paste("/storage/home/gzx103/scratch/S3V2norm_compare/merged/", mk, "_commonpkfdr01_nb.cbg.txt", sep='')
allpk_file=paste("/storage/home/gzx103/scratch/S3V2norm_compare/merged/", mk, "_commonpkfdr01_nb.allpk.txt", sep='')
ct_file=paste(mk, ".ct.list.txt", sep='')
### get cpk & cbg
cpk = scan(cpk_file)!=0
cbg = scan(cbg_file)!=0
ct_list = read.table(ct_file, header=F)
### get info
for (i in 1:dim(ct_list)[1]){
	ct = as.character(ct_list[i,1])
	print(ct)
	### RAW
	f1_file_raw = paste(ct, 'rep1.', mk,'.meanrc.txt.bedgraph', sep='')
	f2_file_raw = paste(ct, 'rep2.', mk,'.meanrc.txt.bedgraph', sep='')
	RAW_tmp_vec = get_pk_bg_non0mean_non0var_mse_r2(f1_file_raw, f2_file_raw, cpk, cbg)
	RAW_mat = rbind(RAW_mat, RAW_tmp_vec)
	### TS
	f1_file_TS = paste(ct, 'rep1.', mk,'.TS.bedgraph', sep='')
	f2_file_TS = paste(ct, 'rep2.', mk,'.TS.bedgraph', sep='')
	TS_tmp_vec = get_pk_bg_non0mean_non0var_mse_r2(f1_file_TS, f2_file_TS, cpk, cbg)
	TS_mat = rbind(TS_mat, TS_tmp_vec)
	### MA
	f1_file_MA = paste(ct, 'rep1.', mk,'.MA.bedgraph', sep='')
	f2_file_MA = paste(ct, 'rep2.', mk,'.MA.bedgraph', sep='')
	MA_tmp_vec = get_pk_bg_non0mean_non0var_mse_r2(f1_file_MA, f2_file_MA, cpk, cbg)
	MA_mat = rbind(MA_mat, MA_tmp_vec)
	### QT
	f1_file_QT = paste(ct, 'rep1.', mk,'.QT.bedgraph', sep='')
	f2_file_QT = paste(ct, 'rep2.', mk,'.QT.bedgraph', sep='')
	QT_tmp_vec = get_pk_bg_non0mean_non0var_mse_r2(f1_file_QT, f2_file_QT, cpk, cbg)
	QT_mat = rbind(QT_mat, QT_tmp_vec)
	### S3
	f1_file_S3 = paste(ct, 'rep1.', mk,'.S3.bedgraph.s3norm.bedgraph', sep='')
	f2_file_S3 = paste(ct, 'rep2.', mk,'.S3.bedgraph.s3norm.bedgraph', sep='')
	S3_tmp_vec = get_pk_bg_non0mean_non0var_mse_r2(f1_file_S3, f2_file_S3, cpk, cbg)
	S3_mat = rbind(S3_mat, S3_tmp_vec)
	### S3V2
	f1_file_S3V2 = paste(ct, 'rep1.', mk,'.S3V2.bedgraph', sep='')
	f2_file_S3V2 = paste(ct, 'rep2.', mk,'.S3V2.bedgraph', sep='')
	S3V2_tmp_vec = get_pk_bg_non0mean_non0var_mse_r2(f1_file_S3V2, f2_file_S3V2, cpk, cbg)
	S3V2_mat = rbind(S3V2_mat, S3V2_tmp_vec)
}
colnames_list = c('MPR1', 'MPR2', 'VPR1', 'VPR2', 'PE', 'MBR1', 'MBR2', 'VBR1', 'VBR2', 'BE', 'AE', 'R2P', 'R2B', 'R2')
rownames_list = ct_list[,1]
colnames(RAW_mat) = colnames_list
colnames(TS_mat) = colnames_list
colnames(MA_mat) = colnames_list
colnames(QT_mat) = colnames_list
colnames(S3_mat) = colnames_list
colnames(S3V2_mat) = colnames_list
rownames(RAW_mat) = rownames_list
rownames(TS_mat) = rownames_list
rownames(MA_mat) = rownames_list
rownames(QT_mat) = rownames_list
rownames(S3_mat) = rownames_list
rownames(S3V2_mat) = rownames_list

print(RAW_mat)
print(TS_mat)
print(MA_mat)
print(QT_mat)
print(S3_mat)
print(S3V2_mat)

write.table(RAW_mat, paste(mk, '.', 'RAW_mat.txt', sep=''))
write.table(TS_mat, paste(mk, '.', 'TS_mat.txt', sep=''))
write.table(MA_mat, paste(mk, '.', 'MA_mat.txt', sep=''))
write.table(QT_mat, paste(mk, '.', 'QT_mat.txt', sep=''))
write.table(S3_mat, paste(mk, '.', 'S3_mat.txt', sep=''))
write.table(S3V2_mat, paste(mk, '.', 'S3V2_mat.txt', sep=''))
}



















