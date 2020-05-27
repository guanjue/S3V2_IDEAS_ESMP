### r
library(data.table)
library(ggplot2)
library(RColorBrewer)

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

### get NB model
get_true_NB_prob_size = function(x){
        m=mean(x[x>0]);
        m2=mean(x[x>0]^2);
        p0 = length(which(x==0)) / length(x);
        p = m/(m2-m^2 * (1-p0));
        s = m * (1 - p0) * p /(1-p);
        rt=c(p,s,p0);

        for(i in 1:100){
                op = p;
                os = s;
                p0=p^s;
                #print(p0)
                p=m/(m2-m^2*(1-p0));
                if (p<0.001){
                        p = 0.001
                }
                if (p>=0.999){
                        p = 0.999
                }
                s=m * (1 - p0) * p / (1-p);
                #rt=rbind(rt,c(p,s,p0));
                rt = c(p,s,p0)
                if(abs(op-p)<0.00001 & abs(os-s)<0.00001) break;
        }
        #print('change best_p0: ')
        #print(p0)
        return(rt);
}

### get p-value
get_pval = function(N, l, sig_0_size, sig_0_prob, num_0){
        if (N != 0){
                pval_new = pnbinom(N-1, sig_0_size, sig_0_prob, lower.tail=FALSE) / pnbinom(0, sig_0_size, sig_0_prob, lower.tail=FALSE) * (l-num_0)/l
        } else {
                pval_new = 1.0
        }
        return(pval_new)
}

### get fdr p-value vector
get_p = function(d){
        ds = d
        ds_notop = ds[ds<=quantile(ds, 0.99)]
        ds_notop_probT_sizeT = get_true_NB_prob_size(ds_notop)
        print(ds_notop_probT_sizeT)
        ds_obs_0_num = sum(ds == 0)
        bin_num = length(ds)
        ### get NB para
        ds_p0 = ds_notop_probT_sizeT[3]
        ds_size = ds_notop_probT_sizeT[2]
        ds_prob = ds_notop_probT_sizeT[1]
        ### set limit for prob
        if (ds_prob<0.001){
                ds_prob = 0.001
        }
        if (ds_prob>=0.999){
                ds_prob = 0.999
        }
        ### get p
        ds_nb_pval = apply(cbind(ds), MARGIN=1, function(x) get_pval(x[1], bin_num, ds_size, ds_prob, ds_obs_0_num) )
        ds_nb_pval[ds_nb_pval>1] = 1
        ### get fdr
        ds_nb_pval_fdr0 = p.adjust(ds_nb_pval, 'fdr')
        return(ds_nb_pval_fdr0)
}

get_pk_bg_non0mean_non0var_mse_r2 = function(f1_file,f2_file, cpk_used, cbg_used){
	### read data
	f1od = as.data.frame(fread(f1_file))
	f2od = as.data.frame(fread(f2_file))
	f1 = f1od[f1od[,1]=='chr1',4]
	f2 = f2od[f2od[,1]=='chr1',4]
	f1fdr = get_p(f1)
	f2fdr = get_p(f2)
	### cpk
#	f1cpk = f1[cpk_used]
#	f2cpk = f2[cpk_used]
	pk_binary1 = f1fdr<0.1
	f1pk = f1[pk_binary1]
	pk_binary2 = f2fdr<0.1
	f2pk = f2[pk_binary2]
	f1pknon0 = log2(f1pk[f1pk!=0])
	f2pknon0 = log2(f2pk[f2pk!=0])
	f1cpk = f1[pk_binary1 & pk_binary2]
	f2cpk = f2[pk_binary1 & pk_binary2]
	r2cpk = max(r2(log2(f1cpk+1),log2(f2cpk+1)),r2(log2(f2cpk+1),log2(f1cpk+1)))
	### cbg
#	f1cbg = f1[cbg_used]
#	f2cbg = f2[cbg_used]
	bg_binary1 = f1fdr>=0.1
	f1bg = f1[bg_binary1]
	bg_binary2 = f2fdr>=0.1
	f2bg = f2[bg_binary2]
	f1bgnon0 = log2(f1bg[f1bg!=0])
	f2bgnon0 = log2(f2bg[f2bg!=0])
	f1cbg = f1[bg_binary1 | bg_binary2]
	f2cbg = f2[bg_binary1 | bg_binary2]
	r2cbg = max(r2(log2(f1cbg+1),log2(f2cbg+1)),r2(log2(f2cbg+1),log2(f1cbg+1)))
	### MSE
	cpk_tmp = c(mean(f1pknon0),mean(f2pknon0), var(f1pknon0),var(f2pknon0), mse(log2(f1cpk+1),log2(f2cpk+1)))
	cbg_tmp = c(mean(f1bgnon0),mean(f2bgnon0), var(f1bgnon0),var(f2bgnon0), mse(log2(f1cbg+1),log2(f2cbg+1)))
	non0mean_non0var_mse_tmp = c( cpk_tmp, cbg_tmp, mse(log2(f1+1),log2(f2+1)), r2cpk, r2cbg, max(r2(log2(f1+1),log2(f2+1)),r2(log2(f2+1),log2(f1+1))) )
	return(non0mean_non0var_mse_tmp)
}


mk_list = c('atac', 'ctcf', 'h3k4me3', 'h3k27ac')

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
print(RAW_mat)
print(TS_mat)
print(MA_mat)
print(QT_mat)
print(S3_mat)
print(S3V2_mat)
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
write.table(RAW_mat, paste(mk, '.', 'RAW_mat.dpk.txt', sep=''))
write.table(TS_mat, paste(mk, '.', 'TS_mat.dpk.txt', sep=''))
write.table(MA_mat, paste(mk, '.', 'MA_mat.dpk.txt', sep=''))
write.table(QT_mat, paste(mk, '.', 'QT_mat.dpk.txt', sep=''))
write.table(S3_mat, paste(mk, '.', 'S3_mat.dpk.txt', sep=''))
write.table(S3V2_mat, paste(mk, '.', 'S3V2_mat.dpk.txt', sep=''))
}

