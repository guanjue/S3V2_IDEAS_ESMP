library(data.table)

### get parameters
args = commandArgs(trailingOnly=TRUE)
file_list_file = args[1]
output = args[2]
thresh = as.numeric(args[3])
method = args[4]
prob_pk = as.numeric(args[5])
prob_bg = as.numeric(args[6])

### get z
get_z_pre = function(x){
        x_notop = x[x<=quantile(x[x>0], 0.95)]
        xz = (x - mean(x_notop)) / sd(x_notop)
        return(xz)
}

get_z = function(x){
	x_notop = x[x<=quantile(x[x>0], 0.95)]
        xz = (x - mean(x_notop)) / sd(x_notop)
        return(xz)
}

### get fdr
get_fdr = function(x){
        z = get_z(x)
        zp = pnorm(-abs(z))
        zpfdr = p.adjust(zp)
        return(zpfdr)
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
get_nb_fdrp = function(d){
	ds = d
	ds_notop = ds[ds<=quantile(ds, 0.99)]
	ds_notop_probT_sizeT = get_true_NB_prob_size(ds_notop)
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


### read input
file_list = read.table(file_list_file, header=F)
set.seed(2020)
used_id = sample(13000000, 100000)
#used_id = 1:2000000
d00 = as.data.frame(fread(toString(file_list[1,1])))
common_pk = matrix(0, nrow=dim(d00)[1],ncol=dim(file_list)[1])
for (i in 1:dim(file_list)[1]){
	print(file_list[i,1])
	d10 = as.data.frame(fread(toString(file_list[i,1])))
	sig_tmp = d10[,4]
	#sig_tmp[sig_tmp<1]=0
	if (method =='z'){
	sig_tmp_fdr = get_fdr(sig_tmp)
	} else if (method =='nb'){
	sig_tmp = sig_tmp/min(sig_tmp[sig_tmp>0])
	sig_tmp_fdr = get_nb_fdrp(sig_tmp)
	}
	print(summary(sig_tmp_fdr))
	sig_tmp_fdr_pk = sig_tmp_fdr<thresh
	print(sum(sig_tmp_fdr_pk))
	common_pk[,i] = sig_tmp_fdr_pk
}

common_pk = common_pk*1
print(summary(common_pk))

cpk_num_thresh = round(prob_pk * dim(common_pk)[2])
common_pk_binary = apply(common_pk, 1, function(x) sum(x!=0))
common_pk_binary1 = common_pk_binary >= cpk_num_thresh
common_pk_binary = common_pk_binary1*1

cbg_num_thresh = round(prob_bg * dim(common_pk)[2])
common_bg_binary = apply(common_pk, 1, function(x) sum(x==0))
common_bg_binary1 = common_bg_binary >= cbg_num_thresh
common_bg_binary = common_bg_binary1*1

write.table(common_pk_binary, paste(output,'.cpk.txt', sep=''), quote=F, col.names=F, row.names=F, sep='\t')
write.table(common_bg_binary, paste(output,'.cbg.txt', sep=''), quote=F, col.names=F, row.names=F, sep='\t')


