### get parameters
args = commandArgs(trailingOnly=TRUE)

input_target = args[1]
input_ref = args[2]
output_target = args[3]
threshold = as.numeric(args[4])
exp_win = as.numeric(args[5])
plot_hs = args[6]
for_ref = args[7]

rank_lim = as.numeric(args[8])
upperlim = as.numeric(args[9])
lowerlim = as.numeric(args[10])

p_method = args[11]

cpk_file = args[12]
cbg_file = args[13]

allpk_file = args[14]

#input_target = 'PBMC_rep1.H3K27ac.s3norm.0_100.chr16.bedgraph'
#input_ref = 'ERY_S002S3.H3K27ac.s3norm.0_100.chr16.bedgraph'
#output_target = 'PBMC_rep1.H3K27ac.s3norm.0_100.chr16.SSnorm.bedgraph'
#threshold = 0.01
#exp_win = 5
#plot_hs = 'F'

#time Rscript S2norm.R PBMC_rep1.H3K27ac.s3norm.0_100.chr16.bedgraph ERY_S002S3.H3K27ac.s3norm.0_100.chr16.bedgraph PBMC_rep1.H3K27ac.s3norm.0_100.chr16.SSnorm.bedgraph 0.1 5 F


### get NB model
get_true_NB_prob_size = function(x){
	m=mean(x[x>0]);
	m2=mean(x[x>0]^2);
	p0 = length(which(x==0)) / length(x);
	p = m/(m2-m^2 * (1-p0));
	s = m * (1 - p0) * p /(1-p);
	rt=c(p,s,p0);
	print(summary(x))
	print(summary(x[x>0]))
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
get_p_r1 = function(d){
	ds = d
	ds[ds<1] = 0
	ds_notop = ds[ds<=quantile(ds[ds>0], 0.99)]
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

get_p_z = function(d, notop_p){
        d_notop = d[d<=quantile(d, notop_p)]
        dz = (d - mean(d_notop))/sd(d_notop)
        dzp = pnorm(-(dz))
        dzpfdr = p.adjust(dzp,'fdr')
        return(dzpfdr)
}


get_p_r2 = function(d, dr1){
	ds_notop = dr1#[dr1<=quantile(dr1, 0.99)]
	ds_notop_probT_sizeT = get_true_NB_prob_size(ds_notop)
	ds_obs_0_num = sum(d == 0)
	bin_num = length(d)
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
	ds_nb_pval = apply(cbind(d), MARGIN=1, function(x) get_pval(x[1], bin_num, ds_size, ds_prob, ds_obs_0_num) )
	ds_nb_pval[ds_nb_pval>1] = 1
	### get fdr
	ds_nb_pval_fdr0 = p.adjust(ds_nb_pval, 'fdr')
	return(ds_nb_pval_fdr0)
}

getsf_pk = function(a0, b0){
        #print('getsf:')
        a = log2(a0[(a0>0)&(b0>0)])
        b = log2(b0[(a0>0)&(b0>0)])
        B = sd(b)/sd(a)
        A = mean(b) - mean(a)/sd(a)*sd(b)
        sf = c(B, A)
        print(sf)
        return(sf)
}

getsf_pk_LM_A0 = function(a0, b0){
        #print('getsf:')
	a = ((a0[(a0>0)&(b0>0)]))
	b = ((b0[(a0>0)&(b0>0)]))
	expLM_linear <- function(x) {
	sum((b-a^x)^2)
	}
	B = optimize(expLM_linear, lower=0, upper=5)$minimum
	sf = c(B, 0.0)
	print(sf)
	return(sf)
}


getsf_pk_LM_A00 = function(a0, b0){
        #print('getsf:')
        a = (log2(a0[(a0>0)&(b0>0)]))
        b = (log2(b0[(a0>0)&(b0>0)]))
        lm1 = lm(b~a)
        B = lm1$coefficients[1]
        A = 0.0#mean(b) - mean(a)/sd(a)*sd(b)
        sf = c(B, A)
        print(sf)
        return(sf)
}


getsf_pk_LM_A01 = function(a0, b0){
        #print('getsf:')
        a = ((a0[(a0>0)&(b0>0)]))
        b = ((b0[(a0>0)&(b0>0)]))
        lm1 = lm(b~a-1)
        #B = lm1$coefficients[1]
        A = lm1$coefficients[1]#mean(b) - mean(a)/sd(a)*sd(b)
        sf = c(1.0, A)
        print(sf)
        return(sf)
}

getsf_bg_iter = function(a0, b0){
	print('getsf_bg:')
	a = (a0[(a0>0)])
	b = (b0[(b0>0)])
	B = 1
	A = 0
	a_tmp = a
	sd_b = sd(b)
	mean_b = mean(b)
	for (i in 1:50){
	print(i)
	B_tmp = sd_b/sd(a_tmp[a_tmp>0])
	A_tmp = mean_b - mean(a_tmp[a_tmp>0])/sd(a_tmp[a_tmp>0])*sd_b
	B_pre = B
	A_pre = A

	B = B_pre * B_tmp
	A = A_pre * B_tmp + A_tmp
#	a_ttest = a0*B+A
#	print(mean(a_ttest[a_ttest>0]))
	a_tmp = a_tmp * B_tmp + A_tmp

	print(mean(a_tmp[a_tmp>0]))
	print(mean_b)
	if ((abs(log2(B/B_pre))<0.001) & (abs(A-A_pre)<0.001)){
	print(i)
	print(sd_b)
	print(sd(a_tmp[a_tmp>0]))
	break
	}
	print(c(B,A))
	}
	sf = c(B, A)
	print(sf)
	return(sf)
}

getsf_bg = function(a0, b0){
        print('getsf_bg:')
        a = (a0[(a0>0)])
        b = (b0[(b0>0)])
        B = 1
        A = 0
        a_tmp = a
        sd_b = sd(b)
        mean_b = mean(b)
        B_tmp = sd_b/sd(a_tmp[a_tmp>0])
        A_tmp = mean_b - mean(a_tmp[a_tmp>0])/sd(a_tmp[a_tmp>0])*sd_b
        B_pre = B
        A_pre = A
        B = B_pre * B_tmp
        A = A_pre * B_tmp + A_tmp
        a_tmp = a_tmp * B_tmp + A_tmp
        print(mean(a_tmp[a_tmp>0]))
        print(mean_b)
        print(sd_b)
        print(sd(a_tmp[a_tmp>0]))
        print(c(B,A))
        sf = c(B, A)
        print(sf)
        return(sf)
}


getsf = function(a0, b0){
	#print('getsf:')
#	a = log2(a0[(a0>0)&(b0>0)])
#	b = log2(b0[(a0>0)&(b0>0)])
	a = log2(a0[(a0>0)])
	b = log2(b0[(b0>0)])
	#a = a[a<quantile(a,0.99)]
	#b = b[b<quantile(b,0.99)]
	B = sd(b)/sd(a)
	A = mean(b) - mean(a)/sd(a)*sd(b)
	sf = c(B, A)
	print(sf)
	return(sf)
}


get_local_bg_sig_each_row = function(x, exp_win, d_lim, xsig, xbinary){
	### only change pk
	if (xbinary!=0){
		sig_tmp = x[(2*exp_win+1):length(x)]
		binary_tmp = x[1:(2*exp_win)]
		### only change pk
			if (sum(binary_tmp)==0){
				sig_bg_tmp = (max(sig_tmp) )#, xsig))
			} else if (prod(binary_tmp)!=0){
				sig_bg_tmp = (d_lim)#, xsig)
			} else{
				sig_bg_tmp = max(sig_tmp[binary_tmp==0])#, xsig)
				if(!is.finite(sig_bg_tmp)){print(x)}
			}
	} else {
		sig_bg_tmp = xsig #max(max(sig_tmp[binary_tmp==0]), xsig)
	}
	return(sig_bg_tmp)
}


### get local bg signal
get_local_bg_sig = function(exp_win, d_sig_all, d_pkb, d_lim){
	### get win sig pkb mat
	d_exp_sig_up = c()
	d_exp_sig_down = c()
	d_exp_sig_up_pkb = c()
	d_exp_sig_down_pkb = c()
	for (i in 1:exp_win){
		d_exp_sig_down = cbind( d_exp_sig_down, c(d_sig_all[(1+i):length(d_sig_all)], rep(0,i)) )
		d_exp_sig_up = cbind( d_exp_sig_up, c(rep(0,i), d_sig_all[1:(length(d_sig_all)-i)]) )
		d_exp_sig_down_pkb = cbind( d_exp_sig_down_pkb, c(d_pkb[(1+i):length(d_pkb)], rep(0,i)) )
		d_exp_sig_up_pkb = cbind( d_exp_sig_up_pkb, c(rep(0,i), d_pkb[1:(length(d_pkb)-i)]) )
	}
	### merge pkb & sig mat
	d_exp_sig = cbind(d_exp_sig_up, d_exp_sig_down)
	d_exp_pkb = cbind(d_exp_sig_up_pkb, d_exp_sig_down_pkb)
	d_exp_pkb_sig = cbind(d_sig_all, d_pkb, d_exp_pkb, d_exp_sig)
	rm(d_exp_sig)
	rm(d_exp_pkb)
	rm(d_exp_sig_up)
	rm(d_exp_sig_down)
	rm(d_exp_sig_up_pkb)
	rm(d_exp_sig_down_pkb)
	rm(d_sig_all)
	rm(d_pkb)
	### get local bg signal
	### the matrix first 1:(2*exp_win) is the binary info, the (2*exp_win+1):length(x) is the signal
#	print(head(d_exp_pkb_sig))
#	print(head(d_pkb))
#	print(sum(d_pkb))
#	print(length(d_pkb))
	d_exp_pkb_sig_bg_sig = apply(d_exp_pkb_sig, 1, function(x) get_local_bg_sig_each_row(x[3:length(x)], exp_win, d_lim, x[1], x[2]))
	### replace all pk region by gloabal bg
	d_exp_pkb_sig_bg_sig[!is.finite(d_exp_pkb_sig_bg_sig)] = d_lim
	return(d_exp_pkb_sig_bg_sig)
}

getmode = function(v) {
   uniqv = unique(v)
   uniqv[which.max(tabulate(match(v, uniqv)))]
}
### read signal
d10 = read.table(input_ref, header=F, sep='\t')
d20 = read.table(input_target, header=F, sep='\t')
d1 = d10[,4]
d2 = d20[,4]

rm(d10)
#d1_3rd_qt = quantile(d1[d1>0],0.75,type=1)
#if (d1_3rd_qt<1){
#        d1 = d1/d1_3rd_qt
#}


if (cpk_file=='F'){
	d2_quantile = quantile(d2[d2>0],0.75,type=1)
	if (d2_quantile<1){
		d2 = d2/d2_quantile
	}
} else{
	cpk_id = (scan(cpk_file)!=0)
	d2_cpk = d2[cpk_id]
	d2_cpk_non0_min = min(d2_cpk[d2_cpk>0])
	if (d2_cpk_non0_min<1){
		d2 = d2/d2_cpk_non0_min
	}	
}

#d1 = d1/(min(d1[d1>0]))
#d2 = d2/(min(d2[d2>0]))

d1[d1>upperlim] = upperlim
d1[d1<lowerlim] = lowerlim
d2[d2>upperlim] = upperlim
d2[d2<lowerlim] = lowerlim
### get FDR NB p-value vector
if (p_method=='neglog10p'){
	d1s_nb_pval_out = p.adjust(10^(-d1), 'fdr')
	d2s_nb_pval_out = p.adjust(10^(-d2), 'fdr')
} else if (p_method=='rc') {
	d1s_nb_pval_out = get_p_r1(d1)
	d2s_nb_pval_out = get_p_r1(d2)
	#d1r1 = d1[(d1s_nb_pval_out>threshold)]
	#d2r1 = d2[(d2s_nb_pval_out>threshold)]
	#d1s_nb_pval_out = get_p_r2(d1, d1r1)
	#d2s_nb_pval_out = get_p_r2(d2, d2r1)
}else if (p_method=='z') {
        d1s_nb_pval_out = get_p_z(d1, 1)
        d2s_nb_pval_out = get_p_z(d2, 1)
}


### get pk bg binary vec
d1s_nb_pval_out_binary_bg = (d1s_nb_pval_out>=threshold)
d2s_nb_pval_out_binary_bg = (d2s_nb_pval_out>=threshold)
d1s_nb_pval_out_binary_pk = (d1s_nb_pval_out<threshold)
d2s_nb_pval_out_binary_pk = (d2s_nb_pval_out<threshold)


### get pk sf
d12_bgb = d1s_nb_pval_out_binary_bg | d2s_nb_pval_out_binary_bg
d1bg = d1[d12_bgb] 
d2bg = d2[d12_bgb] 
d1s_lim = mean(d1bg)
d2s_lim = mean(d2bg)

d12_allpkb = d1s_nb_pval_out_binary_pk | d2s_nb_pval_out_binary_pk

### get local bg
d1_pkb = d1s_nb_pval_out_binary_pk
d2_pkb = d2s_nb_pval_out_binary_pk
d12_pkb = d1s_nb_pval_out_binary_pk | d2s_nb_pval_out_binary_pk
print('ref peak num')
print(sum(d1s_nb_pval_out_binary_pk))
print('tar peak num')
print(sum(d2s_nb_pval_out_binary_pk))

### get all pk
if (allpk_file!='F'){
#all_pk = as.data.frame(fread(allpk_file))
#used_id = all_pk[,4]!=0
all_pk = scan(allpk_file)
print(summary(all_pk))
allpk_used_id = all_pk!=0
}else{
#d2pk_sig_z = (d2-mean(d2))/sd(d2)
#d2pk_sig_zp = 2*pnorm(-(d2pk_sig_z))
#d2pk_sig_zp_fdr = p.adjust(d2pk_sig_zp, 'fdr')
#potential_pk = d2pk_sig_zp_fdr<threshold
allpk_used_id = d2s_nb_pval_out_binary_pk#(potential_pk)|(d12_pkb)
}

d1_sig_all = d1
d1_pkb = allpk_used_id#d1s_nb_pval_out_binary_pk|d2s_nb_pval_out_binary_pk
if (length(d1_sig_all)<100000){
	d1_exp_pkb_sig_bg_sig = get_local_bg_sig(exp_win, d1_sig_all, allpk_used_id, d1s_lim)
} else{
	d1_exp_pkb_sig_bg_sig = rep(0, length(d1_sig_all))
	split_range = cbind(seq(1,length(d1_sig_all), by=100000), c(seq(1,length(d1_sig_all), by=100000)[-1],length(d1_sig_all)) )
	for (i in 1:dim(split_range)[1]){
		used_id_i = split_range[i,1]:split_range[i,2]
		d1_exp_pkb_sig_bg_sig[used_id_i] = get_local_bg_sig(exp_win, d1_sig_all[used_id_i], allpk_used_id[used_id_i], d1s_lim)
	}

}

d2_sig_all = d2
d2_pkb = allpk_used_id#d1s_nb_pval_out_binary_pk|d2s_nb_pval_out_binary_pk
if (length(d1_sig_all)<100000){
	d2_exp_pkb_sig_bg_sig = get_local_bg_sig(exp_win, d2_sig_all, allpk_used_id, d2s_lim)
} else{
	d2_exp_pkb_sig_bg_sig = rep(0, length(d2_sig_all))
	split_range = cbind(seq(1,length(d2_sig_all), by=100000), c(seq(1,length(d2_sig_all), by=100000)[-1],length(d2_sig_all)) )
	for (i in 1:dim(split_range)[1]){
		used_id_i = split_range[i,1]:split_range[i,2]
		d2_exp_pkb_sig_bg_sig[used_id_i] = get_local_bg_sig(exp_win, d2_sig_all[used_id_i], allpk_used_id[used_id_i], d2s_lim)
	}

}


if (for_ref == 'T'){
	### get pk sf
	d1pk_sig = d1_sig_all - d1_exp_pkb_sig_bg_sig
	d2pk_sig = d2_sig_all - d2_exp_pkb_sig_bg_sig
	d1pk_sig[d1pk_sig<0]=0
	d2pk_sig[d2pk_sig<0]=0
	print('PK SF')
	pksf = getsf_pk(d2pk_sig, d1pk_sig)	
	### get bg sf
        d1bg_sig = d1_exp_pkb_sig_bg_sig
        d2bg_sig = d2_exp_pkb_sig_bg_sig
        print('BG SF')
	bgsf_bg = getsf_bg(d2bg_sig, d1bg_sig)
#	bgsf = getsf(d2bg_sig[cbg_id], d1_sig_all[cbg_id])
#	bgsf = getsf(d2bg_sig, d1bg_sig)
} else{
	### cpk
	if (cpk_file=='F'){
	cpk_id = (d1s_nb_pval_out_binary_pk) & (d2s_nb_pval_out_binary_pk)
	} else{
	cpk_id = (scan(cpk_file)!=0)
	}
	if (cbg_file=='F'){
	cbg_id = (d1s_nb_pval_out_binary_bg) & (d2s_nb_pval_out_binary_bg)
	} else{
	cbg_id = (scan(cbg_file)!=0)
	}
	### get pk sf
	d1pk_sig = d1_sig_all - d1_exp_pkb_sig_bg_sig
	d2pk_sig = d2_sig_all - d2_exp_pkb_sig_bg_sig
	### get shoulder
	print('ref non bg num:')
	print(sum(d1pk_sig>0))
	print('tar non bg num:')
	print(sum(d2pk_sig>0))
	d1pk_sig[d1pk_sig<0]=0
	d2pk_sig[d2pk_sig<0]=0
	print('PK SF')
	print(sum(cpk_id))
	print(summary(d2pk_sig[cpk_id]))
	print(summary(d1[cpk_id]))
	if (sum(cpk_id)>(rank_lim*length(cpk_id))){
		pksf = getsf_pk_LM_A0(d2pk_sig[cpk_id], d1[cpk_id])
	} else {
		cpk_id = (d1s_nb_pval_out_binary_pk) & (d2s_nb_pval_out_binary_pk)
		pksf = getsf_pk_LM_A0(d2pk_sig[cpk_id], d1[cpk_id])
	}
	### get bg sf
	d1bg_sig = d1_exp_pkb_sig_bg_sig
	d2bg_sig = d2_exp_pkb_sig_bg_sig
	d2pk_sig_norm0 = (d2pk_sig^pksf[1]) * (2^pksf[2])
	print('BG SF')
	bgsf_bg = getsf_bg(d2bg_sig[cbg_id], d1[cbg_id])
#	bgsf_pk = getsf_pk(d2bg_sig[cbg_id], d1bg_sig[cbg_id])
}

### norm bg
#d2bg_sig_norm = (d2bg_sig^bgsf_bg[1]) * (2^bgsf_bg[2])
d2bg_sig_norm = (d2bg_sig*bgsf_bg[1]) + (bgsf_bg[2])
d2bg_sig_norm[d2bg_sig==0] = 0

#print(mean(d2bg_sig_norm[cbg_id][d2bg_sig_norm[cbg_id]>0]))
#print(mean(d1[cbg_id][d1[cbg_id]>0]))

#print('summary(d2bg_sig_norm)')
#print(mean(d1[cbg_id][d1[cbg_id]>0]))
#print(mean(d1bg_sig[cbg_id][d1bg_sig[cbg_id]>0]))
#print(sd(d1[cbg_id][d1[cbg_id]>0]))
#print(sd(d1bg_sig[cbg_id][d1bg_sig[cbg_id]>0]))

#print(summary(d2bg_sig_norm))
#print(sum(d2bg_sig_norm<0))
#print(sd(d2bg_sig_norm[cbg_id][d2bg_sig_norm[cbg_id]>0]))
#print(sd(d1[cbg_id][d1[cbg_id]>0]))
#print(mean(d2bg_sig_norm[cbg_id][d2bg_sig_norm[cbg_id]>0]))
#print(mean(d1[cbg_id][d1[cbg_id]>0]))




### norm pk
d2pk_sig_norm = d2pk_sig
print('used pk num')
print(sum(allpk_used_id))
print(summary(d2pk_sig_norm))
d2pk_sig_norm[allpk_used_id] = (d2pk_sig_norm[allpk_used_id]^pksf[1]) * (2^pksf[2])
#d2pk_sig_norm[d2s_nb_pval_out_binary_pk] = (d2pk_sig_norm[d2s_nb_pval_out_binary_pk]^pksf[1]) * (2^pksf[2])
#d2pk_sig_norm[used_id] = (d2pk_sig_norm[used_id]*pksf[1]) * (pksf[2])
#d2pk_sig_norm[d2pk_sig<0] = 0
print(sum(d2pk_sig<0))

### pk bg modify
#d2bg_sig_norm[used_id] = (d2bg_sig[used_id]^bgsf_pk[1]) * (2^bgsf_pk[2])
#d2bg_sig_norm[used_id] = (d2bg_sig[used_id]*bgsf_pk[1]) + (bgsf_pk[2])

print('signal head')
print(head(d2_sig_all))
print(head(d2pk_sig))
print(head(d2pk_sig_norm))

print(head(d2bg_sig))
print(head(d2bg_sig_norm))

print(summary(d2bg_sig_norm))
print(summary(d2pk_sig_norm))

### add pk & bg signal
d2_sig_norm = d2bg_sig_norm + d2pk_sig_norm

### remove 0s
d2_sig_norm[d2_sig_norm>upperlim] = upperlim
d2_sig_norm[d2_sig_norm<lowerlim] = lowerlim

print(summary(d2_sig_norm))

#d2_sig_norm_nb_pval_out = get_p_r1(d2_sig_norm)

### write output
write.table(cbind(d20[,1:3], d2_sig_norm), output_target, sep='\t', quote=F, col.names=F, row.names=F)
write.table(c(pksf[1], (2^pksf[2]), bgsf_bg[1], bgsf_bg[2]), paste(output_target, '.info.txt', sep=''), sep='\t', quote=F, col.names=F, row.names=F)

rm(d20)
#write.table(cbind(d20[,1:3], d2pk_sig_norm), paste(output_target, '.pk.sig.bedgraph', sep=''), sep='\t', quote=F, col.names=F, row.names=F)
#write.table(cbind(d20[,1:3], d2bg_sig_norm), paste(output_target, '.bg.sig.bedgraph', sep=''), sep='\t', quote=F, col.names=F, row.names=F)
#write.table(cbind(d20[,1:3], d2pk_sig), paste(output_target, '.pkraw.sig.bedgraph', sep=''), sep='\t', quote=F, col.names=F, row.names=F)
#write.table(cbind(d20[,1:3], d2bg_sig), paste(output_target, '.bgraw.sig.bedgraph', sep=''), sep='\t', quote=F, col.names=F, row.names=F)


#print('check!')
#print('mean')
#print(log2(mean(d1_sig_all[cbg_id][d1_sig_all[cbg_id]>0])))
#print(log2(mean(d2_sig_norm[cbg_id][d2_sig_norm[cbg_id]>0])))
#print(mean(log2(d1_sig_all[cpk_id])))
#print(mean(log2(d2_sig_norm[cpk_id])))
#print(mean(log2(d2_sig_norm[d2_sig_norm_nb_pval_out<threshold])))
#print(mean(d2pk_sig_norm[cbg_id]))
#print(mean(log2(d1_exp_pkb_sig_bg_sig[d1_exp_pkb_sig_bg_sig>0])))
#print(mean(log2(d2_exp_pkb_sig_bg_sig[d2_exp_pkb_sig_bg_sig>0])))
#print('sd')
#print(sd((d1_sig_all[cbg_id][d1_sig_all[cbg_id]>0])))
#print(sd((d2_sig_norm[cbg_id][d2_sig_norm[cbg_id]>0])))
#print(sd((d2bg_sig_norm[cbg_id][d2bg_sig_norm[cbg_id]>0])))
#print(sd((d2bg_sig_norm[d2_sig_norm_nb_pval_out>=threshold][d2bg_sig_norm[d2_sig_norm_nb_pval_out>=threshold]>0])))
#print(sd(log2(d1_sig_all[cpk_id])))
#print(sd(log2(d2_sig_norm[cpk_id])))
#print(sd((d1_exp_pkb_sig_bg_sig[d1_exp_pkb_sig_bg_sig>0])))
#print(sd((d2_exp_pkb_sig_bg_sig[d2_exp_pkb_sig_bg_sig>0])))


### plot heatscatter
if (plot_hs == 'T'){
	set.seed(2019)
	used_id = sample(dim(d1)[1], 100000)
	library(LSD)
	png(paste(output_target, '.png', sep=''), width=1000)
	par(mfrow=c(1,2))
	heatscatter(d1[used_id], d2[used_id], log='xy', xlim=c(0.1, 1000), ylim=c(0.1, 1000))
	abline(0,1)
	heatscatter(d1[used_id], d2_sig_norm[used_id], log='xy', xlim=c(0.1, 1000), ylim=c(0.1, 1000))
	abline(0,1)
	dev.off()
}





