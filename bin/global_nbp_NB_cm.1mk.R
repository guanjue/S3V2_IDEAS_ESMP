args = commandArgs(trailingOnly=TRUE)


###### get NB model prob and size
get_true_NB_prob_size = function(x){
	m=mean(x[x>0]);
	m2=mean(x[x>0]^2);
	x[x<1]=0
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

###### get NB model p-value
get_pval = function(N, l, sig_0_size, sig_0_prob, num_0){
	if (N != 0){
		pval_new = pnbinom(N-1, sig_0_size, sig_0_prob, lower.tail=FALSE) / pnbinom(0, sig_0_size, sig_0_prob, lower.tail=FALSE) * (l-num_0)/l
	} else {
		pval_new = 1.0
	}
	return(pval_new)
}

file_list_file = args[1] #'H3K36me3.file_list.txt'
average_sig_file = args[2]
cbg_file = args[3]

file_list = read.table(file_list_file, header=F)

######## read average signal
print('read average signal')
mk = unlist(strsplit(file_list_file, split='\\.'))[1]
file_tmp = average_sig_file

AVE = read.table(file_tmp, header=F, sep='\t')
bed = AVE[,1:3]
AVEmat = AVE[,4]
### cbg


### all cbg
#cbg = (cbg!=0)
######### get global NB bg model
AVEmat_cbg = as.numeric(AVEmat)
top_sigs = AVEmat_cbg[AVEmat_cbg>10]
scale_down = 200/mean(tail(sort(top_sigs),100))

print(scale_down)
print('Summary signal')
###### get NB model prob and size and p0
AVEmat_cbg = round(AVEmat_cbg)

print(summary(AVEmat_cbg))

min_non0 = min(AVEmat_cbg[AVEmat_cbg>0])
print(min_non0)
AVEmat_cbg = (AVEmat_cbg-min_non0)*scale_down+min_non0
AVEmat_cbg[AVEmat_cbg<0] = 0
AVEmat_cbg[AVEmat_cbg>200] = 200
### get output


#########
### get output sigi
# for average signal
print(summary(AVEmat_cbg))

#use_pois = FALSE
#if ((sum(IP_nb_pval<=(1e-16))<(length(IP_nb_pval)*(1e-05))) | (sum(IP_nb_pval<=(1e-16))>(length(IP_nb_pval)*0.2))){
	use_pois = TRUE
	inflate_thresh = 2
	AVEmat_cbg_non0_num = sum(AVEmat_cbg>inflate_thresh)
	pois_mean0_non0 = mean(AVEmat_cbg[AVEmat_cbg>inflate_thresh])
	pois_mean_all = mean(AVEmat_cbg[AVEmat_cbg>inflate_thresh])
	print(pois_mean_all)
	print('fit pois')
	for (i in 1:100){
		exp_inflate = sum(dpois(0:inflate_thresh, pois_mean_all))
		#print(exp_inflate)
		AVEmat_cbg_all_num = AVEmat_cbg_non0_num / (1-exp_inflate)
		AVEmat_cbg_othersInflate_binnum_sig = 0
		for (k in 1:inflate_thresh){
			exp_inflatek = dpois(k, pois_mean_all)
			AVEmat_cbg_othersInflate_binnum_sig = AVEmat_cbg_othersInflate_binnum_sig+AVEmat_cbg_all_num*exp_inflatek*k
		}
		pois_mean_all_new = (pois_mean0_non0*AVEmat_cbg_non0_num+AVEmat_cbg_othersInflate_binnum_sig) / AVEmat_cbg_all_num
		if (abs(pois_mean_all-pois_mean_all_new)>0.001){
			pois_mean_all = pois_mean_all_new
			print(pois_mean_all)
		} else{
			break
		}
	}
	print(pois_mean_all)
	IP_nb_pval = ppois(AVEmat*scale_down, pois_mean_all, lower.tail=F)
#}
print('use_pois:')
print(use_pois)
IP_nb_pval[IP_nb_pval<=1e-323] = 1e-323
IP_nb_pval[IP_nb_pval>1] = 1.0
IP_neglog10_nb_pval = -log10(IP_nb_pval)
IP_neglog10_nb_pval[IP_neglog10_nb_pval<0] = 0

        print('summary negative log10 NB p-value:')
        print(summary(IP_neglog10_nb_pval))
        print(sum(IP_neglog10_nb_pval>16)/length(IP_neglog10_nb_pval)*dim(bed)[1])
        print(sum(IP_neglog10_nb_pval>10)/length(IP_neglog10_nb_pval)*dim(bed)[1])
        print(sum(IP_neglog10_nb_pval>5)/length(IP_neglog10_nb_pval)*dim(bed)[1])
        print(sum(IP_neglog10_nb_pval>2)/length(IP_neglog10_nb_pval)*dim(bed)[1])
        print(sum(IP_neglog10_nb_pval>1)/length(IP_neglog10_nb_pval)*dim(bed)[1])

neglog10_nb_pval_bedgraph = cbind(bed, IP_neglog10_nb_pval)
### write output
output_file_tmp = paste(average_sig_file, '.NBP.bedgraph', sep='')
write.table(neglog10_nb_pval_bedgraph, output_file_tmp, quote=FALSE, col.names=FALSE, row.names=FALSE, sep='\t')


# for each sample
for (i in 1:dim(file_list)[1]){
	print(i)
	### get IP
	file_tmp = toString(file_list[i,1])
	print(file_tmp)
	IP_tmp0 = read.table(file_tmp, header=F, sep='\t')
	if (i==1){
		bed = IP_tmp0[,1:3]
	}

        IP_tmp = IP_tmp0[,4]
        min_non0 = min(IP_tmp[IP_tmp>0])
        IP_tmp = (IP_tmp-min_non0)*scale_down+min_non0
	IP_tmp[IP_tmp<0] = 0
	print(summary(IP_tmp))
	### get CTRL
	file_tmp1 = toString(file_list[i,2])
	print(file_tmp1)	
	CTRL_tmp = read.table(file_tmp1, header=F, sep='\t')[,4]
	### get both
	CTRL_tmp_mean = mean(CTRL_tmp)
	if (use_pois){
		CTRL_tmp_adj = (CTRL_tmp+1)/(CTRL_tmp_mean+1)
	} else{
		CTRL_tmp_adj = (CTRL_tmp+1)/(CTRL_tmp_mean+1)*AVEmat_cbg_size
	}
	IP_CTRL_tmp = cbind(IP_tmp, CTRL_tmp_adj)

	rm(IP_tmp0)
	rm(CTRL_tmp)
	### get negative binomial p-value 
	if (use_pois){
		IP_nb_pval = ppois(IP_CTRL_tmp[,1], IP_CTRL_tmp[,2]*pois_mean_all, lower.tail=F)
	} else{
		IP_nb_pval = pnbinom(IP_CTRL_tmp[,1]-1, IP_CTRL_tmp[,2], AVEmat_cbg_prob, lower.tail=FALSE)
	}

	### remove extrame p-value
	IP_nb_pval[IP_nb_pval<=1e-323] = 1e-323
	IP_nb_pval[IP_nb_pval>1] = 1.0

	IP_neglog10_nb_pval = -log10(IP_nb_pval)
	IP_neglog10_nb_pval[IP_neglog10_nb_pval<0] = 0
	neglog10_nb_pval_bedgraph = cbind(bed, round(IP_neglog10_nb_pval, 3))
	print(toString(file_list[i,1]))
	print('summary negative log10 NB p-value:')
	print(summary(IP_neglog10_nb_pval))
	print(sum(IP_neglog10_nb_pval>16)/length(IP_neglog10_nb_pval)*dim(bed)[1])
	print(sum(IP_neglog10_nb_pval>10)/length(IP_neglog10_nb_pval)*dim(bed)[1])
	print(sum(IP_neglog10_nb_pval>5)/length(IP_neglog10_nb_pval)*dim(bed)[1])
	print(sum(IP_neglog10_nb_pval>2)/length(IP_neglog10_nb_pval)*dim(bed)[1])
	print(sum(IP_neglog10_nb_pval>1)/length(IP_neglog10_nb_pval)*dim(bed)[1])
	### write output
	output_file_tmp = paste(toString(file_list[i,1]), '.NBP.bedgraph', sep='')
	write.table(neglog10_nb_pval_bedgraph, output_file_tmp, quote=FALSE, col.names=FALSE, row.names=FALSE, sep='\t')
}



