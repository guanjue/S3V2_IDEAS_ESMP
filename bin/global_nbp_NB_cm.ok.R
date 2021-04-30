args = commandArgs(trailingOnly=TRUE)


###### get NB model prob and size
get_true_NB_prob_size = function(x, siglim){
	x_pass_lim = x[x>=siglim]
	m=mean(x_pass_lim)
	v=var(x_pass_lim)
        print('old m & v')
        print(c(m, v))
	passlim_num = sum(x>=siglim)
	for (i in 1:1){
		if (v<m){v=m+0.1}
		m_pre = m
		v_pre = v
		p = m/v
		s = m^2/(v-m)
		exp_siglim_p = dnbinom(0:siglim, s, p)
		exp_total_num = passlim_num/(1-sum(exp_siglim_p))
		### get new mean
		siglim_n_sum_for_m = 0
		for (j in 0:siglim){
			siglim_n_sum_for_m = siglim_n_sum_for_m + j* exp_siglim_p[j+1]*exp_total_num
		}
		m = (siglim_n_sum_for_m + m * passlim_num) / exp_total_num

		### get new var
		siglim_n_sum_for_v = 0
		for (j in 0:siglim){
			siglim_n_sum_for_v = siglim_n_sum_for_v + (j-m)^2 * exp_siglim_p[j+1]*exp_total_num
		}
		v = (siglim_n_sum_for_v + sum((x_pass_lim-m)^2)) / (exp_total_num-1)
		if (v<m){v=m+0.1}
		print(c(i,m,v))
		if(abs(m-m_pre)<0.001 & abs(v-v_pre)<0.00001) {break}
	}
	p = m/v
	s = m^2/(v-m)
	print('new m & v')
	print(c(m, v))
	return(c(p,s))
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
#file_tmp = paste(mk, '.average_sig.bedgraph.S3V2ref.bedgraph', sep='')
#file_tmp = paste(mk, '.average_sig.bedgraph.S3.bedgraph', sep='')
file_tmp = average_sig_file

AVE = read.table(average_sig_file, header=F, sep='\t')
bed = AVE[,1:3]
#AVEmat = matrix(0, nrow=dim(bed)[1],ncol=dim(mk_list)[1])#[used_row,]
#AVE0mat = AVE0[,4]
AVEmat = AVE[,4]
### cbg


### all cbg
#cbg = (cbg!=0)
######### get global NB bg model
### for S3 AVE
#AVE0mat[AVE0mat<1] = 0
AVEmat_cbg = as.numeric(AVEmat)
AVEmat_cbg[AVEmat_cbg<1] = 0
top_sigs = AVEmat_cbg[AVEmat_cbg>10]

print('sum non0s')
print(sum(AVEmat_cbg>0))

top_mean = mean(tail(sort(top_sigs),100))
scale_down = 200/top_mean
print(scale_down)


if (scale_down>1){
print('scale_down = 1')
scale_down = 1
}


NB_thresh = top_mean*scale_down*0.01
print('NB_thresh')
print(NB_thresh)
if (NB_thresh<1){
NB_thresh = 1
}
print(NB_thresh)

AVEmat_cbg = round(AVEmat_cbg)

min_non0 = min(AVEmat_cbg[AVEmat_cbg>0])
print(min_non0)
AVEmat_cbg0 = AVEmat_cbg
AVEmat_cbg = (AVEmat_cbg-min_non0)*scale_down+min_non0
AVEmat_cbg[AVEmat_cbg0==0] = 0
AVEmat_cbg[AVEmat_cbg<0] = 0
AVEmat_cbg[AVEmat_cbg>200] = 200

print(summary(AVEmat_cbg))
print('sum non0s')
print(sum(AVEmat_cbg>0))
#if (max(AVEmat_cbg)<1){

top_rm_thresh = quantile(AVEmat_cbg[AVEmat_cbg>0],0.95)
print('top_rm_thresh')
print(top_rm_thresh)
if (top_rm_thresh<=NB_thresh){
	top_rm_thresh = NB_thresh+1
}

AVEmat_cbg = AVEmat_cbg[AVEmat_cbg<top_rm_thresh]
#}
print(summary(AVEmat_cbg))
###### get NB model prob and size and p0
AVEmat_cbg_NBmodel = get_true_NB_prob_size(AVEmat_cbg, NB_thresh)

print('AVEmat_cbg_NBmodel:')
print(AVEmat_cbg_NBmodel)
AVEmat_cbg_size = AVEmat_cbg_NBmodel[2]
AVEmat_cbg_prob = AVEmat_cbg_NBmodel[1]
### set limit for prob
if (AVEmat_cbg_prob<0.001){
        AVEmat_cbg_prob = 0.001
}
if (AVEmat_cbg_prob>=0.999){
        AVEmat_cbg_prob = 0.999
}


### get output
bin_num = length(AVEmat)
obs_0_num = round(sum(AVEmat_cbg<1)/length(AVEmat_cbg)*length(AVEmat))

print(bin_num)
print(obs_0_num)

get_p_z = function(d, mean_i, sd_i){
        dz = (d - mean_i)/sd_i
#        dzp = ppois(d,mean_i, lower.tail=F)#pnorm(-(dz))
	dzp = pnorm(-(dz))
        return(dzp)
}

#########
### get output sigi
# for average signal
min_non0 = min(AVEmat[AVEmat>0])
IP_CTRL_tmp = cbind((AVEmat-min_non0)*scale_down+min_non0, rep(AVEmat_cbg_size,length(AVEmat)))
IP_nb_pval = pnbinom(IP_CTRL_tmp[,1], IP_CTRL_tmp[,2], AVEmat_cbg_prob, lower.tail=FALSE)


use_pois = FALSE

IP_nb_pval[IP_nb_pval<=1e-323] = 1e-323
IP_nb_pval[IP_nb_pval>1] = 1.0
IP_neglog10_nb_pval = -log10(IP_nb_pval)
IP_neglog10_nb_pval[IP_neglog10_nb_pval<0] = 0
neglog10_nb_pval_bedgraph = cbind(bed, IP_neglog10_nb_pval)
### write output
output_file_tmp = paste(average_sig_file, '.NBP.bedgraph', sep='')
write.table(neglog10_nb_pval_bedgraph, output_file_tmp, quote=FALSE, col.names=FALSE, row.names=FALSE, sep='\t')

    print('summary negative log10 NB p-value:')
    print(summary(IP_neglog10_nb_pval))
        print(sum(IP_neglog10_nb_pval>16)/length(IP_neglog10_nb_pval)*dim(bed)[1])
        print(sum(IP_neglog10_nb_pval>10)/length(IP_neglog10_nb_pval)*dim(bed)[1])
        print(sum(IP_neglog10_nb_pval>5)/length(IP_neglog10_nb_pval)*dim(bed)[1])
        print(sum(IP_neglog10_nb_pval>2)/length(IP_neglog10_nb_pval)*dim(bed)[1])
        print(sum(IP_neglog10_nb_pval>1)/length(IP_neglog10_nb_pval)*dim(bed)[1])

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

	print(summary(IP_tmp))
	### get CTRL
	file_tmp1 = toString(file_list[i,2])
	print(file_tmp1)	
	CTRL_tmp = read.table(file_tmp1, header=F, sep='\t')[,4]
	### get both
#	obs_0_num = sum(IP_tmp==0)
	CTRL_tmp_mean = mean(CTRL_tmp)
    if (use_pois){
    	CTRL_tmp_adj = (CTRL_tmp+1)/(CTRL_tmp_mean+1)
	} else{
	    CTRL_tmp_adj = (CTRL_tmp+1)/(CTRL_tmp_mean+1)*AVEmat_cbg_size
    }
    IP_CTRL_tmp = cbind(IP_tmp, CTRL_tmp_adj)

	IP_CTRL_tmp = cbind(IP_tmp, CTRL_tmp_adj)

	rm(IP_tmp0)
	rm(CTRL_tmp)
        ### get negative binomial p-value 
#	IP_nb_pval = apply(IP_CTRL_tmp, MARGIN=1, function(x) get_pval(x[1], bin_num, x[2], AVEmat_cbg_prob, obs_0_num) )
	if (use_pois){
		IP_nb_pval = ppois(IP_CTRL_tmp[,1], IP_CTRL_tmp[,2]*pois_mean0, lower.tail=F)
	} else{
		IP_nb_pval = pnbinom(IP_CTRL_tmp[,1], IP_CTRL_tmp[,2], AVEmat_cbg_prob, lower.tail=FALSE) #/ pnbinom(0, IP_CTRL_tmp[,2], AVEmat_cbg_prob, lower.tail=FALSE) * (bin_num-obs_0_num)/bin_num
	}

    ### remove extrame p-value
    IP_nb_pval[IP_nb_pval<=1e-323] = 1e-323
    IP_nb_pval[IP_nb_pval>1] = 1.0

    IP_neglog10_nb_pval = -log10(IP_nb_pval)
	IP_neglog10_nb_pval[IP_neglog10_nb_pval<0] = 0
    neglog10_nb_pval_bedgraph = cbind(bed, IP_neglog10_nb_pval)
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



