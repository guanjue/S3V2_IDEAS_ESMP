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
average_sig_file0 = args[2]
average_sig_file = args[3]
cbg_file = args[4]

file_list = read.table(file_list_file, header=F)

######## read average signal
print('read average signal')
mk = unlist(strsplit(file_list_file, split='\\.'))[1]
#file_tmp = paste(mk, '.average_sig.bedgraph.S3V2ref.bedgraph', sep='')
#file_tmp = paste(mk, '.average_sig.bedgraph.S3.bedgraph', sep='')
file_tmp = average_sig_file

AVE0 = read.table(average_sig_file0, header=F, sep='\t')
AVE = read.table(average_sig_file, header=F, sep='\t')
bed = AVE[,1:3]
#AVEmat = matrix(0, nrow=dim(bed)[1],ncol=dim(mk_list)[1])#[used_row,]
AVE0mat = AVE0[,4]
AVEmat = AVE[,4]
### cbg
#cbg_file = paste(mk, '_commonpkfdr01_z.cbg.txt', sep='')
cbg = scan(cbg_file)!=0


### all cbg
#cbg = (cbg!=0)
######### get global NB bg model
print(sum(AVE0mat<1))
print(summary(AVE0mat[AVE0mat<1]))
print(sum(AVE0mat<1))
### for S3 AVE
AVE0mat[AVE0mat<1] = 0
AVE0mat_cbg = as.numeric(AVE0mat[cbg])
if (max(AVE0mat_cbg)<1){
AVE0mat_cbg = as.numeric(AVE0mat)
AVE0mat_cbg = AVE0mat_cbg[AVE0mat_cbg<quantile(AVE0mat_cbg[AVE0mat_cbg>0],0.95)]
}
### for S3V2 AVE
AVEmat[AVEmat<1] = 0
AVEmat_cbg = as.numeric(AVEmat[cbg])
if (max(AVEmat_cbg)<1){
AVEmat_cbg = as.numeric(AVEmat)
AVEmat_cbg = AVEmat_cbg[AVEmat_cbg<quantile(AVEmat_cbg[AVEmat_cbg>0],0.95)]
}

###### get NB model prob and size and p0
AVEmat_cbg_NBmodel = get_true_NB_prob_size(AVE0mat_cbg)

print('AVEmat_cbg_NBmodel:')
print(AVEmat_cbg_NBmodel)
AVEmat_cbg_p0 = AVEmat_cbg_NBmodel[3]
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
IP_CTRL_tmp = cbind(AVEmat, rep(AVEmat_cbg_size,length(AVEmat)))
IP_nb_pval = pnbinom(IP_CTRL_tmp[,1]-1, IP_CTRL_tmp[,2], AVEmat_cbg_prob, lower.tail=FALSE)
IP_nb_pval[IP_nb_pval<=1e-323] = 1e-323
IP_nb_pval[IP_nb_pval>1] = 1.0
IP_neglog10_nb_pval = -log10(IP_nb_pval)
IP_neglog10_nb_pval[IP_neglog10_nb_pval<0] = 0
neglog10_nb_pval_bedgraph = cbind(bed, IP_neglog10_nb_pval)
### write output
output_file_tmp = paste(average_sig_file, '.NBP.bedgraph', sep='')
fwrite(neglog10_nb_pval_bedgraph, output_file_tmp, quote=FALSE, col.names=FALSE, row.names=FALSE, sep='\t')


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
	print(summary(IP_tmp))
	### get CTRL
	file_tmp1 = toString(file_list[i,2])
	print(file_tmp1)	
	CTRL_tmp = read.table(file_tmp1, header=F, sep='\t')[,4]
	### get both
#	obs_0_num = sum(IP_tmp==0)
	CTRL_tmp_mean = mean(CTRL_tmp)
        CTRL_tmp_adj = (CTRL_tmp+1)/(CTRL_tmp_mean+1)*AVEmat_cbg_size
        IP_CTRL_tmp = cbind(IP_tmp, CTRL_tmp_adj)

        rm(IP_tmp0)
	rm(CTRL_tmp)
        ### get negative binomial p-value 
#	IP_nb_pval = apply(IP_CTRL_tmp, MARGIN=1, function(x) get_pval(x[1], bin_num, x[2], AVEmat_cbg_prob, obs_0_num) )

	IP_nb_pval = pnbinom(IP_CTRL_tmp[,1]-1, IP_CTRL_tmp[,2], AVEmat_cbg_prob, lower.tail=FALSE) #/ pnbinom(0, IP_CTRL_tmp[,2], AVEmat_cbg_prob, lower.tail=FALSE) * (bin_num-obs_0_num)/bin_num


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
        fwrite(neglog10_nb_pval_bedgraph, output_file_tmp, quote=FALSE, col.names=FALSE, row.names=FALSE, sep='\t')
}


