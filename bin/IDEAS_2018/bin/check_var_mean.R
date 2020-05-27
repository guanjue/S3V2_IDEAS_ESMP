args = commandArgs(trailingOnly=TRUE)
state_file = args[1]
full_mk_ct = args[2]
mk_list = args[3]
signal_folder = args[4]
lowlim = as.numeric(args[5])
upperlim = as.numeric(args[6])
check_num = as.numeric(args[7])


#state_file = 'ideas_bp23.state'
#full_mk_ct = 'B_B15_50'
#mk_list = '../list_files/mark_list.txt'
#signal_folder = '../s3norm_0_100_sig/'
#upperlim = 100
#lowlim = 0
#check_num = 100000

### read ideas state
d = read.table(state_file, comment.char='~', header=T)

sample_id = sample(dim(d)[1], check_num)
d = d[sample_id,]
### get cell type
dct = d[,colnames(d)==full_mk_ct]

### get states
state_num = unique(dct)

### read mk
mk_list = read.table(mk_list, header=F)

### read signal
signal_mat = c()
for (mk in t(mk_list)){
	sig_file = paste(signal_folder, full_mk_ct, '.', mk, '.s3norm.', lowlim, '_', upperlim, '.txt', sep='')
	print(sig_file)
	sig_tmp = scan(sig_file)
	sig_tmp = sig_tmp[sample_id]
	signal_mat = cbind(signal_mat, sig_tmp)
}

### get state mean & variance
signal_mat_info = c()

for (i in state_num[1:28]){
	print(i)
	signal_mat_i = signal_mat[dct==i,]
	signal_mat_i_colvar = apply(signal_mat_i, 2, var)
	signal_mat_i_colmean = apply(signal_mat_i, 2, mean)
	signal_mat_i_info = sum(signal_mat_i_colvar/(signal_mat_i_colmean)^2)
	signal_mat_info = c(signal_mat_info, signal_mat_i_info)
}

output_mat = cbind(state_num, signal_mat_info)[order(state_num),]

write.table(output_mat, 'state_var_mean.txt', col.names=F, row.names=F, sep='\t', quote=F)

output_mat = read.table('state_var_mean.txt', header=F)

pdf('var_vs_meansq.pdf', width=14)
plot(output_mat[,1], output_mat[,2])
lines(output_mat[,1], output_mat[,2])
dev.off()


colvar = read.table('state_var.txt', header=F)
colmean = read.table('state_mean.txt', header=F)

info_1 = (colvar)^0.5/(colmean)

pdf('var_vs_meansq_mk.pdf', width=14)
plot(0:(dim(info_1)[1]-1), info_1[,1], ylim = c(min(info_1), max(info_1)))
for (i in 1:dim(info_1)[2]){
points(0:(dim(info_1)[1]-1), info_1[,i])
lines(0:(dim(info_1)[1]-1), info_1[,i])
}
dev.off()





