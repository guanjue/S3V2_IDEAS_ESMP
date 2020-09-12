### get parameters
args = commandArgs(trailingOnly=TRUE)

input = args[1]
output = args[2]
win = as.numeric(args[3])

### read info
d = scan(input)
bin_num = length(d)
start_n = 1+win
end_n = bin_num-win
### get body
sig_mat = d[start_n:end_n]
for (i in 1:win){
sig_vec_tmp1 = d[(start_n-i):(end_n-i)]
sig_vec_tmp2 = d[(start_n+i):(end_n+i)]
sig_mat = cbind(sig_mat, sig_vec_tmp1, sig_vec_tmp2)
}
sig_smooth = rowMeans(sig_mat)

### get head & tail smooth
sig_smooth_head = c()
sig_smooth_tail = c()
for (i in 1:win){
sig_head_tmp = mean(d[1:(win+i)])
sig_smooth_head = c(sig_smooth_head, sig_head_tmp)
sig_tail_tmp = mean(d[(bin_num-win-i+1):bin_num])
sig_smooth_tail = c(sig_tail_tmp, sig_smooth_tail)
}

### merge head body tail
sig_smooth_all = c(sig_smooth_head, sig_smooth, sig_smooth_tail)

write.table(sig_smooth_all, output, quote=F, col.names=F, row.names=F, sep='\t')