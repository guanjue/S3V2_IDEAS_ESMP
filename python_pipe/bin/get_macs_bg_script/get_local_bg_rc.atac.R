args = commandArgs(trailingOnly=TRUE)
input_200bp = args[1]
input_rc_5_10kb = args[2]
output = args[3]

### get common mean
input_200bp_sig = scan(input_200bp)
wg_mean_input_sig = mean(input_200bp_sig)

input_rc_5kb_10kb_sig = read.table(input_rc_5_10kb, header=F)
wg_mean_input_rc_sig = rep(wg_mean_input_sig, dim(input_rc_5kb_10kb_sig)[1])

input_sig_mat = cbind(wg_mean_input_rc_sig, input_rc_5kb_10kb_sig[,4:5])

input_sig_mat_max = apply(input_sig_mat, 1, max)
input_sig_mat_max[is.na(input_sig_mat_max)] = wg_mean_input_sig
print(dim(input_sig_mat_max))

write.table(cbind(input_rc_5kb_10kb_sig[,1:3], input_sig_mat_max), output, quote=FALSE, col.names=FALSE, row.names=FALSE, sep='\t')

