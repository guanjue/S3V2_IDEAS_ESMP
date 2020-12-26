args = commandArgs(trailingOnly=TRUE)
input_neglog10 = args[1]
output = args[2]

d = read.table(input_neglog10, header=F, sep='\t')
ds = d[,4]
ds[ds>323] = 323
dp = 10^(-ds)
dp[dp<=1e-323] = 1e-323
dp[dp>1] = 1.0
dpfdr = -log10(p.adjust(dp, 'fdr'))

fwrite(cbind(d[,1:3],dpfdr), output, quote=F, sep='\t', row.names=F, col.names=F)
rm(d)

