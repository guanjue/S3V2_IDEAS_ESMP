args = commandArgs(trailingOnly=TRUE)
input_neglog10 = args[1]
output = args[2]

library(data.table)

d = as.data.frame(fread(input_neglog10))
dp = 10^(-d[,4])
dpfdr = -log10(p.adjust(dp, 'fdr'))

fwrite(cbind(d[,1:3],dpfdr), output, quote=F, sep='\t', row.names=F, col.names=F)
rm(d)

