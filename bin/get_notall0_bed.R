args = commandArgs(trailingOnly=TRUE)

input_state_file = args[1] #'H3K36me3.file_list.txt'
output_nonall0_bed_file = args[2]

#input_state_file = 'S3V2_IDEAS_mm10_r3_withHg38Mm10prior.state'
#output_nonall0_bed_file = 'S3V2_IDEAS_mm10_r3_withHg38Mm10prior.nonall0.bed'

d = read.table(input_state_file)
dbed = d[,2:4]
ds = d[,5:(dim(d)[2]-1)]
dss = rowSums(ds)
dbed_nonall0 = dbed[dss!=0,]

write.table(dbed_nonall0, output_nonall0_bed_file, quote=F, sep='\t', col.names=F, row.names=F)

