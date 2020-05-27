args = commandArgs(trailingOnly=TRUE)
input_parafile = args[1]
input_statefile = args[2]
scale_use = as.logical(args[3])
source_file = args[4]

### get bin number & cell number
d = read.table(input_statefile)
bin_num = dim(d)[1]
cell_num = dim(d)[2]-5

### get previous para file
d_para = read.table(input_parafile, comment.char = "", header=T)
state_num = dim(d_para)[1]

### get state bin number
n_vec = c()
for (i in 0:(state_num-1)){
	n = sum(d[,5:(cell_num+5-1)]==i)
	n_vec = c(n_vec, n)
	print(i)
	print(n)
	print(n/(bin_num*cell_num))
}

### get rescaled para file
d_para_m = c()
for (i in 1:dim(d_para)[1]){
d_para_m = rbind(d_para_m, d_para[i,]/d_para[i,1]*n_vec[i])
}

### write modified para file
output_modified_parafile = paste(input_parafile, '.modified.para', sep='')
write.table(d_para_m, output_modified_parafile, sep=' ', col.names=T, row.names=F, quote=F)

### get heatmap
source(source_file)

pdf(paste(output_modified_parafile, '.pdf', sep=''))
createHeatmap(output_modified_parafile, scale=scale_use)
dev.off()

png(paste(output_modified_parafile, '.png', sep=''))
createHeatmap(output_modified_parafile, scale=scale_use)
dev.off()



#time Rscript get_modified_para_heatmap.R run_IDEAS.para run_IDEAS.state F /storage/home/g/gzx103/group/software/IDEAS/IDEAS_2018/bin/createGenomeTracks.R