args = commandArgs(trailingOnly=TRUE)
input_parafile = args[1]
scale_use = as.logical(args[2])
source_file = args[3]

source(source_file)

pdf(paste(input_parafile, '.pdf', sep=''))
createHeatmap(input_parafile, scale=scale_use)
dev.off()

png(paste(input_parafile, '.png', sep=''))
createHeatmap(input_parafile, scale=scale_use)
dev.off()
