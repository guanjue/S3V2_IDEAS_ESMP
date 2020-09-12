
source('../bin/createGenomeTracks.R')
chr_list = c(1:22, 'X', 'Y')

para_mat = read.table(paste('ideas_bp19_0_16_rerun.chr', chr_list[1], '.para', sep=''))
for (chr in chr_list[-1]){
	print(chr)
	para_mat = para_mat+read.table(paste('ideas_bp19_0_16_rerun.chr', chr_list[1], '.para', sep=''))
}

para_mat = para_mat/length(chr_list)
write.table(para_mat, 'ideas_bp19_0_16_rerun.allchr.para', quote=F, sep=' ', col.names=F, row.names=F)


mk_list = c('H3K4me3', 'H3K4me1', 'H3K27me3', 'H3K27ac', 'H3K9me3', 'H3K36me3', 'ATAC')

for (mk in mk_list){
print(mk)
d=scan(paste(mk,'.fisherp.s3norm.ref.txt',sep=''))
print(summary(d[d>0]))
}

for (chr in chr_list){
	print(chr)
	pdf(paste('ideas_bp19_0_16_rerun_modified.chr', chr, '.para.pdf', sep=''))
	createHeatmap(paste('ideas_bp19_0_16_rerun_modified.chr', chr, '.para', sep=''))
	dev.off()
}


cp ideas_bp19_0_16_rerun.chr1.state ideas_bp19_0_16_rerun.allchr.state
for i in {2..22}
do
	echo $i
	tail -n+2 'ideas_bp19_0_16_rerun.chr'$i'.state' >> ideas_bp19_0_16_rerun.allchr.state
done

tail -n+2 ideas_bp19_0_16_rerun.chrX.state >> ideas_bp19_0_16_rerun.allchr.state
tail -n+2 ideas_bp19_0_16_rerun.chrY.state >> ideas_bp19_0_16_rerun.allchr.state

time Rscript ../bin/get_modified_para_heatmap.R ideas_bp19_0_16_rerun.allchr.para ideas_bp19_0_16_rerun.allchr.state F ../bin/createGenomeTracks.R




plot_heat = function(name, para){
pdf(name)
source('../bin/createGenomeTracks.R')
createHeatmap(para)
dev.off()
}

plot_heat('ideas_bp19_1_16.para0.pdf', 'ideas_bp19_2_16.para0')

