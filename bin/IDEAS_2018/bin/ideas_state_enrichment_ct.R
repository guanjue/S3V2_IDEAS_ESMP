####################################################
### get parameters
args = commandArgs(trailingOnly=TRUE)
state_file = args[1]
output_name = args[2]
add_num = as.numeric(args[3])
cluster_state = as.logical(args[4])
cluster_ct = as.logical(args[5])


library(pheatmap)

#state_file = 'test.state'
#output_name = 'test.state.enrich.pdf'
#add_num = 100

#d =read.table('pknorm_2_16lim_ref1mo_0424_lesshet.state', header=T, sep=' ', comment.char = "~")
d =read.table(state_file, header=T, sep=' ', comment.char = "~")

### get state
d_state = d[,-c(1:4, dim(d)[2])]
allbins = dim(d_state)[1]

### get uniq state
d_state_uniq = unique(d_state[,1])
d_state_uniq = d_state_uniq[order(d_state_uniq)]

### get cell type
ct = colnames(d_state)
ct_num = length(ct)

s_enrich_all = c()
for (s in d_state_uniq){
	print(s)
	s_enrich_s = c()
	s_all = sum(d_state==s)
	for (ci in 1:ct_num){
		s_ct_tmp = sum(d_state[,ci] == s)
		s_ct_tmp_enrich = (s_ct_tmp+add_num) / (s_all/ct_num+add_num)
		s_enrich_s = c(s_enrich_s, s_ct_tmp_enrich)
	}
	s_enrich_all = cbind(s_enrich_all, s_enrich_s)
}

rownames(s_enrich_all) = ct
colnames(s_enrich_all) = d_state_uniq



my_colorbar=colorRampPalette(c('white', 'blue'))(n = 50)

pdf(output_name)
pheatmap(s_enrich_all,cluster_rows=cluster_ct,cluster_cols=cluster_state,color=my_colorbar,show_rownames=TRUE,show_colnames=TRUE,annotation_names_row=FALSE,annotation_names_col=TRUE)
dev.off()


