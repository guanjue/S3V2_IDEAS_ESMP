library(pheatmap)

#IDEAS_folder = '/storage/home/gzx103/scratch/S3V2norm_compare/human_vision/'
#IDEAS_output_name = 'hg38bp0402'
#output_name='.var_mean.pdf'

args = commandArgs(trailingOnly=TRUE)
IDEAS_folder = args[1]
IDEAS_output_name = args[2]
output_name = args[3]
IDEAS_output_folder_tail = args[4]


print(IDEAS_output_folder_tail)

d = read.table(paste(IDEAS_folder, IDEAS_output_name, '.input', sep=''), header=F, sep=' ')

uniq_mk = unique(d[,2])
uniq_mk_num = length(uniq_mk)

ct_num = table(d[,1])
print(ct_num)
fullset_ct = rownames(ct_num)[ct_num==uniq_mk_num]

print(paste(IDEAS_folder, IDEAS_output_name, IDEAS_output_folder_tail, IDEAS_output_name, '.state', sep=''))
library(data.table)
state = as.data.frame(fread(paste(IDEAS_folder, IDEAS_output_name, IDEAS_output_folder_tail, IDEAS_output_name, '.state', sep='')))
state_tmp = state[,5]
state_num = length(unique(state_tmp))

state_var_all = matrix(0, nrow=state_num, ncol=uniq_mk_num)
state_mean_all = matrix(0, nrow=state_num, ncol=uniq_mk_num)
state_var_mean_all = matrix(0, nrow=state_num, ncol=uniq_mk_num)

for (i in 1:length(fullset_ct)){
#for (i in 1:2){

#
### get file list
print(fullset_ct[i])
files_tmp = d[d[,1]==fullset_ct[i],3]
mk_tmp = d[d[,1]==fullset_ct[i],2]
### get signal mat
sigmat_tmp = c()
for (mk_i in uniq_mk){
print(mk_i)
sigmat_tmp_mk = scan(toString(files_tmp[mk_tmp == mk_i]))
sigmat_tmp = cbind(sigmat_tmp, sigmat_tmp_mk)
}
### get ct state
state_tmp = state[,colnames(state) == fullset_ct[i]]
print(head(state_tmp))
state_var = c()
state_mean = c()
for (j in 0:(state_num-1)){
print(j)
sigmat_tmp_j = sigmat_tmp[state_tmp==j,]
#print(dim(sigmat_tmp_j))
sigmat_tmp_j_colmean = colMeans(sigmat_tmp_j)
#print(sigmat_tmp_j_colmean)
sigmat_tmp_j_colvar = apply(sigmat_tmp_j, 2, var)
#sigmat_tmp_j_colvar = apply(sigmat_tmp_j, 2, sd)
#print(sigmat_tmp_j_colvar)
#print(sigmat_tmp_j_colvar)
#print(sigmat_tmp_j_colmean)
state_var = rbind(state_var, (sigmat_tmp_j_colvar+1))
state_mean = rbind(state_mean, (sigmat_tmp_j_colmean+1))
}
colnames(state_var) = uniq_mk
rownames(state_var) = 0:(state_num-1)
colnames(state_mean) = uniq_mk
rownames(state_mean) = 0:(state_num-1)
state_var_mean_i = state_var/state_mean
print(state_mean)
print(state_var)
print(state_var/state_mean)
state_var_all = state_var_all + state_var
state_mean_all = state_mean_all + state_mean
state_var_mean_all = state_var_mean_all+state_var_mean_i
#pdf(paste(fullset_ct[i], output_name, sep=''))
#my_colorbar=colorRampPalette(c('white', 'blue'))(n = 128)
#col_breaks = c(seq(0, 2000,length=33))
#pheatmap(state_var_mean, color=my_colorbar, cluster_cols = FALSE,cluster_rows=FALSE,show_rownames=TRUE,show_colnames=TRUE)
#dev.off()
}

used_order_col = order(colnames(state_var_all))
state_var_all = state_var_all[,used_order_col]
state_mean_all = state_mean_all[,used_order_col]

state_var_mean_all = (state_var_all)/(state_mean_all)

state_var_mean_all = state_var_mean_all/length(fullset_ct)

write.table(state_var_mean_all, paste(IDEAS_folder, IDEAS_output_name, IDEAS_output_folder_tail, 'all', output_name, '.txt', sep=''), quote=F)
write.table(state_var_all, paste(IDEAS_folder, IDEAS_output_name, IDEAS_output_folder_tail, 'all_var', output_name, '.txt', sep=''), quote=F)
write.table(state_mean_all, paste(IDEAS_folder, IDEAS_output_name, IDEAS_output_folder_tail, 'all_mean', output_name, '.txt', sep=''), quote=F)

#state_var_mean_all = cbind(rowMeans(state_var_mean_all), rowMeans(state_var_mean_all))
pdf(paste(IDEAS_folder, IDEAS_output_name, IDEAS_output_folder_tail, 'all', output_name, sep=''))
my_colorbar=colorRampPalette(c('white', 'blue'))(n = 128)
col_breaks = c(seq(0, 2000,length=33))
pheatmap(state_var_mean_all, color=my_colorbar, cluster_cols = FALSE,cluster_rows=FALSE,show_rownames=TRUE,show_colnames=TRUE)
dev.off()



