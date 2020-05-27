args = commandArgs(trailingOnly=TRUE)

ctlist_file = args[1]
input_mat = args[2]
output = args[3]

### get ct list
ct_list = read.table(ctlist_file, header=T)

### get ct vector
ce_vec = c()
for (i in 1:dim(ct_list)[1]){
ce_vec[i] = unlist(strsplit(toString(ct_list[i,1]), "_"))[1]
}

### get merged ct
ct_num = table(ce_vec)

if (sum(ct_num>1)>0){
ct_merge = rownames(ct_num)[ct_num>1]

### read peak od
d = read.table(input_mat, header=F)
dnew = c()
rm_coln = c()

check_rep = function(x){
x_binary = x>0
pk01 = (sum(x_binary)>=2)*1
return(pk01)
}

for (i in 1:length(ct_merge)){
  used_coln = which(ce_vec==ct_merge[i])+3
  newcol = apply(as.matrix(d[,used_coln]), 1, function(x) check_rep(x))
  dnew = cbind(dnew, newcol)
  rm_coln = c(rm_coln, used_coln)
}

drm = d[,-rm_coln]
drepmerge = cbind(drm, dnew)


dc = drepmerge[,-c(1:3)]
dcb = dc>0
dcs = rowSums(dcb)

dbedr = d[dcs>=1,c(1:3)]
write.table(dbedr, output, quote=F, col.names=F, row.names=F, sep='\t')
} else{
d = read.table(input_mat, header=F)
write.table(d[,1:3], output, quote=F, col.names=F, row.names=F, sep='\t')
}

