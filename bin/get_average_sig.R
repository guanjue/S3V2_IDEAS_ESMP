args = commandArgs(trailingOnly=TRUE)

file_list_file = args[1]
output = args[2]

### read input
file_list = read.table(file_list_file, header=F)

set.seed(2020)
if (dim(file_list)[1]>50){
used_id_sample = sample(dim(file_list)[1], 50)
file_list = file_list[used_id_sample,]
}

common_pk = c()

print('read the first one')
d10 = read.table(toString(file_list[1,1]), header=F, sep='\t')
bed = d10[,1:3]
sum_sig = d10[,4]
non0mean = mean(sum_sig[sum_sig>0])
rm(d10)

if (non0mean<0.01){
	print('Is the data average counts data?')
	sum_sig = sum_sig/non0mean
}

print('read other files')
notused_n = 0
if (dim(file_list)[1]>1){
for (i in 2:dim(file_list)[1]){
	print(i)
        print(file_list[i,1])
	if (file.exists(toString(file_list[i,1]))){
        	d10 = read.table(toString(file_list[i,1]), header=F, sep='\t')[,4]
		if (is.na(mean(d10)) || (max(d10)==0)){
        		print('!!!Something wrong with the bigWig to signal step!!!')
        		notused_n = notused_n+1
        	next}
		print(summary(d10[d10>0]))
		non0mean = mean(d10[d10>0])
		if (non0mean<0.01){
			d10 = d10/non0mean
		}
        	sum_sig = sum_sig + d10
		rm(d10)
	}
}
}

average_sig = sum_sig / (dim(file_list)[1]-notused_n)

write.table(cbind(bed, average_sig), output, col.names=F, row.names=F, quote=F, sep='\t')

