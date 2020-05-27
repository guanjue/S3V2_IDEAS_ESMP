library(data.table)

args = commandArgs(trailingOnly=TRUE)

file_list_file = args[1]
output = args[2]

### read input
file_list = read.table(file_list_file, header=F)

common_pk = c()

print('read the first one')
d10 = as.data.frame(fread(toString(file_list[1,1])))
bed = d10[,1:3]
sum_sig = d10[,4]
non0mean = mean(sum_sig[sum_sig>0])
rm(d10)

if (non0mean<0.01){
	print('Is the data average counts data?')
	sum_sig = sum_sig/non0mean
}

print('read other files')
if (dim(file_list)[1]>1){
for (i in 2:dim(file_list)[1]){
	print(i)
        print(file_list[i,1])
        d10 = as.data.frame(fread(toString(file_list[i,1])))[,4]
	print(summary(d10[d10>0]))
	non0mean = mean(d10[d10>0])
	if (non0mean<0.01){
		d10 = d10/non0mean
	}
        sum_sig = sum_sig + d10
	rm(d10)
}
}

average_sig = sum_sig / dim(file_list)[1]

write.table(cbind(bed, average_sig), output, col.names=F, row.names=F, quote=F, sep='\t')

