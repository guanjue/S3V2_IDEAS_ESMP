
args = commandArgs(trailingOnly=TRUE)
filelist_file = args[1]

files = read.table(filelist_file, head=F)
cpk = scan(paste(toString(files[1,4]), '_commonpkfdr01_z.cpk.txt', sep=''))
cpk_b = cpk!=0

for (i in 1:dim(files)[1]){
f = toString(files[i,1])
print(f)
d=read.table(f, header=F, sep='\t')[,4]
print(summary(d[cpk_b]))
print(sum(d>10))
rm(d)
}
