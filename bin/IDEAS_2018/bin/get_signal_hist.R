file_list = read.table('signal_list.txt', header=F)

d_mat = c()
for (file in file_list){
	print(file)
	d = scan(file)
	d_mat = cbind(d_mat, d)
}

d_mat_m = d_mat
d_mat_m[d_mat>16] = 16

pdf('../all.sig.hist.pdf', height=15)
par(mfrow=c(2,1))
hist.data = hist(d_mat_m, plot=F)
hist.data$counts = (hist.data$counts)
plot(hist.data, ylab='(Frequency)')
abline(v=2)
box()
hist.data = hist(d_mat_m[d_mat_m>0], plot=F)
hist.data$counts = (hist.data$counts)
plot(hist.data, ylab='(Frequency)')
abline(v=2)
box()
dev.off()

