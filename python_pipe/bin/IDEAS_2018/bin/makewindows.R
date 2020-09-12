makewindows<-function(genome, blacklist, foutpref, wsz=200)
{
	fout = paste(foutpref, ".bed", sep="");
	str = paste("bedtools makewindows -g", genome, "-w", wsz, "-i winnum >", fout);
	system(str);
	if(length(blacklist)>0)
	{	fout1 = paste(fout,"tmp",sep="");
		for(i in blacklist)
		{	str = paste("bedtools subtract -A -a", fout,"-b", i, ">", fout1);
			system(str);
			system(paste("mv", fout1, fout));
		}
	}
	library("data.table");
	x=fread(fout);
	chr=unique(as.matrix(x[,1]));
	pos=as.matrix(x[,2:3]);
	y=NULL;
	for(i in chr)
	{	t=which(x[,1]==i);
		y=rbind(y,c(i,min(t)-1,max(t)));
	}

	finv=paste(foutpref,".inv",sep="");
	write.table(y,finv,quote=F,row.names=F,col.names=F);	
}
