### get parameters
args = commandArgs(trailingOnly=TRUE)

parafile1 = args[1]
parafile2 = args[2]
outputfile = args[3]
source_createGenomeTracks = args[4]

print(parafile1)
print(parafile2)
print(outputfile)
print(source_createGenomeTracks)
#source('~/group/projects/vision/createGenomeTracks.R')
source(source_createGenomeTracks)


#createHeatmap<-function(parafile, statecolor = NULL, markcolor = NULL, cols=c("white","dark blue"), show=TRUE,fout=NULL, sortstate=TRUE,scale=FALSE)
#{
#pdf('vision_with_without_Meth.pdf', height=14, width=7)
pdf(outputfile, height=14, width=7)
#parafile1='~/group/projects/vision/withMETH/check_order/run_IDEAS.para'
#parafile2='~/group/projects/vision/cell7_noMETH/run_IDEAS.para'
statecolor = NULL
markcolor = NULL
show=TRUE
fout=NULL
sortstate=TRUE
cols=c("white","dark blue")
scale=F
	x1=read.table(parafile1,comment="!",header=T);
	x2=read.table(parafile2,comment="!",header=T);
	x = rbind(x1, x2)
	k=dim(x)[2];
	l=dim(x)[1];
	l1=dim(x1)[1];
	l2=dim(x2)[1];

	p=(sqrt(9+8*(k-1))-3)/2;
	m=as.matrix(x[,1+1:p]/x[,1]);
	colnames(m) = colnames(x)[1+1:p];
	m0 = m[,-6]
	marks=colnames(m);
	rn1 = paste('1r:',1:l1-1," (",round(x1[,1]/sum(x1[,1])*10000)/100,"%)",sep="");
	rn2 = paste('2r:',1:l2-1," (",round(x2[,1]/sum(x2[,1])*10000)/100,"%)",sep="");
	#rownames(m)=paste(1:l-1," (",round(x[,1]/sum(x[,1])*10000)/100,"%)",sep="");
	rownames(m)=c(rn1, rn2)
if(sortstate)
{
o=hclust(dist(m0),method="ward.D2")$order;
pdf(paste(outputfile, '.tree.pdf', sep=''), width=14)
plot(hclust(dist(m0),method="ward.D2"), cex = 0.6, hang = -1)
dev.off()

m=m[o,];
if(length(statecolor) != 0)
{	statecolor=statecolor[o,];
}
}
type = c(rep(1, length(rn1)), rep(2, length(rn2)))
orn1 = o[type==1]
orn2 = o[type==2]
om=m;
if(scale)
{	m = t((t(m) - apply(m,2,min))/(apply(m,2,max)-apply(m,2,min)+1e-10));
}
	if(length(fout)!=0)
	{	pdf(fout);	}	
	par(mar=c(6,1,1,6));
	rg=range(m);
	colors=0:100/100*(rg[2]-rg[1])+rg[1];
        my_palette=colorRampPalette(cols)(n=100);
	defpalette=palette(my_palette);
if(show)
{
	plot(NA,NA,xlim=c(0,p+0.7),ylim=c(0,l),xaxt="n",yaxt="n",xlab=NA,ylab=NA,frame.plot=F);
	axis(1,at=1:p-0.5,labels=colnames(m),las=2);
	axis(4,at=1:l-0.5,labels=rownames(m),las=2);
	rect(rep(1:p-1,l),rep(1:l-1,each=p),rep(1:p,l),rep(1:l,each=p),col=round((t(m)-rg[1])/(rg[2]-rg[1])*100));#,border=NA);
}
#if(scale)
#{	m = om;
#}
	if(length(statecolor)==0)
	{	if(length(markcolor)==0)
		{	markcolor=t(col2rgb(terrain.colors(ceiling(p))[1:p]));	
			for(i in 1:length(marks))
			{	if(regexpr("h3k4me3",tolower(marks[i]))>0)
				{	markcolor[i,]=c(255,0,0);	}
				if(regexpr("h3k4me2",tolower(marks[i]))>0)
				{	markcolor[i,]=c(250,100,0);	}
				if(regexpr("h3k4me1",tolower(marks[i]))>0)
				{	markcolor[i,]=c(250,250,0);	}
				if(regexpr("h3k36me3",tolower(marks[i]))>0)
				{	markcolor[i,]=c(0,150,0);	}
				if(regexpr("h2a",tolower(marks[i]))>0)
				{	markcolor[i,]=c(0,150,150);	}
				if(regexpr("dnase",tolower(marks[i]))>0)
				{	markcolor[i,]=c(0,200,200);	}
				if(regexpr("atac",tolower(marks[i]))>0)
				{	markcolor[i,]=c(200,50,150);	}
				if(regexpr("dnase",tolower(marks[i]))>0)
				{	markcolor[i,]=c(200,50,150);	}
				if(regexpr("h3k9ac",tolower(marks[i]))>0)
				{	markcolor[i,]=c(250,0,200);	}
				if(regexpr("h3k9me3",tolower(marks[i]))>0)
				{	markcolor[i,]=c(100,100,100);	}
				if(regexpr("h3k27ac",tolower(marks[i]))>0)
				{	markcolor[i,]=c(250,150,0);	}
				if(regexpr("h3k27me3",tolower(marks[i]))>0)
				{	markcolor[i,]=c(0,0,225);	}
				if(regexpr("h3k79me2",tolower(marks[i]))>0)
				{	markcolor[i,]=c(200,0,200);	}
				if(regexpr("h4k20me1",tolower(marks[i]))>0)
				{	markcolor[i,]=c(50,200,50);	}
				if(regexpr("ctcf",tolower(marks[i]))>0)
				{	markcolor[i,]=c(200,0,250);	}
				if(regexpr("wgbs",tolower(marks[i]))>0)
				{	markcolor[i,]=c(30,144,255);	}
			}
		}
		statecolor=array(stateColor(m,markcolor),dim=c(dim(m)[1],2));
	}
	if(show)
	{	rect(rep(p+0.2,l),1:l-0.8,rep(p+0.8,l),1:l-0.2,col=statecolor[,2]);
	}
	if(sortstate)	statecolor[o,]=statecolor;
	palette(defpalette);
	if(length(fout)!=0)
	{	dev.off();	}

dev.off()

#	return(statecolor);
#}

