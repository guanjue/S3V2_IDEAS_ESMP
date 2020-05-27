#markcol=rbind(c(1,0,0),c(1,1,0),c(0,1,1),c(0,0,1));
stateColor<-function(statemean, markcolor=NULL)
{	
	if(length(markcolor)==0)
	{	markcolor=rep("",dim(statemean)[2]);
		markcolor[order(apply(statemean,2,sd),decreasing=T)]=hsv((1:dim(statemean)[2]-1)/dim(statemean)[2],1,1)
		markcolor=t(col2rgb(markcolor));
	}

	rg=apply(statemean,1,range);
	mm=NULL;
	for(i in 1:dim(statemean)[1])
	{	mm=rbind(mm,(statemean[i,]-rg[1,i]+1e-10)/(rg[2,i]-rg[1,i]+1e-10));
	}
	mm = mm^5; 
	if(dim(mm)[2]>1) mm = mm / (apply(mm, 1, sum)+1e-10);
	mycol=mm%*%markcolor;
	s=apply(statemean,1,max);
	s=(s-min(s))/(max(s)-min(s)+1e-10);
#s=s^0.5;
	
mycol=round(255-(255-mycol)*s/0.5);
mycol[mycol<0]=0;
rt=paste(mycol[,1],mycol[,2],mycol[,3],sep=",");
h=t(apply(mycol,1,function(x){rgb2hsv(x[1],x[2],x[3])}));
h=apply(h,1,function(x){hsv(x[1],x[2],x[3])});
rt=cbind(rt,h);
return(rt);

	h=t(apply(mycol,1,function(x){rgb2hsv(x[1],x[2],x[3])}));
	h[,2]=h[,2]*s;
	#h[,3]=1-(1-h[,3])*s^0.5;
	h=apply(h,1,function(x){hsv(x[1],x[2],x[3])});
	rt=cbind(apply(t(col2rgb(h)),1,function(x){paste(x,collapse=",")}),h);
	
	return(rt);
}

createTrack<-function(statefiles, genomefile, outpref, statecolor, header, statename=NULL)
{	message("Reading state file: ", appendLF=FALSE);
	library("data.table");
	genomesz = read.table(genomefile);
	g=NULL;
	for(i in statefiles)
	{	message(paste(i," ",sep=""),appendLF=F);
		#tg=as.matrix(read.table(i, comment="!", header=(length(header)==0)));
		tg=as.matrix(fread(i));
		t=NULL;
		for(j in 1:dim(genomesz)[1])
		{	t=c(t,which(tg[,2]==as.character(genomesz[j,1]) & as.numeric(tg[,4])>as.numeric(genomesz[j,2])));
		}	
		if(length(t)>0) { tg = tg[-t,]; }
		t=which(is.na(match(tg[,2], genomesz[,1]))==T);
		if(length(t)>0) { tg = tg[-t,]; }	
print(c(dim(g),dim(tg)));
		g=rbind(g,tg);
	}
	message("Done");
	uchr = sort(unique(as.character(g[,2])));
	g1=NULL;
	for(i in uchr)
	{	t=which(as.character(g[,2])==i);
		g1=rbind(g1,g[t[order(as.numeric(g[t,3]))],]);
	}
	g=NULL;

	chr=as.character(g1[,2]);
	posst=as.numeric(g1[,3]);
	posed=as.numeric(g1[,4]);
	state=as.matrix(g1[,5:(dim(g1)[2]-1)]);
	if(length(statename)==0) statename=0:max(state);
	L=dim(g1)[1];
	n=dim(state)[2];
	if(length(header) > 0) colnames(g1) = header;
	cells=as.character(colnames(g1)[5:(dim(g1)[2]-1)]);
	g1=NULL;
	message("Generating tracks");
	options(scipen=999);

	tt = which(chr[2:L]!=chr[2:L-1]);
	tt = c(tt,which(posst[2:L]!=posed[2:L-1]));
	tt = sort(unique(tt));

	for(i in 1:n)
	{	tstate = state[,i];
		#print(c(i,L,length(tstate),length(chr),length(posst),length(posed)));

		t=c(tt,which(tstate[2:L]!=tstate[2:L-1]));
		t=sort(unique(t));
		t0=c(0,t)+1;
		t=c(t,L);
		np=cbind(chr[t],posst[t0],posed[t],tstate[t]);

		#print("make track");
		x = cbind(np[,1:3],statename[as.integer(np[,4])+1],1000,".",np[,2:3],statecolor[as.numeric(np[,4])+1]);
		write.table(as.matrix(x),paste(outpref,i,"bed1",sep="."),quote=F,row.names=F,col.names=F);
print(x[1,]);
		#x = apply(x,1,function(x){paste(x,collapse="\t")});
		#write.table(x,paste(outpref,i,"bed",sep="."),quote=F,row.names=F,col.names=F);
		#system(paste("sort-bed ", outpref, ".", i, ".bed > ", outpref, ".", i, ".bed1", sep=""));
		system(paste("bedToBigBed ", outpref, ".", i, ".bed1 ", genomefile, " ", outpref, ".",i,".bb",sep=""));
		#system(paste("rm ", paste(outpref, i,"bed",sep=".")));
		system(paste("rm ", paste(outpref, i,"bed1",sep=".")));
	}
	return(cells);
}

#cellinfo: shortid match with those in states, long id, cell description, text color
run<-function(statefiles, hubid, genomeid, genomefile, statecolor, targetURL="http://bx.psu.edu/~yuzhang/tmp/", trackfolder=NULL, hubname = NULL, cellinfo = NULL, header=NULL,statename=NULL)
{
#	if(file.exists("bedToBigBed") == FALSE) { message("Cannot find bedToBigBed in the folder where you run this function."); return(0); }
#	if(file.exists("sort-bed") == FALSE) { message("Cannot find sort-bed in the folder where you run this function."); return(0); }
	if(length(hubname) == 0) hubname = hubid;
	if(length(trackfolder) == 0) trackfolder = paste("tracks_", hubid, "/", sep="");
	if(substring(trackfolder, nchar(trackfolder)) != "/") trackfolder = paste(trackfolder, "/", sep="");
	dir.create(trackfolder,showWarnings=FALSE);

	cells=createTrack(statefiles, genomefile, paste(trackfolder, hubid, sep=""), statecolor, header=header,statename=statename);
	if(length(cellinfo) == 0)
	{	cellinfo = cbind(cells, cells, cells, "#000000");
		cellinfo = array(cellinfo, dim=c(length(cells),4));
	}
	cellinfo = as.matrix(cellinfo);

	trackDb=NULL;
	for(i in 1:length(cells))
	{	ii=which(cells[i] == cellinfo[,1]);
		if(length(ii)==0) next;
		ii=ii[1];
		trackDb=c(trackDb, paste("track bigBed", i, sep=""));
		trackDb=c(trackDb, paste("priority",ii));
		trackDb=c(trackDb, "type bigBed 9 .");
		trackDb=c(trackDb, "itemRgb on");
		trackDb=c(trackDb, "maxItems 100000");
		trackDb=c(trackDb, paste("bigDataUrl ", targetURL, hubid,".",i,".bb",sep=""));
		trackDb=c(trackDb, paste("shortLabel", cellinfo[ii,2]));
		trackDb=c(trackDb, paste("longLabel", paste(hubname, cellinfo[ii,3])));
		trackDb=c(trackDb, paste("color",paste(c(col2rgb(cellinfo[ii,4])),collapse=",")));
		trackDb=c(trackDb, "visibility dense");
		trackDb=c(trackDb, "");	
	}

	write.table(trackDb, paste(trackfolder, "trackDb_", hubid, ".txt",sep=""),quote=F,row.names=F,col.names=F);

	write.table(c(paste("genome", genomeid), paste("trackDb trackDb_", hubid, ".txt", sep="")), paste(trackfolder, "genomes_", hubid, ".txt",sep=""), quote=F,row.names=F,col.names=F);

	write.table(c(paste("hub", hubid), paste("shortLabel", hubid), paste("longLabel", hubname), paste("genomesFile genomes_", hubid, ".txt", sep=""), "email yzz2 at psu.edu"), paste(trackfolder, "hub_", hubid, ".txt",sep=""), quote=F,row.names=F,col.names=F);

	return(1);
}

createHeatmap<-function(parafile, statecolor = NULL, markcolor = NULL, cols=c("white","dark blue"), show=TRUE,fout=NULL, sortstate=TRUE,scale=FALSE)
{	x=read.table(parafile,comment="!",header=T);
	k=dim(x)[2];
	l=dim(x)[1];
	p=(sqrt(9+8*(k-1))-3)/2;
	x[,2:(p+1)]=x[,2:(p+1)][,order(colnames(x)[2:(p+1)])]
	colnames(x)[2:(p+1)] = colnames(x)[2:(p+1)][order(colnames(x)[2:(p+1)])]
	m=as.matrix(x[,1+1:p]/x[,1]);
	colnames(m) = colnames(x)[1+1:p];
	m = m[,order(colnames(m))];
	marks=colnames(m);
	rownames(m)=paste(1:l-1," (",round(x[,1]/sum(x[,1])*10000)/100,"%)",sep="");

if(sortstate)
{
o=hclust(dist(m),method="ward.D2")$order;
m=m[o,];
if(length(statecolor) != 0)
{	statecolor=statecolor[o,];
}
}
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
	return(statecolor);
}
