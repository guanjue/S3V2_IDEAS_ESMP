#draw cell type clusters and epigenetic states
showPop<-function(fname, suffix="", indlist = NULL, inv = NULL, repn = 1)
{	
	border="light grey";

	g = read.table(paste(fname, suffix, ".g", sep="")); g=as.matrix(g[,-(1:3)]);
	N = (dim(g)[2] - 1) / repn;
	L = dim(g)[1];
	if(length(indlist) > 0)
	{	N = max(indlist); }
	else { indlist = 1:N; }
	if(length(inv) > 0)
	{	inv=inv[inv<=L & inv > 0];
		L = max(inv);	
	}
	else { inv = 1:L; }
	
	x = readPop(paste(fname, suffix, sep=""), indlist, inv);
	left = x$l; right = x$r; ss = x$s; index = x$i; len = x$z;

	oo = 1:length(indlist);

r=repn;
a=NULL;for(i in 1:length(indlist)) a=c(a,(indlist[oo[i]]-1)*r+1:r);
g=g[inv,c(a,dim(g)[2])];

	o = 1:length(index);
	allele = 0.5;

	c=unique(as.numeric(ss));
	t=max(ss)+1;
	m=rep(-1,t);
	for(i in 1:length(c)) { m[c[i]+1]=i; }
	par(mfrow=c(1,1));
	par(mar=c(2,2,2,2));
	palette(rainbow(max(6, length(c))));
	#if(allele==TRUE) palette(rainbow(max(6, max(length(c),max(g)+1))));
	Osz = max(allele, (length(o)+allele*dim(g)[2]) / 40);
	offset=max(Osz/2,allele/2);
	plot(-10,-10000,xlim=c(min(inv),max(inv)+1), ylim=c(0.5-allele*(dim(g)[2]-1)-Osz-offset*2,length(o)+0.3));
	ll = bb = rr = tt = cc = 0;
	for(i in 1:length(o))
	{	t = index[o[i]]:(index[o[i]]+len[o[i]]-1);
		ll=c(ll,left[t]);
		bb=c(bb,rep(i-0.5,length(t)));
		rr=c(rr,right[t]);
		tt=c(tt,rep(i+0.51,length(t)));
		cc=c(cc,m[ss[t]+1]);
	}
	rect(ll,bb,rr,tt,border=border,col=cc);

	if(allele>0)
	{	
colmap=rep(-1,max(g)+1);
j=1;
for(i in unique(c(as.matrix(g)))) {colmap[i+1]=j;j=j+1;}
#palette(rainbow(j));
nn=dim(g)[2];
xxx=rep(inv, nn);
yyy=ccc=ttt=NULL;
for(i in 1:(nn-1)) { yyy=c(yyy,rep(0.5-allele*(nn)-offset+i*allele,length(inv))); ccc=c(ccc,colmap[g[,i]+1]); }
yyy=c(yyy,rep(0.5-allele*(nn)-offset*2+allele-Osz,length(inv)));ccc=c(ccc,g[,nn]+1);
ttt=yyy+allele*0.95;
ttt[length(ttt)-1:length(inv)+1]=ttt[length(ttt)-1:length(inv)+1]+Osz-allele*0.95;
rect(xxx,yyy,xxx+0.95,ttt,col=ccc,border=FALSE);#pch=20,cex=0.5);

	}
}

readPara<-function(fprefix, xsz, ysz, suffix="")
{	x=read.table(paste(fprefix,suffix,".para",sep=""));
	n=dim(x)[1];
	m=array(0,dim=c(n*ysz, xsz));
	v=array(0,dim=c(ysz*n,ysz));

	yvsz = ysz * (ysz + 1) / 2;
	xvsz = xsz * (xsz - 1) / 2;

	for(i in 1:n)
	{	kkk=(i-1)*ysz;
		yty = diag(1, ysz);
		l=1;
		for(j in 1:ysz)
		{	for(k in 1:j)
			{	yty[j,k] = yty[j,k] + x[i,1+ysz+l];
				if(k!=j) { yty[k,j] = yty[k,j] + x[i,1+ysz+l]; }
				l=l+1;
			}
		}
		xtx = diag(1, xsz);
		xtx[1,1] = xtx[1,1] + x[i,1];
		if(xsz > 1)
		{	for(j in 1:(xsz-1)) 
			{	xtx[1,1+j] = xtx[1,1+j] + x[i,1+ysz+yvsz+j];
				xtx[1+j,1] = xtx[1+j,1] + x[i,1+ysz+yvsz+j];
			}
			tl=1;
			for(j in 2:xsz)
			{	for(k in 2:j)
				{	xtx[j,k] = xtx[j,k] + x[i,1+ysz+yvsz+xsz-1+tl];
					if(j!=k) xtx[k,j] = xtx[k,j] + x[i,1+ysz+yvsz+xsz-1+tl];
					tl=tl+1;
				}
			}
		}
		xty = array(0,dim=c(xsz,ysz));
		xty[1,] = as.matrix(x[i,1+1:ysz]);
		if(xsz > 1)
		{	l = 1;
			for(j in 1:(xsz-1))
			{	for(k in 1:ysz)
				{	xty[j+1,k] = x[i,1+ysz+yvsz+xsz-1+xvsz+l];
					l = l + 1;
				}
			}
		}
		m[kkk+1:ysz,] = t(solve(xtx)%*%xty);
		yty = (yty - t(xty)%*%solve(xtx)%*%xty) / (x[i,1]+1);
		v[kkk+1:ysz,] = yty;
	}
	para=NULL;
	para$p=x[,1]/sum(x[,1]);
	para$m=m;#beta matrix ysz by xsz
	para$v=v;

	return(para);
}

readPop<-function(fname, indlist, inv)
{
	N = max(indlist);
	L = max(inv) + 1;
	left = right = ss = index = len = NULL;
	str=paste(fname,".pop",sep="");
	for(i in 1:N)
	{	
		x = read.table(str, nrows=1,skip=(i-1));
		if(length(which(indlist == i)) > 0)
		{	tl = tr = ts = NULL;
			llll=length(x);
			for(j in 1:(llll-1))
			{	#p = regexpr(":",x[1,j+1])[1];
				#tl[j] = as.integer(substr(x[1,j+1],1,p-1));
				#ts[j] = as.integer(substr(x[1,j+1],p+1,nchar(as.character(x[1,j+1]))));
				tt = as.numeric(strsplit(as.character(x[1,j+1]), "[^0-9]+")[[1]]);
				tl[j]=tt[1];ts[j]=tt[2];
				if(j > 1) tr[j - 1] = tl[j];
				if(tl[j] > max(inv)) { break; }
			}
			tr[j] = L;
			z = which(tl <= min(inv));
			if(length(z) > 0)
			{	z = max(z);
				l = length(tl);
				tl=tl[z:l]; tl[1] = min(inv);
				tr=tr[z:l];
				ts=ts[z:l];
			}
			z = which(tr >= L);
			if(length(z) > 0)
			{	z = min(z);
				tl=tl[1:z];
				tr=tr[1:z]; tr[z] = L;
				ts=ts[1:z];
			}
			index = c(index, length(left) + 1);
			len = c(len, length(tl));
			left = c(left, tl);
			right = c(right, tr);
			ss = c(ss, ts);
		}
	}
	return(list(l=left, r=right, s=ss, i=index, z=len));
}
