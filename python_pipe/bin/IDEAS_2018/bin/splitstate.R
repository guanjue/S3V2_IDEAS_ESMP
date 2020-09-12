splitstate<-function(x,g,O,fout)
{	k=max(g)+1;
	d0 = apply(x,2,sd) + 1e-5;
	ng = rep(-1,length(g));
	S = 0;
	for(i in 1:k-1)
	{	t = which(g==i);
		if(length(t) == 0) next;
		r=NULL;
		if(length(t)>1)
		{	d = apply(array(x[t,],dim=c(length(t),dim(x)[2])),2,sd)+1e-5;
			tt = sample(length(t),size=min(10000,length(t)));
			r = which(d / d0 > max(1.5,quantile(d/d0,prob=1-3/length(d))));
			print(c(i,(d/d0)[r]));
		}
		if(length(r) == 0 | length(tt) < 2^length(r)*10) 
		{	ng[t] = S;
			S = S + 1;
			next;
		}
		m = cutree(hclust(dist(x[t[tt],r]),method="ward.D"),k=2^length(r));
		cc = NULL;
		for(j in unique(m))
		{	cc = rbind(cc,apply(array(x[t[tt][m==j],],dim=c(length(which(m==j)),dim(x)[2])),2,mean));
		}
		print(c(length(cc),dim(cc)));
		mm = kmeans(x[t,],centers=cc)$cluster;
		ng[t] = S + mm - 1;
		S = S + max(mm);
	}
	para = makepara(x,ng);
	write.table(para,paste(fout,".parasplit",sep=""),quote=F,row.names=F);

	prof = NULL;
	for(i in 0:max(O))
	{	t = which(O == i);
		tt = rep(1, max(ng)+1);
		if(length(t)>0)
		{	tt = tt + tabulate(ng[t]+1,nbins=max(ng)+1);
		}
		prof = rbind(prof, c(length(t), tt / sum(tt)));
	}
	write.table(prof, paste(fout, ".profilesplit",sep=""), quote=F,row.names=F,col.names=F);
	
	return(ng);
}

makepara<-function(x,g)
{	k=max(g)+1;
	p=dim(x)[2];
	para=NULL;
	for(i in 1:k-1)
	{	print(i);
		t=which(g==i);
		n=length(t);
		if(n<10) next;
		m=apply(x[t,],2,sum);
		v=NULL;
		for(j in 1:p)
		{	for(jj in 1:j)
			{	v = c(v, sum(x[t,j]*x[t,jj]));
			}
		}
		para = rbind(para, c(n, m, v));
	}
	mark=colnames(x);
	lab=c("#count",mark);
	for(j in 1:p)
	{	for(jj in 1:j)
		{	lab = c(lab, paste(mark[j],mark[jj],sep="*"));
		}
	}	
	colnames(para) = lab;
	return(para);
}

runsplit<-function(finput, fbed, fpref, log2value)
{	x=read.table(finput);
	mark=unique(as.matrix(x[,2]));
	t = table(as.matrix(x[,1]));
	tt = which(t == length(mark));
	if(length(tt) > 0)
	{	cell = names(t[tt]);
		library("data.table");
		g = fread(paste(fpref,".state",sep=""));
		pid = as.matrix(g[,1]);
		g = as.matrix(g[,-(1:4)]);
		bid = fread(fbed);
		bid = as.matrix(bid[,4]);
		rg = range(match(pid, bid));	
		data = state = O = NULL;
		for(i in cell)
		{	cm = which(colnames(g)==i);
			if(length(cm)==0) next;
			cm = cm[1];
			td = NULL;
			for(j in mark)
			{	print(c(i,j));
				t = which(x[,1] == i & x[,2] == j)[1];
				td = cbind(td, as.matrix(read.table(as.matrix(x[t,3]),skip=rg[1]-1,nrows=rg[2]-rg[1]+1)));
			}
			data = rbind(data, td);
			state = c(state, g[,cm]);
			O = c(O, g[,dim(g)[2]]);
			if(length(state) > 1000000) break;
		}
		if(log2value > 0) data = log2(data + log2value);

		colnames(data) = mark;
		ns = splitstate(data, state, O, fpref);
	} else {
		system(paste("cp ", fpref, ".para ", fpref, ".parasplit", sep=""));	
		system(paste("cp ", fpref, ".profile ", fpref, ".profilesplit", sep=""));	
	}
}
