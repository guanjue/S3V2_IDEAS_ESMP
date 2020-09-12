binning<-function(z, input)
{	cell=unique(input[,1]);
	mark=unique(input[,2]);
	k = length(mark);
	myset = t(sapply(1:(2^k)-1,function(x){as.integer(intToBits(x))[1:k]}));

	cn=rep(0,2^k);
	rt=array(-1,dim=c(dim(z)[1],length(cell)));
	for(i in order(table(input[,1]),decreasing=T))
	{	print(cell[i]);
		m=match(input[input[,1]==cell[i],2],mark);
		tk=length(m);
		tset = sapply(1:(2^tk)-1,function(x){as.integer(intToBits(x))[1:tk]});
		if(tk > 1) tset = t(tset);
		t1=apply(array(myset[,m],dim=c(dim(myset)[1],length(m))),1,function(x){paste(x,collapse="")});
		t2=apply(array(tset,dim=c(length(tset)/tk,tk)),1,function(x){paste(x,collapse="")});
		mm=rep(-1,2^k);
		for(j in 1:length(t2))
		{	mm[which(is.na(match(t1,t2[j]))==F)] = j-1;
		}
			
		t=which(input[,1]==cell[i]);
		u=rep(0,dim(z)[1]);
		for(j in tk:1)
		{	u=u*2+z[,t[j]];
		}
		
		trt=rep(-1,dim(z)[1]);
		for(j in 1:length(t2))
		{	tj=which(u==j-1);
			tm=which(mm==j-1);
			tt=tm[sample(length(tm),size=length(tj),replace=T,prob=cn[tm]+0.5)];
			trt[tj] = tt-1;
			cn = cn + tabulate(tt,nbins=length(cn));
		}
		rt[,i]=trt;
	}
	RT=NULL;
	colnames(rt)=cell;
	RT$state=rt;
	RT$cn=cn;
	return(RT);		
}

createPara<-function(data, input, state, cn)
{	t=which(cn/sum(cn)>=1e-4);
	map=rep(-1,length(cn));
	map[t[order(cn[t],decreasing=T)]]=1:length(t)-1;
	
	cell=unique(input[,1]);
	mark=unique(input[,2]);
	l=length(mark);
	l2=l*(l+1)/2;
	n=length(t);
	para=array(0, dim=c(n,1+l+l2));
	for(i in 1:dim(state)[2])
	{	cat(i, ":");
		state[,i] = map[state[,i]+1];		
		m=which(input[,1]==cell[i]);
		m=m[match(mark,input[m,2])];
		tt=NULL;
		for(jj in 1:l)for(kk in 1:jj){tt=c(tt,(kk-1)*l+jj);}
		for(j in 1:n)
		{	cat(j,",");
			t=which(state[,i]==j-1);
			if(length(t)>0)
			{	para[j,1]=para[j,1]+length(t);
				v=array(data[t,m],dim=c(length(t),length(m)));
				para[j,1+1:l]=para[j,1+1:l]+apply(v,2,sum);
				v=t(v)%*%v;
				para[j,1+l+1:l2] = para[j,1+l+1:l2] + v[tt];
			}
		}
		if(length(which(state[,i]<0))>0)
			state[state[,i]<0,i]=sample(n,replace=T,size=length(which(state[,i]<0)),prob=tabulate(state[,i]+2,nbins=n+1)[-1])-1;
		print("");
	}
	colnames(state) = cell;
	mark2=NULL;
	for(i in 1:l)for(j in 1:i) mark2=c(mark2,paste(mark[i],mark[j],sep="*"));
	colnames(para) = c("#Count",mark,mark2);
	rt$state=state;
	rt$para=para;
	return(rt);
}
