combineState<-function(parafiles, method="ward.D", mycut=0.9, pcut=1., fout=NULL)
{	X=K=I=NULL;
	library("data.table");
	myheader=p=NULL;
	for(i in 1:length(parafiles))
	{	x=fread(parafiles[i]);t=max(which(is.na(x[1,])==F));x=as.matrix(x[,1:t])
x=x[which(x[,1]>=sum(x[,1])/1e4 & x[,1]>10),];
		#x=read.table(parafiles[i]);
		if(i==1)
		{	myheader=colnames(x);
			p=sqrt(9/4-2*(1-length(myheader)))-3/2;
		}
		m=match(myheader[1:p+1],colnames(x)[1:p+1]);
		v=NULL;for(ii in 1:p)for(jj in 1:ii)
		{a=max(m[ii],m[jj]);b=min(m[ii],m[jj]);v=c(v,a*(a+1)/2+b-a);}
		X=rbind(X, array(as.matrix(x[,c(1,1+m,1+p+v)]),dim=c(length(x)/(1+p+length(v)),1+p+length(v))));
		K=c(K,dim(x)[1]);
		I=c(I,rep(i,dim(x)[1]));
	}
	N=length(parafiles);
	
	p = (sqrt(1 + dim(X)[2] * 8) - 3) / 2;
	omycut=mycut;
	mycut=round(length(parafiles)*mycut);

	M=array(X[,1:p+1]/X[,1], dim=c(dim(X)[1],p));
	V=array(0,dim=c(dim(X)[1]*p,p));
	for(i in 1:dim(X)[1])
	{	t = (i - 1) * p;
		l = 1;
		for(j in 1:p)
		{	for(k in 1:j)
			{	V[t + j, k] = V[t + k, j] = X[i,1+p+l] / X[i,1] - M[i,j] * M[i,k];
				l=l+1;
			}
		}
		V[t+1:p,]=t(solve(chol(V[t+1:p,]+diag(1e-1,p))));
#diag(V[t+1:p,])=sqrt(diag(V[t+1:p,]));
	}


	D=array(0,dim=rep(dim(X)[1],2));
	for(i in 2:dim(X)[1])
		for(j in 1:(i-1))
		{	D[i,j]=D[j,i]=sqrt((sum((V[(i-1)*p+1:p,]%*%(M[i,]-M[j,]))^2) + sum((V[(j-1)*p+1:p,]%*%(M[i,]-M[j,]))^2)));
			#D[i,j]=D[j,i]=sqrt(sum((M[i,]-M[j,])^2/(diag(V[(i-1)*p+1:p,])^2+diag(V[(j-1)*p+1:p,])^2)));
			#D[i,j]=D[j,i]=sqrt(sum((M[i,]-M[j,])^2));#/(diag(V[(i-1)*p+1:p,])^2+diag(V[(j-1)*p+1:p,])^2)));
		}

	MM = NULL;
	kk=NULL;
	for(i in 1:N)
	{	t=1:K[i];
		if(i > 1) t = t + sum(K[1:(i-1)]);
		t=(1:dim(D)[1])[-t];
		h=hclust(as.dist(D[t,t]),method=method);
		k=-1;
		tM = NULL;
		for(j in min(K):(min(length(t),max(K)*2)))
		{	m=cutree(h,k=j);
			tt=NULL;
			for(l in 1:j)
			{	tt[l] = length(unique(I[t[which(m==l)]]));
			}
			tk=length(which(tt>=mycut));
			if(tk > k)
			{	k = tk;
				tM = makepara(1:j,I[t],m,mycut,X[t,]);
			} else if(tk < k) { break; }
		}
		kk[i]=k;
		MM = rbind(MM, cbind(i, tM));
	}

	mysel = median(kk);
	h=hclust(as.dist(D),method=method);
	rt=rep(0,max(K)*2);
	k = -1;
	for(i in min(K):min(dim(D)[1],max(K)*2))
	{	m=cutree(h,k=i);
		tt = NULL;
		for(j in 1:i)
		{	tt[j] = length(unique(I[which(m==j)]));
		}
		tk=length(which(tt>=mycut));
		if(tk==mysel | tk < k) break;
		k = tk;
		rt[i] = length(which(tt>=mycut));
	}
	mysel = max(k,tk);
	
	m=cutree(h,k=mysel);
	nn = NULL;
	for(i in 1:mysel)
	{	t=which(m==i);
		nn[i] = sum(X[t,1]);
	}
	oo=order(nn,decreasing=T);
	rt = makepara(oo,I,m,mycut,X);
	onstate = max(rt[,1])+1;
        ooo=NULL;
	for(i in oo)
	{	t=which(m==i);
        	if(length(unique(I[t]))>=mycut) ooo=c(ooo,i);
	}

	d=NULL;
	for(i in 1:N)
	{	d=rbind(d,compareTwo(rt,MM[MM[,1]==i,-1])[1:onstate]);
	}
	#boxplot(d);
	dd=array(cutree(hclust(dist(c(d))),k=2),dim=dim(d));
	kk=table(c(dd));
	kk=which(as.integer(kk)==max(as.integer(kk)))[1];
	pp=apply(dd,2,function(x){length(which(x!=kk))/length(x)});
	pp0=apply(d,2,function(x){length(which(x>0.5))/length(x)});
	#print(rbind(pp,pp0));
	pp[pp0<pp]=pp0[pp0<pp];
	t=which(pp>pcut);
	if(length(t)>0)
	{	j=0;
		nrt = NULL;
		for(i in (1:onstate-1)[-t])
		{	nrt = rbind(nrt, cbind(j, rt[rt[,1]==i,-1]));
			j=j+1;
		}
		rt=nrt;
		ooo=ooo[-t];
	}
	message(paste("nstate =", onstate, "-->", max(rt[,1])+1));

	nrt=NULL;
	for(i in 0:max(rt[,1]))
	{	t=which(rt[,1]==i);
		nrt=rbind(nrt, apply(array(rt[t,],dim=c(length(t),dim(rt)[2])),2,sum)[-1]);
	}
	rt=nrt;
	colnames(rt)=myheader;

	if(length(fout)>0)
	{	#l=as.matrix(read.table(parafiles[1],comment="!",nrows=1));
		#l[1]=gsub("#","",l[1]);
		#l=c("#state",l);
		#colnames(rt)=l;
		write.table(rt,fout,quote=F,row.names=F);
	}

	O=NULL;
	#if(FALSE)
	{	Ip = Xp = NULL;
		k = 0;
		for(i in 1:length(parafiles))
		{	str=gsub(".para",".profile",parafiles[i]);
			p=as.matrix(read.table(str));
			u=array(0,dim=c(dim(p)[1],length(ooo)));
			for(j in 1:length(ooo))
			{	t=which(m[k+1:K[i]] == ooo[j]);
				u[,j] = apply(array(p[,1+t],dim=c(dim(p)[1],length(t))),1,sum);
			}
			k = k + K[i];
			u=u / (apply(u,1,sum)+1e-10);
			Xp = rbind(Xp, cbind(p[,1],u));
			Ip = c(Ip, rep(i,dim(u)[1]));
		}
		#hp=hclust(dist((log(Xp[,-1] + min(1e-3,min(Xp[,-1][Xp[,-1]>0]))))),method="ward");
		hp=hclust(dist(((Xp[,-1] + min(1e-3,min(Xp[,-1][Xp[,-1]>0]))))),method=method);
		ocut = min(mycut / 2, length(parafiles) / 2);
		t=range(as.integer(table(Ip)));
		Kp = NULL;
		for(i in t[1]:(t[2]*2))
		{	m=cutree(hp,k=i);
			tt=table(Ip,m);
			ll=apply(tt,2,function(x){length(which(x>0))});
			Kp = c(Kp, length(which(ll>=ocut)));
		}
		oN=(t[1]:(t[2]*2))[which(Kp==max(Kp))[1]];
		m=cutree(hp,k=oN);
		tt=table(Ip,m);
		ll=apply(tt,2,function(x){length(which(x>0))});
		tt=which(ll>=ocut);
		
		for(i in tt)
		{	t=which(m==i);
			O=rbind(O,c(sum(Xp[t,1]),apply(array(Xp[t,-1]*Xp[t,1],dim=c(length(t),dim(Xp)[2]-1)),2,sum)/sum(Xp[t,1])));
		}
		
		if(length(fout) > 0)
		{	if(as.integer(regexpr(".para",fout)) > 0) 
			{	tfout = gsub(".para",".profile",fout);
			} else
			{	tfout = paste(fout, ".profile", sep="");
			}
			write.table(O,tfout, quote=F,row.names=F,col.names=F);
		}
	}
	
	nrt=NULL;
	nrt$para = rt;
	nrt$profile = O;

	return(nrt);
}
