#---------------------PARAMETERS--------------------------------------
ideaspath="./bin/";
walltime="24:00:00";
library="export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/aci/sw/quantumwise/2014.2/lib/";
library="export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/aci/sw/mathematica/11.0.1/SystemFiles/Libraries/Linux-x86-64/"
folder=getwd();
host="aci-b-test";

#----------------FUNCTION DEFINITION----------------------------------
runideas<-function(args, out, head, inv=NULL)
{	str=paste("ideas", paste(args, collapse=" "), "-thread", threadn, "-o", out);
	if(length(inv)==2)
	{	str=paste(str, "-inv", format(inv[1],scientific=F), format(inv[2],scientific=F)); }
	head = c(head, paste("#PBS -o ", out, ".script.out
", sep=""));
	writeLines(c(head, str),paste(out,".script",sep=""));

	if(file.exists(paste(out, ".para", sep="")))
	{	system(paste("rm ", out, ".para", sep=""));	}
	if(file.exists(paste(out, ".state", sep="")))
	{	system(paste("rm ", out, ".state", sep=""));	}
	if(file.exists(paste(out, ".cluster", sep="")))
	{	system(paste("rm ", out, ".cluster", sep=""));	}
	if(file.exists(paste(out, ".profile", sep="")))
	{	system(paste("rm ", out, ".profile", sep=""));	}

	print(str);
	#rt=system(paste("ssh ", host, " \"cd ",folder, "; qsub -A open ", out, ".script\"", sep=""));
	rt=system(paste("sh ", out, ".script", sep=""));
	system(paste("rm ", out, ".script", sep=""));
	if(rt!=0) stop("IDEAS quite unexpectedly");
}

defheader<-function(threadn=10)
{	head=paste("#PBS -l nodes=1:ppn=",threadn,"
#PBS -l walltime=",walltime,"
#PBS -j oe
cd ",
folder,
"
",library,"
PATH=$PATH:",ideaspath, "
",sep="");

	return(head);
}
makepara<-function(myorder,id,mem,mycut,para)
{       rt=NULL;
        j=0;
        for(i in myorder)
        {       t=which(mem==i);
                if(length(unique(id[t]))>=mycut)
                {       rt=rbind(rt,cbind(j,array(para[t,],dim=c(length(t),dim(para)[2]))));
                        j=j+1;
                }
        }

        return(rt);
}

getmean<-function(n)
{
        N=NULL;
        for(i in sort(unique(n[,1])))
        {       t=which(n[,1]==i);
                N=rbind(N,apply(array(n[t,],dim=c(length(t),dim(n)[2])),2,sum)[-1]);
        }
        NN=N[,-1]/N[,1];

        return(array(NN,dim=c(length(NN)/(dim(n)[2]-2),dim(n)[2]-2)));
}

compareTwo<-function(n, m, myplot=0)
{       NN=getmean(n);
        MM=getmean(m);

	p=(-3+sqrt(9+8*(dim(n)[2]-2)))/2;
        dd=NULL;
        for(i in 1:dim(NN)[1])
        {       dd[i]=min(apply(array(MM[,1:p],dim=c(dim(MM)[1],p)),1,function(x){sqrt(sum((x-NN[i,1:p])^2))}));
        }
        for(i in 1:dim(MM)[1])
        {       dd[i+dim(NN)[1]]=min(apply(array(NN[,1:p],dim=c(dim(NN)[1],p)),1,function(x){sqrt(sum((x-MM[i,1:p])^2))}));
        }

        if(myplot>0)
        {       M=rbind(NN[,1:myplot],MM[,1:myplot]);
                rownames(M)=c(1:dim(NN)[1],paste("R",1:dim(MM)[1],sep=""));
                source("/storage/home/yzz2/work/IDEAS/ideas_source/heatmap3.R");
                source("/storage/home/yzz2/work/IDEAS/ideas_source/draw.R");
                heatmapwrap(M,cols=c("black","blue","cyan","green","yellow","red","pink"));
        }
        return(dd);
}

combineState<-function(parafiles, method="ward.D", mycut=0.9, pcut=1., fout=NULL)
{	X=K=I=NULL;
	library("data.table");
	myheader=p=NULL;

	### remove para only has 1 state
	parafiles_no1 = c()
	for(i in 1:length(parafiles))
	{
		x=read.table(parafiles[i], header=T, comment.char='~')
		if (dim(x)[1]>1){
			parafiles_no1 = c(parafiles_no1, parafiles[i])
		}
	}
	parafiles = parafiles_no1

	for(i in 1:length(parafiles))
	{	x=fread(parafiles[i]);t=max(which(is.na(x[1,])==F));x=as.matrix(x[,1:t])
		x=as.matrix(rbind(x))
		x=as.matrix(x[which(x[,1]>=sum(x[,1])/1e4 & x[,1]>10),]);
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

waitFiles<-function(target)
{
	k=1;
	sel=rep(0,length(target));
	while(min(sel)==0)
	{	for(i in 1:length(target))
		{	sel[i] = as.integer(file.exists(paste(target[i],".para",sep="")));
			if(file.exists(paste(target[i], ".script.out", sep="")))
			{	x=tolower(readLines(paste(target[i], ".script.out", sep="")));
				if(max(as.integer(regexpr("abort",x)))>0) sel[i]=111;
				if(max(as.integer(regexpr("error",x)))>0) sel[i]=111;
				if(max(as.integer(regexpr("bad",x)))>0) sel[i]=111;
				if(max(as.integer(regexpr("unknown",x)))>0) sel[i]=111;
				if(max(as.integer(regexpr("cannot",x)))>0) sel[i]=111;
			}
		}
		if(min(sel)>0) break;
		if(k%%30==0) { cat("|"); }
		else if(k%%10==0) { cat("!"); }
		else if(k%%5==0) { cat(":"); }
		else { cat("."); }
		if(k%%60==0) { cat("\n"); }
		Sys.sleep(60);
		k=k+1;
		if(k > 4320) {
			message("ERROR: IDEAS didn't finish within 3 days\n");
			break;
		}
	}
	cat("\n");
	for(i in 1:length(target))
	{	if(sel[i] == 111) next;
		if(file.exists(paste(i,".script.out", sep="")))
		{	system(paste("rm ", i,".script.out", sep=""));	}
	}
	t=which(sel==111);
	if(length(t)>0) message(paste("IDEAS failed on these files:", target[t]));

	return(sel);
}

#--------------------RUN---------------------------------------------
args<-commandArgs(trailingOnly=TRUE);
#message(paste("args = ", paste(args,collapse=" "), collapse=""));
if(length(args) == 0)
{	message("Input arguments needed");
	quit();
}

splitstateflag=FALSE;
t=which(args == "-splitstate");
if(length(t)>0) splitstateflag = TRUE;

out = args[1];
t=which(args == "-o");
if(length(t)>0)
{	out = args[t[length(t)] + 1];
	args = args[-c(t,t+1)];
}

threadn = 10;
t=which(args == "-thread");
if(length(t)>0)
{	threadn = as.integer(args[t[length(t)]+1]);
	args = args[-c(t,t+1)];
}

type=0;
if(length(args)>1 & substring(args[2],1,1)!='-') 
{	type = 1;
	ncells=system(paste("awk -F ' |\t' '{print $1}'", args[1]), intern=T);
	ncells=length(unique(ncells));
} else
{	ncells=system(paste("head -n 1",args[1]), intern=T);
	ncells=unlist(strsplit(ncells,' |\t'));
	if(tolower(substr(ncells[3],1,3))=="pos") 
	{	ncells=ncells[-(1:3)];
	} else
	{	ncells=ncells[-(1:2)];
	}
	ncells=length(unique(unlist(strsplit(ncells,"\\."))[1:length(ncells)*2-1]));
}

inv=NULL;
t=which(args == "-inv");
if(length(t)>0)
{	inv = c(as.integer(args[t[length(t)]+1]), as.integer(args[t[length(t)]+2]));
	args = args[-c(t,t+1,t+2)];
}else
{	if(type == 1)
	{	a=system(paste("wc -l", args[2]), intern=TRUE)
	}else 
	{	a=system(paste("wc -l", args[1]), intern=TRUE)
	}
	inv=c(0, as.integer(strsplit(a," ")[[1]][1]));
}

randstart = c(0, max(50000, 1e6/ncells));
t=which(args == "-randstart");
if(length(t)>0)
{	randstart[1] = as.integer(args[t[length(t)]+1]);
	randstart[2] = as.integer(args[t[length(t)]+2]);
	args = args[-c(t, t+1, t+2)];
}
randstart[2] = min(randstart[2], inv[2] - inv[1]);
print(randstart);

split = NULL;
t=which(args == "-split");
if(length(t)>0)
{	split = args[t[length(t)]+1];
	args = args[-c(t,t+1)];
}

burnin = 20; mcmc = 1;
t=which(args == "-sample");
if(length(t)>0)
{	burnin = as.integer(args[t[length(t)]+1]);
	mcmc = as.integer(args[t[length(t)]+2]);
	args = args[-c(t,t+1,t+2)];
}
args = c(args, "-sample", burnin, mcmc);

head=defheader(threadn);

listfile = args[1];
if(type==1)
{	y=as.matrix(read.table(args[1]));
	ny=y;
	for(i in 1:dim(y)[1])
	{	a=unlist(strsplit(y[i,3],"\\."));
		if(a[length(a)]=="gz")
		{	ny[i,3] = paste(out, ".tmpdata.", i, sep="");	
			system(paste("gunzip -cf", y[i,3], ">", ny[i,3]));
		}
	}
	if(length(which(y[,3]!=ny[,3]))>0)
	{	args[1] = paste(out, ".newinput", sep="");
		write.table(ny, args[1], quote=F, row.names=F,col.names=F);
	}
}

if(randstart[1]>1)
{	
	targs = args;
	#tburnin = 20; tmcmc = 1;
	t=which(targs == "-impute");
	if(length(t)>0) 
	{	targs[t+1]="none";
	} else
	{	targs = c(targs, "-impute none");
	}
#targs=c(targs, "-C 50");
#targs = c(targs, "-K", "1");
mycut=0.5;

	#if(length(t)>0)
	#{	targs = targs[-c(t,t+1,t+2)];
	#}
	#targs = c(targs, "-sample", tburnin, tmcmc);
	set.seed(2019);
	for(i in 1:randstart[1])
	{	tout = paste(out, ".tmp.",i, sep=""); 
		tinv = round(runif(1, inv[1], inv[2]-randstart[2]));
		tinv = c(tinv, tinv + randstart[2]);
		if(type==1)
		{	tf=as.matrix(read.table(args[1]));
			tt=table(tf[,1]);
			spb = (tt / max(tt))^5;
			tn=min(length(tt), min(min(10,min(which(cumsum(spb[order(spb,decreasing=T)])/sum(spb)>0.8))), ceiling(1e6/randstart[2])));
			nn=sample(length(tt),size=tn,prob=spb);
			mm=which(is.na(match(tf[,1],names(tt)[nn]))==F);
			while(length(unique(tf[mm,2])) < length(unique(tf[,2])))
			{	tt = tt[-nn];
				spb = spb[-nn];
				nn=sample(length(tt),size=1,prob=spb);
				mm = sort(c(mm, which(is.na(match(tf[,1], names(tt)[nn]))==F)));
			}
			targs[1] = paste(args[1],".",i,sep="");	
			write.table(tf[mm,],targs[1],quote=F,row.names=F,col.names=F);
		}
		runideas(targs, tout, head, tinv);
		if(type==1)
		{	system(paste("rm",targs[1]));
		}
		Sys.sleep(1);
	}
	
	rt = waitFiles(paste(out, ".tmp.", 1:randstart[1], sep=""));
	if(length(which(rt==111))==randstart[1]) stop("IDEAS failed on all training data");

	if(splitstateflag)
	{	source("bin/splitstate.R");
		log2value = 0;
		t = which(targs == "-log2");
		if(length(t)>0) log2value = as.numeric(targs[t[1]+1]);
		for(i in 1:randstart[1])
		{	runsplit(targs[1], targs[2], paste(out,".tmp.",i,sep=""), log2value)
		}
		tpara = combineState(paste(out,".tmp.",(1:randstart[1])[rt!=111],".parasplit",sep=""), mycut=1-(1-mycut)/2);#,fout=paste(out,".mixpara",sep=""));
	} else
	{	tpara = combineState(paste(out,".tmp.",(1:randstart[1])[rt!=111],".para",sep=""), mycut=mycut);#,fout=paste(out,".mixpara",sep=""));
	}
	para = tpara$para;
	write.table(tpara$profile, paste(out, ".profile0", sep=""), quote=F,row.names=F,col.names=F);

	para = apply(para, 1, function(x){paste(x,collapse=" ")});
	para = c(readLines(paste(out, ".tmp.1.para",sep=""),n=1), para);
	outp = paste(out, ".para0", sep="");
	writeLines(para, outp);

	tt=which(args=="-prevrun");
	if(length(tt)>0) args=args[-c(tt,tt+1)];
	tt=which(args=="-otherpara");
	if(length(tt)>0) args=args[-c(tt,tt+1)];
	args = c(args, "-otherpara", outp, paste(out, ".profile0", sep="")); 
	if(length(which(args == "-G")) == 0) 
	{ args = c(args, "-G", length(para)-1); 
	} else
	{	tt = which(args == "-G");
		args[tt + 1] = length(para)-1;
	}
	tt = which(args == '-C');
	if(length(tt)>0) args = args[-c(tt,tt+1)];
	#args = c(args, "-mixpara", paste(out, ".mixpara",sep=""));

	for(i in 1:randstart[1])
	{	for(j in c("state","cluster","profile"))
		{	if(file.exists(paste(out, ".tmp.",i,".",j, sep="")))
			{	system(paste("rm ", out, ".tmp.",i,".", j, sep=""));	
			}
		}
	}
}

if(length(split) > 0)
{	tinv = as.matrix(read.table(split));
	inv = tinv[which(as.integer(tinv[,2]) + 1 < inv[2]),];	
	inv = array(inv, dim=c(length(inv)/3,3));
}else
{	inv = array(c("",inv), dim=c(1, 3));
}

target = NULL;
for(i in 1:dim(inv)[1])
{	if(dim(inv)[1]==1) { target[i] = out; }
	else { target[i] = paste(out, inv[i,1], sep="."); }
	runideas(args, target[i], head, inv[i,-1]);
}

rt = waitFiles(target);

if(listfile != args[1])
{	y=as.matrix(read.table(listfile));
	ny=as.matrix(read.table(args[1]));
	for(i in 1:dim(y)[1])
	{	if(y[i,3]!=ny[i,3])
		{	system(paste("rm",ny[i,3]));
		}
	}
	system(paste("rm",args[1]));
}

