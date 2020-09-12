#include "genomicTensor.h"
#define STDFORMAT 1 

extern bool hpass, splitmerge, gzip;
extern int bST, bED; //only works for list data
extern vector<string> imputelist;

unsigned int timeA, timeB;
bool adjustB = true; //false: no adjust for B when update hyperparameter; true: adjust
bool useSoftState = false; //use tmpSpace to calculate marginalized emission probability integrating out position state
bool stateinteraction = false;//individual state x position state, or independence

typedef struct MySortType {
        int index;
        double score, weight; 
} MySortType;

bool operator<(const MySortType &a, const MySortType &b)
{       return a.score < b.score;
} 

typedef struct MySortType2 {
        int id1, id2;
	string str;
} MySortType2;

bool operator<(const MySortType2 &a, const MySortType2 &b)
{       return a.str < b.str;
}
//things learned so far
//1) position specific Vh helps pops to converge
//2) taking imputation maximums helps gauss to converge
//3) hyper-parameters of priors should not change, although priors can be updated, otherwise priors may go to 0
//4) including position specific classes may need stronger priors, otherwise it will be hard to converge in small samples
//5) priors for emitMatrix should depend on distribution of gauss components, otherwise may be hard to converge in small samples

clock_t Time1, Time2, Time3;

genomicTensor::genomicTensor(unsigned int rs):tensorHMMbase(rs)
{
timeA=timeB=0;
Time1=Time2=Time3=0;
		
	TLP = NULL;
	dataYP = dataXP = NULL;
	emitMatrix = NULL;
	gClass = NULL;
	posClassSZ = gClassSZ = 0;
	emitK = 0;
	maxG = 100;
	maxPosSZ = 100;
	independentInd = false;
	log2 = -1.;
	norm = false;

	xmsz = ymsz = NULL;
	xmap = ymap = NULL;

	fixstateN = false;
	nthread = 1;
}

genomicTensor::~genomicTensor()
{	if(dataYP != NULL) 
	{	for(int i = 0; i < mydata.totalN; i++) if(dataYP[i] != NULL) delete dataYP[i]; 
		delete dataYP;
	}
	if(dataXP != NULL) 
	{	for(int i = 0; i < mydata.totalN; i++) if(dataXP[i] != NULL) delete dataXP[i];
		delete dataXP;
	}
	if(emitMatrix != NULL) delete emitMatrix;
	if(gClass != NULL) delete gClass;

	if(ymsz != NULL) delete ymsz;
	if(xmsz != NULL) delete xmsz;
	if(xmap != NULL) { for(int i = 0; i < mydata.totalN; i++) delete xmap[i]; delete xmap; }
	if(ymap != NULL) { for(int i = 0; i < mydata.totalN; i++) delete ymap[i]; delete ymap; }
}

void genomicTensor::_getAsz(MYDATA const &mydata)
{
	int i;
        for(i = 0; i < mydata.L; i++) mydata.asz[i] = mylik.clustersz;
}

/////////////////////////////////////////////////////////////////////////////////////
void genomicTensor::run(char const *finput, char const *fbed, char const *fcov, char const *fcovbed, int burnin, int mcmc, int thread, int maxHapK, int maxGG, int maxPos, double AA, double recr, double heteroVh, bool samplemaximum, char const *foutput, bool sqc, char const *fparam, bool add2, double merr, double mxerr, bool ind, bool indind, int startK, int fixC, vector<vector<int> > const &removelist, bool nb, char const *fmixpara, char const *parafile, char const *statefile, char const *clusterfile, char const *para0file, char const *state0file, char const *cluster0file, char const *profile0file)
{

	int i, j;
	fixstateN = (bool)(fixC > 0);
    nthread = thread;
	
	tpara0.param0 = NULL;
	tpara0.nh0 = NULL;
	tpara0.recombN0 = NULL;
	tpara0.Q0 = NULL;
	tpara0.clustersz0 = tpara0.popn0 = tpara0.unitsz0 = tpara0.totalN0 = 0;

	mydata.ploidity = 1;
	mydata.L = -1;
	initPara(tpara);
	tpara.priorW = 10.;
	tpara.maxHapK = maxHapK;
	if(maxGG > 0) maxG = maxGG;
	maxG = max(maxG, startK);
	if(maxPos > 0) maxPosSZ = maxPos;
	tpara.burnin = burnin;
	tpara.mcmc = mcmc;
	tpara.changeAlpha = 0;
	tpara.maximum = 0;//abc, recomb, state, allele
	tpara.samplemaximum = samplemaximum;//true;
	tpara.defmut = 0.1;
	tpara.A = 1.;
	tpara.recrate = recr;
	tpara.recombcut = 0;
	tpara.heteroVh = heteroVh; 
	tpara.samplefull = false; //fullsample is for haplotypes only, sample genome data by gClass

	mylik.minerr = merr;//min(0.01, overallsd / 10.);
	mylik.maxerr = mxerr;
	mylik.indc = ind;
	independentInd = indind;
	mylik.nb = nb;
	overallmean = overallmeanx = 0;
	overallsd = overallsdx = 1.;

	if(fparam != NULL) 
	{
		parseParameter(fparam);
		bool flag = false;
		for(i = (int)mylik.modelparameter.size() - 1; i > maxG; i--) 
			if((int)mylik.modelparameter[i].size() > 0)
			{	printf("cluster %d ignored as it is beyond the maximum # of clusters allowed (%d).\n", i-1, maxG);	
				flag = true;
			}
		if(flag) mylik.modelparameter.resize(maxG + 1);
	}

	readInput(mydata, finput, fbed, fcov, fcovbed, add2); //need revision for fix state and fix alleles
	removeData(removelist);

	if(fmixpara != NULL)
	{	ifstream f(fmixpara);
		string tmp;
		getline(f, tmp);
		int k = 0;
		while(getline(f, tmp))
		{	istringstream buf(tmp);
			istream_iterator<string> beg(buf), end;
			vector<string> tokens(beg, end);
			k = max(k, atoi(tokens[0].c_str()));
		}
		maxG = startK = mylik.clustersz = k + 1;
		splitmerge = false;
	}

/////////////////
	//tpara.probAlpha = mydata.totalN/mydata.ploidity;//0.1;//min(1., max(1e-1, 10./ (double)mydata.totalN));
	tpara.probAlpha = tpara.trueindN;
	tpara.probAlpha = max(1., tpara.probAlpha / 5.);
	if(AA > 0) tpara.probAlpha = AA;//tpara.probAlpha * AA; //changed 02/07/2015

//	tpara.A = tpara.probAlpha; //added on Jul 11th

	if(parafile != NULL)
	{	FILE *ff = fopen(parafile, "r");
		char tmppara[10000];
		startK = 0;
		fgets(tmppara, 10000, ff);
		while(fgets(tmppara, 10000, ff) != NULL) startK++;
		fclose(ff);
	}
	if(para0file != NULL)
	{	ifstream f(para0file);
		string tmp;
		getline(f, tmp);
		istringstream buf(tmp);
		istream_iterator<string> beg(buf), end;
	        vector<string> tokens(beg, end);

		int ymsz0 = 0, xmsz0 = 0;
		int i, j, k;
		for(i = 1; i < (int)tokens.size(); i++)
		{	k = (int)tokens[i].find("*");
			if(k >= 0) break;
			else 
			{	for(j = 0; j < maxymsz; j++)
					if(mydata.fyinfo[j] == tokens[i]) break;
				if(j >= maxymsz)
					mydata.fyinfo.push_back(tokens[i]);
			}
		}
		ymsz0 = i - 1;
		j = 1 + ymsz0 + ymsz0 * (ymsz0 + 1) / 2;
		for(i = j; i < (int)tokens.size(); i++)
		{	k = (int)tokens[i].find("*");
			if(k >= 0) break;
			else
			{	for(j = 0; j < maxxmsz; j++)
					if(mydata.fxinfo[j] == tokens[i]) break;
				if(j >= maxxmsz)
					mydata.fxinfo.push_back(tokens[i]);
			}
		}
		xmsz0 = i - j;

		maxymsz = (int)mydata.fyinfo.size();
		maxxmsz = (int)mydata.fxinfo.size();

		k = 0;
		while(getline(f, tmp)) k++;
		f.close();
		if(startK < k) startK = k;
		mylik.loadGaussPrior0(para0file, mydata.fyinfo, mydata.fxinfo);
	}

	dataNorm();
	printf("N=%d Marks=%d, Cov=%d, Mean=%f, Sd=%f\n", mydata.totalN, maxymsz, maxxmsz, overallmean, overallsd);fflush(stdout);

	if(startK<=0) startK=20;
	tpara.A = (double)startK / 2.;
	if(fixstateN) { maxG = fixC; }
	mylik.clustersz = startK = min(startK, maxG); 
	posSZ = max(1,min(startK / 2, maxPosSZ));
	if(profile0file != NULL)
	{	FILE *f = fopen(profile0file, "r");
		char tmp[10000];
		while(fgets(tmp, 10000, f) != NULL)
		{	vector<double> row;
			row.push_back(atof(&tmp[0]));
			for(i = 1; i < (int)strlen(tmp) - 1; i++)
			{	if(tmp[i] == ' ')
				{	row.push_back(atof(&tmp[i + 1]));
					i++;
				}
			}
			emit0.push_back(row);
		}	
		fclose(f);	
		maxPosSZ = posSZ = (int)emit0.size();
	}

	mylik.initializePara(startK, dataYP, dataXP, mydata.totalN, mydata.L, maxymsz, ymsz, ymap, maxxmsz, xmsz, xmap, tpara.priorW, parafile);

	double posRate = 0.1;
	basePop = new int[mydata.L];
	baseR = new bool[mydata.L];
	double tm[posSZ], tv[posSZ];
	for(i = 0; i < posSZ; i++) 
	{	tm[i] = 0;
		if((int)emit0.size() > 0) tm[i] = emit0[i][0];
	}
	_updateVh(tm, posSZ, tpara.A, tv);
	for(i = 1; i < posSZ; i++) tv[i] += tv[i - 1];
	for(i = 0; i < mydata.L; i++) 
	{	baseR[i] = false;
		if(i == 0 || gsl_ran_flat(gammar, 0, 1.) < posRate)
		{	double un = gsl_ran_flat(gammar, 0, tv[posSZ-1]);
			for(j = 0; j < posSZ; j++) if(un <= tv[j]) break;
			basePop[i] = j;
			baseR[i] = true;
		}
		else basePop[i] = basePop[i - 1];
	}

	vector<bool> myinit(mydata.indN, (statefile == NULL));
	if(state0file != NULL && cluster0file != NULL) loadOtherResult(state0file, cluster0file, posRate);
	if(statefile != NULL) loadPreviousRun(statefile, clusterfile, posRate);

	if(fmixpara != NULL)
	{	mylik.preComputeLP(fmixpara, mydata.fyinfo, mydata.fxinfo, dataYP, dataXP, mydata.totalN, mydata.L);
		maxG = startK = mylik.clustersz;
		splitmerge = false;
	}

	tmpSpace = tmpSpace2 = NULL;//new double[100 * mydata.L * maxG];
	mylik.updateLambda(tpara.priorW, true, false);
//printf("<%f %f %f | %f,%f,%f>\n", mylik.gaussprior[0], mylik.gaussprior[mylik.gausssz], mylik.gaussprior[mylik.gausssz*2],mylik.gaussprior0[0], mylik.gaussprior0[mylik.gausssz], mylik.gaussprior0[mylik.gausssz*2]);fflush(stdout);
/*
for(i=0;i<mylik.clustersz;i++)
{	for(j=0;j<mylik.maxymsz+1;j++) printf("%f,",mylik.gaussprior[i*mylik.gausssz+j]);
	printf("\n");
}
for(i = 0; i < mydata.totalN; i++)
{	float *dpy = dataYP[i];
	double tlp[mylik.clustersz];
	for(j = 0; j < mydata.L; j++, dpy += ymsz[i])
	{	mylik.computeLP(dpy, NULL, i, tpara.priorW, tlp);
		int k, l = 0;
		for(k = 1; k < mylik.clustersz; k++)
			if(tlp[l] < tlp[k]) l = k;
		mydata.data[i * mydata.L + j] = l;
if(i == 0 && j < 10) 
{	for(k=0;k<mylik.clustersz; k++) printf("%f ", tlp[k]);
	printf("\n");
}
	}
}
for(i=0;i<mylik.clustersz;i++)
{	for(j=0;j<mylik.maxymsz+1;j++) printf("%f,",mylik.gaussprior[i*mylik.gausssz+j]);
	printf("\n");
}
FILE *fff = fopen("tmp.txt","w");
for(i = 0; i < mydata.L; i++)
{	for(j = 0; j < mydata.totalN; j++)
	{	fprintf(fff, "%d ", (int)mydata.data[j * mydata.L + i]);
		mydata.data[j*mydata.L+i]=-1;
	}
	fprintf(fff, "\n");
}
fclose(fff);
exit(0);
*/
	inferStructure(mydata, tpara, foutput, sqc, myinit);
	if(tmpSpace != NULL) delete tmpSpace;
	if(tmpSpace2 != NULL) delete tmpSpace2;
	
	//outputResult(foutput);
//printf("updatePara=%f, updateBase=%f\n", (double)timeA/1000000., (double)timeB/1000000.);fflush(stdout);
//printf("baseForward=%f, baseBackward=%f\n", (double)Time1/1000000., (double)Time2/1000000.);fflush(stdout);

//printf("Time3=%f\n", (double)Time3/1e6);fflush(stdout);
	
	if(mydata.data != NULL) delete mydata.data;
	if(mydata.indIndex != NULL) delete mydata.indIndex;
	if(tpara.indWeight != NULL) delete tpara.indWeight;
	
	delete basePop;
	delete baseR;

	if(tpara0.param0 != NULL) delete tpara0.param0;
	if(tpara0.nh0 != NULL) delete tpara0.nh0;
	if(tpara0.recombN0 != NULL) delete tpara0.recombN0;
	if(tpara0.Q0 != NULL) delete tpara0.Q0;
}

/////////////////////////////////////////////////////////////////////////////////////
void genomicTensor::_updatePara_remove(PARAUNIT &para, float data, float b, int id, int pos)
{	//float *pdy = dataYP[id] + pos * ymsz[id], *pdx = NULL;
	//if(maxxmsz > 0) pdx = dataXP[id] + pos * xmsz[id];
/*	if(pos==0) 
	{	pdy = dataYP[id];
		ystepsz = ymsz[id];
		if(maxxmsz > 0) pdx = dataXP[id];
		xstepsz = xmsz[id];
	}
*/
	int i;
	for(i = 0; i < (int)para.allele.size(); i+=2)
		if(para.allele[i] == data) break;
/*if(i >= (int)para.allele.size()-1)
{	printf("data=%d, sz=%d, maxK=%d, id=%d, pos=%d ", (int)data, (int)para.allele.size(), tpara.maxK, id, pos);
	for(i = 0; i < (int)para.allele.size(); i++) printf("%f,",para.allele[i]);
	printf("\n");
	for(i = 0; i < (int)tpara.param[pos].size(); i++) 
	{	printf("%d: ", i);fflush(stdout);
		for(int j = 0; j < (int)tpara.param[pos][i].allele.size(); j++)
			printf("%f,",tpara.param[pos][i].allele[j]);
		printf("\n");
	}
	exit(0);
}*/
	para.lapseN -= pos;
	para.popN--;
//if(para.allele[i+1]<(int)b)
//{	printf("%d,%d", para.allele[i+1], (int)b);fflush(stdout);
// 	exit(0);
//}
	para.allele[i+1] -= b;
	if(para.allele[i+1] < NUMPRECISION) para.allele.erase(para.allele.begin() + i, para.allele.begin() + i + 2);
//printf(";");fflush(stdout);
	//mylik.removePara((int)data, pdy, pdx, id);

//	pdy += ystepsz;
//	pdx += xstepsz;
}

void genomicTensor::_updatePara_add(PARAUNIT &para, float data, float b, int id, int pos)
{	//float *pdy = dataYP[id] + pos * ymsz[id], *pdx = NULL;  
	//if(maxxmsz > 0) pdx = dataXP[id] + pos * xmsz[id];
/*	if(pos==0) 
	{	pdy = dataYP[id];
		ystepsz = ymsz[id];
		if(maxxmsz > 0) pdx = dataXP[id];
		xstepsz = xmsz[id];
	}
*/	
	int i;
	para.lapseN += pos;
	para.popN ++;
	for(i = 0; i < (int)para.allele.size(); i+=2)
		if(para.allele[i] == data) break;
	if(i >= (int)para.allele.size()) 
	{	para.allele.resize(i + 2, 0);
		para.allele[i] = data;
	}
	para.allele[i+1] += b;
	//mylik.addPara((int)data, pdy, pdx, id);

//	pdy += ystepsz;
//	pdx += xstepsz;
}


//need to update pointers, bst and bed are breaks index, if they exist, otherwise they are positions
void genomicTensor::_sumforwardOne(int id, vector<int> const &breaks, int bst, int bed, TENSORPARA &tpara, int ploidity, float *&dp, float *&dp2, double const *&rp, float const *&nhp, vector<vector<PARAUNIT> >::const_iterator &ipm, int const *&aszp, double *&lpp, double *&lmm, double *&lss, int SSZ, int MK, vector<vector<int> > const &ss, vector<vector<int> > const &mlist, double *Vh, double *Vh2)
{	int i, j;
	double f, t, maxlp = MINUSINFINITE;

double aaa = tpara.heteroVh;
double NNN = aaa * (double)(tpara.trueindN * mydata.ploidity + tpara0.totalN0) + 1.;

	int ided = mydata.indIndex[id+1], idst = mydata.indIndex[id];
//clock_t cst, ced;
//cst=clock();
	if(*dp < 0 && bst == 0)
	{	if(TLP != NULL) { printf("!!!TLP!=NULL\n"),exit(0); }
		TLP = new double[SSZ * mydata.L];
		int nstep = max(1000,(int)ceil((double)mydata.L / (double)nthread));
		int nsz = max(1,(int)ceil((double)mydata.L / (double)nstep));
		

    #pragma omp parallel num_threads(nthread)
	{
        //std::cout << "Thread number: " << omp_get_thread_num() << endl;
		# pragma omp for
		for(int nn = 0; nn < nsz; nn++)
		{	
			int posst = nn * nstep, posed = min(mydata.L, posst + nstep);	
			double *tlp = TLP + SSZ * posst;
			for(int j = posst; j < posed; j++)
			{	for(int K = 0; K < SSZ; K++, tlp++)
				{	double prop[*aszp + 1];
					_getProp(tpara, id, j, K, *aszp, ipm + j, nhp + tpara.maxK * j, basePop + j, prop);
					*tlp = 0;
					for(int i = idst; i < ided; i++)
					{	float *dpy = dataYP[i] + j * ymsz[i], *dpx;
						if(maxxmsz > 0) dpx = dataXP[i] + j * xmsz[i];
						double dlp[*aszp];
						if(mylik.preLP != NULL)
						{	float *prelpp = mylik.preLP + mylik.clustersz * (i * mydata.L + j);
							for(int ii = 0; ii < mylik.clustersz; ii++, prelpp++)
								dlp[ii] = (double)(*prelpp);
						}
						else
						{	mylik.computeLP(dpy, dpx, i, tpara.priorW, dlp);	
						}
						int k = 0;
						dlp[0] += fast_log(prop[0]);
						for(int l = 1; l < *aszp; l++)
						{	dlp[l] += fast_log(prop[l]);
							if(dlp[k] < dlp[l]) k = l;
						}
						*tlp += dlp[k] - fast_log(prop[*aszp]);
						if(i < ided - 1) { prop[k]++; prop[*aszp]++; }
					}
				}
			}
		}
	}
	}
//ced=clock();
//Time3+=ced-cst;

	(*lss) = 0;
	if(bst <= 0)
	{	for(j = 0; j < SSZ; j++, lpp++)
        	{       (*lpp) = 0;
			f = _getLp(id, breaks, bst, false, ploidity, tpara, dp, dp2, &ss[j][0], &ss[j][1], ipm, nhp, aszp);

			lmm[j] = f;
			//t = fast_log(tpara.Vh[ss[j][0]]);
			int state = ss[j][0];
			//double tnnn = (double)(*(nhp + state));
			//if(j < tpara0.popn0) tnnn += (double)(*(tpara0.nh0 + j));
			t = fast_log((Vh[state] + aaa*(double)(state < tpara.maxK)*(double)(*(nhp+state))/*tnnn*/) / NNN);
			(*lpp) = t + f;
			if(independentInd && j != id) { (*lpp) = MINUSINFINITE; }
			if(j == 0 || (*lpp) > maxlp) maxlp = (*lpp);
        	}
		lpp -= SSZ;
		for(j = 0; j < SSZ; j++, lpp++)
			(*lss) += exp((*lpp) - maxlp);
		(*lss) = fast_log(*lss) + maxlp;
	}
	else
	{	double recp = *rp, nrecp = 1. - recp;
		if(tpara.singlecall) { recp = 1.-1e-10; nrecp = 1e-10; }//recp2 = 1.-1e-10;nrecp2 = 1e-10; }
		recp/=NNN;
		
	        double *tlpp = lpp - SSZ;
		if(TLP==NULL)
		{
			int i, j;

			double a = tpara.probAlpha, sss = 1.;
			double *ppp;
			int step, sst = 0, o = *(basePop+bst);

			if(emitMatrix != NULL) 
			{	if(useSoftState && tmpSpace != NULL)
				{	if(stateinteraction) 
					{	ppp = tmpSpace + ((emitK + 1) * bst) * (mylik.clustersz + 1);
						sst = mylik.clustersz + 1;
					}
					else ppp = tmpSpace + bst * (mylik.clustersz + 1);
					sss = ppp[mylik.clustersz];
				}
				else	
				{	if(stateinteraction) 
					{	ppp = emitMatrix + o * mylik.clustersz;
						sst = posSZ * mylik.clustersz;
					}
					else ppp = emitMatrix + (emitK * posSZ + o) * mylik.clustersz;
					sss = 1.;
				}
				step = 1;
			}
			else mylik.getStatePrior(ppp, step, tpara.A + (double)fixstateN * 1000000., tpara.priorW);

			sss*=a;
			double *statep = new double[SSZ * (mylik.clustersz + 1)];
			float *ddp = dp;
			
			for(i = idst; i < ided; i++, ddp += mydata.L)
			{	double *statepd = &statep[(int)(*ddp)], *statepa = &statep[mylik.clustersz], *pppd = ppp + ((int)(*ddp)) * step;
				//double aa = a * tpara.indWeight[i], ssss = sss * tpara.indWeight[i];
				for(j = 0; j < SSZ; j++, statepd += mylik.clustersz + 1, statepa += mylik.clustersz + 1)
				{	double tnnn = 0;
                    			if(j < tpara.maxK) tnnn = (double)nhp[j];
					if(j < tpara0.popn0) 
					{	tnnn += (double)tpara0.nh0[tpara0.popn0 * bst + j];
					}
					if(j < tpara.maxK && tnnn > NUMPRECISION)
					{	double taaa = 0;
						for(int l = 0; l < (int)(*ipm)[j].allele.size(); l+=2)
							if((*ipm)[j].allele[l] == (unsigned short)(*ddp))
							{	taaa = (*ipm)[j].allele[l+1]; break; }
						if(j < tpara0.popn0 && ((int)(*ddp)) < tpara0.clustersz0 - 1) taaa += (double)(*(tpara0.param0 + bst * tpara0.unitsz0 + j * tpara0.clustersz0 + 1 + ((int)(*ddp))));
						*statepd = taaa + a*(*pppd);
						*statepa = tnnn + sss;
/*if(id==16 && bst==953)
{	printf("(%d: %f,%f %f,%f)\n",j, taaa, aa*(*pppd), tnnn, ssss);
	for(int ii = 0; ii < mydata.indN; ii++)
	{	printf("%d:%d ", (int)tpara.pop[ii*mydata.L+953], (int)mydata.data[ii*mydata.L+953]);
	}
	printf("\n");
}*/
					}
					else
					{	*statepd = a*(*pppd);
						*statepa = sss;
					}
					if(j < SSZ - 1) pppd += sst;
				}
			}

			double sslp[SSZ];
			for(i = 0; i < SSZ; i++) sslp[i] = 1.;
			ddp = dp;
			bool *tlapse = tpara.lapse + mydata.L * idst + bst;
			for(i = idst; i < ided; i++, ddp += mydata.L, tlapse += mydata.L)
			{	double *statepd = &statep[(int)(*ddp)], *statepa = &statep[mylik.clustersz];
				for(j = 0; j < SSZ; j++, statepd += mylik.clustersz + 1, statepa += mylik.clustersz+1)
				{	//double lapsep = 0.95 * (double)(i > 0);
					//if(j < tpara.maxK) lapsep = ((*ipm)[j].lapseN + lapsep * 10.) / ((*ipm)[j].popN + 10.);
					if(!hpass || !(*tlapse)) sslp[j] *= ((*statepd / (*statepa)));// * (1. - lapsep) + lapsep * (int)(*tlapse));
					if(i < ided - 1)
					{	(*statepd) += tpara.indWeight[i];
						(*statepa) += tpara.indWeight[i];
					}
				}
			}
			delete [] statep;

			for(j = 0; j < SSZ; j++, lpp++, tlpp++)
			{	//double tnnn = (double)(*(nhp + j));	
				//if(j < tpara0.popn0) tnnn += (double)(*(tpara0.nh0 + tpara0.popn0 * bst + j));
				double ttt;
        		    	if(j<tpara.maxK)
                    		{   //t = fast_log((Vh[j] + aaa*(double)(*(nhp+j))/*tnnn*/) * recp + exp((*tlpp) - (*(lss-1))) * nrecp) + (*(lss-1));
ttt=(Vh[j]+aaa*(double)(*(nhp+j)))*recp;

//t+=log(min(1.,*(nhp+j)+NUMPRECISION));
                    		}
                    		else
                    		{   //t = fast_log((Vh[j]) * recp + exp((*tlpp) - (*(lss-1))) * nrecp) + (*(lss-1));
ttt=(Vh[j]*recp);
                    		}
				t=fast_log(ttt)+fast_log(1.+exp((*tlpp)-(*(lss-1)))*nrecp/ttt) + (*(lss-1));
//if(id==16 && bst==953)
//{	printf("j=%d Vh=%f, aaa=%f, nh=%f, recp=%f, nrecp=%f, tlpp=%f, t=%f, sslp=%f\n", j, Vh[j], aaa, nhp[j], recp, nrecp, *tlpp, t, sslp[j]);
//}
//if(id==16 && bst == 9750) printf("%d:%f(%f,%f|%f,%f)%f, ", j, t, Vh[j]+aaa*(double)(j<tpara.maxK)*(double)(nhp[j]), exp((*tlpp)-(*(lss-1))),recp,nrecp,sslp[j]);fflush(stdout);

				lmm[j] = fast_log(sslp[j]);
	        	        (*lpp) = t + lmm[j];//fast_log(sslp[j]);
				if(independentInd && j != id) { (*lpp) = MINUSINFINITE;  }
				if(j == 0 || maxlp < (*lpp)) maxlp = (*lpp);
			}
//if(id==16 && bst == 9750) 
//{	for(int ii = 0; ii < mydata.indN; ii++) printf("%d;",tpara.pop[ii*mydata.L+9750]);
//	printf("\n");fflush(stdout);
//}
		}
		else
		{
//clock_t cst, ced;
//cst=clock();
			for(j = 0; j < SSZ; j++, lpp++, tlpp++)
			{	
				f = TLP[bst * SSZ + j];//_getLp(id, breaks, bst, false, ploidity, tpara, dp, dp2, &ss[j][0], &ss[j][1], ipm, nhp, aszp);
				lmm[j] = f;

  	      	        	int state = ss[j][0];
				//double tnnn = (double)(*(nhp + state));	
				//if(state < tpara0.popn0) tnnn += (double)(*(nh0pp + state));
double ttt;
                        if(state < tpara.maxK)
                        {  // t = fast_log((Vh[state] + aaa*(double)(*(nhp+state))/*tnnn*/) * recp + exp((*tlpp) - (*(lss-1))) * nrecp) + (*(lss-1));
ttt=(Vh[state]+aaa*(double)(*(nhp+state)))*recp;
//t+=log(min(1.,*(nhp+state)+NUMPRECISION));
                        }
                        else    
			{  //t = fast_log((Vh[state]) * recp + exp((*tlpp) - (*(lss-1))) * nrecp) + (*(lss-1));
ttt=Vh[state]*recp;
			}
				t=fast_log(ttt)+fast_log(1.+exp((*tlpp)-(*(lss-1)))*nrecp/ttt) + (*(lss-1));
//printf("<%d:%f>", j, f);
	        	        (*lpp) = t + f;
				if(independentInd && j != id) { (*lpp) = MINUSINFINITE;  }
				if(j == 0 || maxlp < (*lpp)) maxlp = (*lpp);
			}
//ced=clock();
//Time3+=ced-cst;
//printf("\n");
		}
		tlpp = lpp - SSZ;
//if(id==16 ||(bst==mydata.L-1||bst==0)) printf("%d:%d ",id, bst),fflush(stdout);
		for(j = 0; j < SSZ; j++, tlpp++)
		{	(*lss) += exp((*tlpp) - maxlp);
//if(id==16 ||(bst==mydata.L-1||bst==0)) printf("[%d,%f] ",j,*tlpp),fflush(stdout);
		}
//if(id==16 ||(bst==mydata.L-1||bst==0)) printf("\n"),fflush(stdout);
		(*lss) = fast_log(*lss) + maxlp;
//if(isnan(*lss))exit(0);
	}

	if(false && (bst==0 || bst == mydata.L - 1))
	{	double *tlpp = lpp-SSZ;
		printf("%d: ", id);
		for(i = 0; i < SSZ; i++, tlpp++) printf("%f(%f), ", *tlpp, fast_log(Vh[i]));
		printf("| %f\n", *lss);
	}

	if(bst == mydata.L - 1 && TLP != NULL) 
	{
		delete []TLP; TLP = NULL; 
	}

	lmm += MK;
	lss++;
	ipm++;
	nhp += tpara.maxK;
	aszp++;
	rp++;
	dp++;
}

int genomicTensor::_tracebackOne(int id, vector<int> const &breaks, int bst, int bed, TENSORPARA &tpara, int ploidity, int *&pp, int *&pp2, bool *&rr, bool *&rr2, double const *&rp, float const *&nhp, vector<vector<PARAUNIT> >::const_iterator &ipm, double *&lpp, double *&lmm, double *&lss, int SSZ, int MK, vector<vector<int> > const &ss, vector<vector<int> > const &mlist, double *Vh, double *Vh2)
{
	pp--;
	nhp -= tpara.maxK;
	ipm--;
	lpp -= SSZ;
	lmm -= MK;
	lss --;
	rr--;
	rp--; 
	return 0;
}

double genomicTensor::_getLp(int id, vector<int> const &breaks, int bid, bool flip, int ploidity, TENSORPARA const &tpara, float const *dp, float const *dp2, int const *pp, int const *pp2, vector<vector<PARAUNIT> >::const_iterator ipm, float const *nhp, int const *aszp)
{
	int i, asz = *aszp;
	int ided = mydata.indIndex[id+1], idst = mydata.indIndex[id];
	double prop[asz + 1];

	_getProp(tpara, id, bid, (*pp), asz, ipm, nhp, basePop + bid, prop);
	////////////////////////
	
	double rt = 0;//, flp[asz], *tlp = TLP + asz * bid;
	for(i = idst; i < ided; i++, dp += mydata.L)//, tlp += asz * mydata.L)
	{	
		int k = (int)(*dp);
		rt+=fast_log(prop[k]/prop[asz]);//prop[k]-prop[asz];
		if(i<ided-1)
		{	prop[k]+=tpara.indWeight[i];
			prop[asz]+=tpara.indWeight[i];
		}
	}

	return(rt);
}

double genomicTensor::_imputeOne(int id, vector<int> const &breaks, int bst, int bed, int flip, int ploidity, TENSORPARA const &tpara, float *dp, float *dp2, int const *pp, int const *pp2, vector<vector<PARAUNIT> >::const_iterator ipm, float const *nhp, int const *aszp)
{
	int i, j, asz = *aszp;
	int ided = mydata.indIndex[id+1], idst = mydata.indIndex[id];

	double prop[asz + 1];
//	if(tmpM == NULL || gID > 50)
		_getProp(tpara, id, bst, (*pp), asz, ipm, nhp, basePop + bst, prop);
//	else	_getProp(tpara, id, bst, -1, asz, ipm, nhp, basePop + bst, prop);
/*if(hpass)//encourage smoothness of state assignments across genome 
{	if(bst > 0) prop[(int)(*(dp-1))]+=tpara.probAlpha;
	if(bed < mydata.L) prop[(int)(*(dp+1))]+=tpara.probAlpha;
}*/
	////////////////////////////
	
	double a = tpara.priorW;
	double RT = 0;
	for(i = idst; i < ided; i++, dp += mydata.L)
	{	double tlp[asz], maxlp;
		int *np = NULL;
		float *dpy, *dpx = NULL;
		dpy = dataYP[i] + bst * ymsz[i];
		if(maxxmsz > 0) dpx = dataXP[i] + bst * xmsz[i];

		int o = (int)(*dp);//, oo;
		if(mylik.preLP != NULL)
		{	float *prelpp = mylik.preLP + mylik.clustersz * (i * mydata.L + bst);
			for(int ii = 0; ii < mylik.clustersz; ii++, prelpp++)
				tlp[ii] = (double)(*prelpp);
		}
		else
		{	if(o >= 0) np = mylik.computeLP_subset(dpy, dpx, i, o, a, tlp);//???will Neighbor works here?
			else mylik.computeLP(dpy, dpx, i, a, tlp);
		}
double otlp[asz];
for(j=0;j<mylik.clustersz;j++) otlp[j]=tlp[j];

		if(hpass)
		{	int step = (ided - idst) * mylik.clustersz * 2;
			double *ppointer = tpointer + (i-idst)*mylik.clustersz*2 + step * bst;
			for(j = 0; j < mylik.clustersz; j++) 
			{	ppointer[j] = MINUSINFINITE;
				ppointer[mylik.clustersz + j] = prop[j] / prop[asz]; 
			}
/*if(gID>=50 && bst==99979) 
{	for(int jj = 0; jj < mylik.clustersz; jj++)
		printf("<%d:%f> ", jj, ppointer[mylik.clustersz+jj]);
	printf("\n");fflush(stdout);
}*/
			if(o >= 0 && mylik.preLP == NULL)
			{	for(j = 0; j < np[0]; j++) ppointer[np[j + 1]] = tlp[j];
				j = o;
			}
			else
			{	int m = 0;
				for(j = 0; j < mylik.clustersz; j++) 
				{	ppointer[j] = tlp[j];	
					if(tlp[j] > tlp[m]) m = j;
				}
				j = m;
			}
		}	
		else 
		{	if(o >= 0 && mylik.preLP == NULL)
			{	if((tpara.maximum & 1) != 0)
				{	j = 0;
					tlp[0] += fast_log(prop[np[1]]);// / prop[asz]);
					for(int k = 1; k < np[0]; k++)
					{	tlp[k] += fast_log(prop[np[k+1]]);// / prop[asz]);
						if(tlp[k] > tlp[j]) j = k;
					}
					RT += tlp[j];// - tlp[oo];
					j = np[j+1];
				}
				else
				{	maxlp = tlp[0];
					for(j = 1; j < np[0]; j++)
						maxlp = max(maxlp, tlp[j]);
					for(j = 0; j < np[0]; j++)
					{	tlp[j] = exp(tlp[j]-maxlp) * prop[np[j+1]];
					}
					j = _sample(tlp, np[0], 0);
					RT += maxlp;
					j = np[j+1];
				}
			}
			else
			{	if((tpara.maximum & 1) != 0)
				{	j = 0;
					tlp[0] += fast_log(prop[0]);// / prop[asz]);
					for(int k = 1; k < asz; k++)
					{	tlp[k] += fast_log(prop[k]);// / prop[asz]);
						if(tlp[k] > tlp[j]) j = k;
					}
					RT += tlp[j];// - tlp[o];
				}
				else
				{	maxlp = tlp[0];
					for(j = 1; j < asz; j++)
					{	maxlp = max(maxlp, tlp[j]);
					}
					for(j = 0; j < asz; j++)
					{	tlp[j] = exp(tlp[j]-maxlp) * prop[j];///prop[asz];
					}
					j = _sample(tlp, asz, 0);//tpara.maximum & 1);
					RT += maxlp;//fast_log(tlp[j]/tlp[o]);//maxlp;
				}
			}
			if(j != o) *dp = (float)j;
/*if(id==16 && bst == 862 && o>=0)
{	for(int ii = 0; ii < np[0]; ii++) printf("%d:%f+%f=%f\n", np[ii+1], fast_log(prop[np[ii+1]]),otlp[ii],tlp[ii]); printf("| %d\n", j);fflush(stdout);
	for(int jj = 0; jj < 5; jj++) printf("(%f)",dataYP[id][bst*ymsz[id]+jj]);
	printf("\n");
}*/
		}
	
		if(i < ided - 1)
		{	prop[j]+=1.;
			prop[asz]+=1.;
		}	
	}

	return(RT);
}

double genomicTensor::_updateLapseData(int id, double *P, int posst, int posed, float *nlapse)
{	int i, j, k = 0, idst = mydata.indIndex[id / mydata.ploidity], ided = mydata.indIndex[id / mydata.ploidity + 1];;
	int tk, step = mylik.clustersz * (ided - idst) * 2;
	double RT = 0;
	int ln = 0;
	for(i = idst; i < ided; i++)
	{	float *dp = mydata.data + i * mydata.L + posed;
		int *popp = tpara.pop + id * mydata.L + posed;
		vector<vector<PARAUNIT> >::iterator ipm = tpara.param.begin() + posed;
		float *dpy, *dpx = NULL;
		dpy = dataYP[i] + (posed - 1) * ymsz[i];
		if(maxxmsz > 0) dpx = dataXP[i] + (posed - 1) * xmsz[i];
		double *ttpointer = P + ((posed - 1) * (ided - idst) + i - idst) * mylik.clustersz * 2;
		double *npointer = NULL;
		int o = 0;
		bool *tlapse = tpara.lapse + i * mydata.L + posed;
	//	float *od = new float[posed - posst], prevd, *odp = od;
	//	for(j = posst; j < posed; j++, dp++, odp++, tlapse++) 
	//	{	if(!(*tlapse)) prevd = *dp;
	//		*odp = prevd;
	//	}
	//	odp--;
		for(j = posed - 1; j >= posst; j--, dp--, /*odp--, */tlapse--, popp--, ipm--)
		{	if(j == posed - 1)
			{	k = _sample(ttpointer, mylik.clustersz, tpara.maximum & 1);
			}
			else
			{	double lapsep = 0.95;//check if next is a lapse
				int tg = (int)(*(dp+1));
				lapsep = nlapse[tg*2+1] / nlapse[tg*2];
				if(*popp < tpara.maxK) lapsep = ((*ipm)[*popp].lapseN + lapsep * 10.) / ((*ipm)[*popp].popN + 10.);
lapsep *= tpara.indWeight[i];
//lapsep=min(lapsep,0.1);//min(0.95,pow((double)gID / (double)tpara.burnin, 2.)));
				int l;
				if(tpara.maximum & 1)
				{	l = (int)(ttpointer[k] * lapsep > npointer[k] * (1. - lapsep));
				}
				else
				{	double un = gsl_ran_flat(gammar, 0, ttpointer[k] * lapsep + npointer[k] * (1. - lapsep));
					l = (int)(un <= ttpointer[k] * lapsep);
				}
				if(l == 1)
				{	k = o;
					tk = _sample(npointer/* + mylik.clustersz*/, mylik.clustersz, tpara.maximum & 1);
//if(tk==3)
//{	for(int ii = 0; ii < mylik.clustersz; ii++) printf("%d:%f ", ii, npointer[ii]);
//	printf(" | l=1, %d\n", j);fflush(stdout);
//}
//added 10/29/2016, otherwise too many states, because skipped "true" state will not depend on data
				}	
				else
				{	k = _sample(ttpointer, mylik.clustersz, tpara.maximum & 1);
					tk = o;
/*if(gID > 50 && k==3)
{	for(int ii = 0; ii < mylik.clustersz; ii++) printf("%d:%f(%f) ", ii, ttpointer[ii], ttpointer[ii+mylik.clustersz]);
	printf(" | l=0, %d\n", j);fflush(stdout);
for(int ii = 0; ii < mydata.totalN; ii++) printf("%d;",(int)mydata.data[ii*mydata.L+j]);printf(" posSZ=%d,%d\n", posSZ, (int)basePop[j]);
exit(0);
}*/
				}
				ln += (int)(o != tk);
				*dp = (float)tk;
				*tlapse = (l==1);
			}

		//	if(*odp != k)
		//	{	mylik.removePara((int)(*odp), dpy, dpx, i);
		//		mylik.addPara(k, dpy, dpx, i);
		//		RT += fast_log(ttpointer[k] / (ttpointer[(int)(*odp)] + 1e-100));
		//	}

			o = k;
			npointer = ttpointer;
			ttpointer -= step;
			dpy -= ymsz[i];
			if(maxxmsz) dpx -= xmsz[i];
		}
		*dp = (float)o;
		*tlapse = false;
		//delete od;
	}	
	//printf("%d:%d", idst, ln);fflush(stdout);
	return(RT);
}

void genomicTensor::_getProp(TENSORPARA const &tpara, int id, int pos, int state, int asz, vector<vector<PARAUNIT> >::const_iterator ipm, float const *nhp, int const *op, double *prop)
{
	int i;
	
	////////////////////////////
	double a = tpara.probAlpha, ss = 1.;
	double *ppp;
	int step;

	if(state<0)
	{	_getProp_soft(tpara, id, pos, asz, ipm, nhp, op, prop);
		return;
	}

	if(emitMatrix != NULL) 
	{	if(useSoftState && tmpSpace != NULL)
		{	if(stateinteraction) ppp = tmpSpace + ((emitK + 1) * pos + min(state, emitK)) * (mylik.clustersz + 1);
			else ppp = tmpSpace + pos * (mylik.clustersz + 1);// + emitK/*min(tpara.maxK, state)*/ * mylik.clustersz;
			ss = ppp[mylik.clustersz];
		}
		else	
		{	if(stateinteraction) ppp = emitMatrix + (min(state, emitK) * posSZ + (*op)) * mylik.clustersz;
			else ppp = emitMatrix + (emitK * posSZ + (*op)) * mylik.clustersz;
			ss = 1.;
		}
		step = 1;
	}
	else mylik.getStatePrior(ppp, step, tpara.A + (double)fixstateN * 1000000., tpara.priorW);

	ss*=a;
	for(i = 0; i < asz; i++, ppp += step) prop[i] = a * (*ppp);
	prop[asz] = ss;
int ns=state;
//for(int ns = 0; ns <= tpara.maxK; ns++)
{
	double tnnn = 0;
       	double ww = 1.;//(double)(ns==state) * 0.9 + (1.-0.1)/(double)tpara.maxK;
    	if(ns < tpara.maxK)
	{	tnnn = (double)(*(nhp + ns));
        	if(ns < tpara0.popn0) tnnn += (double)(*(tpara0.nh0 + pos * tpara0.popn0 + ns));
	}
    	if(ns < tpara.maxK && tnnn > NUMPRECISION)
	{	for(i = 0; i < (int)(*ipm)[ns].allele.size(); i+=2)
			prop[(int)(*ipm)[ns].allele[i]] += (*ipm)[ns].allele[i+1] * ww;
		if(ns < tpara0.popn0)
		{	unsigned short *tppp0 = tpara0.param0 + pos * tpara0.unitsz0 + ns * tpara0.clustersz0;
			for(i = 0; i < min(tpara0.clustersz0 - 1, asz); i++) 
				prop[i] += (double)tppp0[i + 1] * ww;

		}
		prop[asz] += tnnn * ww;
	}
}
}

void genomicTensor::_getProp_soft(TENSORPARA const &tpara, int id, int pos, int asz, vector<vector<PARAUNIT> >::const_iterator ipm, float const *nhp, int const *op, double *prop)
{
	int i;
	
	////////////////////////////
	double a = tpara.probAlpha, ss = 1.;
	double *ppp;
	int step;
	int MK = min(tpara.maxK + 1, tpara.maxHapK);
	double *lmm = tmpM + (mydata.L + 1) * MK + pos * MK;
	for(i = 0; i < asz + 1; i++) prop[i] = 0;
	double tprop[asz + 1];

	for(int state = 0; state < MK; state++)
	{	if(emitMatrix != NULL) 
		{	if(useSoftState && tmpSpace != NULL)
			{	if(stateinteraction) ppp = tmpSpace + ((emitK + 1) * pos + min(state, emitK)) * (mylik.clustersz + 1);
				else ppp = tmpSpace + pos * (mylik.clustersz + 1);// + emitK/*min(tpara.maxK, state)*/ * mylik.clustersz;
				ss = ppp[mylik.clustersz];
			}
			else	
			{	if(stateinteraction) ppp = emitMatrix + (min(state, emitK) * posSZ + (*op)) * mylik.clustersz;
				else ppp = emitMatrix + (emitK * posSZ + (*op)) * mylik.clustersz;
				ss = 1.;
			}
			step = 1;
		}
		else mylik.getStatePrior(ppp, step, tpara.A + (double)fixstateN * 1000000., tpara.priorW);

		ss*=a;
		double tnnn = 0;
	    if(state < tpara.maxK)
	    {   tnnn = (double)(*(nhp + state));
	        if(state < tpara0.popn0) tnnn += (double)(*(tpara0.nh0 + pos * tpara0.popn0 + state));
		}
	    if(state < tpara.maxK && tnnn > NUMPRECISION)
		{	for(i = 0; i < asz; i++, ppp += step) 
				tprop[i] = a * (*ppp);
			for(i = 0; i < (int)(*ipm)[state].allele.size(); i+=2)
				tprop[(int)(*ipm)[state].allele[i]] += (*ipm)[state].allele[i+1];
			if(state < tpara0.popn0)
			{	unsigned short *tppp0 = tpara0.param0 + pos * tpara0.unitsz0 + state * tpara0.clustersz0;
				for(i = 0; i < min(tpara0.clustersz0 - 1, asz); i++) 
					tprop[i] += (double)tppp0[i + 1];

			}
			tprop[asz] = tnnn + ss;
		}
		else
		{	
			for(i = 0; i < asz; i++, ppp += step) tprop[i] = a * (*ppp);
			tprop[asz] = ss;
		}
		
		for(i = 0; i < asz + 1; i++) prop[i] += tprop[i] / tprop[asz] * lmm[state];
	}
}

void genomicTensor::_updatePrior(vector<bool> const &init, bool oneround)
{	int j, k;
/*if(gID == 0)
{
FILE *fff = fopen("tmp.txt","w");
for(i = 0; i < mydata.L; i++)
{	for(j = 0; j < mydata.totalN; j++)
	{	fprintf(fff, "%d ", (int)mydata.data[j * mydata.L + i]);
	}
	fprintf(fff, "\n");
}
fclose(fff);
}
*/

//printf("<%f %f %f | %f,%f,%f>\n", mylik.gaussprior[0], mylik.gaussprior[mylik.gausssz], mylik.gaussprior[mylik.gausssz*2],mylik.gaussprior0[0], mylik.gaussprior0[mylik.gausssz], mylik.gaussprior0[mylik.gausssz*2]);fflush(stdout);

clock_t cst, ced;
cst=clock();
	if(oneround) 
	{	tpara.A = 1.;
	
		printf(" C:%d", mylik.clustersz);

if(splitmerge)
{
for(int i = 0; i < 10; i++)
{if(gID >= 0 && gID > tpara.burnin*1/4  && gID < tpara.burnin) 
{		_splitmergeCluster(0);
}
if(gID >= 0 && gID <= tpara.burnin *1/ 2 && gsl_ran_flat(gammar, 0, 1.)<0.2)//gID >= tpara.burnin*3/5 && gID < tpara.burnin) 
{		_splitmergeCluster(1);
}
}
}

		if(mylik.preLP==NULL && gID >= 0 && gID <= tpara.burnin) _rearrangeState(1+1+(int)(tpara.A+0.5), maxG);
		printf("->%d ", mylik.clustersz);fflush(stdout);
	}
ced=clock();
Time3+=ced-cst;

	if(oneround && gID < tpara.burnin)
	{	double b = 1. + pow(1. - min(1., (double)gID / (double)tpara.burnin), 0.5) * 4.;
		int i, j;
		vector<int> tymsz(mydata.totalN);
		for(i = 0; i < mydata.totalN; i++) tymsz[i] = ymsz[i];
		sort(tymsz.begin(), tymsz.end());
		int medianymsz = tymsz[mydata.totalN -1];/// 2];
		tpara.trueindN = 0;
		for(i = 0; i < mydata.indN; i++)
		{	tpara.idw[i] = 0;
			for(j = mydata.indIndex[i]; j < mydata.indIndex[i + 1]; j++)
			{	tpara.indWeight[j] = pow((double)ymsz[j] / (double)medianymsz, b);
				tpara.indWeight[j] = max((float)1e-5, tpara.indWeight[j]);
				tpara.indWeight[j] = min((float)1e+5, tpara.indWeight[j]);
				tpara.idw[i] += tpara.indWeight[j];
//printf("%d\t%f\t%f\n", j, tpara.indWeight[j], b);fflush(stdout);
			}
			tpara.trueindN += tpara.idw[i];
		}
		_refreshParameter();
	}

//	double omaxerr = mylik.maxerr;
//	if(gID < tpara.burnin / 2) mylik.maxerr = mylik.minerr + (double)gID / (double)tpara.burnin * 2.;
	mylik.updateLambda(tpara.priorW, true, gID>0);
//	mylik.maxerr = omaxerr;
	//update emission matrix	
	if(oneround)
	{	
/*
if(gID >= tpara.burnin)
{
float tmpy[mydata.L * (maxymsz - ymsz[0]) + 1], tmpx[mydata.L * (maxxmsz - xmsz[0]) + 1];
mylik.imputeMissingTrack(0, mydata.L, dataYP[0], dataXP[0], tmpy, tmpx, mydata.data);
float *py = &tmpy[0];
FILE *f = fopen("2.txt","a");
for(i=0;i<mydata.L;i++)
{	for(j=0;j<maxymsz-ymsz[0];j++,py++)
		fprintf(f, "%f ", *py);
	fprintf(f, "\n");
}
fclose(f);
}
*/
//if(gID < tpara.burnin) maxPosSZ = gID / 2 + 1; //this is to avoid huge memory usage at the begining iterations when data is too large
//else maxPosSZ = 100;
		if(gID >= 0 && gID < tpara.burnin * 4 / 5 && tpara.samplemaximum) tpara.maximum = 3 * (int)(gsl_ran_flat(gammar, 0, 1.) < min(0.5,(double)gID / (tpara.burnin)));	
//tpara.probAlpha = max(1, tpara.burnin - gID);
clock_t cst, ced;
cst=clock();

		int K = min(tpara.maxK, tpara.maxHapK);
		int pc2 = posSZ * mylik.clustersz, esz = (K + 1) * pc2;
		double *nemitMatrix = new double[esz];
		int n = mydata.totalN;
		int idnn = n * mydata.ploidity;
		///////////////////////
		int *op = basePop; 
		float *nhp = tpara.nh;
		for(j = 0; j < esz; j++) nemitMatrix[j] = 0;
		double *ep;

//double *logp = new double[(K+1)*posSZ*mylik.clustersz], *logpN = new double[(K+1)*posSZ];
//for(j = 0; j < (K+1)*posSZ*mylik.clustersz; j++) logp[j] = 0;
//for(j = 0; j < (K+1)*posSZ; j++) logpN[j] = 0;

		double *priorpp;
		int priorstep;
		mylik.getStatePrior(priorpp, priorstep, tpara.A + (double)fixstateN * 1000000., tpara.priorW); 

		for(j = 0; j < mydata.L; j++, /*ipm++, */op++, nhp += tpara.maxK)
		{	ep = nemitMatrix + (*op) * mylik.clustersz;

			int tpop[idnn];
			int *tp = tpara.pop + j;
			int tttg[idnn];
			float *pdata = mydata.data + j;
			int jj = 1;
			for(k = 0; k < idnn; k++, pdata += mydata.L)
			{	if(k >= mydata.indIndex[jj]) { tp += mydata.L; jj++; }
				tpop[k] = *tp;
				tttg[k] = (int)(*pdata);
			}
			double PP[(K + 1) * mylik.clustersz];
			_getMLE_P(tpop, n, (*op), emitMatrix, emitK, K, PP, tttg, priorpp, priorstep);

			double sumn = 0;
			//double sumk[mylik.clustersz], nnn = 0;
			//for(k = 0; k < mylik.clustersz; k++) sumk[k] = 0;
//double lpn=0;
			
			double *pPP = PP;
			bool tttsel[K];
			for(k = 0; k < K; k++, ep += pc2)
			{	double tnnn = (double)nhp[k];
				if(k < tpara0.popn0) tnnn += (double)tpara0.nh0[j * tpara0.popn0 + k];
				if(tnnn > NUMPRECISION)
				{	if(stateinteraction)
					{ for(int l = 0; l < mylik.clustersz; l++, pPP++)
					{	ep[l] += tnnn*(*pPP+1e-10);
					}
					}
					else pPP += mylik.clustersz;
					sumn += tnnn;
					tttsel[k] = true;
				}
				else 
				{	pPP += mylik.clustersz;
					tttsel[k] = false;
				}
			}
			pPP=PP + K * mylik.clustersz;
			for(k = 0; k < mylik.clustersz; k++, pPP++)
				ep[k] += sumn * (*pPP+1e-10);
			/*for(int l = 0; l < K; l++)
				if(tttsel[l])//nhp[l] > 0)
				{	for(k = 0; k < mylik.clustersz; k++, pPP++)
					{	ep[k] += sumn * (*pPP+1e-10);
					}
				}
				else pPP+=mylik.clustersz;
			*/
		/*
			for(k = 0; k < mylik.clustersz; k++)
			{	for(int l = 0; l < K; l++) 
				{	ep[k] += sumn*(double)(nhp[l]>0)*(PP[l*mylik.clustersz+k]+1e-10);
//logp[K*pc2+(*op)*mylik.clustersz+k] += fast_log(PP[l*mylik.clustersz+k] + 1e-100);
				}
			}	
		*/
//logpN[K*posSZ+(*op)]+=lpn;
		}

		ep = nemitMatrix;
		double totalep[mylik.clustersz], tttp[mylik.clustersz], tp[mylik.clustersz + 1];
		for(j = 0; j < mylik.clustersz; j++)
		{	totalep[j] = 0;
			for(k = 0; k < posSZ; k++) totalep[j] += ep[K*pc2+k*mylik.clustersz+j];
			totalep[j] *= tpara.probAlpha / ((double)adjustB*(double)(tpara.trueindN - 1.) + 1. + tpara.probAlpha);
		}
		_updateVh(totalep, mylik.clustersz, tpara.A, tttp, false);
		for(j = 0; j < (K+1)*posSZ; j++, ep += mylik.clustersz)
		{	for(k = 0; k < mylik.clustersz; k++) 
			{	ep[k] *= tpara.probAlpha / ((double)adjustB*(double)(tpara.trueindN - 1.) + 1. + tpara.probAlpha);
				ep[k] += 100.*tttp[k];
			}
			if(j >= K * posSZ && (int)emit0.size() > j - K * posSZ)
			{	for(k = 0; k < mylik.clustersz; k++)
					ep[k] += emit0[j-K*posSZ][k + 1] * emit0[j-K*posSZ][0];
			}
			_updateVh(ep, mylik.clustersz, tpara.A, tp, false);
			for(k = 0; k < mylik.clustersz; k++) ep[k] = tp[k];
		}
		///////////////////////	

		if(emitMatrix != NULL) delete emitMatrix;
		emitMatrix = nemitMatrix;
		emitK = K;
ced=clock();
timeA+=(unsigned int)(ced-cst);

		//update basePop
cst=clock();
		if(gID >= 0) _updateBase();
ced=clock();
timeB+=(unsigned int)(ced-cst);

		//sample state members
		if(false && gID >= tpara.burnin)
		{	int *pc;
			if(gClass == NULL) 
			{	gClassSZ = mylik.clustersz + 5;
				posClassSZ = posSZ + 5;
				unsigned int sz = gClassSZ * n * mydata.L + posClassSZ * mydata.L;
				gClass = new int[sz];
				pc = gClass;
				for(unsigned int j = 0; j < sz; j++, pc++) *pc = 0;
				maxG = min(maxG, gClassSZ);
				maxPosSZ = min(maxPosSZ, posClassSZ);
			}
			pc = gClass;
			float *pd = mydata.data;
			bool *tlapse = tpara.lapse;
			for(j = 0; j < n; j++)
			{	float  prevd = *pd;
				for(k = 0; k < mydata.L; k++, pc += gClassSZ, pd++, tlapse++)
				{	if(!(*tlapse)) prevd = *pd;	
					(*(pc + (int)prevd)) += 1;
				}
			}
			int *po = basePop;
			for(k = 0; k < mydata.L; k++, pc += posClassSZ, po++)
			{	(*(pc + *po)) = (*(pc + *po)) + 1;
			}
		}
	}

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void genomicTensor::readInput(MYDATA &mydata, char const *input, char const *fbed, char const *fcov, char const *fcovbed, bool add2)
{
	int i, j, k;
	
	if(fbed != NULL) _loadData_list(dataYP, input, fbed, mydata.indinfo, mydata.fyinfo, mydata.snpinfo, mydata.indN, mydata.L, mydata.totalN, maxymsz, ymsz, ymap, mydata.indIndex);
	else _loadData(dataYP, input, mydata.indinfo, mydata.fyinfo, mydata.snpinfo, mydata.indN, mydata.L, mydata.totalN, maxymsz, ymsz, ymap, mydata.indIndex);
	
	tpara.indWeight = new float[mydata.totalN + mydata.indN];
	tpara.idw = tpara.indWeight + mydata.totalN;
	tpara.trueindN = 0;
	for(i = 0; i < mydata.indN; i++)
	{	tpara.idw[i] = 0;
		for(j = mydata.indIndex[i]; j < mydata.indIndex[i + 1]; j++)
		{	tpara.indWeight[j] = pow((double)ymsz[j] / (double)maxymsz, 2.);
			tpara.indWeight[j] = max((float)1e-5, tpara.indWeight[j]);
			tpara.idw[i] += tpara.indWeight[j];
		}
		tpara.trueindN += tpara.idw[i];
	}

	maxxmsz = 0; 
	if(fcov != NULL)
	{	vector<string> tmpindinfo;
		vector<SNPINFO> tmpsnpinfo;
		int N, L, tN, *tmpindIndex = NULL;
		if(fcovbed != NULL) _loadData_list(dataXP, fcov, fcovbed, tmpindinfo, mydata.fxinfo, tmpsnpinfo, N, L, tN, maxxmsz, xmsz, xmap, tmpindIndex);
		else _loadData(dataXP, fcov, tmpindinfo, mydata.fxinfo, tmpsnpinfo, N, L, tN, maxxmsz, xmsz, xmap, tmpindIndex);
		if(add2 && maxxmsz > 0)
		{	for(i = 0; i < mydata.totalN; i++)
			{	float *tdata = new float[L * xmsz[i] * 2];
				for(j = 0; j < xmsz[i]; j++)
				{	double m = 0;
					for(k = 0; k < L; k++) m += dataXP[i][k * xmsz[i] + j];
					m /= (double)L;
					for(k = 0; k < L; k++) 
					{	tdata[k * xmsz[i] * 2 + j * 2] = dataXP[i][k * xmsz[i] + j];
						tdata[k * xmsz[i] * 2 + j * 2 + 1] = pow((double)dataXP[i][k * xmsz[i] + j] - m, 2.);
					}
				}
				delete dataXP[i];
				dataXP[i] = tdata;
				int *tmap = new int[xmsz[i] * 2];
				for(j = 0; j < xmsz[i]; j++) 
				{	tmap[j * 2] = xmap[i][j] * 2;
					tmap[j * 2 + 1] = xmap[i][j] * 2 + 1;
				}
				delete xmap[i];
				xmap[i] = tmap;
				xmsz[i] *= 2;
			}
		}
		if(tmpindIndex != NULL) tmpindIndex = NULL;
	}
	else
	{	xmsz = new int[mydata.totalN];
		xmap = new int*[mydata.totalN];
		for(i = 0; i < mydata.totalN; i++)
		{	xmsz[i] = 0;
			xmap[i] = NULL;
		} 
		dataXP = new float*[mydata.totalN];
		for(i = 0; i < mydata.totalN; i++) dataXP[i] = NULL;
	}

	mydata.data = new float[mydata.totalN * mydata.L];
	float *dp = mydata.data;
	for(i = 0; i < mydata.totalN * mydata.L; i++) *dp++ = -1.;
	mydata.indN *= mydata.ploidity;//requires preset of ploidity

}

void genomicTensor::outputResult(char const *fname)
{
	int i, j, k;
	char str[200];
	int oclustersz = mylik.clustersz;
	mylik.clustersz = 0;
	sprintf(str, "%s.state", fname);
	FILE *f = fopen(str, "w");
	fprintf(f, "#ID CHR POSst POSed ");
	for(i = 0; i < mydata.totalN; i++) fprintf(f, "%s ", mydata.indinfo[i].c_str());
	fprintf(f, "PosClass\n");
	
	vector<int> baseN(posSZ, 0);
	for(i = 0; i < mydata.L; i++)
	{	fprintf(f, "%s %s %d %d ", mydata.snpinfo[i].snpid.c_str(), mydata.snpinfo[i].chr.c_str(), mydata.snpinfo[i].posst, mydata.snpinfo[i].posed);
		for(j = 0; j < mydata.totalN; j++) 
		{	//fprintf(f, "%d ", (int)mydata.data[j * mydata.L + i]);
			int mm = mydata.data[j * mydata.L + i];
			if(gClass != NULL)
			{	int *gp = gClass + (j * mydata.L + i) * gClassSZ;
			 	for(k = 1; k < gClassSZ; k++) if(gp[mm] < gp[k]) mm = k;
				//fprintf(f, "%d(%d)%d ", mm, gp[mm], (int)mydata.data[j*mydata.L+i]);
			}
			fprintf(f, "%d ", mm);
			mydata.data[j * mydata.L + i] = (float)mm;
			mylik.clustersz = max(mylik.clustersz, mm);
		}
		//fprintf(f, "%d\n", basePop[i]);
		int mm = (int)basePop[i];
		if(gClass != NULL)
		{	int *gp = gClass + mydata.totalN * gClassSZ * mydata.L + i * posClassSZ;
			for(k = 1; k < posClassSZ; k++) if(gp[mm] < gp[k]) mm = k;
			//fprintf(f, "%d(%d)%d\n", mm, gp[mm], basePop[i]);
		}
		fprintf(f, "%d\n", mm);
		baseN[mm]++;
	}
	fclose(f);
	mylik.clustersz++;

	if(outputproportion && gClass != NULL)
	{	sprintf(str, "%s.stateprop", fname);
		f = fopen(str, "w");
		for(i = 0; i < mydata.L; i++)
		{	fprintf(f, "%d ", i);
			for(j = 0; j < mydata.totalN; j++)
			{	int *gp = gClass + (j * mydata.L + i) * gClassSZ;
				for(k = 0; k < gClassSZ; k++) fprintf(f, "%d ", gp[k]);
			}
			fprintf(f, "\n");
		}
		fclose(f);
	}

	_refreshParameter();
	sprintf(str, "%s.para", fname);
	mylik.outputParameter(str, tpara.priorW, mydata.fyinfo, mydata.fxinfo);

	sprintf(str, "%s.profile", fname);
	f=fopen(str, "w");
	for(i = 0; i < posSZ; i++)
	{	fprintf(f, "%d ", baseN[i]);
		for(j = 0; j < mylik.clustersz; j++)
		{	fprintf(f, "%f ", emitMatrix[(emitK * posSZ + i) * oclustersz + j]);
		}
		fprintf(f, "\n");
	}
	fclose(f);


/*
for(int id = 0; id < mydata.totalN; id++)
{	if(maxymsz - ymsz[id] == 0) continue;
	for(i = 0; i < (int)imputelist.size(); i++)
	{	string myid = mydata.indinfo[id];
		for(j = myid.size() - 1; j >=0; j--) if(myid[j] == '.') break;
		if(j >= 0) myid = myid.substr(0, j);
		if(myid == imputelist[i]) break;
	}
	if(i >= (int)imputelist.size() && i > 0) continue;

	float *tmpdata = new float[mydata.L * (maxymsz + maxxmsz - ymsz[id] - xmsz[id]) + 1], *tmpy = tmpdata, *tmpx = tmpy + mydata.L * (maxymsz - ymsz[id]);
	mylik.imputeMissingTrack(id, 0, mydata.L, dataYP[id], dataXP[id], tmpy, tmpx, mydata.data + id * mydata.L);//tlp);
//	delete tlp;

	char str[100];
	sprintf(str, "%s.impute.%s", fname, mydata.indinfo[id].c_str());
	f = fopen(str,"w");
	for(i = 0; i < maxymsz; i++)
	{	for(j = 0; j < ymsz[id]; j++) if(i == ymap[id][j]) break;
		if(j >= ymsz[id]) fprintf(f, "%s ", mydata.fyinfo[i].c_str());
	}
	for(i = 0; i < maxxmsz; i++)
	{	for(j = 0; j < xmsz[id]; j++) if(i == xmap[id][j]) break;
		if(j >= xmsz[id]) fprintf(f, "%s ", mydata.fxinfo[i].c_str());
	}
	fprintf(f, "\n");
	float *py = tmpy, *px = tmpx;
	for(i = 0; i < mydata.L; i++)
	{	for(j = 0; j < maxymsz - ymsz[id]; j++, py++) fprintf(f, "%5.3f ", *py);
		for(j = 0; j < maxxmsz - xmsz[id]; j++, px++) fprintf(f, "%5.3f ", *px);
		fprintf(f, "\n");
	}
	fclose(f);
	delete []tmpdata;
	if(gzip)
	{	sprintf(str, "gzip -f %s.impute.%s", fname, mydata.indinfo[id].c_str());
		system(str);
	}
}
	
	mylik.imputeMissing2(mydata, dataYP, dataXP, tpara.priorW, fname);
*/	_imputeMissingData(fname, -1);

	double *varstatus = NULL;
/*
	if(varstatus == NULL) varstatus = new double[mydata.L * 3];
	vector<int> response(mydata.indN);
	for(i = 0; i < mydata.indN; i++) response[i] = i;
	mylik.updateLambda(tpara.priorW * 1e-3);
	_testVariability(response, varstatus);
*/
	if(varstatus != NULL)
	{	sprintf(str, "%s.var", fname);
		f = fopen(str, "w");
		for(i = 0; i < mydata.L; i++)
		{	fprintf(f, "%f %f %f\n", varstatus[i], varstatus[i + mydata.L], varstatus[i + 2 * mydata.L]);	
		}
		fclose(f);
	}
}

double genomicTensor::_lgammaM(int q, double x)
{
	int i;
	double rt = (double)q * (double)(q - 1) / 4. * fast_log(PI);
	/*for(i = 0; i < q; i++)
	{	rt += gsl_sf_lngamma(x);
		x -= 0.5;
	}*/
	double s = gsl_sf_lngamma(x);
	for(i = 0; i < q; i += 2)
	{	rt += s;
		if(i + 2 < q) s -= fast_log(x - (double)i / 2. - 1.);
	}
	s = gsl_sf_lngamma(x - 0.5);
	for(i = 1; i < q; i += 2)
	{	rt += s;
		if(i + 2 < q) s -= fast_log(x - (double)i / 2. - 1.);
	}
	return rt;
}

void genomicTensor::_updateBase()
{
	int pz1 = posSZ + 1;
        int i, j, k, p2 = pz1 * pz1;
        double mn[posSZ], mp[pz1], tn[p2], trans[p2], a = (tpara.probAlpha);//changed 05/26/2015 //1.;//tpara.probAlpha;//changed 02/07/2015
	double posRate[pz1], rr[posSZ], rn[posSZ];

        for(i = 0; i < posSZ; i++) mn[i] = rr[i] = rn[i] = 0;
        for(i = 0; i < p2; i++) tn[i] = 0;//((double)((int)(i/posSZ)==(i%posSZ)) + 1./(double)posSZ) / 2.;
        int *op = basePop;
	bool *rp = baseR;
	mn[*op]++; op++; rp++;
        for(i = 1; i < mydata.L; i++, op++, rp++)
       	{	mn[*op]++;
		if(*rp)
		{	tn[(*(op-1)) * pz1 + (*op)]++;
       		}
		if(*op == *(op-1) && !(*rp))
		{	rr[*(op-1)]++;
		}
		rn[*(op-1)]++;
	}
	_updateVh(mn, posSZ, tpara.A, mp, false);
	if(posSZ >= maxPosSZ) mp[posSZ] = 0;
	for(i = 0; i < posSZ; i++) 
	{	posRate[i] = 1. - (rr[i]+tpara.defmut*tpara.priorW * 1.) / ((double)rn[i] + tpara.priorW * 1.);
	}
	posRate[posSZ] = tpara.defmut;

	double *np = tn, *tp = trans;
	for(i = 0; i < pz1; i++, np += pz1, tp += pz1)
	{	_updateVh(np, posSZ, tpara.A, tp);
		if(posSZ >= maxPosSZ) tp[posSZ] = 0; 
	}

	int nn = mydata.totalN + 1;
	for(i = 0; i < mydata.totalN; i++) nn += (int)(max(1., (double)tpara.indWeight[i]) + 0.999999);
	//printf("nn=%d\n",nn);fflush(stdout);

	int pc2 = posSZ * mylik.clustersz;
	double *A = new double[(emitK + 1) * pc2 * nn], *app = A;
	double *epp = emitMatrix;
	double *priorpp;
	int priorstep;
	mylik.getStatePrior(priorpp, priorstep, tpara.A, tpara.priorW); 
	for(i = 0; i < emitK + 1; i++)
	{	for(j = 0; j < pc2; j++, epp++)
		{	double f = a * priorpp[(j % mylik.clustersz) * priorstep];
			if(emitMatrix != NULL) f = a * (*epp);
			(*app) = 0;//gsl_sf_lngamma(f);
			app++;
			for(k = 1; k < nn; k++, app++)
			{	(*app) = *(app-1) + fast_log((double)k - 1. + f);
			}
		}
	}
	double B[nn];
	B[0] = 0;//gsl_sf_lngamma(a);
	for(k = 1; k < nn; k++) B[k] = B[k-1] + fast_log((double)k - 1. + a);

	double A0[mylik.clustersz * nn];
	for(i = 0; i < mylik.clustersz; i++) 
	{	double f = a * priorpp[i * priorstep];
		A0[i * nn] = 0;//gsl_sf_lngamma(f);
		for(j = 1; j < nn; j++)
			A0[i * nn + j] = A0[i * nn + j - 1] + fast_log((double)j - 1. + f);
	}

	int nstep = max(1000,(int)ceil((double)mydata.L / (double)nthread));
	vector<int> region(1,0);
	j = 0;
	for(i = 1; i < mydata.L; i++) if(mydata.snpinfo[i].chr != mydata.snpinfo[j].chr) { region.push_back(i); j = i; }
	region.push_back(i);
	vector<int> bkregion = region;
	region.clear();
	region.push_back(0);
	for(i = 1; i < (int)bkregion.size(); i++)
	{	int tnsz = max(1, (int)ceil((double)(bkregion[i] - bkregion[i - 1]) / nstep));
		if(tnsz > 1)
		{	int tnstep = (int)ceil((double)(bkregion[i] - bkregion[i - 1]) / (double)tnsz);
			int ok = bkregion[i - 1];
			for(j = 1; j < tnsz; j++)
			{	k = bkregion[i - 1] + tnstep * j;
				int lk, rk;
				for(lk = k; lk > max(ok, k - tnstep / 2); lk--)
					if(basePop[lk] != basePop[lk - 1]) break;
				for(rk = k + 1; rk < min(bkregion[i], k + tnstep / 2); rk++)
					if(basePop[rk] != basePop[rk - 1]) break;
				if(k - lk < rk - k && lk > max(ok, k - tnstep / 2)) k = lk;
				else if(rk - k < k - lk && rk < min(bkregion[i], k + tnstep / 2)) k = rk;
				if(k > region[(int)region.size() - 1]) 
				{	region.push_back(k);
					ok = k;
				}
			}
		}
		region.push_back(bkregion[i]);
	}
	//for(i = 0; i < (int)region.size(); i++) printf("%d: %d\n", i, region[i]); fflush(stdout);

	int tspstep = 1;
	if(stateinteraction) tspstep = emitK + 1;
	if(useSoftState)
	{	if(tmpSpace != NULL) delete tmpSpace;
		tmpSpace = new double[tspstep * mydata.L * (mylik.clustersz + 1)];
	}
#pragma omp parallel num_threads(nthread)
{

# pragma omp for

	for(int r = 1; r < (int)region.size(); r++)
	{	int tL = region[r] - region[r - 1];
		int i, j, k, l;
	
		double *E = new double[2*pz1 * (region[r] - region[r - 1] + 1)], *M = E + pz1 * (region[r] - region[r - 1] + 1);
      		double *P = new double[(pz1 + 1) * (region[r] - region[r - 1] + 1)];
        	vector<vector<PARAUNIT> >::iterator ipm = tpara.param.begin() + region[r - 1];
        	float *nhp = tpara.nh + tpara.maxK * region[r - 1];
		int *aszp = mydata.asz + region[r - 1];

	        //forward
clock_t cst, ced;
cst=clock();
		double *pE = E, *maxP = P + pz1 * (region[r] - region[r - 1] + 1); 
	        for(i = 0; i < tL; i++, P += pz1, ipm++, nhp += tpara.maxK, maxP++, aszp++)
        	{       double prob[pz1];
	                if(i==0)
        	        {       for(j = 0; j < pz1; j++) prob[j] = (mp[j]);
                	}
	                else
	                {       double *tP = P - pz1;
				double mmm = *(maxP-1);
                	        for(j = 0; j < pz1; j++) prob[j] = 0;
	                        double *tr = trans;
        	                for(j = 0; j < pz1; j++, tP++)
                	        {       double pr = posRate[j], npr = 1. - posRate[j], ttp = exp(*tP-mmm), ttpr = ttp * pr;
					for(k = 0; k < pz1; k++, tr++)
                                	        prob[k] += ttpr * (*tr);
	                        	prob[j] += ttp * npr;
				}
                	}
			for(j = 0; j < pz1; j++) prob[j] = fast_log(prob[j]);

	                int asz = *aszp, K = emitK, jz=0, jzstep = mylik.clustersz * nn, Kz = K*pc2*nn;
        	        double emit[pz1], maxe = MINUSINFINITE;
			int numbers[K * asz], kid[K * asz], astep[K * asz], kn = 0, kzstep = pc2 * nn, k0 = 0, kh[K];

			for(k = 0; k < K; k++)
			{	double tnnn = (double)nhp[k];
				if(k < tpara0.popn0) tnnn += (double)tpara0.nh0[i * tpara0.popn0 + k];
				if(tnnn > NUMPRECISION)
				{	float taaa[asz];
					for(l = 0; l < asz; l++) taaa[l] = 0;
					for(l = 0; l < (int)(*ipm)[k].allele.size(); l+=2)
					{	taaa[(int)(*ipm)[k].allele[l]] += (*ipm)[k].allele[l + 1];	
					}
					if(k < tpara0.popn0)
					{	
						unsigned short *tppp0 = tpara0.param0 + i * tpara0.unitsz0 + k * tpara0.clustersz0;
						for(l = 0; l < min(tpara0.clustersz0 - 1, asz); l++)
							taaa[l] += (float)tppp0[l + 1];	
					}
					for(l = 0; l < asz; l++)
						if(taaa[l] > NUMPRECISION) { numbers[kn] = (int)taaa[l]; kid[kn] = k; astep[kn] = l * nn; kn++; }
					kh[k0++] = tnnn;
				}
			}

			for(j = 0; j < pz1; j++, jz+=jzstep)
			{	double eee = 0, *aaa;
				aaa = A + Kz + jz;
				if(j == posSZ) aaa = A0;
				for(k = 0; k < kn; k++)
				{	
					if(stateinteraction && (k == 0 || kid[k] != kid[k-1])) aaa = A + kid[k] * kzstep + jz;
					eee += aaa[astep[k] + numbers[k]];
				}
				for(k = 0; k < k0; k++)	eee -= B[kh[k]];
				emit[j] = eee;
       	         		if(j == 0  || maxe < eee) maxe = eee;
			}

			for(j = 0; j < pz1; j++) 
			{	P[j] = prob[j] + emit[j];//exp(emit[j] - maxe);
				*pE++ = emit[j];
				if(j == 0 || *maxP < P[j]) *maxP = P[j];
			}
        	}

ced=clock();Time1+=ced-cst;cst=ced;

		double *tmpP = NULL, emitW[posSZ + 1], gW[mylik.clustersz];
		int tspstep = 1;
		if(useSoftState)
		{	tmpP = tmpSpace + tspstep * (region[r] - 1) * (mylik.clustersz + 1);

			for(i = 0; i < mylik.clustersz; i++)
			{	gW[i] = 0;
			}
			for(i = 0; i < posSZ; i++)
			{	emitW[i] = 0;
				for(j = 0; j < mylik.clustersz; j++) emitW[i] += emitMatrix[emitK * pc2 + i * mylik.clustersz + j] * exp(gW[j] / (double)maxymsz);
			}
			emitW[posSZ] = 1;
		}

	        //backward
		double *pM = M + pz1 * (tL - 1);
		int oo = -1;
		int *op = basePop + region[r] - 1;
		bool *rp = baseR + region[r] - 1;
		for(i = tL - 1; i >= 0; i--, op--, rp--)
		{	P = P - pz1;
			maxP--;
			pE -= pz1;
			double ss = 0;
			double prob[pz1], ttt[2];
			if(i == tL - 1)
			{	for(j = 0; j < pz1; j++) prob[j] = exp(P[j] - *maxP);
				if(useSoftState) { for(j = 0; j < pz1; j++) { pM[j] = P[j], P[j] = pE[j]; /*ss+= pM[j];*/ } }
				ss = *maxP;
			}
			else
			{	double ms;
				ms = ttt[0] = P[oo] + fast_log(1. - posRate[oo]);
				int jj = oo;
				for(j = 0; j < pz1; j++, jj += pz1) 
				{	prob[j] = P[j] + fast_log(posRate[j] * trans[jj]);
					if(ms < prob[j]) ms = prob[j];
				}
				for(j = 0; j < pz1; j++) prob[j] = exp(prob[j] - ms);
				ttt[0] = exp(ttt[0] - ms);
				prob[oo] += ttt[0];
				ttt[1] = prob[oo];
				if(useSoftState)
				{
					int jj=0;
					for(j = 0; j < pz1; j++,jj+=pz1)
					{	double s = 0, ms = 0, f[pz1], pr = posRate[j], npr = 1. - posRate[j];
						for(k = 0; k < pz1; k++)
						{	f[k] = fast_log(((double)(j==k)*npr+pr*trans[jj + k])) + P[pz1+k];
							if(k == 0 || ms < f[k]) ms = f[k];
						}
						for(k = 0; k < pz1; k++)
						{	s += exp(f[k] - ms);
						}
						if(s > 0 && s < 100000);
						else
						{	printf("i=%d: ", i);
							for(int k = 0; k < pz1; k++) printf("%f(%f),",P[pz1+k], fast_log(trans[j*pz1+k]));printf(" max=%f\n", *(maxP+1));
							exit(0);
						}
						s = fast_log(s);
						pM[j] = P[j] + s;

						P[j] = pE[j] + s;
						if(j == 0 || ss < pM[j]) ss = pM[j];
					}
				}
			}
			if(useSoftState)
			{
				double *tsp = tmpP;
				int kkkk=0;
				for(int kkk = emitK + 1 - tspstep; kkk < emitK + 1; kkk++,kkkk+=pc2)
				{	double sumtsp = 0, *otsp = tsp, sw = 0, swn = 0, f;
					for(j = 0; j < mylik.clustersz; j++, tsp++)
					{	*tsp = 0;
						int tkk=j;
						for(int k = 0; k < pz1; k++,tkk+=mylik.clustersz)
						{	f = exp(pM[k] - ss);
							if(k < posSZ) (*tsp) += f * emitMatrix[kkkk+tkk]; 
							else (*tsp) += f * priorpp[j * priorstep];
							sw += f * emitW[k];
							swn += f;
						}
						sumtsp += *tsp;
						if(*tsp > 0 && *tsp < 100000);
						else
						{	for(int k = 0; k < pz1; k++) printf("%f(%f), ", pM[k], fast_log(emitMatrix[kkk*pc2+k*mylik.clustersz+j])); printf("| %f, tsp=%f\n", ss, *tsp);
							for(int k = 0; k < pz1; k++) printf("%d:%f, ", k, pM[k]);
							exit(0);
						}
					}
					sw = max(0.1, 0. + sw / swn);
					sumtsp /= sw;
					for(j = 0; j < mylik.clustersz; j++) otsp[j] /= sumtsp;
					*tsp++ = sw;
				}
				tmpP -= tspstep * (mylik.clustersz + 1);
			}

			int ooo = oo;
			oo = (*op) = _sample(prob, pz1, tpara.maximum & 2);
			(*rp) = true;
			if(oo == ooo)
			{	if(gsl_ran_flat(gammar, 0, ttt[1]) <= ttt[0]) (*rp) = false;
			}
			if(oo == posSZ)
			{	if(!(*rp) && i < tL - 1) (*op) = *(op+1);
				else if((tpara.maximum & 1) == 0)
                			(*op) = min(maxPosSZ - 1, posSZ + (int)(- fast_log(gsl_ran_flat(gammar, 0, 1.)) / fast_log((1. + tpara.A) / tpara.A)));
			}
			pM -= pz1;
		}
		
		delete P;
		delete E;
Time2+=clock()-cst;
	}
}

	////collapse basePop, update posSZ;
	vector<bool> on;
	op = basePop;
	for(i = 0; i < mydata.L; i++, op++) 
	{	if(*op >= (int)on.size()) on.resize((int)(*op)+1,false);	
		on[*op] = true;
	}
	vector<int> omap(max(posSZ, (int)on.size()), -1);
	int newpossz = 0;
	for(i = 0; i < (int)on.size(); i++) 
		if(on[i]) omap[i] = newpossz++;
	op = basePop;
	for(i = 0; i < mydata.L; i++, op++) *op = omap[*op];
	double *nemitMatrix = new double[(emitK + 1) * newpossz * mylik.clustersz];
	double *ep = emitMatrix, *nep = nemitMatrix;
	for(i = 0; i < emitK + 1; i++)
	{	int l = 0;
		for(j = 0; j < posSZ; j++, ep += mylik.clustersz)
			if(j < (int)on.size() && on[j])
			{	for(k = 0; k < mylik.clustersz; k++, nep++) *nep = ep[k];
				l++;
			}
		for(j = l; j < newpossz; j++)
			for(k = 0; k < mylik.clustersz; k++, nep++) *nep = priorpp[k * priorstep];
	}
	int oz = posSZ;
	posSZ = newpossz;
	printf("O=%d(%d) ", posSZ, tpara.maximum & 1); 
	delete emitMatrix;
	emitMatrix = nemitMatrix;
	
	delete []A;
}

void genomicTensor::_getMLE_P(int *tpop, int N, int o, double *E, int oK, int K, double *PP, int *tttg, double *priorpp, int priorstep)
{	
	int i, j;
	int pc2 = posSZ * mylik.clustersz;

//	double *prop[K + 1], *ip = PP;
//	int step[K + 1];
	double *ip = PP, *prop = priorpp;
	int step = priorstep;
	if(oK != 0) { prop = E + oK * pc2 + o * mylik.clustersz; step = 1; }
	double nn[K + 1];
	i = 0;
	if(!stateinteraction) { i = K; ip = PP + K * mylik.clustersz; }
	for(i = i; i < K + 1; i++)
	{	//step[i] = 1;
		//if(oK == 0) { prop[i] = priorpp; step[i] = priorstep; }
		//else prop[i] = E + oK * pc2 + o * mylik.clustersz;
		for(j = 0; j < mylik.clustersz; j++, ip++) 
		{	*ip = 1e-10/(double)mylik.clustersz+(double)(gID > tpara.burnin / 2) * sqrt(tpara.probAlpha) * (*(prop + j * step));//changed 05/26/2015 
			//*ip = 1e-10/(double)mylik.clustersz+tpara.probAlpha * (*(prop+j*step));
								//(double)(gID*0 > tpara.burnin / 2) * tpara.probAlpha * (*(prop[i] + j * step[i]));//changed 02/07/2015
//if(gID > tpara.burnin/2)			*ip += 100. * (double)mydata.totalN * tppp[j];
		}
		nn[i] = //tpara.probAlpha;
			1e-10 + (double)(gID > tpara.burnin / 2) * sqrt(tpara.probAlpha);//changed 05/26/2015
				//(double)(gID*0 > tpara.burnin / 2) * tpara.probAlpha;//changed 02/07/2015
//if(gID > tpara.burnin/2)		nn[i] += 100. * (double)mydata.totalN;
	}
	
	ip = PP + K * mylik.clustersz;
	for(i = 0; i < N; i++)
	{	int ddd = tttg[i];//not compatible with ploidity
		if(stateinteraction)
		{	 *(PP+tpop[i]*mylik.clustersz+ddd) = *(PP+tpop[i]*mylik.clustersz+ddd) + tpara.indWeight[i];
                	nn[tpop[i]]+= tpara.indWeight[i];
		}
		ip[ddd]+=tpara.indWeight[i];
		nn[K]+=tpara.indWeight[i];
	}
	
	i=0;
	ip = PP;
	if(!stateinteraction) { i = K; ip = PP + K * mylik.clustersz; }
	for(i = i; i < K + 1; i++)
	{	for(j = 0; j < mylik.clustersz; j++, ip++) (*ip) /= nn[i];
	}
}

void genomicTensor::_getMLE_HyperP(double *logp, double *ak, int sz)
{
	int i;
	double oak[sz];
	for(i = 0; i < sz; i++) oak[i] = ak[i] * tpara.probAlpha;

	double psiC = gsl_sf_psi(tpara.probAlpha);
	double tsiC = gsl_sf_psi_1(tpara.probAlpha);
	bool a_inner;
	int cn = 0;
	do
	{	a_inner = true;
		double dF[sz], q[sz], b = 0, iq = 0;
		for(i = 0; i < sz; i++)
		{	dF[i] = (psiC - gsl_sf_psi(oak[i]) + logp[i]);
			q[i] = -gsl_sf_psi_1(oak[i]);
			b += dF[i] / q[i];
			iq += 1./q[i];
		}
		double s = 0;
		for(i = 0; i < sz; i++)
		{	ak[i] = max(0., oak[i] - (dF[i] - b / (1./tsiC+iq)) / q[i]);
			s += ak[i];
		}
		double e2 = 0;
		for(i = 0; i < sz; i++)
		{	ak[i] /= s / tpara.probAlpha;
			e2 += (oak[i]-ak[i])*(oak[i]-ak[i])/oak[i];
			oak[i] = ak[i];
		}
		cn++;
		if(e2 < 0.01 || cn >= 10) a_inner = false;
	} while(a_inner);
}

void genomicTensor::_rearrangeState(int add, int maxG)
{	int i, j, k;
	vector<char> sel(mylik.clustersz, 0);
//add=3;//added on Jul 11th

	int n = mydata.totalN * mydata.ploidity;
	float *dp = mydata.data;	
	for(i = 0; i < n; i++)
	{	for(j = 0; j < mydata.L; j++, dp++)
		{	
			sel[(int)(*dp)] = 1;
		}
	}
	double cn[mylik.clustersz];
	mylik.getStateCount(cn);
	for(i = 0; i < mylik.clustersz; i++) sel[i] |= (int)(cn[i] > 0);
	for(i = 0; i < tpara0.clustersz0 - 1; i++) sel[i] = 2;
	
	for(i = 1; i < (int)mylik.modelparameter.size(); i++)
		if((int)mylik.modelparameter[i].size() > 0) sel[i-1] = 2;

	int map[mylik.clustersz], remap[mylik.clustersz], remain[mylik.clustersz], newclustersz = 0;
	j = k = 0;
	for(i = 0; i < mylik.clustersz; i++) 
	{	if(sel[i] == 2) { map[i] = i; remap[i] = i; newclustersz = max(newclustersz, i + 1); }
		else if(sel[i] == 1) { while(sel[j] == 2) j++;  map[i] = j; remap[j] = i; j++; newclustersz = max(newclustersz, j); }
		else { map[i] = -1; remain[k++] = i; }
	}
	k = 0;
	for(i = j; i < mylik.clustersz; i++) { if(sel[i] != 2) remap[i] = remain[k++]; }
	int trueclustersz = max(newclustersz, min(newclustersz + add, maxG));
	
	if(newclustersz < mylik.clustersz || trueclustersz != mylik.clustersz)
	{	

//clock_t cst, ced;
//cst=clock();
		vector<vector<PARAUNIT> >::iterator ip = tpara.param.begin();
		for(i = 0; i < mydata.L; i++, ip++)//, inp++)
		{	
			vector<PARAUNIT>::iterator ipp = (*ip).begin();
			for(j = 0; j < tpara.maxK; j++, ipp++)//, inpp++)
			{	int kk = 0;
				for(k = 0; k < (int)(*ipp).allele.size(); k+=2)
				{	if(map[(int)(*ipp).allele[k]] >= 0 && (*ipp).allele[k + 1] > NUMPRECISION)
					{	(*ipp).allele[kk] = map[(int)(*ipp).allele[k]];
						(*ipp).allele[kk+1] = (*ipp).allele[k + 1];
						kk += 2;
					}
				}
				if(kk < k) (*ipp).allele.resize(kk);
			}
			mydata.asz[i] = trueclustersz;
		}
//ced=clock();
//Time3+=ced-cst;

		dp = mydata.data;
		for(i = 0; i < n; i++)
			for(j = 0; j < mydata.L; j++, dp++)
			{	if((int)(*dp) >= mylik.clustersz) printf("!!!%d,%d,%d\n",(int)(*dp), i, j),exit(0);
				(*dp) = (float)map[(int)(*dp)];
			}
	
		int oclustersz = mylik.clustersz;
		mylik.rearrangeParameter(remap, trueclustersz, dataYP, dataXP, mydata.totalN, mydata.L, tpara.priorW);
	//	_getAsz(mydata);

		if(emitMatrix != NULL)
		{	n = (emitK + 1) * posSZ;
			double *ep = emitMatrix;
			double *newemitMatrix = new double[n * trueclustersz], *nep = newemitMatrix;
			double tppp = 1. / (1. + tpara.A);;
			for(i = 0; i < n; i++, ep += oclustersz, nep += trueclustersz)
			{	double msum = 0;
				for(j = 0; j < trueclustersz; j++) 
				{	if(j < oclustersz) { nep[j] = ep[remap[j]]; msum += nep[j]; }
					else nep[j] = max(1e-100,(1. - msum)) * (1. - tppp) * pow(tppp, (double)(j - oclustersz));
				}
			}
			delete emitMatrix;
			emitMatrix = newemitMatrix;
		}
	}
}

void genomicTensor::_loadData(float **&datamatrix, char const *input, vector<string> &sid, vector<string> &fid, vector<SNPINFO> &snpinfo, int &indN, int &L, int &totalN, int &maxmsz, int *&msz, int **&map, int *&indIndex)
{
	ifstream f(input);
	string tmp;
	int i, j, k, n, l;
	bool oldformat = false;
	
	if(f.fail()) printf("Cannot open data file %s\n", input), exit(0);
	for(L = 0; getline(f, tmp); L++) ;
	L--;
	f.close();

	f.open(input);
	getline(f, tmp);
	l = (int)tmp.size(); l--; while(l>0 && (int)tmp[l] < 32) l--; l++;
	while(l > 0 && (tmp[l-1] == ' ' || tmp[l-1] == '\t')) l--;
	for(i = 0; i < l; i++) if(tmp[i] == ' ' || tmp[i] == '\t') break;
	for(i = i + 1; i < l; i++) if(tmp[i] == ' ' || tmp[i] == '\t') break;
	if((tmp[i + 1] == 'P' || tmp[i + 1] == 'p') && (tmp[i + 2] == 'o' || tmp[i + 2] == 'O') && (tmp[i + 3] == 's' || tmp[i + 3] == 'S'))
	{	for(i = i + 1; i < l; i++) if(tmp[i] == ' ' || tmp[i] == '\t') break;
	}
	else oldformat = true;
	vector<string> sampleid, factorid, ids;
	while(i < l)
	{	for(j = i + 1; j < l; j++) if(tmp[j] == ' ' || tmp[j] == '\t') break;
		tmp[j] = 0;
		ids.push_back(&tmp[i + 1]);
		for(k = i + 1; k < j; k++) if(tmp[k] == '.') break;
		if(k < j) tmp[k] = 0;
		sampleid.push_back(&tmp[i + 1]);
		factorid.push_back(&tmp[k + (int)(k < j)]);
		if(k < j) tmp[k] = '.';
		if(j < l) tmp[j] = ' ';
		i = j;
	}
	n = (int)sampleid.size();
	sid.clear(); fid.clear();
	for(i = 0; i < n; i++)
	{	for(j = 0; j < (int)sid.size(); j++) if(sid[j] == sampleid[i]) break;
		if(j >= (int)sid.size()) sid.push_back(sampleid[i]);
	}	
	for(i = 0; i < (int)factorid.size(); i++)
	{	for(j = 0; j < (int)fid.size(); j++) if(fid[j] == factorid[i]) break;
		if(j >= (int)fid.size()) fid.push_back(factorid[i]);
	}	
	indN = (int)sid.size();
	maxmsz = (int)fid.size();
	vector<vector<int> > repn(indN, vector<int>(maxmsz, 0)), trepn = repn;
	for(i = 0; i < n; i++)
	{	for(j = 0; j < (int)sid.size(); j++) if(sid[j] == sampleid[i]) break;
		for(k = 0; k < (int)fid.size(); k++) if(fid[k] == factorid[i]) break;
		repn[j][k]++;
	}
	indIndex = new int[indN + 1];
	indIndex[0] = 0;
	for(i = 0; i < indN; i++)
	{	k = 0;
		for(j = 0; j < maxmsz; j++) k = max(k, repn[i][j]);
		indIndex[i + 1] = indIndex[i] + k;
	}
	
	totalN = indIndex[indN];
	msz = new int[totalN];
	map = new int*[totalN];
	for(i = 0; i < indN; i++)
	{	for(j = indIndex[i]; j < indIndex[i + 1]; j++)
		{	map[j] = new int[maxmsz];
			k = msz[j] = 0;
			for(l = 0; l < maxmsz; l++) 
				if(repn[i][l] > j - indIndex[i]) 
				{	msz[j]++;
					map[j][k++] = l;
				}
		}
	}
	
//	vector<int> tn(totalN, 0);
	vector<vector<int> > samplemap;
	for(i = 0; i < n; i++)
	{	int ss, ff, tt;
		for(ss = 0; ss < indN; ss++) if(sampleid[i] == sid[ss]) break;
		for(ff = 0; ff < maxmsz; ff++) if(factorid[i] == fid[ff]) break;
		vector<int> row(2, 0);
		row[0] = indIndex[ss] + trepn[ss][ff];
		for(tt = 0; tt < msz[row[0]]; tt++) if(map[row[0]][tt] == ff) break;
		row[1] = tt;//tn[indIndex[ss] + trepn[ss][ff]];
		samplemap.push_back(row); 
//		tn[indIndex[ss] + trepn[ss][ff]]++;
		trepn[ss][ff]++;
	}
	vector<string> newsid(totalN);
	for(i = 0; i < indN; i++)
	{	if(indIndex[i + 1] - indIndex[i] > 1)
		{	char tstr[1000];
			for(j = indIndex[i]; j < indIndex[i + 1]; j++)
			{	sprintf(tstr, "%s.%d", sid[i].c_str(), j - indIndex[i]); 
				newsid[j] = tstr;
			}
		}
		else newsid[indIndex[i]] = sid[i];
	}
	sid = newsid;
/*
for(i = 0; i < (int)samplemap.size(); i++)
{	printf("%d: %d,%d\n", i, samplemap[i][0], samplemap[i][1]);
}
printf("maxymsz=%d\n", maxmsz);
for(i=0;i<(int)fid.size();i++) printf("%s ", fid[i].c_str());printf("\n");
for(i=0;i<totalN;i++)
{	printf("%d: %d ", i, msz[i]);
	for(j = 0; j < msz[i]; j++) printf("%d,",map[i][j]);
	printf("\n");
}
for(i=0;i<indN;i++)
{	printf("%d: ", i);
	for(j=0;j<maxmsz;j++) printf("%d,",repn[i][j]);
	printf("\n");
}
*/

	snpinfo.clear();
	float row[n];
	datamatrix = new float*[totalN];
	float *dp[totalN];
	if(bST < 0)
	{	bST = 0; bED = L;
	}
	else L = min(L, bED - bST);
	for(i = 0; i < totalN; i++) dp[i] = datamatrix[i] = new float[(bED-bST) * msz[i]];
	for(i = 0; getline(f, tmp); i++)
	{	if(i<bST) continue;
		if(i>=bED) break;
		l = (int)tmp.size();
		SNPINFO tsnp;
		char nd[100];
		sprintf(nd, "%d", i);
		tsnp.snpid = nd;
		/*if(tmp[0]=='c')
		{	tsnp.chr = atoi(&tmp[3]);
			if(tmp[3]=='X') tsnp.chr=23;
			else if(tmp[3]=='Y') tsnp.chr=24;
			else if(tmp[3] == 'M') tsnp.chr=25;
		}
		else tsnp.chr = atoi(&tmp[0]);
		*/
		for(j = 1; j < l; j++) if(tmp[j] == ' ' || tmp[j] == '\t') break;
		tsnp.chr = tmp.substr(0,j);

		tsnp.posst = atoi(&tmp[j+1]);
		if(!oldformat) { for(j = j + 1; j < l; j++) if(tmp[j] == ' ' || tmp[j] == '\t') break; }
		tsnp.posed = atoi(&tmp[j+1]);
		//tsnp.qual = -1;
		//tsnp.label = 1;
		snpinfo.push_back(tsnp);
		for(j = j + 1; j < l; j++) if(tmp[j] == ' ' || tmp[j] == '\t') break;
				
		k = 0;
		while(j < l - 1)
		{	double value = atof(&tmp[j + 1]);
			//if(log2 > 0) value = fast_log(value + log2) / fast_log(2.);
			row[k++] = value; 
			for(j = j + 1; j < l - 1; j++) if(tmp[j] == ' ' || tmp[j] == '\t') break;
		}
		for(j = 0; j < n; j++) 
		{	//mean += (double)row[j]; sd += (double)(row[j] * row[j]); 
			dp[samplemap[j][0]][samplemap[j][1]] = row[j];
		}
//		datan += n;
		for(j = 0; j < totalN; j++) dp[j] += msz[j];
	}
	f.close();
//	mean /= ((double)datan + 1e-10);
//	sd = sqrt(max(1e-5, sd / ((double)datan + 1e-10) - mean * mean));

}

//solve psi(x)=y
double genomicTensor::_solveDigamma(double x0, double y, double precision, int maxitern)
{	int i = 0;
	double x = x0, err;
	while(i < maxitern)
	{	x = x0 - (gsl_sf_psi(x0) - y) / gsl_sf_psi_1(x0);
		err = fabs((x - x0) / x0);
		if(err < precision) break;
		x0 = x;
		i++;
	}
	return(x);
}

//mle of mean and precision parameters
void genomicTensor::_getMLE_HyperP1(double *logp, double *ak, double &M, int sz, double N, double precision, int maxitern)
{	int i, cn = 0;
	double newak[sz], newM;
	if(N<=0) return;
	while(cn < maxitern)
	{	double s, err, err2, y = 0;
		for(i = 0; i < sz; i++)
		{	
			y += ak[i] * (logp[i] - gsl_sf_psi(M * ak[i]));
		}
		s = 0;
		for(i = 0; i < sz; i++)
		{	
			newak[i] = _solveDigamma(ak[i], logp[i] - y, precision, maxitern);
			s += newak[i];
		}
		err = y = 0;
		for(i = 0; i < sz; i++) 
		{	newak[i] /= s;
			err += (ak[i] - newak[i]) * (ak[i] - newak[i]) / ak[i];
			ak[i] = newak[i];
			y += ak[i] * (logp[i] - gsl_sf_psi(M * ak[i]));
		}
		y = 0;
		for(i = 0; i < sz; i++)
			y += ak[i] * ak[i] * gsl_sf_psi_1(M * ak[i]);	
		newM = M;//min(100.,max(1.,M - ds / ds2));
		err2 = (M - newM) * (M - newM) / M;
		M = newM;
		if(err < precision && err2 < precision) break;
		cn++;
	}
}

void genomicTensor::parseParameter(char const *parafile)
{
	FILE *f = fopen(parafile, "r");
	if(f != NULL)
	{	char tmp[1000];
		while(fgets(tmp, 1000, f) != NULL)
		{	int i, j, k, l = (int)strlen(tmp) - 1;
			double f;
			if(l < 3) continue;
			string str(tmp);
			for(i = 0; i < l; i++) if(tmp[i] == '=') break;
			if(str.substr(0,6) == "burnin")
			{	j = atoi(&tmp[i + 1]);
				if(j > 0) tpara.burnin = j; 
				printf("burnin = %d\n", tpara.burnin);
			}
			else if(str.substr(0,4) == "mcmc")
			{	j = atoi(&tmp[i+1]);
				if(j > 0) tpara.mcmc = j;
				printf("mcmc = %d\n", tpara.mcmc);
			}
			else if(str.substr(0,13) == "samplemaximum")
			{	tpara.samplemaximum = (bool)atoi(&tmp[i+1]);
				printf("samplemaximum = %d\n", (int)tpara.samplemaximum);
			}
			else if(str.substr(0,6) == "maxN_S")
			{	tpara.maxHapK = atoi(&tmp[i+1]);
				printf("maxHapK = %d\n", tpara.maxHapK);
			}
			else if(str.substr(0,6) == "maxN_O")
			{	j = atoi(&tmp[i + 1]);
				if(j > 0) maxPosSZ = j;
				printf("maxPosSZ = %d\n", maxPosSZ);
			}
			else if(str.substr(0,6) == "maxN_K")
			{	j = atoi(&tmp[i+1]);
				if(j > 0) maxG = j;
				printf("maxG = %d\n", maxG);
			}
			else if(str.substr(0,4) == "log2")
			{	log2 = atof(&tmp[i+1]);
				printf("log2 = %f\n", log2);
			}
			else if(str.substr(0,1) == "A")
			{	f = atof(&tmp[i+1]);
				if(f > 0) tpara.A = f;
				printf("A = %f\n", tpara.A);
			}
			else if(str.substr(0,7) == "recrate")
			{	f = atof(&tmp[i+1]);
				if(f >= 0) tpara.recrate = f;
				printf("recrate = %f\n", tpara.recrate);
			}
			else if(str.substr(0,5) == "minsd")
			{	f = atof(&tmp[i+1]);
				if(f >= 0) mylik.minerr = f;
				printf("minerr = %f\n", mylik.minerr);
			}
			else if(str.substr(0,7) == "fixcoef")
			{	j = atoi(&tmp[i + 1]);
				for(k = i + 2; k < l; k++)
				{	if(tmp[k] == ' ' || tmp[k] == '\t') break;
				}
				vector<double> row;
				while(k < l)
				{	f = atof(&tmp[k + 1]);
					for(k = k + 1; k < l; k++) if(tmp[k] == ':') break;
					double w = atof(&tmp[k + 1]);
					for(k = k + 1; k < l; k++) if(tmp[k] == ' ' || tmp[k] == '\t') break;
					row.push_back(f);
					row.push_back(w);
				}
				if((int)mylik.modelparameter.size() <= max(-1,j)+1) mylik.modelparameter.resize(max(-1,j)+2, vector<double>());
				mylik.modelparameter[max(-1,j)+1] = row;
				printf("%d(%d): ", j, (int)mylik.modelparameter.size());
				for(k = 0; k < (int)row.size(); k+=2) printf("%f:%f ", row[k], row[k+1]);
				printf("\n");
			}
			
		}
		fclose(f);
	}
}

void genomicTensor::_splitmergeCluster(int type)
{	
	double *oemt = emitMatrix;
	emitMatrix = NULL;
//	double lp0 = _logP(mydata, tpara);
//	float *newdata = new float[mydata.totalN * mydata.L];
	int i, j, mi, mj, tid;
//for(i=0;i<mydata.totalN;i++) printf("%d ",(int)mydata.data[i*mydata.L+166]);printf("\n");
	double dlp = mylik.splitmergeCluster(type, dataYP, dataXP, mydata.totalN, mydata.L, mydata.data, tpara.lapse, tpara.priorW, mi, mj, tid, tpara.indWeight);
//	double odlp = dlp;
	
//if(type==0) printf("mi=%d,mj=%d,tid=%d,dlp=%f ", mi, mj, tid, dlp),fflush(stdout);
	if(mi>=0)
	{
		if(type == 1) mj = tid;
		int estep = 0;
		double tepp[mylik.clustersz], *epp = &tepp[0];
		if(emitMatrix != NULL) 
		{	epp = emitMatrix + emitK * posSZ * mylik.clustersz;
			estep = mylik.clustersz;
		}
		else
		{	mylik.getStateCount(&tepp[0]);
			double ss = 1e-100;
			for(i = 0; i < mylik.clustersz; i++) ss += tepp[i];
			for(i = 0; i < mylik.clustersz; i++) tepp[i] /= ss;
		}
		int *op = basePop;
		double a = tpara.probAlpha * 1.;
		vector<vector<PARAUNIT> >::iterator ipm = tpara.param.begin();
		double N1 = 1., N2 = 1.;//, sign = 2. * (double)(type == 0) - 1.;
		for(i = 0; i < mydata.L; i++, ipm++, op++)
		{	double *eee = epp + estep * (*op);
			for(j = 0; j < tpara.maxK; j++)
			{	int n1 = 0, n2 = 0, l;
				//if((int)(*ipm)[j].size()%2 !=0) printf("!!!%d",(int)(*ipm)[j].size()),exit(0);
				for(l = 0; l < (int)(*ipm)[j].allele.size(); l+=2)
				{	if((*ipm)[j].allele[l] == mi) { n1 = (*ipm)[j].allele[l+1]; }
					if((*ipm)[j].allele[l] == mj) { n2 = (*ipm)[j].allele[l+1]; }
				}
				if(n1 > NUMPRECISION && n2 > NUMPRECISION)
				{	
					dlp += gsl_sf_lngamma(n1 + n2 + a * (*(eee + mi) + *(eee + mj)));
					dlp -= (gsl_sf_lngamma(n1 + a * *(eee + mi)) + gsl_sf_lngamma(n2 + a * *(eee + mj)));
					dlp -= gsl_sf_lngamma(a * (*(eee + mi) + *(eee + mj)));
					dlp += (gsl_sf_lngamma(a * *(eee + mi)) + gsl_sf_lngamma(a * *(eee + mj)));
				}
				N1+=n1;N2+=n2;
			}
		}
//dlp += gsl_sf_lngamma(N1+N2)-gsl_sf_lngamma(N1)-gsl_sf_lngamma(N2);

		double un = gsl_ran_flat(gammar, 0, 1.);

		if(type == 0)
		{	
			if(fast_log(un / (1. - un)) <= dlp) //to merge
			{	float *dp = mydata.data;
				for(i = 0; i < mydata.totalN; i++)
					for(j = 0; j < mydata.L; j++, dp++)
					{	if((int)(*dp) == mi || (int)(*dp) == mj) (*dp) = (float)(tid);
					}
			}
			_refreshParameter();
		}
		else
		{	if(fast_log(un / (1. - un)) <= dlp) //fail to split
			{	float *dp = mydata.data;
				for(i = 0; i < mydata.totalN; i++)
					for(j = 0; j < mydata.L; j++, dp++)
					{	if((int)(*dp) == tid) (*dp) = (float)(mi);
					}
			}
			_refreshParameter();
		}
	}
	emitMatrix = oemt;
}

void genomicTensor::_getLpVar(float *ypp, float *xpp, int nst, int ned, int *tymsz, int *txmsz, vector<vector<double> > const &props, double priorW, int asz, vector<double> &rt)
{
	int i, j, k, psz = (int)props.size();
	int n = ned - nst;

	TLP = new double[n * asz];
	int tlpzz=0;
	for(i = 0; i < n; i++, ypp += tymsz[i], xpp += txmsz[i], tlpzz += asz)
		mylik.computeLP(ypp, xpp, nst + i, priorW, &TLP[tlpzz]);

	rt.clear();
	rt.resize(psz, 0);
	double *tlp = TLP, maxlp = MINUSINFINITE, flp[asz];
	for(i = 0; i < n; i++, tlp += asz)
	{	for(k = 0; k < psz; k++)
		{	for(j = 0; j < asz; j++)
			{	flp[j] = tlp[j] + props[k][j];
				if(j == 0 || maxlp < flp[j]) { maxlp = flp[j]; }
			}
			double s = 0;
			for(j = 0; j < asz; j++)
				s += exp(flp[j] - maxlp);
			rt[k] += fast_log(s) + maxlp;
		}
	}

	delete []TLP; TLP = NULL;
}

void genomicTensor::_testVariability(vector<int> const &response, double *vprob)
{       int i, j, k;
	vector<MySortType> groupid(mydata.indN);
	for(i = 0; i < mydata.indN; i++)
	{	groupid[i].index = i;
		groupid[i].score = (double)response[i];
	}
	sort(groupid.begin(), groupid.end());
	int rN = 0;
	for(i = 0; i < mydata.indN; i++)
		if(i == mydata.indN - 1 || groupid[i].score != groupid[i+1].score) rN++;
	
	vector<MySortType2> ggg;
	for(j = 0; j < mydata.L; j++)
	{	
		vector<vector<int> > tggg(rN);
		float *dp = mydata.data + j;
		int rid = 0;
		for(int ii = 0; ii < mydata.indN; ii++)
		{	i = groupid[ii].index;
			int idst = mydata.indIndex[i], ided = mydata.indIndex[i + 1];
			for(k = 0; k < ided - idst; k++, dp += mydata.L)
			{	tggg[rid].push_back((int)(*dp));	
			}
			if(ii == mydata.indN - 1 || groupid[ii].score != groupid[ii + 1].score) rid++;
		}
		for(i = 0; i < rN; i++)
		{	
			sort(tggg[i].begin(), tggg[i].end());
			MySortType2 my;
			my.id1 = j; my.id2 = i;
			for(k = 0; k < (int)tggg[i].size(); k++)
			{	//ostringstream convert;
				char tmpstr[100];
				//convert<<tggg[i][k];
				sprintf(tmpstr, "%d,", tggg[i][k]);
				my.str += tmpstr;//convert.str() + ",";
			}
			ggg.push_back(my);
		}
	}

	vector<vector<int> > glabel(mydata.L, vector<int>(rN, -1));
	sort(ggg.begin(), ggg.end());
	vector<double> nnn;
	j = 0;
	int ok = 0;
	int pN = (int)ggg.size();
	double psuedo = 0, cn = psuedo+1, totalcn = 0;//, cutvalue = 0;
	for(k = 1; k < (int)ggg.size(); k++) 
	{	if(ggg[k].str != ggg[j].str || k == (int)ggg.size() - 1)
		{	if(true)//cn > cutvalue)
			{	nnn.push_back(cn + (int)(k==pN-1));
				totalcn += cn;
				for(int l = ok; l < k + (int)(k==pN-1); l++) glabel[ggg[l].id1][ggg[l].id2] = j;
				j++;
			}
	//		else printf("cut=%d\n", (int)cn);
			ggg[j] = ggg[k];
			cn = psuedo;
			ok = k;
		}
		cn++;
	}
	ggg.resize(j);
	pN = j;
	
	vector<vector<double> > prop0(pN, vector<double>(mylik.clustersz, 0));
	for(i = 0; i < pN; i++)
	{	char str[100];
		sprintf(str, "%s", ggg[i].str.c_str());
		k = (int)strlen(str) - 1;
		j = 0;
		int c = 0;
		do
		{	prop0[i][atoi(&str[j])]++;	
			c++;
			for(j = j; j < k; j++) if(str[j] == ',') break;
			j++;
		} while(j < k);
//		for(j = 0; j < mylik.clustersz; j++) prop0[i][j] = fast_log((prop0[i][j] + 1e-100) / (double)c);
	}

vector<vector<double> > score;
_pairEnrich(prop0, nnn, score);
printf("\n{\n");
for(i=0;i<(int)score.size();i++)
{	for(j=0;j<(int)score[i].size();j++)
		printf("%3.3f ", score[i][j]);
	printf("\n");
}
printf("}\n");
for(i = 0; i < (int)nnn.size(); i++)
{	printf("%d: ", (int)nnn[i]);
	for(j = 0; j < mylik.clustersz; j++) printf("%d ", (int)prop0[i][j]);
	printf("\n");
}

	printf("pN=%d(%d): ", (int)nnn.size(), (int)prop0.size());for(i=0;i<(int)nnn.size();i++) printf("%d,",(int)nnn[i]);printf("\n");
	for(i = 0; i < pN; i++)
	{	printf("%d: ", i);
		for(j = 0; j < mylik.clustersz; j++) printf("%d,",(int)prop0[i][j]); printf("\n");
	}
	printf("\n");
//prop0.erase(prop0.begin() + 4);
//prop0.erase(prop0.begin() + 2);
//prop0.erase(prop0.begin() + 0);
//nnn.erase(nnn.begin() + 4);
//nnn.erase(nnn.begin() + 2);
//nnn.erase(nnn.begin() + 0);
//	_propErrorCorrect(prop0, nnn, response, glabel, 1);
	printf("pN=%d(%d): ", (int)nnn.size(), (int)prop0.size());for(i=0;i<(int)nnn.size();i++) printf("%d,",(int)nnn[i]);printf("\n");
	pN = (int)prop0.size();
	for(i = 0; i < pN; i++)
	{	printf("%d: ", i);
		for(j = 0; j < mylik.clustersz; j++) printf("%d,",(int)prop0[i][j]); printf("\n");
	}
	printf("\n");

/*	vector<int> P(mylik.clustersz), R, X;
	vector<vector<int> > clique;
	for(i = 0; i < mylik.clustersz; i++) P[i] = i;
	BronKerbosch1(score, P, R, X, clique, -6);

	vector<vector<int> > gcluster;
	vector<int> gmap(mylik.clustersz, -1);
	while((int)clique.size() > 0)
	{	vector<MySortType> cliqueW((int)clique.size());
		for(i = 0; i < (int)cliqueW.size(); i++)
		{	cliqueW[i].index = i;
			cliqueW[i].score = (double)clique[i].size();
			for(j = 0; j < (int)clique[i].size() - 1; j++)
				for(k = j + 1; k < (int)clique[i].size(); k++)
					cliqueW[i].score *= 1. - 1. / (1. + exp(score[clique[i][j]][clique[i][k]]));
		}
		sort(cliqueW.begin(), cliqueW.end());
		k = cliqueW[(int)cliqueW.size() - 1].index;
		vector<int> t = clique[k];
		for(i = 0; i < (int)t.size(); i++) gmap[t[i]] = (int)gcluster.size();
		gcluster.push_back(t);
		clique.erase(clique.begin() + k);
		for(i = (int)clique.size() - 1; i >= 0; i--)
		{	for(j = (int)clique[i].size() - 1; j >= 0; j--)
			{	for(k = 0; k < (int)t.size(); k++)
					if(t[k] == clique[i][j]) break;
				if(k < (int)t.size()) clique[i].erase(clique[i].begin() + j);
			}
			if((int)clique[i].size() == 0) clique.erase(clique.begin() + i);
		}
	}
	printf("{\n");
	for(i = 0; i < (int)gcluster.size(); i++)	
	{	for(j = 0; j < (int)gcluster[i].size(); j++)
			printf("%d,",gcluster[i][j]);
		printf("\n");
	}
	printf("}\n");

	vector<int> lab(pN, -1);
	for(i = 0; i < pN; i++)
	{	for(j = 0; j < mylik.clustersz; j++)
			if(prop0[i][j] > 0) 
			{	if(lab[i] < 0) lab[i] = gmap[j];
				else if(lab[i] != gmap[j]) { lab[i] = -1; break; }
			}
	}
	vector<double> ttnn;
	vector<vector<double> > ttpp;
	vector<int> pmap(mylik.clustersz, -1);
	for(i = 0; i < pN; i++)
	{	if(lab[i] < 0)
		{	ttnn.push_back(nnn[i]);
			ttpp.push_back(prop0[i]);
		}
		else 
		{	if(pmap[lab[i]] < 0)
			{	pmap[lab[i]] = (int)ttpp.size();
				ttnn.push_back(0);
				ttpp.push_back(vector<double>(mylik.clustersz, 0));
			}
			j = pmap[lab[i]];
			for(k = 0; k < mylik.clustersz; k++) ttpp[j][k] = (ttpp[j][k] * ttnn[j] + prop0[i][k] * nnn[i]) / (ttnn[j] + nnn[i]);
			ttnn[j] += nnn[i];
		}
	}
	prop0 = ttpp;
	nnn = ttnn;
	printf("pN=%d-->%d ", pN, (int)prop0.size());
	pN = (int)prop0.size();
*/

vector<double> tn = nnn;
sort(tn.begin(), tn.end());
for(i = 0; i < (int)tn.size(); i++) printf("%d,", (int)tn[i]);
printf("\n");
	vector<vector<double> > newprop;
	vector<double> newnnn;
	_groupComposition(prop0, nnn, newprop, newnnn, 1., 1.);
	printf("pN=%d --> %d\n", pN, (int)newprop.size());
tn = newnnn;
sort(tn.begin(), tn.end());
for(i = 0; i < (int)tn.size(); i++) printf("%d,", (int)tn[i]);
printf("\n");
	prop0 = newprop;
	nnn = newnnn;
	pN = (int)prop0.size();
	for(i = 0; i < pN; i++)
	{	printf("%d: ", i);
		for(j = 0; j < mylik.clustersz; j++) printf("%d,",(int)prop0[i][j]); printf("\n");
	}
	printf("\n");

if(false)
{
vector<vector<double> > newprop;
vector<double> newnnn;
for(i=0;i<pN;i++)for(j=0;j<mylik.clustersz;j++) prop0[i][j]=exp(prop0[i][j]);
_groupComposition1(prop0, nnn, newprop, newnnn);
printf("pN=%d-->",pN);
prop0 = newprop;
nnn = newnnn;
pN = (int)prop0.size();
for(i=0;i<pN;i++)for(j=0;j<mylik.clustersz;j++) prop0[i][j]=fast_log(prop0[i][j]+1e-100);
printf("%d\n", pN);
}
/*
FILE *fff = fopen("t.txt","w");
for(i = 0; i < pN; i++)
{	fprintf(fff,"%f ", i, nnn[i]);
	for(j = 0; j < mylik.clustersz; j++) fprintf(fff,"%f ", prop0[i][j]);
	fprintf(fff,"\n");
}
fclose(fff);
*/
	vector<vector<double> > cprop0 = prop0;
	for(i = 0; i < pN; i++)
	{	for(j = 0; j < mylik.clustersz; j++)
		{	cprop0[i][j] = exp(prop0[i][j]);
			if(j > 0) cprop0[i][j] += cprop0[i][j - 1];
		}
	}
	vector<double> ncumsum(pN, 0);
	for(i = 0; i < pN; i++) 
	{	ncumsum[i] = nnn[i] / totalcn;
		nnn[i] = fast_log(ncumsum[i]);
		ncumsum[i] = pow((exp(nnn[i])*totalcn),0.9);
		if(i > 0) ncumsum[i] += ncumsum[i - 1];
	}


/*
{
FILE *ff = fopen("t.txt", "w");
for(i = 0; i < (int)prop0.size(); i++)
{	for(j = 0; j < (int)prop0[i].size(); j++)
		fprintf(ff, "%f ", prop0[i][j]);
	fprintf(ff, "%d\n", (int)(exp(nnn[i]) * (double)totalcn));
}
fclose(ff);
}
*/
//int ccccc = 0;

	vector<MySortType> pp(mydata.L), ppt(mydata.L * 10);
        for(i = 0; i < (int)ppt.size(); i++)
        {       
//double uuuuu = gsl_ran_flat(gammar, 0, ncumsum[pN - 1]);
//for(k = 0; k < pN; k++) if(uuuuu <= ncumsum[k]) break;
//int kkkkk = k;
		int kkk = -1;
		bool renew = true;
                vector<double> lp0(pN, 0), tlp1 = lp0, wwlp0 = lp0, wwlp1 = lp0, twlp1 = lp0;
		double wlp1 = 0, wlp0 = 0, dwlp1 = 0;
                double lp1 = 0, maxlp, s, un;
                float yp[maxymsz * mydata.totalN], xp[maxxmsz * mydata.totalN + 1], *ypp, *xpp;
		ypp = &yp[0]; xpp = &xp[0];
		for(int jj = 0; jj < mydata.indN; jj++)
                {       j = groupid[jj].index;
			vector<double> rt;
			int l, ll, lll, idst = mydata.indIndex[j], ided = mydata.indIndex[j + 1];
			if(i < mydata.L)
			{	lll = 0;
				for(k = idst; k < ided; k++)
				{	ll = i * ymsz[k];
					for(l = 0; l < ymsz[k]; l++)
						ypp[lll++] = dataYP[k][ll++];
				}
				if(maxxmsz > 0)
				{	lll = 0;
					for(k = idst; k < ided; k++)
					{	ll = i * xmsz[k];
						for(l = 0; l < xmsz[k]; l++)
						{	xpp[lll++] = dataXP[k][ll++];
						}
					}
				}
				_getLpVar(ypp, xpp, idst, ided, ymsz + idst, xmsz + idst, prop0, tpara.priorW, mylik.clustersz, rt);
        	                for(k = 0; k < pN; k++) 
				{	lp0[k] += rt[k];
					tlp1[k] += rt[k];
				}
				if(jj ==  mydata.indN - 1 || groupid[jj].score != groupid[jj+1].score)
				{	maxlp = tlp1[0] + nnn[0]; s = 0;
					for(k = 1; k < pN; k++) maxlp = max(maxlp, tlp1[k] + nnn[k]);
					for(k = 0; k < pN; k++) 
					{	s += exp(tlp1[k] + nnn[k] - maxlp); 
						tlp1[k] = 0;
					}
					lp1 += fast_log(s) + maxlp;
				}
			}
			
			if(jj==0 || groupid[jj].score != groupid[jj-1].score)
			{	if(renew)
				{	un = gsl_ran_flat(gammar, 0, ncumsum[pN - 1]);
					for(kkk = 0; kkk < pN; kkk++) if(un <= ncumsum[kkk]) break;
				}
				un = gsl_ran_flat(gammar, 0., 1.);
				if(un <= 1. - pow(1. - 0.5, 1. / (double)rN)) renew = true;
				else renew = false;
			}
			k = kkk;
//k=kkkkk;
			int sel[ided - idst];
			int llly = 0, lllx = 0;
			for(ll = 0; ll < ided - idst; ll++)
			{	un = gsl_ran_flat(gammar, 0, cprop0[k][mylik.clustersz - 1]);	
				for(l = 0; l < mylik.clustersz; l++) if(un <= cprop0[k][l]) break;
				mylik.simData(1, idst + ll, l, ypp + llly);
				sel[ll] = l;
				for(k = 0; k < xmsz[idst + ll]; k++) xpp[lllx + k] = 0;
				llly += ymsz[idst + ll];
				lllx += xmsz[idst + ll];
			}
			_getLpVar(ypp, xpp, idst, ided, ymsz + idst, xmsz + idst, prop0, tpara.priorW, mylik.clustersz, rt);
                        for(k = 0; k < pN; k++) 
			{	wwlp0[k] += rt[k];
				wwlp1[k] += rt[k];
				twlp1[k] += rt[k];
			}
                	
			if(jj == mydata.indN - 1 || groupid[jj].score != groupid[jj+1].score)
			{
				double tmaxlp = twlp1[0] + fast_log(ncumsum[0] / ncumsum[pN-1]), ts = 0; 
				maxlp = wwlp1[0] + nnn[0]; s = 0;
				for(k = 1; k < pN; k++)
				{	tmaxlp = max(tmaxlp, twlp1[k] + fast_log((ncumsum[k]-ncumsum[k-1])/ncumsum[pN-1]));
					maxlp = max(maxlp, wwlp1[k] + nnn[k]);
				}
				for(k = 0; k < pN; k++) 
				{	if(k==0) ts += exp(twlp1[k] + fast_log(ncumsum[k]/ncumsum[pN-1]) - tmaxlp);
					else ts += exp(twlp1[k] + fast_log((ncumsum[k]-ncumsum[k-1])/ncumsum[pN-1]) - tmaxlp);
					if(jj == mydata.indN - 1 || renew) twlp1[k] = 0;
					s += exp(wwlp1[k] + nnn[k] - maxlp);
					wwlp1[k] = 0;
				}
				wlp1 += fast_log(s) + maxlp;
				if(jj == mydata.indN - 1 || renew) dwlp1 += fast_log(ts) + tmaxlp;
			}

			ypp += llly;
			xpp += lllx;
		}       
		double maxlp0, s0, tmax;
		if(i < mydata.L)
 		{	maxlp0 = lp0[0] + nnn[0] * (double)rN; s0 = 0;
			for(j = 1; j < pN; j++) maxlp0 = max(maxlp0, lp0[j] + nnn[j] * (double)rN);
			for(j = 0; j < pN; j++) s0 += exp(lp0[j] + nnn[j] * (double)rN - maxlp0);
			tmax = max(lp1, maxlp0);
			lp1 = fast_log(max(exp(lp1 - tmax) - s0 * exp(maxlp0 - tmax),1e-100)) + tmax;

			maxlp = lp0[0] + nnn[0];
            s = 0;
			for(j = 1; j < pN; j++) maxlp = max(maxlp, lp0[j] + nnn[j]);
			for(j = 0; j < pN; j++) s += exp(lp0[j] + nnn[j] - maxlp);
			s += s0 * exp(maxlp0-maxlp);
			s = fast_log(s) + maxlp;
			//maxlp = max(maxlp, lp1);maxlp=max(maxlp,s);
			vprob[i] = lp1 - s;//exp(lp1 - maxlp) / (exp(lp1 - maxlp) + exp(s - maxlp));
			pp[i].index = i;
			pp[i].score = lp1-s;
		}

 		maxlp0 = wwlp0[0] + nnn[0] * (double)rN; s0 = 0;
		for(j = 1; j < pN; j++) maxlp0 = max(maxlp0, wwlp0[j] + nnn[j] * (double)rN);
		for(j = 0; j < pN; j++) s0 += exp(wwlp0[j] + nnn[j] * (double)rN - maxlp0);
		tmax = max(wlp1, maxlp0);
		wlp1 = fast_log(max(exp(wlp1 - tmax) - s0 * exp(maxlp0 - tmax),1e-100)) + tmax;

        maxlp = wwlp0[0] + nnn[0];
        s = 0;
		for(j = 1; j < pN; j++) maxlp = max(maxlp, wwlp0[j] + nnn[j]);
		for(j = 0; j < pN; j++) s += exp(wwlp0[j] + nnn[j] - maxlp);
		wlp0 = fast_log(s) + maxlp;
		s += s0 * exp(maxlp0-maxlp);
		s = fast_log(s) + maxlp;
		//maxlp = max(maxlp, wlp1);maxlp = max(maxlp, s);

		ppt[i].index = i;
		ppt[i].score = //fabs(yp[0]-yp[1]);
				wlp1-s;
				//exp(wlp1 - maxlp) / (exp(wlp1 - maxlp) + exp(s - maxlp));
		ppt[i].weight = exp(wlp0 - dwlp1);//k;
//printf("%d %f %f | %f %f %f | %f\n", i, yp[0], yp[1], wlp0, dwlp1, ppt[i].weight, ppt[i].score);
//printf("c %f %d\n",vprob[i], i);fflush(stdout);

//if(ppt[i].score >= 0.9) ccccc++;
//ppt[i].weight = 1.;

//pN -= rN;
	}

	sort(pp.begin(), pp.end());
	sort(ppt.begin(), ppt.end());
	reverse(pp.begin(), pp.end());	
	reverse(ppt.begin(), ppt.end());
	double ss = 0, ff, off = 1.;
	j = 0;
	for(i = 0; i < mydata.L; i++)
	{	while(j < (int)ppt.size() && ppt[j].score >= pp[i].score) { ss += ppt[j].weight; j++; }
		vprob[pp[i].index + mydata.L] = min(1., ss / (double)ppt.size());
	}
	for(i = mydata.L - 1; i >= 0; i--)
	{	ff = vprob[pp[i].index + mydata.L] * (double)(mydata.L) / (double)(i + 1);
		off = vprob[pp[i].index + mydata.L*2] = min(off, ff);
	}
	
}

void genomicTensor::_groupComposition(vector<vector<double> > &ggg, vector<double> &nnn, vector<vector<double> > &prop, vector<double> &newnnn, double a, double b)
{	int i, j, k, l, N = (int)ggg.size(), L = (int)ggg[0].size();

	vector<vector<double> > score;
	_pairEnrich(ggg, nnn, score);

	double buff  = 3.;
	vector<vector<double> > sim(N, vector<double>(N, -100000000.));
	for(i = 0; i < N; i++)
	{	for(j = i; j < N; j++)
		{	double f = 0, n = 0;
			for(k = 0; k < L; k++)
				if(ggg[i][k] > 0)
				{	for(l = 0; l < L; l++)
						if(ggg[j][l] > 0)
						{	f += min(0.,score[k][l]);//fast_log(1. - 1. / (1. + exp(score[k][l])));
							n++;
						}
				}
			sim[i][j] = sim[j][i] = f + buff;
		}
	}
for(i=0;i<N;i++)
{	for(j=0;j<N;j++)printf("%3.3f ",sim[i][j]);printf("\n");
}
/*
for(i = N - 1; i >= 0; i--)
	if(sim[i][i] < 0)
	{	ggg.erase(ggg.begin() + i);
		nnn.erase(nnn.begin() + i);
		sim.erase(sim.begin() + i);
		for(j = 0; j < (int)sim.size(); j++) sim[j].erase(sim[j].begin() + i);
	}
printf("N=%d-->%d\n", N, (int)ggg.size());
N = (int)ggg.size();
*/
a=1e-10;b=1e-10;
	double A = a * (double)L, B = b * N;
	double *GGG = new double[N * L * 2 + N * 2], *NNN = GGG + N * L, *PROP = NNN + N, *NEWNNN = PROP + N * L;
	double *tggg, *tprop;
	double totaln = 0;
	tggg = GGG; tprop = PROP;
	for(i = 0; i < N; i++)
	{	for(j = 0; j < L; j++, tggg++, tprop++)
		{	(*tggg) = ggg[i][j];
			NNN[i] = nnn[i];
			(*tprop) = (*tggg) * NNN[i]; 
		}
		totaln += NNN[i];
	}
	int gN = 0;
	int cn[N], map[N], nn[N];
	newnnn.resize(N);
	tggg = GGG;
	for(i = 0; i < N; i++)
	{	k = 0;
		for(j = 0; j < L; j++, tggg++) k += (int)(*tggg);
		NEWNNN[i] = NNN[i] * k;
		cn[i] = (int)NEWNNN[i];
		nn[i] = (int)NNN[i];
		map[i] = i;
		if(nnn[i] > 0) gN++;
	}
	
	double *gLP = new double[N * N * 2], *TP1 = gLP + N * N, *gp, *gm;
	double sumLP[N], F0[N];
	tggg = GGG;
	gm = gLP;
	for(l = 0; l < N; l++, tggg += L)
	{	double sum = 0, f0 = 0;
		for(k = 0; k < L; k++)
		{	f0 -= gsl_sf_lngamma(tggg[k]+1);
			sum += tggg[k];
		}
		f0 += gsl_sf_lngamma(sum+1);
		F0[l] = f0;
		for(k = 0; k < N; k++, gm++) *gm = f0;
		for(k = 0; k < L; k++)
		{	gp = gLP + l * N;
			tprop = PROP + k;
			for(j = 0; j < N; j++, gp++, tprop += L)
				(*gp) += fast_log((*tprop + a) / (cn[j] + A)) * tggg[k];
		}
	}
	gp = gLP;
	for(k = 0; k < N; k++, gp += N)
	{	double s = 0, maxlp = *gp + fast_log((nn[0] + b) / (totaln + B));
		for(l = 1; l < N; l++) maxlp = max(maxlp, gp[l] + fast_log((nn[l] + b) / (totaln + B)));
		for(l = 0; l < N; l++) s += exp(gp[l] + fast_log((nn[l] + b) / (totaln + B)) - maxlp);
		sumLP[k] = fast_log(s) + maxlp;
	}

	int burnin = 0, mcmc = 50;
	for(i = 0; i < burnin + mcmc; i++)
	{	bool flag = false;
		printf("%d ", i);fflush(stdout);
		for(j = 0; j < N; j++) printf("%d:%d ", j,nn[j]);printf("\n");
		tggg = GGG;
		for(j = 0; j < N; j++, tggg += L)
		{	printf("%d:%d ", j,map[j]);fflush(stdout);
			int mj = map[j];
			tprop = PROP + mj * L;
			for(k = 0; k < L; k++, tprop++)
				(*tprop) -= tggg[k] * NNN[j];
			cn[mj] -= (int)NEWNNN[j];
			nn[mj] -= (int)NNN[j];
			double lp[N], lp0 = MINUSINFINITE;
			double maxlp = (double)MINUSINFINITE;
			
			double *tm = TP1 + mj, *tp, *ttgg = GGG;
			for(k = 0; k < N; k++, tm += N)
			{	*tm = F0[k];
				tprop = PROP + mj * L;
				for(l = 0; l < L; l++, ttgg++, tprop++)
					if(*ttgg > 0) (*tm) += fast_log((*tprop + a) / (cn[mj] + A)) * (*ttgg);
			}
			for(k = 0; k < N; k++)
			{	lp[k] = (double)MINUSINFINITE;
				if(cn[k] > 0 || (cn[k] == 0 && lp0 <= MINUSINFINITE))
				{	lp[k] = 0;
			/*		for(l = 0; l < L; l++)
					{	if(GGG[j*N+l] > 0.01)
							lp[k] += gsl_sf_lngamma(PROP[k*N+l] + GGG[j*N+l] * NNN[j] + a) - gsl_sf_lngamma(PROP[k*N+l] + a);
					}
					lp[k] -= gsl_sf_lngamma((double)cn[k] + NEWNNN[j] + A) - gsl_sf_lngamma((double)cn[k] + A); 
			*/		
					if(k != mj)
					{	tp = TP1 + k; ttgg = GGG;
						for(int kk = 0; kk < N; kk++, tp += N)
						{	*tp = F0[kk];
							tprop = PROP + k * L;
							for(l = 0; l < L; l++, ttgg++, tprop++)
								if(*ttgg > 0) (*tp) += fast_log((*tprop + tggg[l] * NNN[j] + a) / (cn[k] + NEWNNN[j] + A)) * (*ttgg);
						}
					}
					gp = gLP + k; tp = TP1 + k;
					gm = gLP + mj; tm = TP1 + mj;
					for(l = 0; l < N; l++, gp += N, gm += N, tp += N, tm += N)
					{	double newsumLP = sumLP[l], tmax = max(newsumLP, *tp);
						if(k != mj)
						{	newsumLP = fast_log(max(exp(newsumLP - tmax) - exp(*gp + fast_log((nn[k] + b) / (totaln + B)) - tmax) + exp(*tp + fast_log((nn[k] + NNN[j] + b) / (totaln + B)) - tmax), 1e-100)) + tmax;
							tmax = max(newsumLP, *tm);
							newsumLP = fast_log(max(exp(newsumLP - tmax) - exp(*gm + fast_log((nn[mj] + NNN[j] + b) / (totaln + B)) - tmax) + exp(*tm + fast_log((nn[mj] + b) / (totaln + B)) - tmax), 1e-100)) + tmax;
						}
						lp[k] += newsumLP * NNN[l];
					}
			
//printf("%d(%d).%d:(%d) %f->",j,map[j],k,cn[k],lp[k]);double dd = lp[k];
				double bb=totaln / 1.;	
					lp[k] += gsl_sf_lngamma((double)nn[k] + nnn[j] + bb) - gsl_sf_lngamma((double)nn[k] + bb); 
			//	int remove = (int)(nn[mj] == 0), add = (int)(nn[k] == 0);
			//	lp[k]+=(double)(remove-add)*10.*(double)(gN+add-remove > 10);
//printf("%f (%f)| ",lp[k],lp[k]-dd);
				for(int tt = 0; tt < N; tt++) 
				{	//if(tt != j && map[tt] == k && sim[tt][j] < 0) { lp[k] += sim[tt][j] * 10.; }
				//	if((map[tt] == k || tt == j) && sim[tt][tt] < 0 && nn[k] <= nnn[tt]) lp[k] += sim[tt][tt];//1000000; 
				}
				//if(sim[j][j] < 0 && nn[k]  < nnn[j]) lp[k] -= 1000000.;

					maxlp = max(maxlp, lp[k]);
					if(cn[k] == 0) lp0 = lp[k];
				}
			}
			if(i < burnin)
			{	for(k = 0; k < N; k++)
				{	lp[k] = exp(lp[k] - maxlp);
				}
			}
			k = _sample(lp, N, (int)(i>=burnin));
//printf("%d\n",k);fflush(stdout);
			if(mj != k && i >= burnin) //printf("%d %d:%d->%d %f %f, %d %d\n",i,j,map[j],k, lp[map[j]], lp[k], nn[map[j]], nn[k]),
				flag = true;
			if(mj != k)
			{	gp = gLP + k; tp = TP1 + k;
				gm = gLP + mj; tm = TP1 + mj;
				for(l = 0; l < N; l++, gp += N, gm += N, tp += N, tm += N)
				{	double newsumLP = sumLP[l], tmax = max(newsumLP, *tp);
					newsumLP = fast_log(max(exp(newsumLP - tmax) - exp(*gp + fast_log((nn[k] + b) / (totaln + B)) - tmax) + exp(*tp + fast_log((nn[k] + NNN[j] + b) / (totaln + B)) - tmax), 1e-100)) + tmax;
					tmax = max(newsumLP, *tm);
					newsumLP = fast_log(max(exp(newsumLP - tmax) - exp(*gm + fast_log((nn[mj] + NNN[j] + b) / (totaln + B)) - tmax) + exp(*tm + fast_log((nn[mj] + b) / (totaln + B)) - tmax), 1e-100)) + tmax;
					sumLP[l] = newsumLP;
				}
				gp = gLP + k; tp = TP1 + k;
				gm = gLP + mj; tm = TP1 + mj;
				for(l = 0; l < N; l++, gp += N, gm += N, tp += N, tm += N) 
				{	*gp = *tp;
					*gm = *tm;
				}
			}
			if(nn[mj]==0) gN--; if(nn[k]==0) gN++;
			mj = map[j] = k;
			tprop = PROP + mj * L;
			for(k = 0; k < L; k++, tprop++)
				*tprop += tggg[k] * NNN[j];
			cn[mj] += (int)NEWNNN[j];
			nn[mj] += (int)NNN[j];
		}
		if(!flag) break;
	}
	delete []gLP;
	
	vector<double> totalprop(L,10);
	tprop = PROP;
	for(i = 0; i < N; i++)
		for(j = 0; j < L; j++, tprop++)
			totalprop[j] += (*tprop + a);
	double totalcn = 0;
	for(i = 0; i < L; i++) totalcn += totalprop[i];
	
	tprop = PROP;
	for(i = 0; i < N; i++, tprop += L)
	{	for(j = 0; j < L; j++) tprop[j] = fast_log((tprop[j] + a) / (cn[i] + A));
						//fast_log((tprop[j] + a * totalprop[j]/totalcn)/(cn[i]+A));
		NEWNNN[i] = (double)nn[i] + b;
	}

	prop.clear();
	newnnn.clear();
	tprop = PROP;
	for(i = 0; i < N; i++, tprop += L)
	{	if(NEWNNN[i] > b + 10)
		{	vector<double> row(L);
			for(j = 0; j < L; j++) row[j] = tprop[j];
			prop.push_back(row);
			newnnn.push_back(NEWNNN[i]);
		}
	}
	
	delete []GGG;
}

void genomicTensor::_pairEnrich(vector<vector<double> > const &prop, vector<double> const &nnn, vector<vector<double> > &score)
{	int i, j, k;
	int K = (int)prop[0].size();
	double sum = 0;
	vector<vector<double> > A(K, vector<double>(K, 0));
	for(i = 0; i < (int)prop.size(); i++)
	{	sum = 1e-100;
		for(j = 0; j < K; j++) sum += prop[i][j];
		for(j = 0; j < K; j++)
			for(k = j; k < K; k++)
				if(prop[i][j] > 0 && prop[i][k] > 0)
				{	double t = prop[i][j] * (prop[i][k] - (double)(j==k)) / (1. + (double)(j==k)) * nnn[i] / (sum - 1.);
					A[j][k] += t;
					A[k][j] += t;
				}
	}
	sum = 1e-100;
	vector<double> P(K, 0);
	for(i = 0; i < K; i++)
	{	for(j = 0; j < K; j++) P[i] += A[i][j];
		sum += P[i];
	}
	for(i = 0; i < K; i++) P[i] /= sum;
	vector<vector<double> > E(K, vector<double>(K));
	for(i = 0; i < K; i++)
	{	for(j = i; j < K; j++)
			E[i][j] = E[j][i] = P[i] * P[j] * sum;
	}
	score = A;
	for(i = 0; i < K; i++)
	{	for(j = i; j < K; j++)
		{	score[i][j] = score[j][i] = (A[i][j] - E[i][j]) / sqrt(E[i][j] + 1.);
		}
	}
}

void genomicTensor::BronKerbosch1(vector<vector<double> > const &score, vector<int> P, vector<int> R, vector<int> X, vector<vector<int> > &clique, double scorecut)
{	int i;

/*			printf("[P: ");
			for(i = 0; i < (int)P.size(); i++) printf("%d,",P[i]);
			printf("]--");
			printf("[R: ");
			for(i = 0; i < (int)R.size(); i++) printf("%d,",R[i]);
			printf("]--");
			printf("[X: ");
			for(i = 0; i < (int)X.size(); i++) printf("%d,",X[i]);
			printf("]\n");
*/
	if((int)P.size() == 0 && (int)X.size() == 0) 
	{	if((int)R.size() > 0)
		{	clique.push_back(R);
			printf("[ ");
			for(i = 0; i < (int)R.size(); i++) printf("%d,",R[i]);
			printf("]\n");
		}
		return;
	}
	
	while((int)P.size() > 0)
	{	int v = P[0];
		vector<int> nP, nX, nR = R;
		nR.push_back(v);
		for(i = 1; i < (int)P.size(); i++)
		{	if(score[v][P[i]] > scorecut) nP.push_back(P[i]);
		}	
		for(i = 0; i < (int)X.size(); i++)
		{	if(score[v][X[i]] > scorecut) nX.push_back(X[i]);
		}	
		BronKerbosch1(score, nP, nR, nX, clique, scorecut);
		P.erase(P.begin());
		X.push_back(v);
	}
}

void genomicTensor::_groupComposition1(vector<vector<double> > const &ggg, vector<double> const &nnn, vector<vector<double> > &prop, vector<double> &newnnn)
{
	vector<vector<double> > score;
	_pairEnrich(ggg, nnn, score);

	int N = (int)ggg.size(), L = (int)ggg[0].size();
	int i, j, k, l;
	double buff  = 3.;
	vector<vector<double> > sim(N, vector<double>(N, -100000000.));
	for(i = 0; i < N - 1; i++)
	{	for(j = i; j < N; j++)
		{	double f = 0, n = 0;
			for(k = 0; k < L; k++)
				if(ggg[i][k] > 0)
				{	for(l = 0; l < L; l++)
						if(ggg[j][l] > 0)
						{	f += min(0., score[k][l]);//fast_log(1. - 1. / (1. + exp(score[k][l])));
							n++;
						}
				}
			sim[i][j] = sim[j][i] = f + buff;
		}
	}
	
	vector<int> list(N);
	prop = ggg;
	for(i = 0; i < N; i++)
	{	for(j = 0; j < L; j++) prop[i][j] *= nnn[i];
		list[i] = i;
	}
	newnnn = nnn;
	bool a_inner = true;
	do
	{	int tN = (int)list.size();
		vector<MySortType> me(tN * (tN - 1) / 2);
		k = 0;
		for(i = 0; i < tN - 1; i++)
			for(j = i + 1; j < tN; j++)
			{	me[k].score = sim[list[i]][list[j]];
				me[k].index = i * tN + j;
				k++;
			}
		sort(me.begin(), me.end());
		if(me[(int)me.size() - 1].score > 0)
		{	i = me[(int)me.size() - 1].index;
			j = i % tN;
			i = (int)(i / tN);
			for(k = 0; k < L; k++) 
			{	prop[list[i]][k] += prop[list[j]][k];
				prop[list[j]][k] = 0;
			}
			newnnn[list[i]] += newnnn[list[j]];
			newnnn[list[j]] = 0;
			list.erase(list.begin() + j);
			int ii = list[i];
			for(k = 0; k < tN - 1; k++)
				if(k != i)
				{	int jj = list[k];
					double f = 0, n = 0;
					for(int kk = 0; kk < L; kk++)
						if(prop[ii][kk] > 0)
						{	for(l = 0; l < L; l++)
								if(prop[jj][l] > 0)
								{	f += fast_log(1. - 1. / (1. + exp(score[kk][l])));
									n++;
								}
						}	
					sim[ii][jj] = sim[jj][ii] = f + buff;	
				}
		}
		else a_inner = false;
	} while(a_inner);

	j = 0;
	for(i = 0; i < N; i++)
		if(newnnn[i] > 0)
		{	prop[j] = prop[i];
			for(k = 0; k < L; k++) prop[j][k] = ((prop[i][k] + 0) / newnnn[i]);
			newnnn[j] = newnnn[i];
			j++;
		}
	prop.resize(j);
	newnnn.resize(j);
}

void genomicTensor::_propErrorCorrect(vector<vector<double> > &ggg, vector<double> &nnn, vector<int> const &response, vector<vector<int> > &glabel, double a)
{	int i, j, k, l, pN = (int)ggg.size();
	vector<vector<double> > prop0 = ggg;
	for(i = 0; i < pN; i++)
	{	double s = 0;
		for(j = 0; j < mylik.clustersz; j++) s += ggg[i][j] + 1e-100;
		for(j = 0; j < mylik.clustersz; j++) prop0[i][j] = fast_log((ggg[i][j] + 1e-100)/s);
	}
	vector<MySortType> groupid(mydata.indN);
	for(i = 0; i < mydata.indN; i++)
	{	groupid[i].index = i;
		groupid[i].score = (double)response[i];
	}
	sort(groupid.begin(), groupid.end());
	int rN = 0;
	for(i = 0; i < mydata.indN; i++)
		if(i == mydata.indN - 1 || groupid[i].score != groupid[i+1].score) rN++;

	vector<vector<double> > score;
	_pairEnrich(ggg, nnn, score);

	double buff  = 6.;
	vector<vector<double> > sim(pN, vector<double>(pN, -100000000.));
	for(i = 0; i < pN; i++)
	{	for(j = i; j < pN; j++)
		{	double f = 0, n = 0;
			for(k = 0; k < mylik.clustersz; k++)
				if(ggg[i][k] > 0)
				{	for(l = 0; l < mylik.clustersz; l++)
						if(ggg[j][l] > 0)
						{	f += min(0.,score[k][l]);//fast_log(1. - 1. / (1. + exp(score[k][l])));
							n++;
						}
				}
			sim[i][j] = sim[j][i] = f + buff;
		}
	}
for(i=0;i<pN;i++)printf("%f ",sim[i][i]);

        float yp[maxymsz * mydata.totalN], xp[maxxmsz * mydata.totalN + 1];
	int burnin = 0, mcmc = 50;
	for(int r = 0; r < burnin + mcmc; r++)
	{	bool flag = false;
		for(i = 0; i < mydata.L; i++)
		{	vector<double> lp(pN, 0);
			int rid = 0;
			for(int jj = 0; jj < mydata.indN; jj++)
			{	j = groupid[jj].index;
				vector<double> rt;
				int ll, lll, idst = mydata.indIndex[j], ided = mydata.indIndex[j + 1];
				lll = 0;
				for(k = idst; k < ided; k++)
				{	ll = i * ymsz[k];
					for(l = 0; l < ymsz[k]; l++)
						yp[lll++] = dataYP[k][ll++];
				}
				if(maxxmsz > 0)
				{	lll = 0;
					for(k = idst; k < ided; k++)
					{	ll = i * xmsz[k];
						for(l = 0; l < xmsz[k]; l++)
						{	xp[lll++] = dataXP[k][ll++];
						}
					}
				}
				_getLpVar(yp, xp, idst, ided, ymsz + idst, xmsz + idst, prop0, tpara.priorW, mylik.clustersz, rt);
        	                for(k = 0; k < pN; k++) 
				{	lp[k] += rt[k];
				}
				if(jj ==  mydata.indN - 1 || groupid[jj].score != groupid[jj+1].score) 
				{	double tlp[pN], maxlp = MINUSINFINITE;
if(i>=(int)glabel.size() || rid >= (int)glabel[0].size() || glabel[i][rid] < 0 || glabel[i][rid] >= pN)
{	printf("!!!%d,%d,%d,%d,%d,%d\n",i,(int)glabel.size(), rid,(int)glabel[0].size(),glabel[i][rid],pN);fflush(stdout);
}
					for(k = 0; k < pN; k++)
					{	tlp[k] = lp[k] + fast_log(nnn[k] - (double)(k == glabel[i][rid]) + a);
						if(k == 0 || maxlp < lp[k]) maxlp = tlp[k];
					}
					if(r < burnin)
					{	for(k = 0; k < pN; k++) tlp[k] = exp(tlp[k] - maxlp);	
					}
					k = _sample(tlp, pN, (int)(r >= burnin));
int t=glabel[i][rid];if(sim[t][t]>=0) k=glabel[i][rid];
					if(k != glabel[i][rid])
					{	nnn[glabel[i][rid]]--;
						glabel[i][rid] = k;
						nnn[k]++;
						flag = true;
					}
					lp.clear();
					lp.resize(pN, 0);
					rid++;
				}
			}
		}
		if(!flag && r >= burnin) break;
	}

	for(i = (int)nnn.size() - 1; i >= 0; i--)
		if(nnn[i] == 0)
		{	ggg.erase(ggg.begin() + i);	
			nnn.erase(nnn.begin() + i);
		}
}	

void genomicTensor::_refreshParameter()
{	int i, j;
	mylik.clearParameter();
	tpara.param.clear();
	PARAUNIT pu;
	pu.lapseN = pu.popN = 0;
	vector<PARAUNIT> eee(tpara.maxK, pu);//, vector<float>(mylik.clustersz, 0));
	tpara.param.resize(mydata.L, eee);
        float *dp = mydata.data;

	for(int id = 0; id < mydata.indN; id++)
        	for(i = mydata.indIndex[id]; i < mydata.indIndex[id+1]; i++)
      	       	{	vector<vector<PARAUNIT> > ::iterator ipm = tpara.param.begin();
			int *pp = tpara.pop + id * mydata.L;
			float *dpy = dataYP[i], *dpx = NULL;
			if(maxxmsz > 0) dpx = dataXP[i];
			float prevd = *dp;
			bool *tlapse = tpara.lapse + i * mydata.L;
                      	for(j = 0; j < mydata.L; j++, dp++, ipm++, pp++, dpy += ymsz[i], dpx += xmsz[i], tlapse++)
			{	
//if(*pp >= tpara.maxK) printf("!!!%d>=%d, %d %d\n", (int)*pp, tpara.maxK, id, j),exit(0);
				_updatePara_add((*ipm)[(int)(*pp)], (int)(*dp), tpara.indWeight[i], i, (int)(tpara.lapse!=NULL && *tlapse));//j);
				if(tpara.lapse == NULL || !(*tlapse)) prevd = *dp;
				mylik.addPara((int)prevd, dpy, dpx, i, tpara.indWeight[i]);
			}
       		}
	mylik.updateLambda(tpara.priorW);
	
	if(tpara.nh == NULL) tpara.nh = new float[tpara.maxK * mydata.L];
	
	for(i = 0; i < mydata.L; i++) 
		for(j = 0; j < tpara.maxK; j++)
		{	tpara.nh[i * tpara.maxK + j] = 0;
			for(int k = 0; k < (int)tpara.param[i][j].allele.size(); k += 2)
				tpara.nh[i * tpara.maxK + j] += tpara.param[i][j].allele[k + 1];
		}
}

void genomicTensor::removeData(vector<vector<int> > const &list)
{
	int i, j, k;
	for(i = 0; i < (int)list.size(); i++)
		if((int)list[i].size() > 0)
		{	int keep[ymsz[i]];
			for(j = 0; j < ymsz[i]; j++) keep[j] = 0;
			printf("%d: ", i);
			for(j = 0; j < (int)list[i].size(); j++) { keep[list[i][j]] = 1; printf("%d,",list[i][j]); }; 
			printf("\n");
			k = 0; for(j = 0; j < ymsz[i]; j++) if(keep[j] == 0) keep[k++] = j;
			int omsz = ymsz[i];
			ymsz[i] = k;
			for(j = 0; j < k; j++) ymap[i][j] = ymap[i][keep[j]]; 
			float *dp = dataYP[i], *dp1 = dp; 
			for(j = 0; j < mydata.L; j++, dp1 += omsz)
			{	for(k = 0; k < ymsz[i]; k++)
				{	*dp = dp1[keep[k]];
					dp++;
				}
			}
			printf("%d:sz=%d ", i, ymsz[i]);
			for(j = 0; j < ymsz[i]; j++) printf("%d,", ymap[i][j]); printf("\n");
		}
//	for(i = 0; i < mydata.totalN; i++) printf("%f ", dataYP[i][3]); printf("\n");
//	exit(0);
//	for(i=0;i<mydata.totalN;i++) printf("%d:%d ", i, ymsz[i]);
}

void genomicTensor::loadPreviousRun(char const *statefile, char const *clusterfile, double posRate)
{	
	ifstream f(statefile);
	string tmp;
	getline(f, tmp);
	int i, j, maxg = 0, maxp = 0;
	for(i = 0; i < mydata.L; i++)
	{	getline(f, tmp);
		istringstream buf(tmp);
                istream_iterator<string> beg(buf), end;
                vector<string> tokens(beg, end);
		for(j = 0; j < mydata.totalN; j++)
		{	int s = atoi(tokens[j + 4].c_str());
			mydata.data[j * mydata.L + i] = (float)s;
			maxg = max(maxg, s);
		}
		basePop[i] = atoi(tokens[4 + mydata.totalN].c_str());
		maxp = max(maxp, basePop[i]);
		if(i == 0) baseR[i] = 1.;
		else baseR[i] = (basePop[i] != basePop[i - 1] | gsl_ran_flat(gammar, 0, 1.) < posRate);
	}
	f.close();
	mylik.clustersz = max(mylik.clustersz, maxg + 1);
	posSZ = max(posSZ, maxp + 1);
	
	j = mydata.indN * mydata.L;
	tpara.pop = new int[j];
	int *pp = tpara.pop, maxk = 0;
	tpara.rec = new bool[j];
	bool *rr = tpara.rec;
	if(clusterfile != NULL) f.open(clusterfile);
	if(clusterfile == NULL || f.fail())
	{	for(i = 0; i < j; i++, pp++) (*pp) = 0;
		for(i = 0; i < j; i++, rr++) 
		{	(*rr) = false;
                        if((int)mydata.snpinfo.size() > 0 && i % mydata.L > 0)
                        {       if(mydata.snpinfo[i % mydata.L].chr == mydata.snpinfo[i % mydata.L - 1].chr && mydata.snpinfo[i % mydata.L].posst >= mydata.snpinfo[i%mydata.L-1].posst)
                                        (*rr) = (gsl_ran_flat(gammar, 0, 1.) < 0.99 * (1. - exp(- tpara.recrate*1e-6 * (double)max(1, mydata.snpinfo[i%mydata.L].posst- mydata.snpinfo[i%mydata.L - 1].posst))));
                                else (*rr) = true;
                        }
		}
	}
	else
	{	for(i = 0; i < mydata.indN; i++, pp += mydata.L, rr += mydata.L)
		{	getline(f, tmp);
			istringstream buf(tmp);
			istream_iterator<string> beg(buf), end;
			vector<string> tokens(beg, end);
			int k = (int)tokens[1].find(":"), l;
			int pos, lastpos = atoi(tokens[1].c_str());
			int state, laststate = atoi(&tokens[1].c_str()[k + 1]);
			maxk = max(maxk, laststate);	
			for(j = 2; j < (int)tokens.size(); j++)
			{	k = (int)tokens[j].find(":");
				pos = atoi(tokens[j].c_str());
				state = atoi(&tokens[j].c_str()[k + 1]);		
				for(l = lastpos; l < pos; l++) pp[l] = laststate;
				rr[pos] = true;
				laststate = state;
				maxk = max(maxk, laststate);	
				lastpos = pos;
			}
			for(l = lastpos; l < mydata.L; l++) pp[l] = laststate;
		}
		f.close();
	}
	tpara.maxK = max(tpara.maxK, maxk + 1);

	_refreshParameter();
/*for(i=0;i<tpara.param[1194].size();i++)
{	printf("%d: ", i);
	for(j=0;j<(int)tpara.param[1194][i].allele.size(); j++)
		printf("%f ", tpara.param[1194][i].allele[j]);
	printf("\n");
}*/
//	gID = -1;
//	_updatePrior(vector<bool>(), true);
//	gID = 0;
}

void genomicTensor::loadOtherResult(char const *state0file, char const *cluster0file, double posRate)
{
//	if(cluster0file != NULL && state0file != NULL)
	{	ifstream f(state0file);
		if(f.fail()) { printf("Cannot read %s\n", state0file);return; }
		string tmp;
		getline(f, tmp);
		istringstream buf(tmp);
		istream_iterator<string> beg(buf), end;

		vector<string> ind0(beg, end);
		int offset = 4;
		if(ind0[3]!="POSed") offset--;
		ind0.erase(ind0.begin(), ind0.begin() + offset);
		ind0.resize((int)ind0.size() - 1);
		vector<string> id0 = ind0;
		int n0 = 0, i, j, k;
		vector<int> map0((int)id0.size(), -1), poprep0;
		for(i = 0; i < (int)id0.size(); i++)
		{	j = (int)id0[i].find(".");
			if(j >= 0) id0[i] = id0[i].substr(0,j);
			for(j = 0; j < i; j++) if(id0[i] == id0[j]) break;
			if(j >= i)
			{	map0[i] = n0;
				n0++;
				poprep0.resize(n0, 0);
			}
			else map0[i] = map0[j];
			poprep0[map0[i]]++;
		}
		f.close();
		tpara0.totalN0 = (int)id0.size();
		mylik.totalN0 += tpara0.totalN0;

		f.open(cluster0file);
		if(f.fail()) { printf("Cannot read %s\n", cluster0file);return; }
		int *tpop0 = new int[n0 * mydata.L];
		tpara0.recombN0 = new float[mydata.L];
		tpara0.Q0 = new double[n0];
		for(i = 0; i < mydata.L; i++) tpara0.recombN0[i] = 0;
		for(i = 0; i < n0; i++) tpara0.Q0[i] = 0;
		int *pp0 = tpop0, maxk = 0;
		for(i = 0; i < n0; i++, pp0 += mydata.L)
		{	getline(f, tmp);
			istringstream buf(tmp);
			istream_iterator<string> beg(buf), end;
			vector<string> tokens(beg, end);
			int k = (int)tokens[1].find(":"), l;
			int pos, lastpos = atoi(tokens[1].c_str());
			int state, laststate = atoi(&tokens[1].c_str()[k + 1]);		
			maxk = max(maxk, laststate);
			for(j = 2; j < (int)tokens.size(); j++)
			{	k = (int)tokens[j].find(":");
				pos = atoi(tokens[j].c_str());
				state = atoi(&tokens[j].c_str()[k + 1]);		
				if(pos >= mydata.L) break;
				for(l = lastpos; l < pos; l++) pp0[l] = laststate;
				tpara0.recombN0[lastpos]++;
				tpara0.Q0[laststate]++;
				laststate = state;
				lastpos = pos;
				maxk = max(maxk, laststate);
			}
			for(l = lastpos; l < mydata.L; l++) pp0[l] = laststate;
			tpara0.recombN0[lastpos]++;
			tpara0.Q0[laststate]++;
		}
		f.close();
	
		tpara0.popn0 = maxk + 1;
		tpara.maxK = max(tpara.maxK, tpara0.popn0);
		tpara0.nh0 = new float[tpara0.popn0 * mydata.L];
		for(i = 0; i < tpara0.popn0 * mydata.L; i++) tpara0.nh0[i] = 0;
		pp0 = tpop0;
		for(i = 0; i < n0; i++, pp0 += mydata.L)
		{	for(j = 0; j < mydata.L; j++)
				tpara0.nh0[j * tpara0.popn0 + pp0[j]] += poprep0[i];
		}

		f.open(state0file);
		getline(f, tmp);
		int maxg = 0, maxp = 0;
		for(i = 0; i < mydata.L; i++)
		{	getline(f, tmp);
			istringstream buf(tmp);
			istream_iterator<string> beg(buf), end;
			vector<string> tokens(beg, end);
			for(j = offset; j < (int)tokens.size() - 1; j++)
			{	k = atoi(tokens[j].c_str());
				maxg = max(maxg, k);
			}
	
			basePop[i] = atoi(tokens[j].c_str());
			maxp = max(maxp, basePop[i]);
			if(i == 0) baseR[i] = 1.;
			else baseR[i] = (basePop[i] != basePop[i - 1] | gsl_ran_flat(gammar, 0, 1.) < posRate);
		}
		tpara0.clustersz0 = maxg + 2;
		tpara0.unitsz0 = tpara0.popn0 * tpara0.clustersz0;
		mylik.clustersz = max(mylik.clustersz, maxg + 1);
		posSZ = max(posSZ, maxp + 1);
		f.close();	

		f.open(state0file);
		getline(f, tmp);
		tpara0.param0 = new unsigned short[mydata.L * tpara0.unitsz0];
		unsigned short *tparam0 = tpara0.param0;
		for(i = 0; i < mydata.L; i++, tparam0 += tpara0.unitsz0)
		{	getline(f, tmp);
			istringstream buf(tmp);
			istream_iterator<string> beg(buf), end;
			vector<string> tokens(beg, end);
			for(j = 0; j < tpara0.unitsz0; j++) tparam0[j] = 0;
			for(j = offset; j < (int)tokens.size() - 1; j++)
			{	k = atoi(tokens[j].c_str());
				tparam0[tpop0[map0[j - offset] * mydata.L + i] * tpara0.clustersz0 + k + 1]+=1.;
			}
		}
		f.close();

		delete []tpop0;
	}
}

void genomicTensor::_loadData_list(float **&datamatrix, char const *listfile, char const *bedfile, vector<string> &sid, vector<string> &fid, vector<SNPINFO> &snpinfo, int &indN, int &L, int &totalN, int &maxmsz, int *&msz, int **&map, int *&indIndex)
{

	ifstream f(bedfile);
	string tmp;
	int i, j, k, n, l;
	
//do{cout<<'\n'<<"0) Press a key to continue ..."; } while(cin.get() != '\n');
	snpinfo.clear();
	if(f.fail()) printf("Cannot open file %s\n", bedfile), exit(0);
	if(bST >= 0) for(i = 0; i < bST; i++) getline(f, tmp);
	for(L = 0; getline(f, tmp); L++)
	{	if(bED >= 0 && L >= bED - bST) break;
		SNPINFO tsnp;
		istringstream buf(tmp);
                istream_iterator<string> beg(buf), end;
                vector<string> tokens(beg, end);
		if((int)tokens.size() > 3) tsnp.snpid = tokens[3];
		else
		{	char nd[100];
			sprintf(nd, "%d", L);
			tsnp.snpid = nd;
		}
		tsnp.chr = tokens[0];
		tsnp.posst = atoi(tokens[1].c_str());
		tsnp.posed = atoi(tokens[2].c_str());
		//tsnp.qual = -1;
		//tsnp.label = 1;
		snpinfo.push_back(tsnp);
//printf("%d %d %d %s\n", tsnp.chr, tsnp.posst, tsnp.posed, tsnp.snpid.c_str());
	}
	f.close();
//do{cout<<'\n'<<"a) Press a key to continue ..."; } while(cin.get() != '\n');

	vector<string> sampleid, factorid, ids, filenames;
	f.open(listfile);
	if(f.fail()) printf("Cannot open file %s\n", listfile), exit(0);
	for(n = 0; getline(f, tmp); n++)
	{	istringstream buf(tmp);
                istream_iterator<string> beg(buf), end;
                vector<string> tokens(beg, end);
		ids.push_back(tokens[0]+"."+tokens[1]);
		sampleid.push_back(tokens[0]);
		factorid.push_back(tokens[1]);
		filenames.push_back(tokens[2]);
//printf("%s\n", tokens[2].c_str());
	}
	f.close();

	sid.clear(); fid.clear();
	for(i = 0; i < n; i++)
	{	for(j = 0; j < (int)sid.size(); j++) if(sid[j] == sampleid[i]) break;
		if(j >= (int)sid.size()) sid.push_back(sampleid[i]);
	}	
	for(i = 0; i < (int)factorid.size(); i++)
	{	for(j = 0; j < (int)fid.size(); j++) if(fid[j] == factorid[i]) break;
		if(j >= (int)fid.size()) fid.push_back(factorid[i]);
	}	
	indN = (int)sid.size();
	maxmsz = (int)fid.size();
	vector<vector<int> > repn(indN, vector<int>(maxmsz, 0)), trepn = repn;
	for(i = 0; i < n; i++)
	{	for(j = 0; j < (int)sid.size(); j++) if(sid[j] == sampleid[i]) break;
		for(k = 0; k < (int)fid.size(); k++) if(fid[k] == factorid[i]) break;
		repn[j][k]++;
	}
	indIndex = new int[indN + 1];
	indIndex[0] = 0;
	for(i = 0; i < indN; i++)
	{	k = 0;
		for(j = 0; j < maxmsz; j++) k = max(k, repn[i][j]);
		indIndex[i + 1] = indIndex[i] + k;
	}
	
	totalN = indIndex[indN];
	msz = new int[totalN];
	map = new int*[totalN];
	for(i = 0; i < indN; i++)
	{	for(j = indIndex[i]; j < indIndex[i + 1]; j++)
		{	map[j] = new int[maxmsz];
			k = msz[j] = 0;
			for(l = 0; l < maxmsz; l++) 
				if(repn[i][l] > j - indIndex[i]) 
				{	msz[j]++;
					map[j][k++] = l;
				}
		}
	}
	
	vector<vector<int> > samplemap;
	for(i = 0; i < n; i++)
	{	int ss, ff, tt;
		for(ss = 0; ss < indN; ss++) if(sampleid[i] == sid[ss]) break;
		for(ff = 0; ff < maxmsz; ff++) if(factorid[i] == fid[ff]) break;
		vector<int> row(2, 0);
		row[0] = indIndex[ss] + trepn[ss][ff];
		for(tt = 0; tt < msz[row[0]]; tt++) if(map[row[0]][tt] == ff) break;
		row[1] = tt;//tn[indIndex[ss] + trepn[ss][ff]];
		samplemap.push_back(row); 
		trepn[ss][ff]++;
	}
	vector<string> newsid(totalN);
	for(i = 0; i < indN; i++)
	{	if(indIndex[i + 1] - indIndex[i] > 1)
		{	char tstr[1000];
			for(j = indIndex[i]; j < indIndex[i + 1]; j++)
			{	sprintf(tstr, "%s.%d", sid[i].c_str(), j - indIndex[i]); 
				newsid[j] = tstr;
			}
		}
		else newsid[indIndex[i]] = sid[i];
	}
	sid = newsid;
/*
for(i = 0; i < (int)samplemap.size(); i++)
{	printf("%d: %d,%d\n", i, samplemap[i][0], samplemap[i][1]);
}
printf("maxymsz=%d\n", maxmsz);
for(i=0;i<(int)fid.size();i++) printf("%s ", fid[i].c_str());printf("\n");
for(i=0;i<totalN;i++)
{	printf("%d: %d ", i, msz[i]);
	for(j = 0; j < msz[i]; j++) printf("%d,",map[i][j]);
	printf("\n");
}
for(i=0;i<indN;i++)
{	printf("%d: ", i);
	for(j=0;j<maxmsz;j++) printf("%d,",repn[i][j]);
	printf("\n");
}
*/
//double mean = 0, sd = 0;

//do{cout<<'\n'<<"b) Press a key to continue ..."; } while(cin.get() != '\n');

	datamatrix = new float*[totalN];
	for(i = 0; i < totalN; i++) datamatrix[i] = NULL;
	for(i = 0; i < n; i++) 
	{	int mi = samplemap[i][0], mj = samplemap[i][1];
		if(datamatrix[mi] == NULL) 
		{	datamatrix[mi] = new float[L * msz[mi]];
			printf("%d: %d x %d\n", i, L, msz[mi]);fflush(stdout);
		}
		float *dp = datamatrix[mi];
		l = (int)filenames[i].size() - 1;
		bool gz = false;
		int ccc;
		string tfile;
		while(l > 0 && filenames[i][l] != '.') l--;
		if(filenames[i][l] == '.' && filenames[i][l + 1] == 'g' && filenames[i][l + 2] == 'z')
		{	char tmp[100];
			struct stat buffer;
			do{ 
				ccc = (int)(gsl_ran_flat(gammar, 0, 1.) * 100000000.);
				sprintf(tmp, "tmp%d", ccc);
			} while(stat(tmp, &buffer)==0);
			sprintf(tmp, "gunzip -cf %s > tmp%d", filenames[i].c_str(), ccc);
			system(tmp);
			gz = true;
			sprintf(tmp, "tmp%d", ccc);
			tfile = filenames[i];
			filenames[i] = tmp;//filenames[i].substr(0, l);
		}
	
		f.open(filenames[i].c_str());
		if(f.fail()) printf("Failed to open %s\n", filenames[i].c_str());
		if(bST >= 0) for(j = 0; j < bST; j++) getline(f, tmp);
		for(j = 0; j < L; j++, dp += msz[mi])
		{	getline(f, tmp);
			double value = atof(tmp.c_str());
			//if(log2 > 0) value = fast_log(value + log2) / fast_log(2.);	
			dp[mj] = value;
//			mean += value; sd += value * value;
		}
		f.close();
		if(gz)
		{	char tmp[100];
			sprintf(tmp, "rm tmp%d", ccc);//filenames[i].c_str());
			system(tmp);
			filenames[i] = tfile;
		}
	}
/*printf("mean=%f, sd=%f, n=%d, L=%d\n", mean, sd, n, mydata.L);
	int datan = n * mydata.L;
	mean /= ((double)datan + 1e-10);
	sd = sqrt(max(1e-5, sd / ((double)datan + 1e-10) - mean * mean));
	printf("mean=%f, sd=%f\n", mean, sd);
*/
}

void genomicTensor::dataNorm()
{	
	int i, j, k;
	if(norm)
	{	vector<double> mean(maxymsz+maxxmsz, 0), sd=mean, n=mean;
		for(i = 0; i < mydata.totalN; i++)
		{	float *dyp = dataYP[i], *dxp = NULL;
			if(maxxmsz > 0) dxp = dataXP[i];
			for(j = 0; j < mydata.L; j++)
			{	for(k = 0; k < ymsz[i]; k++, dyp++)
				{	mean[ymap[i][k]] += *dyp;
					sd[ymap[i][k]] += (*dyp) * (*dyp);
					n[ymap[i][k]]++;
				}
				for(k = 0; k < xmsz[i]; k++, dxp++)
				{	mean[maxymsz + xmap[i][k]] += *dxp;
					sd[maxymsz + xmap[i][k]] += (*dxp) * (*dxp);
					n[maxymsz + xmap[i][k]]++;
				}
			}
		}
		double gm = 0, gs = 0, gn = 0;
		for(i = 0; i < maxymsz; i++)
		{	gm += mean[i];	gs += sd[i]; gn += n[i];
		}
		for(i = 0; i < maxymsz; i++)
		{	if(log2 > 0)
			{	sd[i] = 1.;
				mean[i] = mean[i] / (n[i] + 1e-100) / (gm + 1e-100) * gn;
			}
			else
			{	sd[i] = sd[i] / (n[i] + 1e-100) - mean[i] / (n[i] + 1e-100) * mean[i] / (n[i] + 1e-100);
				sd[i] /= gs / (gn + 1e-100) - gm / (gn + 1e-100) * gm / (gn + 1e-100);
				sd[i] = sqrt(sd[i]);
				mean[i] = mean[i] / (n[i] + 1e-100) - gm / (gn + 1e-100);
			}
		}
		if(maxxmsz > 0)
		{	gm = 0; gs = 0; gn = 0;
			for(i = 0; i < maxxmsz; i++)
			{	gm += mean[maxymsz + i]; gs += sd[maxymsz + i]; gn += n[maxymsz + i];
			}
			for(i = 0; i < maxxmsz; i++)
			{	if(log2 > 0)
				{	sd[maxymsz + i] = 1.;
					mean[maxymsz + i] = mean[maxymsz + i] / (n[maxymsz + i] + 1e-100) / (gm + 1e-100) * gn;
				}
				else
				{	sd[maxymsz + i] = sd[maxymsz + i] / (n[maxymsz + i] + 1e-100) - mean[maxymsz + i] / (n[maxymsz + i] + 1e-100) * mean[maxymsz + i] / (n[maxymsz + i] + 1e-100);
					sd[maxymsz + i] /= gs / (gn + 1e-100) - gm / (gn + 1e-100) * gm / (gn + 1e-100);
					sd[maxymsz + i] = sqrt(sd[maxymsz + i]);
					mean[maxymsz + i] = mean[maxymsz + i] / (n[maxymsz + i] + 1e-100) - gm / (gn + 1e-100);
				}
			}
		}

		if(log2 > 0)
		{	for(i = 0; i < mydata.totalN; i++)
			{	float *dyp = dataYP[i], *dxp = NULL;
				if(maxxmsz > 0) dxp = dataXP[i];
				for(j = 0; j < mydata.L; j++)
				{	for(k = 0; k < ymsz[i]; k++, dyp++)
						*dyp /= (mean[ymap[i][k]] + 1e-100);
					for(k = 0; k < xmsz[i]; k++, dxp++)
						*dxp /= (mean[maxymsz + xmap[i][k]] + 1e-100);
				}
			}
		}
		else
		{	for(i = 0; i < mydata.totalN; i++)
			{	float *dyp = dataYP[i], *dxp = NULL;
				if(maxxmsz > 0) dxp = dataXP[i];
				for(j = 0; j < mydata.L; j++)
				{	for(k = 0; k < ymsz[i]; k++, dyp++)
						*dyp = (*dyp - mean[ymap[i][k]]) / (sd[ymap[i][k]] + 1e-100);
					for(k = 0; k < xmsz[i]; k++, dxp++)
						*dxp = (*dxp - mean[maxymsz + xmap[i][k]]) / (sd[xmap[i][k]] + 1e-100);
				}
			}
		}
	}
	if(log2 > 0)
	{	double ln2 = log(2.);
		vector<double> mean(maxymsz+maxxmsz, 0), sd=mean, n=mean;
		for(i = 0; i < mydata.totalN; i++)
		{	float *dyp = dataYP[i], *dxp = NULL;
			if(maxxmsz > 0) dxp = dataXP[i];
			for(j = 0; j < mydata.L; j++)
			{	for(k = 0; k < ymsz[i]; k++, dyp++)
				{	(*dyp) = log(*dyp + log2) / ln2;
					if(norm)
					{	mean[ymap[i][k]] += *dyp;
						sd[ymap[i][k]] += (*dyp) * (*dyp);
						n[ymap[i][k]]++;
					}
				}
				for(k = 0; k < xmsz[i]; k++, dxp++)
				{	(*dxp) = log(*dxp + log2) / ln2;
					if(norm)
					{	mean[maxymsz + xmap[i][k]] += *dxp;
						sd[maxymsz + xmap[i][k]] += (*dxp) * (*dxp);
						n[maxymsz + xmap[i][k]]++;
					}
				}
			}
		}
		if(norm)
		{	double gm = 0, gs = 0, gn = 0;
			for(i = 0; i < maxymsz; i++)
			{	gm += mean[i];	gs += sd[i]; gn += n[i];
			}
			for(i = 0; i < maxymsz; i++)
			{	sd[i] = sd[i] / (n[i] + 1e-100) - mean[i] / (n[i] + 1e-100) * mean[i] / (n[i] + 1e-100);
				sd[i] /= gs / (gn + 1e-100) - gm / (gn + 1e-100) * gm / (gn + 1e-100);
				sd[i] = sqrt(sd[i]);
				mean[i] = mean[i] / (n[i] + 1e-100) - gm / (gn + 1e-100);
			}
			if(maxxmsz > 0)
			{	gm = 0; gs = 0; gn = 0;
				for(i = 0; i < maxxmsz; i++)
				{	gm += mean[maxymsz + i]; gs += sd[maxymsz + i]; gn += n[maxymsz + i];
				}
				for(i = 0; i < maxxmsz; i++)
				{	sd[maxymsz + i] = sd[maxymsz + i] / (n[maxymsz + i] + 1e-100) - mean[maxymsz + i] / (n[maxymsz + i] + 1e-100) * mean[maxymsz + i] / (n[maxymsz + i] + 1e-100);
					sd[maxymsz + i] /= gs / (gn + 1e-100) - gm / (gn + 1e-100) * gm / (gn + 1e-100);
					sd[maxymsz + i] = sqrt(sd[maxymsz + i]);
					mean[maxymsz + i] = mean[maxymsz + i] / (n[maxymsz + i] + 1e-100) - gm / (gn + 1e-100);
				}
			}
			for(i = 0; i < mydata.totalN; i++)
			{	float *dyp = dataYP[i], *dxp = NULL;
				if(maxxmsz > 0) dxp = dataXP[i];
				for(j = 0; j < mydata.L; j++)
				{	for(k = 0; k < ymsz[i]; k++, dyp++)
						*dyp = (*dyp - mean[ymap[i][k]]) / (sd[ymap[i][k]] + 1e-100);
					for(k = 0; k < xmsz[i]; k++, dxp++)
						*dxp = (*dxp - mean[maxymsz + xmap[i][k]]) / (sd[xmap[i][k]] + 1e-100);
				}
			}
		}
	}

	mylik.markmean.resize(maxymsz + maxxmsz, 0);
	mylik.marksd.resize(maxymsz + maxxmsz, 0.001);
	vector<double> n(maxymsz + maxxmsz, 1);
	for(i = 0; i < mydata.totalN; i++)
	{	float *dyp = dataYP[i], *dxp = NULL;
		if(maxxmsz > 0) dxp = dataXP[i];
		for(j = 0; j < mydata.L; j++)
		{	for(k = 0; k < ymsz[i]; k++, dyp++)
			{	mylik.markmean[ymap[i][k]] += *dyp;
				mylik.marksd[ymap[i][k]] += (*dyp) * (*dyp);
				n[ymap[i][k]]++;
			}
	            	if(maxxmsz > 0)
        	    	{   	for(k = 0; k < xmsz[i]; k++, dxp++)
	               		{	mylik.markmean[maxymsz + xmap[i][k]] += *dxp;
        	            		mylik.marksd[maxymsz + xmap[i][k]] += (*dxp) * (*dxp);
	                    		n[maxymsz + xmap[i][k]]++;
        	        	}
            		}
		}
	}

	overallmean = 0; overallsd = 0; 
	double gn = 0;
	for(i = 0; i < maxymsz; i++)
	{	overallmean += mylik.markmean[i];	overallsd += mylik.marksd[i]; gn += n[i];
	}
	overallmean /= (gn + 1e-100);
	overallsd = sqrt(overallsd / (gn + 1e-100) - overallmean * overallmean);
	overallmeanx = 0; overallsdx = 0; 
	gn = 0;
	for(i = 0; i < maxxmsz; i++)
	{	overallmeanx += mylik.markmean[maxymsz + i]; overallsdx += mylik.marksd[maxymsz + i]; gn += n[maxymsz + i];
	}
	overallmeanx /= (gn + 1e-100);
	overallsdx = sqrt(overallsdx / (gn + 1e-100) - overallmeanx * overallmeanx);
	for(i = 0; i < maxymsz + maxxmsz; i++)
	{	mylik.markmean[i] /= n[i];
		mylik.marksd[i] = max(0.0001, sqrt(mylik.marksd[i] / n[i] - mylik.markmean[i] * mylik.markmean[i]));
	}
}

void genomicTensor::groupData(float *data, int index, int msz, double step, vector<int> const &inputlist, vector<vector<double> > &values, vector<vector<int> > &outputlist)
{	int i, j, k;
	vector<MySortType> tv((int)inputlist.size());
	for(i = 0; i < (int)inputlist.size(); i++)
	{	tv[i].score = data[inputlist[i] * msz + index];
		tv[i].index = i;
	}
	sort(tv.begin(), tv.end());
	
	values.clear();
	outputlist.clear();
	double stv = (double)((int)floor(tv[0].score / step)) * step;
	j = 0;
	for(i = 1; i < (int)tv.size(); i++) 
	{	if(tv[i].score >= stv + step)
		{	vector<int> row(i - j);
			for(k = j; k < i; k++) row[k - j] = inputlist[tv[k].index];
			outputlist.push_back(row);
			values.push_back(vector<double>(1, stv + step / 2.));
			j = i;
			stv = (double)((int)floor(tv[i].score / step)) * step;
		}
	}
	vector<int> row(i - j);
	for(k = j; k < i; k++) row[k - j] = inputlist[tv[k].index];
	outputlist.push_back(row);
	values.push_back(vector<double>(1, stv + step / 2.));

	if(index + 1 < msz)
	{	vector<vector<double> > tvalue, newvalue;
		vector<vector<int> > tout, newout;
		for(i = 0; i < (int)outputlist.size(); i++)
		{	groupData(data, index + 1, msz, step, outputlist[i], tvalue, tout);
			for(j = 0; j < (int)tvalue.size(); j++)
			{	tvalue[j].insert(tvalue[j].begin(), values[i][0]);
			}
			newvalue.insert(newvalue.end(), tvalue.begin(), tvalue.end());
			newout.insert(newout.end(), tout.begin(), tout.end());
		}
		values = newvalue;
		outputlist = newout;
	}
}

double genomicTensor::_imputeData(MYDATA const &mydata, TENSORPARA const &tpara, int id, vector<int> const &dbreaks, bool init)
{	
//	printf("my impute data\n");fflush(stdout);
	int i, j, k, idst = mydata.indIndex[id / mydata.ploidity], ided = mydata.indIndex[id / mydata.ploidity + 1];;
	int const *pp, *pp2 = NULL, *aszp;
	float *dp, *dp2 = NULL;
	float const *nhp;

	bool ohpass = hpass;
	hpass &= (!init);

    pp = tpara.pop + id * mydata.L;
	dp = mydata.data + mydata.indIndex[id/mydata.ploidity] * mydata.ploidity * mydata.L;
	if(mydata.ploidity == 2) 
	{	pp2 = pp + mydata.L;
		dp2 = dp + mydata.L;
	}
	vector<vector<PARAUNIT> >::const_iterator ipm = tpara.param.begin();
	nhp = tpara.nh;
	aszp = mydata.asz;

	float *odata = new float[mydata.L * (ided - idst)];
	float *nlapse = new float[mylik.clustersz * 2];
//clock_t cst, ced;
//cst=clock();
	if(hpass) 
	{	for(j = 0; j < mylik.clustersz; j++) 
		{	nlapse[j*2] = 100.; nlapse[j*2+1] = 0.95 * nlapse[j*2]; 	}
		int step = (ided - idst) * mylik.clustersz * 2;
		if(tmpSpace2 != NULL) delete tmpSpace2;
		int tsz = (step + ided - idst) * mydata.L;
		tpointer = tmpSpace2 = new double[tsz]; 
		for(j = 0; j < tsz; j++) tmpSpace2[j] = 0;
		dpointer = tmpSpace2 + step * mydata.L;
		bool *tlapse = tpara.lapse + idst * mydata.L;
		float prevd[ided-idst], *ddp = mydata.data + idst * mydata.L;
		for(j = 0; j < mydata.L; j++, ddp++, tlapse++)
			for(int kk = 0; kk < ided - idst; kk++, dpointer++)
			{	if(!(*(tlapse + kk * mydata.L))) prevd[kk] = *(ddp + kk * mydata.L);
				if(gID > tpara.burnin / 3)
				{
				if(j > 0 && *(ddp+kk*mydata.L) == *(ddp+kk*mydata.L-1)) nlapse[(int)prevd[kk]*2+1]+=tpara.indWeight[kk];
				if(j > 0) nlapse[(int)prevd[kk]*2]+=tpara.indWeight[kk];
				}
				*dpointer = prevd[kk];
			}
//for(j=0;j<mylik.clustersz;j++)
//printf("%d: %d %d = %f\n", j, (int)nlapse[j*2],(int)nlapse[j*2+1],nlapse[j*2+1]/nlapse[j*2]);

		dpointer = tmpSpace2 + step * mydata.L;
		
		tlapse = tpara.lapse + idst * mydata.L;
		k = 0;
		for(i = idst; i < ided; i++)
		{	float prevd = dp[k];
			for(j = 0; j < mydata.L; j++, k++, tlapse++)
			{	if(!(*tlapse)) prevd = dp[k];
				odata[k] = prevd;
			}
		}
	}
	else
	{	k = 0;
		for(i = idst; i < ided; i++) for(j = 0; j < mydata.L; j++, k++) odata[k] = dp[k];
	}


	double RT = 0;
	int nstep = max(1000,(int)ceil((double)mydata.L / (double)nthread));
	vector<int> region(1,0);
	j = 0;
	for(i = 1; i < mydata.L; i++) if(mydata.snpinfo[i].chr != mydata.snpinfo[j].chr) { region.push_back(i); j = i; }
	region.push_back(i);
	vector<int> bkregion = region;
	region.clear();
	region.push_back(0);
	bool *tlapse = tpara.lapse + idst * mydata.L;
	for(i = 1; i < (int)bkregion.size(); i++)
	{	int tnsz = max(1, (int)ceil((double)(bkregion[i] - bkregion[i - 1]) / nstep));
		if(tnsz > 1)
		{	int tnstep = (int)ceil((double)(bkregion[i] - bkregion[i - 1]) / (double)tnsz);
			int ok = bkregion[i - 1];
			for(j = 1; j < tnsz; j++)
			{	k = bkregion[i - 1] + tnstep * j;
				if(hpass)
				{	int lk, rk;
					for(lk = k; lk > max(ok, k - tnstep / 2); lk--)
						if(tlapse[lk] == 0) break;
					for(rk = k + 1; rk < min(bkregion[i], k + tnstep / 2); rk++)
						if(tlapse[rk] == 0) break;
					if(k - lk < rk - k && lk > max(ok, k - tnstep / 2)) k = lk;
					else if(rk - k < k - lk && rk < min(bkregion[i], k + tnstep / 2)) k = rk;
				}
				if(k > region[(int)region.size() - 1]) 
				{	region.push_back(k);
					ok = k;
				}
			}
		}
		region.push_back(bkregion[i]);
	}
	//for(i = 0; i < (int)region.size(); i++) printf("%d: %d\n", i, region[i]); fflush(stdout);


# pragma omp parallel num_threads(nthread)
{

# pragma omp for
    
	for(int nn = 1; nn < (int)region.size(); nn++)
	{	

		int posst = region[nn - 1], posed = region[nn];
		int flip = 0;
		int j;
		for(j = posst; j < posed; j++)
		{	if((int)tpara.flips.size() > 0) flip = tpara.flips[id/mydata.ploidity][j*2+1];
			if((int)mydata.fixAllele.size() > 0 && (int)mydata.fixAllele[(id)].size() > 0 && !mydata.fixAllele[(id)][j]) ;
			else 
			{	
				RT += _imputeOne(id, vector<int>(), j, j+1, flip, mydata.ploidity, tpara, dp + j, dp2 + j, pp + j, pp2 + j, ipm + j, nhp + tpara.maxK * j, aszp + j);
			}
		}
		if(hpass)
		{	int step = (ided - idst) * mylik.clustersz * 2;
			double *ttpointer = tpointer + posst * step, *tdpointer = dpointer + posst * (ided - idst);
/*if(gID>50)
{	for(int jj = 0; jj < mylik.clustersz; jj++)
	{	printf("{%d:%f} ", jj, tpointer[mylik.clustersz*2*99979+mylik.clustersz+jj]);
	}
	printf("\n");fflush(stdout);
}*/
			for(int pos = posst; pos < posed; pos++)
			{	for(int idd = idst; idd < ided; idd++)
				{	double lapsep = 0.95 * (double)(pos > 0), maxlp = MINUSINFINITE, tsum = 0, *ppointer = ttpointer - step;
					if(pos > 0)
					{	int tg = (int)(*(dp+pos-1));
						lapsep = nlapse[tg*2+1] / nlapse[tg*2];
					}
					if(pp[pos] < tpara.maxK) lapsep = (ipm[pos][pp[pos]].lapseN + lapsep * 10.) / (ipm[pos][pp[pos]].popN + 10.);
lapsep *= tpara.indWeight[idd];
//lapsep=min(lapsep,0.1);//min(0.95,pow((double)gID / (double)tpara.burnin, 2.)));
/*if(pos>99970) 
{	for(int jj = 0; jj < mylik.clustersz; jj++)
	{	printf("(%d)%f ",jj,ttpointer[jj]);fflush(stdout);
	}
	printf(" |%f| ", lapsep);
}*/
					if(pos > 0)
                    			{	for(j = 0; j < mylik.clustersz; j++)
							ttpointer[mylik.clustersz + j] = ttpointer[mylik.clustersz + j] * (1. - lapsep) + ppointer[j] * lapsep;
					}
					for(j = 0; j < mylik.clustersz; j++) maxlp = max(maxlp, ttpointer[j]);
					for(j = 0; j < mylik.clustersz; j++) 
						if(ttpointer[j] > MINUSINFINITE)
						{	ttpointer[j] = exp(ttpointer[j] - maxlp) * ttpointer[mylik.clustersz + j];
							tsum += ttpointer[j];	
						}
						else ttpointer[j] = 0;
					for(j = 0; j < mylik.clustersz; j++) 
					{	ttpointer[j] /= tsum;
				//		if(pos>99970) printf("[%d]%f ",j,ttpointer[j]);fflush(stdout);
					}
//if(pos>99970) printf("| %d\n", pos);
					ttpointer += mylik.clustersz * 2;
					tdpointer++;
				}
			}
/*if(gID>50)
{	for(int jj = 0; jj < mylik.clustersz; jj++)
	{	printf("[%d:%f] ", jj, tpointer[mylik.clustersz*2*99979+mylik.clustersz+jj]);
	}
	printf("\n");fflush(stdout);
}*/
	

			RT += _updateLapseData(id, tmpSpace2, posst, posed, nlapse);
		}
	}
}

	if(hpass)
	{	delete tmpSpace2;
		tmpSpace2 = NULL;
	}

	k = 0;
	for(i = idst; i < ided; i++)
	{	float *dpy, *dpx = NULL;
		dpy = dataYP[i];
		if(maxxmsz > 0) dpx = dataXP[i];
		if(hpass)
		{	float prevd = dp[k];
			bool *tlapse = tpara.lapse + i * mydata.L;
			for(j = 0; j < mydata.L; j++, k++, dpy += ymsz[i], dpx += xmsz[i], tlapse++)
			{	if(!(*tlapse)) prevd = dp[k];
				if(odata[k] != prevd)
				{	if(odata[k] >= 0) mylik.removePara(odata[k], dpy, dpx, i, tpara.indWeight[i]);
					mylik.addPara(prevd, dpy, dpx, i, tpara.indWeight[i]);
				}
			}
		}
		else
		{	for(j = 0; j < mydata.L; j++, k++, dpy += ymsz[i], dpx += xmsz[i])
				if(odata[k] != dp[k])
				{	if(odata[k] >= 0) mylik.removePara(odata[k], dpy, dpx, i, tpara.indWeight[i]);
					mylik.addPara(dp[k], dpy, dpx, i, tpara.indWeight[i]);
				}
		}
	}	
	delete []odata;
	delete []nlapse;

	hpass = ohpass;

/*
for(int yy = 0; yy < mylik.gausssz; yy++) printf("%f, ", mylik.gausspara[0*mylik.gausssz+yy]);
printf("\n");
double pg[10000], gn[1000];
for(i = 0; i < 10000; i++) pg[i] = 0;
for(i = 0; i < 1000; i++) gn[i] = 0;
if(dp != mydata.data + id * mydata.L) printf("!!!");
dp=mydata.data + id * mydata.L;
float *dpy = dataYP[id];
for(i = 0; i < mydata.L; i++)
{	k=(int)dp[i];
	gn[k]++;
	pg[k*mylik.gausssz]++;
	for(j = 0; j < ymsz[id]; j++) pg[k*mylik.gausssz + 1 + j] += dpy[i * ymsz[id] + j];	
	int kk = 0;
	for(j = 0; j < ymsz[id]; j++) 
		for(int l = 0; l <= j; l++, kk++) 
			pg[k*mylik.gausssz + 1 + ymsz[id] + kk] += dpy[i * ymsz[id] + j] * dpy[i * ymsz[id] + l];	
}
for(int yy = 0; yy < mylik.gausssz; yy++) printf("%f, ", pg[0*mylik.gausssz+yy]);
printf("\n");
for(int yy = 0; yy < mylik.clustersz; yy++) printf("%d,",(int)mylik.gausspara[yy*mylik.gausssz]);
printf("\n");
for(int yy = 0; yy < mylik.clustersz; yy++) printf("%d,",(int)gn[yy]);
printf("\n");
exit(0);
//ced=clock();
//printf("impute of %d = %fsec\n", mydata.L, (double)(ced-cst)/1000000.);
*/
	return(RT);
}

void genomicTensor::_imputeMissingData(char const *fname, int gid)
{
	int i, j, k;
	char str[200];
	
for(int id = 0; id < mydata.totalN; id++)
{	if(maxymsz - ymsz[id] == 0) continue;
	for(i = 0; i < (int)imputelist.size(); i++)
	{	string myid = mydata.indinfo[id];
		for(j = myid.size() - 1; j >=0; j--) if(myid[j] == '.') break;
		if(j >= 0) myid = myid.substr(0, j);
		if(myid == imputelist[i]) break;
	}
	if(i >= (int)imputelist.size() && i > 0) continue;

	float *tmpdata = new float[mydata.L * (maxymsz + maxxmsz - ymsz[id] - xmsz[id]) + 1], *tmpy = tmpdata, *tmpx = tmpy + mydata.L * (maxymsz - ymsz[id]);
	mylik.imputeMissingTrack(id, 0, mydata.L, dataYP[id], dataXP[id], tmpy, tmpx, mydata.data + id * mydata.L);//tlp);

	char str[100];
	if(gid >= 0) sprintf(str, "%s.impute.%s.%d", fname, mydata.indinfo[id].c_str(),gid);
	else sprintf(str, "%s.impute.%s", fname, mydata.indinfo[id].c_str());
	FILE *f = fopen(str,"w");
	for(i = 0; i < maxymsz; i++)
	{	for(j = 0; j < ymsz[id]; j++) if(i == ymap[id][j]) break;
		if(j >= ymsz[id]) fprintf(f, "%s ", mydata.fyinfo[i].c_str());
	}
	for(i = 0; i < maxxmsz; i++)
	{	for(j = 0; j < xmsz[id]; j++) if(i == xmap[id][j]) break;
		if(j >= xmsz[id]) fprintf(f, "%s ", mydata.fxinfo[i].c_str());
	}
	fprintf(f, "\n");
	float *py = tmpy, *px = tmpx;
	for(i = 0; i < mydata.L; i++)
	{	for(j = 0; j < maxymsz - ymsz[id]; j++, py++) fprintf(f, "%5.3f ", *py);
		for(j = 0; j < maxxmsz - xmsz[id]; j++, px++) fprintf(f, "%5.3f ", *px);
		fprintf(f, "\n");
	}
	fclose(f);
	delete []tmpdata;
	if(gzip)
	{	if(gid >= 0) sprintf(str, "gzip -f %s.impute.%s.%d", fname, mydata.indinfo[id].c_str(), gid);
		else sprintf(str, "gzip -f %s.impute.%s", fname, mydata.indinfo[id].c_str());
		system(str);
	}
}
	
	mylik.imputeMissing2(mydata, tpara, dataYP, dataXP, tpara.priorW, fname, gid, tpara.indWeight);
}
