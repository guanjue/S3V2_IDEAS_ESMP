#include "tensorHMMbase.h"


clock_t time1, time2;
extern bool SA;

tensorHMMbase::tensorHMMbase(unsigned int rs)
{	

	gsl_rng_env_setup();
	T = gsl_rng_default;
	gammar = gsl_rng_alloc(T);
	unsigned int rseed = rs;
	if(rs==0) rseed = (unsigned int)time(NULL) * 3 + (unsigned int)clock();
	gsl_rng_set(gammar, rseed);
	outputproportion = false;
	
	tmpM = NULL;
}

tensorHMMbase::~tensorHMMbase(void)
{
	if(tmpM != NULL) delete tmpM;	
	gsl_rng_free(gammar);
}

void tensorHMMbase::inferStructure(MYDATA &mydata, TENSORPARA &tpara, char const *outname, bool outqc, vector<bool> const &myinit)
{	int i, j, k, l;
    //initialize
	//if(tpara.shrinkclass) _shrinkClasses(mydata.totalN, mydata.L, mydata.data);//newchange
	mydata.asz = new int[mydata.L];
	_getAsz(mydata);

//do{cout<<'\n'<<"1) Press a key to continue ..."; } while(cin.get() != '\n');

	printf("N=%d, L=%d, A=%f, ", mydata.indN / mydata.ploidity, mydata.L, tpara.probAlpha);fflush(stdout);
	if(tpara.maxHapK <= 0) tpara.maxHapK = max(tpara0.popn0, min(20, mydata.indN));

	if((int)mydata.fixState.size() > 0)
	{	k = -1;
		for(i = 0; i < (int)mydata.fixState.size(); i++) 
		{	if(mydata.fixState[i].size() > 0) tpara.trueindN -= tpara.indWeight[i];
			for(j = 0; j < (int)mydata.fixState[i].size(); j += 2) 
				k = max(k, mydata.fixState[i][j + 1]);
		}
		tpara.maxHapK = max(tpara.maxHapK, k + 1);
	}

	if(tpara.rPrior == NULL)
	{	tpara.rPrior = new double[mydata.L];
		tpara.rPrior[0] = 0.;
		for(i = 1; i < mydata.L; i++)
		{	tpara.rPrior[i] = 0.01;
			if((int)mydata.snpinfo.size() > 0)
			{	if(mydata.snpinfo[i].chr == mydata.snpinfo[i - 1].chr && mydata.snpinfo[i].posst >= mydata.snpinfo[i-1].posst) 
					tpara.rPrior[i] = 0.99 * (1. - exp(- tpara.recrate*1e-6 * (double)max(1, mydata.snpinfo[i].posst - mydata.snpinfo[i - 1].posst)));
				else tpara.rPrior[i] = 1.;
			}
		}
	}

	if(mydata.ploidity == 2 && mydata.breaks.size() > 0)
        {       tpara.flips.clear();
                tpara.flips.resize(mydata.indN / mydata.ploidity);
        }

	j = mydata.indN * mydata.L;
	if(tpara.pop == NULL)
	{	tpara.pop = new int[j];
		int *pp = tpara.pop;
		for(i = 0; i < j; i++, pp++) (*pp) = 0;
	}
	if(tpara.rec == NULL)
	{	tpara.rec = new bool[j];
		bool *rr = tpara.rec;
		for(i = 0; i < j; i++, rr++) (*rr) = false;
	}
	if(tpara.lapse == NULL)
	{	tpara.lapse = new bool[mydata.totalN * mydata.L];
		bool *ll = tpara.lapse;
		for(i = 0; i < mydata.totalN * mydata.L; i++, ll++) *ll = false;
	}

	tpara.popr = new double[mydata.L];
	for(i = 0; i < mydata.L; i++) tpara.popr[i] = tpara.rPrior[i];

//do{cout<<'\n'<<"2) Press a key to continue ..."; } while(cin.get() != '\n');

	tpara.addsample= (tpara.maxK == 0);
	if(tpara.maxK == 0)
	{	double sum = 0, t1, t2 = 1.;
		while(sum < 0.99)
		{       t1 = max(0.1, gsl_ran_beta(gammar, 1., tpara.A));
			sum += t1 * t2;
			t2 *= (1. - t1);
			tpara.maxK++;
		}
		if(tpara.maxK > min(20, mydata.indN)) tpara.maxK = min(20, mydata.indN);
//		tpara.maxK = 1;
	}
	int maxK0 = tpara.maxK;
	printf("K=%d\n", tpara.maxK);fflush(stdout);

	if(tpara.nh == NULL)
	{	j = tpara.maxK * mydata.L;
		tpara.nh = new float[j];
		for(i = 0; i < j; i++) tpara.nh[i] = 0;
		tpara.param.clear();
		PARAUNIT pu;
		pu.lapseN = pu.popN = 0;
		tpara.param.resize(mydata.L, vector<PARAUNIT>(tpara.maxK, pu));
	//	for(i = 0; i < mydata.L; i++) 
	//		tpara.param[i].resize(tpara.maxK, pu);
	}
//do{cout<<'\n'<<"3) Press a key to continue ..."; } while(cin.get() != '\n');

	//double *Q = new double[1000], *Vh = new double[1000];
	//tpara.Q = Q; tpara.Vh = Vh;
	double Q[1000], Vh[1000];
	tpara.Q = &Q[0]; tpara.Vh = &Vh[0]; 
	for(i = 0; i < tpara.maxK; i++) 
	{	tpara.Q[i] = 0; tpara.Vh[i] = 1. / (double)tpara.maxK;	}
	
	tpara.recombN = new float[mydata.L];
	for(i = 0; i < mydata.L; i++) tpara.recombN[i] = 0;

	Temp = 1.;

	//sampling
	vector<vector<char> > fullsample; //used to store full info about data
	vector<vector<int> > cumsample; //used to store cumulative frequency of data
	vector<vector<int> > drecsample; //used to store double recombination sites
	vector<vector<vector<int> > > popsample; //used to store state info
	double *recombs = new double[mydata.L]; //sample recomb rates
	for(i = 0; i < mydata.L; i++) recombs[i] = 0;
	if(tpara.samplefull)//imputedata)
	{	j = 0;
		for(i = 0; i < mydata.L; i++)
			j += (int)pow((double)mydata.asz[i], (double)mydata.ploidity);
		cumsample.resize(mydata.totalN, vector<int>(j, 0));
	}

	double lp;
	clock_t tst, ted; time1=time2=0;
	imputeT = forwardT = backwardT = 0;

	int tttt = 0,otttt, indNbak = mydata.indN, mMaxK = 0;
	int MAXHAPK = tpara.maxHapK;
	double oprobalpha = tpara.probAlpha;
//do{cout<<'\n'<<"4) Press a key to continue ..."; } while(cin.get() != '\n');

	vector<bool> init(mydata.indN / mydata.ploidity, true);
	if((int)myinit.size() > 0) init = myinit; 

	clock_t sttimer = clock();
    for(i = 0; i < (tpara.burnin + tpara.mcmc); i++)
	{	int curindN = indNbak;
        //emit ideas.on_progressBar_valueChanged(i);

		if(MAXHAPK > 0) tpara.maxHapK = max(MAXHAPK, min(min(20, mydata.indN), min(MAXHAPK+5,MAXHAPK * 2) * (tpara.burnin - i * 5 / 4)/ tpara.burnin));//newchange
		//if(i < tpara.burnin / 2) tpara.maxHapK = maxK0; 

		if(tpara.addsample)
		{	otttt=tttt;
			int ssss = max(10, (int)(sqrt((double)curindN) + 1));//checkthis
			ssss+=ssss%2;
			tttt+=ssss;
			curindN = min(tttt, curindN);
			if(otttt < indNbak) 
			{ //if(otttt>0) init.resize(otttt/mydata.ploidity);
			  //else init.clear();
			  //init.resize(curindN / mydata.ploidity, i==0);
			  if(i > tpara.burnin / 3) i = 0;
			}
		}
//for(j = 0; j < curindN; j++) printf("[%d]",(int)init[j]);printf("|i=%d\n", i);fflush(stdout);

		gID = i;
		if(SA) Temp = 1.+4.*max(0,tpara.burnin/5-gID)*5./(double)tpara.burnin;
		if(tpara.samplemaximum)
		{	if(i >= tpara.burnin * 4 / 5) tpara.maximum = 3;
			//else tpara.maximum = 3 * (i % 2);
			//if(i == tpara.burnin + tpara.maximum - 1) tpara.maximum = 7;
		}

		if(i > 0)
		{	if(tpara.changeAlpha > 0) 
				tpara.probAlpha = oprobalpha + (double)max(0, tpara.burnin / 3 - i) / ((double)tpara.burnin / 3.) * min(5., curindN /20.);
			else if(tpara.changeAlpha < 0) 
				tpara.probAlpha = max(0.1,oprobalpha * min((double)i, (double)tpara.burnin / 2.) / ((double)tpara.burnin / 2.));
		}
		if(i < tpara.burnin) printf("[%d] itern%d : ", curindN, i);
		else printf("[%d] itern%dS: ", curindN, i);
		fflush(stdout);

		//if(i>0)
		{	for(j = 0; j < tpara.maxK; j++) tpara.Q[j] = 0;
			for(j = 0; j < mydata.L; j++) tpara.recombN[j] = 0;
			for(k = 0; k < curindN; k++) //include fixed
				if(!init[k / mydata.ploidity])
				{	for(int l = 0; l < mydata.L; l++)
					{       if(l==0 || tpara.rec[k*mydata.L+l]) 
						{	tpara.Q[tpara.pop[k*mydata.L+l]]+=tpara.idw[k];
							tpara.recombN[l]+=tpara.idw[k];
						}
					}
				}
		}
		float ttrueN = 0;
		for(int ii = 0; ii < curindN; ii++) ttrueN += tpara.idw[ii];
		mydata.indN = curindN;
		mydata.totalN = mydata.indIndex[mydata.indN];
tst=clock();

		lp = _oneRound(mydata, tpara, init, i);
//printf("<");for(j = 0; j < 6; j++) printf("%d ", (int)mydata.data[j*mydata.L+483]);printf(">\n");fflush(stdout);
ted=clock();
time1+=ted-tst;
		_collapseGroups_new(mydata, tpara);

		//_checkCounts(tpara);

		if(i < tpara.burnin)  //split
		{	int sz = max(4, tpara.burnin / 2 / 10 + 1);
			if(i < tpara.burnin / 2 && !tpara.singlecall && i % sz == (sz-1) && (tpara.maxHapK == 0 || tpara.maxHapK >= tpara.maxK * 2) && tpara.maxK <= min(tpara.splitstop, curindN / 10))//09132013, was /2 for the last term
			{	int minstate = tpara.maxK;
				for(j = 0; j < curindN; j += mydata.ploidity)
				{	if((int)mydata.fixState.size() > 0 && mydata.fixState[j].size() > 0) continue;
					//if(gsl_ran_flat(gammar, 0, 1.) < 0.5) continue;
					_removeInd(mydata, tpara, j);
					for(k = 0; k < mydata.ploidity; k++)
					{	int *pp = tpara.pop + (j + k) * mydata.L, ns = *pp;
						bool *rr = tpara.rec + (j + k) * mydata.L;
						for(int l = 0; l < mydata.L; l++, pp++, rr++) 
						{	if(l == 0 || (*rr))
							{	ns = (*pp);
								if(gsl_ran_flat(gammar, 0, 1.) > 0.5) ns += minstate;
							}
							(*pp) = ns;
						}
					}
					_addInd(mydata, tpara, j);
				}
			}
		}

		_updateVh(tpara.Q, tpara.maxK, tpara.A, tpara.Vh, false);
//for(int ii = 0; ii < tpara.maxK; ii++) printf("<%f>",tpara.Q[ii]);printf("\n");fflush(stdout);
//		ted=clock();
//		time1+=ted-tst;

		tst=clock();
		_updatePrior(init, true);
		//lp = _logP(mydata, tpara);
                printf("K=%d ", tpara.maxK/*, lp*/);fflush(stdout);
		ted=clock();
		time2+=ted-tst;	

		mydata.indN = indNbak;

		double opw = tpara.priorW;
		tpara.priorW *= max(1., 0.5*(double)mydata.indN*(double)(tpara.burnin/2-i)/tpara.burnin);
              	int totaln = _updateRecomb(tpara, min(ttrueN, tpara.trueindN), mydata.L);
		tpara.priorW = opw;

		printf("R=%d ", totaln);fflush(stdout);
//if(gID < tpara.burnin)		
		_greedySwitch_new(curindN, mydata.L, mydata.ploidity, tpara, init);

		//_checkCounts(tpara);
	
                //update A 
		//tpara.A = _updateA(tpara.maxK, tpara.Q, tpara.A); printf("A=%f ", tpara.A);

		//sampling
		if(i >= tpara.burnin)
		{	if(mMaxK < tpara.maxK) mMaxK = tpara.maxK;
			//sample data
			if(tpara.samplefull)//imputedata)
			{	int x = 0;
				vector<vector<char> > tmphap;
				if(mydata.ploidity == 2) tmphap.resize(mydata.totalN * mydata.ploidity, vector<char>((mydata.L + 7) / 8, 0));
				float *ddp = mydata.data;
				float prevdata[mydata.totalN];
				bool *tlapse[mydata.totalN];
				for(j = 0; j < mydata.totalN; j++) tlapse[j] = tpara.lapse + j * mydata.L;
				for(j = 0; j < mydata.L; j++, x += mydata.ploidity + 2, ddp++)
				{	int t = j / 8, shift = 7 - j % 8, dddsz=mydata.L*mydata.ploidity;
					float *dddp = ddp;
					for(k = 0; k < mydata.totalN; k++, dddp += dddsz, tlapse[k]++)
					{	if(mydata.ploidity == 1) 
						{	if(!(*tlapse[k])) prevdata[k] = *dddp;
							cumsample[k][x + (int)prevdata[k]]++;
						}
        	                        	else 
						{	
							cumsample[k][x + (int)(*dddp) * mydata.asz[j] + (int)(*(dddp+mydata.L))] ++;
							if(tpara.samplefull && mydata.ploidity == 2)
							{	tmphap[k * 2][t] |= (((int)(*dddp)) << shift);
								tmphap[k * 2 + 1][t] |= (((int)(*(dddp+mydata.L))) << shift);
							}
						}
        	                	}
				}
				if((int)tmphap.size() > 0)
				{	if((int)fullsample.size() >= mydata.totalN * mydata.ploidity * 50) fullsample.erase(fullsample.begin(), fullsample.begin() + mydata.totalN * mydata.ploidity);
					fullsample.insert(fullsample.end(), tmphap.begin(), tmphap.end());
				}
			}

			//sample pops
			int *tp[1];
			vector<vector<int> > sp, ttsp;
			for(j = 0; j < mydata.indN; j++)
			{	tp[0] = tpara.pop + j * mydata.L;
				_LongtoShort(tp, tpara.rec + j * mydata.L, ttsp, 1, mydata.L);
				sp.push_back(ttsp[0]);
			}
			popsample.push_back(sp);
			
			//sample recombs
			int x = 0;
			for(j = 0; j < mydata.L; j++, x += mydata.ploidity + 2) recombs[j] += tpara.popr[j];

			//sample double rec sites
			if(mydata.ploidity == 2)
			{	for(l = 0; l < mydata.indN / mydata.ploidity; l++)
				{	vector<int> tmpsamplebreaks(1, mydata.L - 1);
					bool *r1 = tpara.rec + l * mydata.ploidity * mydata.L + mydata.L - 1, *r2 = r1 + mydata.L;
					for(j = mydata.L - 1; j > 0; j--, r1--, r2--)
						if(*r1 && *r2) tmpsamplebreaks.push_back(j-1);
					drecsample.push_back(tmpsamplebreaks);
				}
			}
		}
//		FILE *ft;
//		if(i==0) ft=fopen("tmp.out", "w");
//		else ft=fopen("tmp.out", "a");
//		fprintf(ft, "%d: CSZ=%d, K=%d, dp=%f, %3.2f\n", i, mydata.asz[0], tpara.maxK, lp, (double)(clock()-sttimer)/1000000.);
//		fclose(ft);
		printf("%3.2f sec\n", (double)(clock()-sttimer)/1000000.);
	}

//	printf("Foward time = %f sec\nBackward time = %f sec\nImpute time = %f sec\n", (double)forwardT / 1000000., (double)backwardT / 1000000., (double)imputeT/1000000.);
	//printf("one round time=%f sec\tupdate para time=%f sec\n", (double)time1/1000000., (double)time2/1000000.);

	printf("Summarize posteriors...");fflush(stdout);	
	for(i = 0; i < mydata.L; i++) recombs[i] /= (double)tpara.mcmc;
	
	char str[200];
	FILE *f;

	if(!tpara.samplefull) 
	{	
		if(mydata.ploidity == 1) 
		{	vector<int> pqc;
			mMaxK = max(mMaxK, tpara.maxK);
		//	tpara.maxK = _matchPopsample(popsample, tpara.pop, mMaxK, pqc);

			if(outputproportion)
			{	sprintf(str, "%s.ppp", outname);	
				f = fopen(str, "w");
				int l = (int)pqc.size() / (mydata.indN * mMaxK + 1);
				for(j = 0; j < l; j++)
				{	for(k = 0; k < mydata.indN * mMaxK + 1; k++) fprintf(f, "%d ", pqc[j * (1 + mydata.indN * mMaxK) + k]);
                    fprintf(f, "\n");
				}
				fclose(f);
			}
		}
		else summarizePosterior(popsample, drecsample, tpara.pop, mydata.ploidity, recombs, mydata.L, mMaxK);////done
/*
vector<vector<int> > tsp0, tsp1;
int *tp[1];
tp[0] = tpara.pop;
_LongtoShort(tp, tpara.rec, tsp0, 1, mydata.L);
tp[0] = tpara.pop+mydata.L;
_LongtoShort(tp, tpara.rec+mydata.L, tsp1, 1, mydata.L);
FILE *f = fopen("tmp.txt", "w");
for(j = 0; j < (int)popsample.size(); j++)
{	fprintf(f, "id%d ", j*2);
	for(i=0;i<(int)popsample[j][0].size();i+=2) fprintf(f, "%d:%d ", popsample[j][0][i], popsample[j][0][i+1]);
	fprintf(f, "\n");
	fprintf(f, "id%d ", j*2+1);
	for(i=0;i<(int)popsample[j][1].size();i+=2) fprintf(f, "%d:%d ", popsample[j][1][i], popsample[j][1][i+1]);
	fprintf(f, "\n");
}
fprintf(f, "id%d ", j*2);
for(i=0;i<(int)tsp0[0].size();i+=2) fprintf(f, "%d:%d ", tsp0[0][i], tsp0[0][i+1]);
fprintf(f, "\n");
fprintf(f, "id%d ", j*2+1);
for(i=0;i<(int)tsp1[0].size();i+=2) fprintf(f, "%d:%d ", tsp1[0][i], tsp1[0][i+1]);
fprintf(f, "\n");
fclose(f);
*/	
/*		vector<int> dbreaks;

	for(i = 0; i < 5; i++)
	{	printf("%d ", i);fflush(stdout);
		for(j = 0; j < mydata.indN; j += mydata.ploidity)
		{	if((int)mydata.breaks.size() == 0) _getDoubleBreaks(vector<int>(), mydata.L, dbreaks);
			else _getDoubleBreaks(mydata.breaks[j/mydata.ploidity], mydata.L, dbreaks);
			_removeInd(mydata, tpara, j);
			if(tpara.imputedata) _imputeData(mydata, tpara, j, dbreaks);
			_addInd(mydata, tpara, j);
			_updatePrior();
		}
	}
*/
	}
	else
	{	vector<vector<char> > haps;
//printf("a");fflush(stdout);
		summarizeHaps(mydata.totalN * mydata.ploidity, mydata.L, mydata.ploidity, fullsample, cumsample, tpara.singlecall, haps);
//printf("d");fflush(stdout);

		if(tpara.imputedata)
		{	for(j = 0; j < mydata.indN; j+=mydata.ploidity) 
			{	_removeInd(mydata, tpara, j);
				int idst = mydata.indIndex[j / mydata.ploidity], ided = mydata.indIndex[j / mydata.ploidity + 1];
				for(l = idst * mydata.ploidity; l < ided * mydata.ploidity; l+=mydata.ploidity)	
				{	for(k = 0; k < mydata.L; k++) 
					{	mydata.data[l*mydata.L+k] = (float)haps[l][k];
						if(mydata.ploidity == 2) mydata.data[(l+1)*mydata.L+k] = (float)haps[l+1][k];
					}
				}
				_addInd(mydata, tpara, j);
			}
//printf("b");fflush(stdout);
		
			char omax = tpara.maximum;
			tpara.maximum = 2;
			vector<vector<bool> > ofixallele = mydata.fixAllele;
			mydata.fixAllele = vector<vector<bool> >(mydata.totalN * mydata.ploidity, vector<bool>(mydata.L, false));
			vector<int> dbreaks;
			_getDoubleBreaks(vector<int>(), mydata.L, dbreaks);
			for(k = 0; k < max(5, (int)(log((double)mydata.indN)/log(2.))); k++)
			{	for(j = 0; j < mydata.indN; j+=mydata.ploidity)
					if((int)mydata.fixState.size() == 0 || (int)mydata.fixState[j].size() == 0)
					{	_removeInd(mydata, tpara, j);
					        _updateStructure(mydata, tpara, j, mydata.ploidity, dbreaks);
						_addInd(mydata, tpara, j);
					}
			}
			mydata.fixAllele = ofixallele;
			tpara.maximum = omax;
		}
	}



	if(outqc) 
	{	sprintf(str, "%s.clusterprop", outname);
		_outputPQC(str, popsample, mydata.indN, mydata.L, mMaxK);
	}
	printf("Done\n");fflush(stdout);

	if(!tpara.singlecall)
	{	sprintf(str, "%s.cluster", outname);
		_outputPop(str, tpara.pop, mydata.indN, mydata.L, tpara.rec);
		//sprintf(str, "%s.pops", outname);
		//_outputPops(str, popsample);
	}
	popsample.clear();

	bool *tlapse = tpara.lapse;
	float *dp = mydata.data;
	for(i = 0; i < mydata.totalN; i++)
	{	int ln = 0;
		float prevd = *dp;
		for(j = 0; j < mydata.L; j++, tlapse++, dp++)
		{	if(!(*tlapse)) prevd = (int)(*dp);
			if(*dp != prevd) /*printf("%d:%d(%d,%d) ", i, j, (int)(*dp), (int)prevd),*/ln++;
			*dp = prevd;
		}
		//printf("%d: %d\n", i, ln);fflush(stdout);
	}
  
	outputResult(outname);

	//release memory
	delete [] tpara.recombN;
        delete [] recombs;
        delete []tpara.popr;
        delete tpara.nh;
	delete tpara.pop;
	delete tpara.rec;
	delete tpara.lapse;
	delete tpara.rPrior;
	delete []mydata.asz;
//	delete Q;
//	delete Vh;
}

double tensorHMMbase::_oneRound(MYDATA &mydata, TENSORPARA &tpara, vector<bool> &init, int iter)
{
	int j, k, l, st = 0, ed = mydata.L;
	bool fixst = false, fixed = false;
	if((int)mydata.fixState.size() > 0)
	{	int rn = 0;
		for(j = 0; j < mydata.indN; j += mydata.ploidity)
			if((int)mydata.fixState[j].size() > 0) rn++;
		int wb = tpara.refweight;
		for(j = 0; j < mydata.indN; j += mydata.ploidity)
			if(init[j / mydata.ploidity] && (int)mydata.fixState[j].size() > 0)
			{	int *pp = tpara.pop + j * mydata.L, x;
				bool *rr = tpara.rec + j * mydata.L;
				for(l = 2; l < (int)mydata.fixState[j].size(); l += 2)
				{	for(x = mydata.fixState[j][l - 2]; x < mydata.fixState[j][l]; x++) 
						if((x > st && x < ed - 1) || (x==st && !fixst) || (x==ed-1 && !fixed)) *(pp + x) = mydata.fixState[j][l - 1];	
					x = mydata.fixState[j][l];
					if((x > st && x < ed) || (x==st && !fixst)) *(rr + x) = true;
				}
				for(x = mydata.fixState[j][l - 2]; x < mydata.L; x++) 
					if((x > st && x < ed - 1) || (x==st && !fixst) || (x==ed-1 && !fixed)) *(pp + x) = mydata.fixState[j][l - 1];
				if(mydata.ploidity == 2)
				{	pp = tpara.pop + (j + 1) * mydata.L;
					rr = tpara.rec + (j + 1) * mydata.L;
					for(l = 2; l < (int)mydata.fixState[j + 1].size(); l += 2)
					{	for(x = mydata.fixState[j + 1][l - 2]; x < mydata.fixState[j + 1][l]; x++) 
							if((x > st && x < ed - 1) || (x==st && !fixst) || (x==ed-1 && !fixed)) *(pp + x) = mydata.fixState[j + 1][l - 1];	
						x = mydata.fixState[j + 1][l];
						if((x > st && x < ed) || (x==st && !fixst)) *(rr + x) = true;
					}
					for(x=mydata.fixState[j+1][l-2];x<mydata.L;x++) 
						if((x > st && x < ed - 1) || (x==st && !fixst) || (x==ed-1 && !fixed)) *(pp + x) = mydata.fixState[j + 1][l - 1];
				}
				_addInd(mydata, tpara, j, wb, st, ed);
				init[j / mydata.ploidity] = false;
			}
	}

	int tsz = mydata.indN / mydata.ploidity;
	vector<int> order, idlist(tsz);
	vector<double> pp(tsz);
	for(j = 0; j < tsz; j++) 
	{ 	pp[j] = tpara.idw[j]; 
		if(j > 0) pp[j] += pp[j - 1]; 
		idlist[j] = j;
	}
//printf("<");
	for(j = 0; j < tsz; j++)
	{	double un = gsl_ran_flat(gammar, 0, pp[(int)pp.size() - 1]);
		for(k = 0; k < (int)pp.size(); k++) if(un <= pp[k]) break;
		order.push_back(idlist[k]);
//printf("%d,",idlist[k]);
		idlist.erase(idlist.begin() + k);
		un = pp[k]; if(k > 0) un -= pp[k - 1];
		pp.erase(pp.begin() + k);
		for(k; k < (int)pp.size(); k++) pp[k] -= un;
	}
//printf(">\n");
	
	int cstep = 5, m, lastprop = -1;
	vector<int> dbreaks;
	if((int)mydata.breaks.size() == 0) _getDoubleBreaks(vector<int>(), mydata.L, dbreaks);

	double RT = 0;
        for(m = 0; m < tsz; m++)
	{	j = order[m] * mydata.ploidity;
		if((int)mydata.breaks.size() > 0) _getDoubleBreaks(mydata.breaks[j/mydata.ploidity], mydata.L, dbreaks);
		
clock_t cst, ced;
		if((int)mydata.fixState.size() == 0 || (int)mydata.fixState[j].size() == 0)
		{
			if(!init[j/mydata.ploidity])
			{	
				_removeInd(mydata, tpara, j, 1, st, ed);
				double opw = tpara.priorW;
				tpara.priorW *=max(1., 0.5*(double)tpara.trueindN*(double)(tpara.burnin/2 - iter)/tpara.burnin);
                		_updateRecomb(tpara, tpara.trueindN, mydata.L);
				tpara.priorW = opw;
			}
            
			_updateVh(tpara.Q, tpara.maxK, tpara.A, tpara.Vh, false);

			if((int)mydata.fixAllele.size() > 0 && (int)mydata.fixAllele[j].size() > 0)
			{	RT+=_updateStructure(mydata, tpara, j, 1, dbreaks); //to be changed, consider possibility of partial fix, originally, it was done by using dtype=0
				if(mydata.ploidity == 2) RT+=_updateStructure(mydata, tpara, j + 1, 1, dbreaks);
			}
			else 
			{	RT+=_updateStructure(mydata, tpara, j, mydata.ploidity, dbreaks);//-------
			}

cst=clock();
			if(tpara.imputedata) RT += _imputeData(mydata, tpara, j, dbreaks, init[j / mydata.ploidity]);
						//to be changed, partial imputation should be allowed, originally, it was done using dtype=0
ced = clock();
imputeT+=ced-cst;
			_addInd(mydata, tpara, j, 1, st,ed);
			init[j / mydata.ploidity] = false;

	//if(gID<tpara.burnin)		
			if(!fixst && !fixed) _collapseGroups(mydata.indN, mydata.L, mydata.ploidity, tpara, init);//don't use collapseGroups_new(), otherwise greedyswitch_interval won't work
	//if(gID<tpara.burnin)	
/*{int i, j;
for(i=0;i<tpara.param[1194].size();i++)
{	printf("%d: ", i);
	for(j=0;j<(int)tpara.param[1194][i].allele.size(); j++)
		printf("%f ", tpara.param[1194][i].allele[j]);
	printf("\n");
}
printf("<%d>", (int)mydata.data[1194]);
exit(0);
}*/
			if(m % cstep == (cstep - 1) && !fixst && !fixed) _greedySwitch_new(mydata.indN, mydata.L, mydata.ploidity, tpara, init);

		}

		_updatePrior(init);

		k = (m + 1) * 100 / tsz;
		/*printf("%d%%", k);
		if(m < tsz - 1)
		{	printf("%c%c", 8, 8);
			if(k >= 10) printf("%c", 8);
		}
		else printf(" ");
		*/
		if((int)(k / 5) > lastprop)
		{	for(int x = lastprop; x < (int)(k/5); x++) printf(".");
			lastprop = (int)(k/5);;
		}
//printf("%d,",tpara.maxK);
		fflush(stdout);
        }
	return RT;
}

double tensorHMMbase::_updateStructure(MYDATA &mydata, TENSORPARA &tpara, int id, int ploidity, vector<int> const &dbreaks)
{	int i, j, SSZ;
	vector<vector<int> > ss, mlist;

	int MK = min(tpara.maxK + 1, tpara.maxHapK);
	if(ploidity == 2) 
	{	SSZ = MK*(MK+1)/2;
		mlist.resize(MK);
		ss.resize(SSZ, vector<int>(2));
		for(i = 0; i < SSZ; i++)
		{	ss[i][0] = (int)((sqrt(1. + 8. * (double)i) - 1.) / 2.);
			ss[i][1] = i - ss[i][0] * (ss[i][0] + 1) / 2;
			mlist[ss[i][0]].push_back(i);
			if(ss[i][1] != ss[i][0]) mlist[ss[i][1]].push_back(i);
		}
	}
	else
	{	SSZ=MK;
		ss.resize(SSZ, vector<int>(1));
		for(i = 0; i < SSZ; i++) ss[i][0] = i;
	}

	int dl = (int)dbreaks.size();

	j = dl * (SSZ + MK + 1);
	if(tmpM != NULL) delete tmpM;
	tmpM = new double[j];
	double *P = tmpM, *p = P;
	for(i = 0; i < j; i++, p++) (*p)=0;
	double Vh2[SSZ], Vh[tpara.maxK + 1];
	for(i = 0; i <= tpara.maxK; i++) Vh[i] = tpara.Vh[i];// * tpara.idw[id];//(double)(tpara.idw[id]==1 || i < tpara.maxK);
	for(i = 0; i < SSZ; i++) { Vh2[i] = tpara.Vh[ss[i][0]]; if(ploidity > 1) Vh2[i] *= tpara.Vh[ss[i][1]]; }

    _forwardSum(mydata, tpara, id, ploidity, ss, mlist, dbreaks, Vh, Vh2, P);
//double sump = 0;
//for(i = 0; i < SSZ; i++) sump += exp(P[(mydata.L-1)*SSZ+i]-P[dl*(SSZ+MK)+mydata.L-1]); 
//printf("%d:%f\n", id, log(sump)+P[dl*(SSZ+MK)+mydata.L-1]); 
    
	double RT = _backwardSample(mydata, tpara, id, ploidity, ss, mlist, dbreaks, Vh, Vh2, P);
/*int rn = 0, rn0=0, rn1=0, flip = tpara.flips[id/2][1], flipi=0;
for(i=2;i<(int)mydata.breaks[id/2].size()-2;i+=2)
{	if(mydata.breaks[id/2][i+1]==2) { flipi+=2; flip = tpara.flips[id/2][flipi+1]; }
	rn += tpara.rec[id*mydata.L+mydata.breaks[id/2][i]];
	if(mydata.breaks[id/2][i+1]==flip) rn0 += tpara.rec[id*mydata.L+mydata.breaks[id/2][i]];
	if(mydata.breaks[id/2][i+1]==1-flip) 
	{	rn1 += tpara.rec[id*mydata.L+mydata.breaks[id/2][i]];
		if(tpara.rec[id*mydata.L+mydata.breaks[id/2][i]]) printf("id=%d, pos=%d, flip=%d, %d->%d\n",id, mydata.breaks[id/2][i],flip, tpara.pop[id*mydata.L+mydata.breaks[id/2][i]-1], tpara.pop[id*mydata.L+mydata.breaks[id/2][i]]),exit(0);
	}
}
printf("R=%d,%d,%d\n", rn, rn0, rn1);
*/
//	delete [] P;
	return RT;
}

void tensorHMMbase::_forwardSum(MYDATA &mydata, TENSORPARA &tpara, int id, int ploidity, vector<vector<int> > const &ss, vector<vector<int> > const &mlist, vector<int> const &dbreaks, double *Vh, double *Vh2, double *P)
{
	int SSZ = (int)ss.size(), MK = (int)mlist.size();
	if(ploidity == 1) MK = SSZ;
	int i, dl = (int)dbreaks.size();
	double *lpp = P, *M = P + dl * SSZ, *lmm = M, *S = M + dl * MK, *lss = S;
	float *dp1 = mydata.data + (mydata.indIndex[id/mydata.ploidity] * mydata.ploidity + id % mydata.ploidity) * mydata.L, *dp2;
//this part is right only if ploidity is within replicates, i.e., each replicate takes ploidity, rather than each ploidity takes replicates.

	double const *rp = tpara.popr;
	if(ploidity > 1) dp2 = dp1 + mydata.L;
	float const *nhp = tpara.nh;
	int const *aszp = mydata.asz;
	vector<vector<PARAUNIT> >::const_iterator ipm = tpara.param.begin();
	int lasti = dbreaks[0];
	vector<int>::const_iterator idb = dbreaks.begin();

	vector<int> tbreaks;
	if((int)mydata.breaks.size() > 0) tbreaks = mydata.breaks[id / mydata.ploidity];
	idb++;
clock_t cst, ced;
cst=clock();
	for(i = 1; i < dl; i++, idb++)
	{
		_sumforwardOne(id, tbreaks, lasti, *idb, tpara, ploidity, dp1, dp2, rp, nhp, ipm, aszp, lpp, lmm, lss, SSZ, MK, ss, mlist, Vh, Vh2);
		lasti = (*idb);
	}
ced=clock();
forwardT += ced - cst;
//printf("forward of %d = %fsec\n", dl, (double)(ced-cst)/1000000.);
}

double tensorHMMbase::_backwardSample(MYDATA &mydata, TENSORPARA &tpara, int id, int ploidity, vector<vector<int> > const &ss, vector<vector<int> > const &mlist, vector<int> const &dbreaks, double *Vh, double *Vh2, double *P)
{	int i, j;
	float const *nhp;
	int *pp, *pp2 = NULL;
        double const *rp;
	bool *rr, *rr2 = NULL;

	pp = tpara.pop + (id * mydata.L +mydata.L - 1);
	rr = tpara.rec + id * mydata.L + mydata.L;
	rp = tpara.popr + mydata.L;
	if(ploidity > 1) 
	{	pp2 = pp + mydata.L;
		rr2 = rr + mydata.L;
	}
        
	int dl = (int)dbreaks.size();
	int SSZ = (int)ss.size(), MK = (int)mlist.size();
	if(ploidity == 1) MK = SSZ;
	double *lpp, *lmm, *sump,  *mp;
	lpp = P + (dl - 2) * SSZ;
	lmm = P + dl * SSZ + (dl - 2) * MK;
	sump = P + dl * (SSZ + MK) + (dl - 2);
	double RT = *sump;
	mp = P + dl * SSZ + (dl - 2) * MK;
	nhp = tpara.nh + tpara.maxK * (mydata.L - 1);

	double cumsum[SSZ];
	vector<vector<PARAUNIT> >::const_iterator pit = tpara.param.begin() + (mydata.L - 1);
	double recp = 0, recp2 = 0, nrecp2 = 0;
	clock_t St = clock(), Et = 0;

	vector<int> flips, tbreaks;
	if((int)mydata.breaks.size() > 0) 
	{	flips.resize((dl - 1) * 2, 0);
		tbreaks = mydata.breaks[id / mydata.ploidity];
	}
	vector<int>::const_iterator idb = dbreaks.end() - 1;	

	double aaa = tpara.heteroVh;
	double NNN = aaa * (double)(tpara.trueindN * mydata.ploidity + tpara0.totalN0);

	int lscode = 0;
	for(i = dl - 1; i > 0; i--, idb--)
	{	int rn = 0;
                if(i < dl - 1)
		{	recp = *rp;
			recp2 = recp * recp;
			nrecp2 = (1. - recp) * (1.-recp);
			if(tpara.singlecall) { recp = recp2 = 1.-1e-10; nrecp2 = 1e-10; }

			double ttt[4];
			if(ploidity == 1)
			{   	if(lscode < tpara.maxK) ttt[1] = (double)(*(nhp + tpara.maxK + lscode));
                		else ttt[1] = 0;
                		if(lscode < tpara0.popn0)
					ttt[1] += (double)(*(tpara0.nh0 + tpara0.popn0 * i + tpara0.popn0 + lscode));
double ottt1=ttt[1];
				ttt[1] = (Vh[lscode] + aaa * (double)(lscode < tpara.maxK)*ttt[1]) / (1. + NNN) * recp;
//if(lscode<tpara.maxK)ttt[1]*=min(1.,ottt1+NUMPRECISION);

                        	if(tpara.logscale) ttt[0] = exp(*(lpp + lscode) - (*sump)) * (1. - recp);
				else ttt[0] = (*(lpp + lscode)) * (1. - recp) / (*sump);
			}
			else
			{	int ls1 = min(tpara.maxK, *(pp+1)), ls2 = min(tpara.maxK, *(pp2 + 1));
				if(tpara.logscale)
				{	ttt[0] = exp(*(lpp + lscode) - (*sump)) * nrecp2 / (1. + (double)(ls1!=ls2));
					ttt[1] = exp((*(mp + ls1) - (*sump)) * Vh[ls2]) / 2. * recp * (1. - recp);
					ttt[2] = exp((*(mp + ls2) - (*sump)) * Vh[ls1]) / 2. * recp * (1. - recp);	
				}
				else
				{	ttt[0] = (*(lpp + lscode)) * nrecp2 / (*sump) / (1. + (double)(ls1!=ls2));
					ttt[1] = ((*(mp + ls1)) * Vh[ls2]) / 2. * recp * (1. - recp) / (*sump);
					ttt[2] = ((*(mp + ls2)) * Vh[ls1]) / 2. * recp * (1. - recp) / (*sump);	
				}
				ttt[3] = recp2 * tpara.Vh[ls1] * tpara.Vh[ls2];
			}
			rn = _sample(ttt, ploidity * 2, tpara.maximum & (1 << 2));
//if(id==15 && i==10058) printf("rn=%d [%f,%f] ", rn, ttt[0], ttt[1]);fflush(stdout);
if(false)//i==1)
{
	int kk;
	printf("id=%d,pos=%d, ttt[0]=%f, ttt[1]=%f, r=%f, lscode=%d ", id, i, log(ttt[0]), log(ttt[1]), recp, lscode);
	for(kk = 0; kk < SSZ; kk++) printf("%f,", lpp[kk]);
    printf(" sump=%f\n", *sump);
	for(kk = 0; kk < SSZ; kk++) printf("Vh[%d]=%f ", kk, Vh[kk]);
    //printf("\n"); // Sachin
//	if(i == 3281) exit(0);
}
                }

		if(i < dl - 1)
		{	(*rr) = false;
			if(ploidity > 1) (*rr2) = false;
		}
                if(i == dl - 1 || rn > 0)
                {       int ls1 = 0, ls2 = 0;
			if(i < dl - 1)
			{	ls1 = min(tpara.maxK, *(pp+1)); 
				if(ploidity > 1) ls2 = min(tpara.maxK, *(pp2+1)); 
				if(rn == 1) 
				{	if(ploidity == 1) (*rr) = true;
					else (*rr2) = true;
				}
				else if(rn == 2) (*rr) = true;
				else if(rn == 3) (*rr) = (*rr2) = true;
			}

			if(ploidity == 1)
			{	
		//	printf("%d: ", id);	
				for(j = 0; j < SSZ; j++)
                        	{       if(tpara.logscale) cumsum[j] = exp(*(lpp + j) - (*sump));
					else cumsum[j] = (*(lpp + j) / (*sump));
		//	printf("%f ", *(lpp+j)-(*sump));
				}
		//	printf("\n");
				j = _sample(cumsum, SSZ, tpara.maximum & (1<<1));
/*if(id ==15 && i == 10058)
{	printf("id=%d:%d j=%d %d | ", id, i, j, (int)(*(idb-1)));
	int kk;
	for(kk = 0; kk < SSZ; kk++) printf("%f(%f) ", log(cumsum[kk]), tpara.nh[(*(idb-1))*tpara.maxK+kk]);
	printf("\n");
}*/
	                        *pp = j;
	                        if(j == tpara.maxK && (tpara.maximum & (1<<1)) == 0)
                        	{	do { (*pp) = tpara.maxK + (int)(- log(gsl_ran_flat(gammar, 0, 1.)/*tpara.Vh[tpara.maxK]) / tpara.Vh[tpara.maxK]*/) / log((1. + tpara.A) / tpara.A));
					} while(tpara.maxHapK > 0 && (*pp) >= tpara.maxHapK && (*pp) > tpara.maxK);
	                        }
			}
			else
			{	int nn = SSZ, ii = 0, ti;
				if(rn > 2 || i == dl - 1)
				{	for(j = 0; j < SSZ; j++)
                        		{       if(tpara.logscale) cumsum[j] = exp(*(lpp + j) - (*sump));
						else cumsum[j] = *(lpp + j) / (*sump);
					}
				}
                        	else 
				{	ii = ls1;
					if(rn == 2) ii = ls2;
					nn = (int)mlist[ii].size();
					for(j = 0; j < nn; j++)
					{	ti = mlist[ii][j];
						if(tpara.logscale) cumsum[j] = exp(*(lpp + ti) - (*sump)) / (1. + (double)(ss[ti][0] != ss[ti][1]));
						else cumsum[j] = *(lpp + ti) / (*sump) / (1. + (double)(ss[ti][0] != ss[ti][1]));
					}
				}
				j = _sample(cumsum, nn, tpara.maximum & (1<<1));
/*if(i==dl-1)
{	printf("\n[i=%d ", i);
	for(k=0;k<nn;k++) printf("%d:%f ", k,cumsum[k]);
	printf("] j=%d, maximum=%d\n", j, (int)tpara.maximum);
}*/
				if(i != dl - 1 && rn <= 2) j = mlist[ii][j];
				
				if(rn > 2 || i == dl-1)
				{	double sp = 0.5;
					if(i!=dl-1 && ls1!=ls2 && (ls1 == ss[j][0] || ls2 == ss[j][1])) sp = 1.;
					else if(i!=dl-1 && ls1!=ls2 && (ls1 == ss[j][1] || ls2== ss[j][0])) sp = 0;
					if(gsl_ran_flat(gammar, 0, 1.) < sp)
					{	*pp = ss[j][0]; *pp2 = ss[j][1];	}
					else {	*pp = ss[j][1]; *pp2 = ss[j][0];	}
				}
				else if(rn == 1)
				{	*pp = *(pp + 1);
					if(ss[j][0] != ls1) *pp2 = ss[j][0];
					else *pp2 = ss[j][1];
				}
				else 
				{	*pp2 = *(pp2 + 1);
					if(ss[j][0] != ls2) *pp = ss[j][0];
					else *pp = ss[j][1];
				}

                                if(*pp == tpara.maxK && (rn > 1 || i == dl - 1) && (tpara.maximum & (1<<1)) == 0)
				{	do { (*pp) = tpara.maxK + (int)(- log(gsl_ran_flat(gammar, 0, 1.)/*tpara.Vh[tpara.maxK]) / tpara.Vh[tpara.maxK]*/) / log((1. + tpara.A) / tpara.A));
                                        } while(tpara.maxHapK > 0 && (*pp) >= tpara.maxHapK && (*pp) > tpara.maxK);
                                }
                                if(*pp2 == tpara.maxK && (rn == 1 || rn > 2 || i == dl - 1) && (tpara.maximum & (1<<1)) == 0)
				{	do { (*pp2) = tpara.maxK + (int)(- log(gsl_ran_flat(gammar, 0, 1.)/*tpara.Vh[tpara.maxK]) / tpara.Vh[tpara.maxK]*/) / log((1. + tpara.A) / tpara.A));
                                        } while(tpara.maxHapK > 0 && (*pp2) >= tpara.maxHapK && (*pp2) > tpara.maxK);
                                }
			}
              	}
                else 
		{	(*pp) = *(pp + 1);
			if(ploidity == 2) *pp2 = *(pp2 + 1);
		}


/*		if(i < dl - 1)
		{	int ii;
			double tmpp, tsum = 0, tttt[MK];
			for(ii = 0; ii < MK; ii++)
			{	if(ii < tpara.maxK) tmpp = (double)(*(nhp + tpara.maxK + ii));
				else tmpp = 0;
                		if(ii < tpara0.popn0)
					tmpp += (double)(*(tpara0.nh0 + tpara0.popn0 * i + tpara0.popn0 + ii));
				tmpp = (Vh[lscode] + aaa * (double)(ii < tpara.maxK)*tmpp) / (1. + NNN) * recp;
				tttt[ii] = exp(lpp[ii + SSZ] - sump[1]);

				tmpp *= tttt[ii];
				tsum += tmpp;
			}
			double t0, tmax, tmaxm; 
			for(ii = 0; ii < MK; ii++)
			{	tmpp = lmm[ii];
				t0 = fast_log((1. - recp) * tttt[ii] + recp * tsum);
				lmm[ii] = lpp[ii] - sump[0] + t0;
				if(ii == 0 || tmaxm < lmm[ii]) tmaxm = lmm[ii];
				lpp[ii] = tmpp + t0;
				if(ii == 0 || tmax < lpp[ii]) tmax = lpp[ii];
			}
			sump[0] = 0;
			double tss = 0;
			for(ii = 0; ii < MK; ii++)
			{	sump[0] += exp(lpp[ii] - tmax);
				lmm[ii] = exp(lmm[ii] - tmaxm);
				tss += lmm[ii];
			}
			for(ii = 0; ii < MK; ii++) lmm[ii] /= tss;
			sump[0] = fast_log(sump[0]) + tmax;
		}
		else
		{	int ii;
			double tmpp, tmax, tmaxm;
			for(ii = 0; ii < MK; ii++)
			{	tmpp = lmm[ii];
				lmm[ii] = lpp[ii];
				if(ii == 0 || tmaxm < lmm[ii]) tmaxm = lmm[ii];
				lpp[ii] = tmpp;
				if(ii == 0 || tmax < lpp[ii]) tmax = lpp[ii];
			}
			sump[0] = 0;
			double tss = 0;
			for(ii = 0; ii < MK; ii++)
			{	sump[0] += exp(lpp[ii] - tmax);
				lmm[ii] = exp(lmm[ii] - tmaxm);
				tss += lmm[ii];
			}
			for(ii = 0; ii < MK; ii++) lmm[ii] /= tss;
			sump[0] = fast_log(sump[0]) + tmax;
		}
*/
/*		double tmax; 
		for(int ii = 0; ii < MK; ii++) 
		{	
			if(ii < tpara.maxK) lmm[ii] = (double)(*(nhp + ii));
			else lmm[ii] = 0.1;
			if(ii==0 || tmax < lmm[ii]) tmax=lmm[ii];
		}
		double tss = 0;
		for(int ii = 0; ii < MK; ii++) { lmm[ii] = exp(lmm[ii] - tmax); tss += lmm[ii]; }
		for(int ii = 0; ii < MK; ii++) { lmm[ii] /= tss; }
		lmm -= MK;
*/
		
		int st = *(idb-1);
		int fff = _tracebackOne(id, tbreaks, st, *idb, tpara, ploidity, pp, pp2, rr, rr2, rp, nhp, pit, lpp, mp, sump, SSZ, MK, ss, mlist, Vh, Vh2);
		if(ploidity == 2 && tbreaks.size() > 0) 
		{	flips[i * 2 - 2] = tbreaks[st]; 
			flips[i * 2 - 1] = fff;
		}

		if(ploidity == 1) lscode = min(MK - 1, *(pp+1));
		else 
		{	int i1 = *(pp+1), i2= *(pp2+1);
			if(i1 < i2) { j = i1; i1 = i2; i2 = j; }
			lscode = i1 * (i1 + 1) / 2 + i2;
		}
	}

	if(ploidity == 2 && tbreaks.size() > 0) tpara.flips[id / mydata.ploidity] = flips;
	Et = clock();
	backwardT += Et - St;

	return RT;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int tensorHMMbase::_shrinkClasses(int indN, int L, float *z)
{	int i, j, k;
	int maxK = 0;
	float *zp = z;
	for(i = 0; i < indN; i++)
		for(j = 0; j < L; j++, zp++)
			if(maxK < (int)(*zp)) maxK = (int)(*zp);
	maxK++;

	int newmaxK = 0;
	vector<int> map(maxK, -1), remap(maxK, -1), cn(maxK, 0);
	for(i = 0; i < L; i++)
	{	k = 0;
		zp = z + i;
		for(j = 0; j < indN; j++, zp += L)
		{	if(map[(int)(*zp)] < 0)
			{	while(remap[k] >= 0) k++;
				map[(int)(*zp)] = k;
				remap[k] = (int)(*zp);
				newmaxK = max(newmaxK, k);
			}
			(*zp) = map[(int)(*zp)];
			cn[(int)(*zp)]++;
		}
		for(j = 0; j < maxK; j++)
		{	if(cn[j] == 0 && remap[j] >= 0)
			{	remap[j] = map[remap[j]] = -1;	}
			cn[j] = 0;
		}
	}

	return newmaxK;
}

//This function is specific to categorical data
//return maxK
//nh: group size per position
//Q: group size overall
//r: recombination probability 
void tensorHMMbase::_removeInd(MYDATA const &mydata, TENSORPARA &tpara, int id, float wb, int st, int ed, bool clean)
{       int i, j, k, idd = id / mydata.ploidity, idst = mydata.indIndex[idd], ided = mydata.indIndex[idd + 1];
	if(st < 0) st = 0;
	if(ed < 0) ed = mydata.L;
        int const *pp = tpara.pop + id * mydata.L + st, *pp2 = pp + (mydata.ploidity - 1) * mydata.L;
	float const *dp = mydata.data + idst*mydata.ploidity * mydata.L + st, *dp2 = dp + (mydata.ploidity - 1) * mydata.L;
	bool const *rr = tpara.rec + id * mydata.L + st, *rr2 = rr + (mydata.ploidity - 1) * mydata.L;
	float *rn = tpara.recombN + st;
        float *nhp = tpara.nh + st * tpara.maxK;
        int K = 0;

	bool *tlapse = tpara.lapse + idst * mydata.L + st;

//	int flipi = 0, curflip = 0;
//	int fsz = 0;
//	if((int)tpara.flips.size() > 0) fsz = (int)tpara.flips[idd].size();

	vector<vector<PARAUNIT> >::iterator pit = tpara.param.begin() + st;
	int ksz = mydata.L * mydata.ploidity, kkk = idst * mydata.ploidity, kkksz = mydata.ploidity;
        for(i = st; i < ed; i++, pp++, pp2++, dp++, dp2++, nhp += tpara.maxK, pit++, rr++, rr2++, rn++, tlapse++)
        { //      while(flipi < fsz && tpara.flips[idd][flipi] < i) flipi += 2;
	//	if(flipi < fsz && i == tpara.flips[idd][flipi])
	//	{	if(curflip != tpara.flips[idd][flipi + 1])
	//		{	//float const *tdp = dp; dp = dp2; dp2 = tdp;
	//			curflip = tpara.flips[idd][flipi + 1];
	//		}
	//	}		

		if(*rr) *rn -= tpara.idw[id];
		int kk = 0, kks = kkk;
		for(k = 0; k < ided - idst; k++, kk += ksz, kks += kkksz) 
		{	float b = tpara.indWeight[idst + k];//wb;//ploidity;
			*(nhp + (*pp)) -= b;
//printf("%d:%d ", id, i);fflush(stdout);
			_updatePara_remove((*pit)[*pp], (int)(*(dp + kk)), (float)b, kks, (int)(tpara.lapse != NULL && *(tlapse+kk)));//i);
		}
		if(i==0 || *rr) tpara.Q[*pp]-=tpara.idw[id];

		if(mydata.ploidity == 2)
		{	 
			if(*rr2) *rn -= tpara.idw[id];
			kk = 0; kks = kkk + 1;
			for(k = 0; k < ided - idst; k++, kk+=ksz, kks+=kkksz) 
			{	float b = tpara.indWeight[idst + k];
				*(nhp + (*pp2)) -= b;
				_updatePara_remove((*pit)[*pp2], (int)(*(dp2 + kk)), (float)b, kks, (int)(tpara.lapse != NULL && *(tlapse+kk)));//i);
			}
			if(i == 0 || (*rr2)) tpara.Q[*pp2]-=tpara.idw[id];
		}
        }

	if(clean)
	{	K = tpara0.popn0;
		nhp = tpara.nh;
		for(i = 0; i < mydata.L; i++, nhp += tpara.maxK)
		{	for(j = tpara.maxK - 1; j >= K; j--) if(*(nhp + j) > NUMPRECISION) break;
			K = j + 1;
		}
		if(K < tpara.maxK)
		{	PARAUNIT pu;
			pu.lapseN = pu.popN = 0;
			float *newnh = new float[K * mydata.L], *newnhp = newnh;
	                nhp = tpara.nh;
        	        for(i = 0; i < mydata.L; i++)
                	{       for(j = 0; j < K; j++, newnhp++, nhp++)
                        	        (*newnhp) = (*nhp);
	                        nhp += tpara.maxK - K;
        	                tpara.param[i].resize(K);//, pu);
                	}
	                delete tpara.nh;
        	        tpara.nh = newnh;

        	        tpara.maxK = K;
		}
        }
}

//This function is specific to categorical data
//return maxK
//nh: group size per position
//Q: group size overall
//r: recombination probability
//asz: number of distinct values per position
void tensorHMMbase::_addInd(MYDATA const &mydata, TENSORPARA &tpara, int id, float wb, int st, int ed)
{       int i, j, k, idd = id / mydata.ploidity, idst = mydata.indIndex[idd], ided = mydata.indIndex[idd + 1];
	if(st < 0) st = 0;
	if(ed < 0) ed = mydata.L;
	int const *pp = tpara.pop + id * mydata.L + st, *pp2 = pp + (mydata.ploidity - 1) * mydata.L;
	float const *dp = mydata.data + idst * mydata.ploidity * mydata.L + st, *dp2 = dp + (mydata.ploidity - 1) * mydata.L;
	bool const *rr = tpara.rec + id * mydata.L + st, *rr2 = rr + (mydata.ploidity - 1) * mydata.L + st;
	float *rn = tpara.recombN + st;
	float *nhp = tpara.nh + st * tpara.maxK;
	int K = tpara.maxK;
	
	int const *tpp = pp, *tpp2 = pp2;
	for(i = st; i < ed; i++, tpp++, tpp2++)
	{       if(*tpp >= K) K = (*tpp) + 1;
		if(*tpp2 >= K) K = (*tpp2) + 1;
	}

//printf("K=%d maxK=%d\n", K, tpara.maxK);fflush(stdout);
	if(K > tpara.maxK)
	{	PARAUNIT pu;
		pu.lapseN = pu.popN = 0;	
		float *newnh = new float[K * mydata.L], *newnhp = newnh; 
		int *aszp = mydata.asz;
		for(i = 0; i < K * mydata.L; i++, newnhp++) *newnhp = 0;
		newnhp = newnh;
		nhp = tpara.nh;
		for(i = 0; i < mydata.L; i++, aszp++)
		{       for(j = 0; j < tpara.maxK; j++, newnhp++, nhp++)
				(*newnhp) = (*nhp);
			newnhp += K - tpara.maxK;
			tpara.param[i].resize(K, pu);//, vector<float>(*aszp, 0));
		}
		delete tpara.nh;
		tpara.nh = newnh;

		for(i = tpara.maxK; i < K + 1; i++) tpara.Q[i] = 0;
		tpara.maxK = K;
		nhp = tpara.nh + st * tpara.maxK;
	}

	bool *tlapse = tpara.lapse + idst * mydata.L + st;

//	int flipi = 0, curflip = 0;
//	int fsz = 0;
//	if((int)tpara.flips.size() > 0) fsz = (int)tpara.flips[idd].size();

	vector<vector<PARAUNIT> >::iterator pit = tpara.param.begin() + st;
	int ksz = mydata.L * mydata.ploidity, kkk = idst * mydata.ploidity, kkksz = mydata.ploidity;
	for(i = st; i < ed; i++, pp++, pp2++, dp++, dp2++, nhp += tpara.maxK, pit++, rr++, rr2++, rn++, tlapse++)
	{ //      while(flipi < fsz && tpara.flips[idd][flipi] < i) flipi += 2;
	//	if(flipi < fsz && i == tpara.flips[idd][flipi])
	//	{	if(curflip != tpara.flips[idd][flipi + 1])
	//		{//	float const *tdp = dp; dp = dp2; dp2 = tdp;
	//			curflip = tpara.flips[idd][flipi + 1];
	//		}
	//	}

		if(*rr) *rn += tpara.idw[id];
		int kk = 0, kks = kkk;
		for(k = 0; k < ided - idst; k++, kk += ksz, kks += kkksz) 
		{	float b = tpara.indWeight[idst + k];
			*(nhp + (*pp)) += b;
			_updatePara_add((*pit)[*pp], (int)(*(dp + kk)), (float)b, kks, (int)(tpara.lapse!=NULL && *(tlapse+kk)));//i);
		}

		if(i == 0 || *rr) tpara.Q[*pp]+=tpara.idw[id];

		if(mydata.ploidity == 2)
		{        
			if(*rr2) *rn += tpara.idw[id];
			kk=0; kks=kkk+1;
			for(k = 0; k < ided - idst; k++, kk+=ksz, kks+=kkksz) 
			{	float b = tpara.indWeight[idst + k];
				*(nhp + (*pp2)) += b;
				_updatePara_add((*pit)[*pp2], (int)(*(dp2 + kk)), (float)b, kks,(int)(tpara.lapse!=NULL && *(tlapse+kk)));//i);
			}
			if(i == 0 || *rr2) tpara.Q[*pp2]+=tpara.idw[id];
		}
	}
}

int tensorHMMbase::_updateRecomb(TENSORPARA &tpara, float trueindN, int L)
{	int i;
	int *list = new int[L], lsz = 0; 
	float *rn = tpara.recombN, *rn0 = tpara0.recombN0;
	double *rp = tpara.rPrior, *rr = tpara.popr, k, t, n;
	t = n = 0;
	double tt = 0;
	for(i = 1; i < L; i++, rn++, rr++, rp++)
	{	k = *rn;//tpara.recombN[i];
		if(rn0 != NULL) { k += *rn0; rn0++; } //if(tpara0.recombN0 != NULL) k += tpara0.recombN0[i];
		tt += k;
		if(k > tpara.recombcut) *rr = (double)(*rp > 0) * ((double)k + (*rp) * tpara.priorW) / ((double)(trueindN + tpara0.totalN0) + tpara.priorW);
		else { list[lsz++] = i; t += k; n += trueindN; }
	}
	double baserate = //gsl_ran_beta(gammar, (double)t + 1. / (double)trueindN, (double)(n - t) + 1.);
			((double)t + 0.0) / ((double)n + 1.);
	for(i = 0; i < lsz; i++) 
		tpara.popr[list[i]] = (double)(tpara.rPrior[list[i]] > 0) * (baserate * (double)trueindN + tpara.rPrior[list[i]] * tpara.priorW) / (double)(trueindN + tpara.priorW);
	delete []list;

	return((int)tt);
}

void tensorHMMbase::_updateVh(double *Q, int K, double A, double *V, bool verbose)
{       int i;
	double totalN = 0, *qq = Q, *vv;

	for(i = 0; i < K; i++, qq++) totalN += *qq;
	for(i = 0; i < tpara0.popn0; i++) totalN += tpara0.Q0[i];

	double vhprod = 1., v;
	qq = Q;
	v = 0;
	if(tpara0.popn0 > 0) v = tpara0.Q0[0];
	v = V[0] = (v + *qq + 1.) / (totalN + 1. + A);
	vv = V + 1;
	for(i = 1; i < K; i++, vv++)
	{       totalN -= *qq;
		if(tpara0.popn0 >= i) totalN -= tpara0.Q0[i - 1];
		vhprod *= 1. - v;
		qq++;
		v = 0;
		if(tpara0.popn0 > i) v = tpara0.Q0[i];
		v = (v + *qq + 1.) / (totalN + 1. + A);
		*vv = v * vhprod;
	}
	vhprod *= 1. - v;
	V[K] = vhprod;

	if(verbose)
	{       double sum = 0, vhprod = 1.;
        printf("\nK=%d,A=%f\n", K, A);
		for(i = 0; i < K + 1; i++)
		{       double v = V[i] / vhprod;
            printf("%d:\t%f (%f, %f)\n", i, V[i], v, Q[i]);
			sum += V[i];
			vhprod *= 1. - v;
		}
        printf("sum=%f\n", sum);
	}
}

double tensorHMMbase::_updateA(int maxK, double const *Q, double A)
{       int i;
	double totalN = 0;
	for(i = 0; i < maxK; i++) totalN += Q[i];
	for(i = 0; i < tpara0.popn0; i++) totalN += tpara0.Q0[i];

	double vhprod = 1., v;
	v = 0;
	if(tpara0.popn0 > 0) v = tpara0.Q0[0];
	v = (v + Q[0] + 1.) / (totalN + 1. + A);
	for(i = 1; i < maxK; i++)
	{       totalN -= Q[i - 1];
		if(tpara0.popn0 >= i) totalN -= tpara0.Q0[i - 1];
		vhprod *= 1. - v;
		v = 0;
		if(tpara0.popn0 > i) v += tpara0.Q0[i];
		v = (v + Q[i] + 1.) / (totalN + 1. + A);
	}
	vhprod *= 1. - v;

	A = gsl_ran_gamma(gammar, 5. + (double)maxK, 1./(1. - log(vhprod)));
	return A;
}

void tensorHMMbase::_collapseGroups(int indN, int L, int ploidity, TENSORPARA &tpara, vector<bool> const &init)
{       int i, j;

	vector<int> sel, rem, map(tpara.maxK, -1);
	j = 0;
	for(i = 0; i < tpara.maxK; i++)
	{	if(tpara.Q[i] > NUMPRECISION || tpara0.popn0 > i)
		{       sel.push_back(i);
			map[i] = j++;
		}
		else 
		{	rem.push_back(i);
		}
	}

	int K = (int)sel.size();
	if(K == tpara.maxK) return;

	for(i = 0; i < indN; i++)
		if(!init[i / ploidity])
		{       int *zp = tpara.pop + i * L;
			for(j = 0; j < L; j++, zp++)
				(*zp) = map[*zp];
		}

	float *newnh = new float[K * L], *newnhp = newnh;
	float *nhp = tpara.nh;
	for(i = 0; i < L; i++)
	{       for(j = 0; j < K; j++, newnhp++)
		*newnhp = *(nhp + sel[j]);
		nhp += tpara.maxK;
	}
	delete tpara.nh;
	tpara.nh = newnh;

	PARAUNIT pu;
	pu.lapseN = pu.popN = 0;
	for(i = 0; i < L; i++)
	{       
		for(j = 0; j < K; j++)
			tpara.param[i][j] = tpara.param[i][sel[j]];
		tpara.param[i].resize(K);//, pu);
	}

	double Qbak[tpara.maxK];
	for(i = 0; i < tpara.maxK; i++) Qbak[i] = tpara.Q[i];
	for(i = 0; i < K; i++)
		tpara.Q[i] = Qbak[sel[i]];
	for(i=i; i < tpara.maxK; i++) tpara.Q[i] = 0;
	tpara.maxK = K;
}

void tensorHMMbase::_collapseGroups_new(MYDATA const &mydata, TENSORPARA &tpara)
{	int i, j, k;

	int newmaxK = 0;
	vector<vector<int> > maps(mydata.L);
	vector<int> map(tpara.maxK, -1), remap = map, cn(tpara.maxK, 0);
	int *zp[mydata.indN];
	for(i = 0; i < mydata.indN; i++) zp[i] = tpara.pop + i * mydata.L;

	for(i = 0; i < mydata.L; i++)
	{	k = 0;
		bool sel[tpara.maxK];
		for(j = 0; j < tpara.maxK; j++) sel[j] = (j < tpara0.popn0);
		for(j = 0; j < mydata.indN; j++) sel[*zp[j]] = true;
		for(j = 0; j < tpara.maxK; j++)
			if(sel[j])
			{	if(map[j] < 0)
				{	while(remap[k] >= 0) k++;
					map[j] = k;
					remap[k] = j;
					newmaxK = max(newmaxK, k);
				}
			}
		for(j = 0; j < mydata.indN; j++)
		{	(*zp[j]) = map[*zp[j]];
			cn[*zp[j]]++;
			zp[j] += 1;
		}
		maps[i] = map;
		for(j = 0; j < tpara.maxK; j++)
		{	if(cn[j] == 0 && remap[j] >= 0)
			{	remap[j] = map[remap[j]] = -1;	}
			cn[j] = 0;
		}
	}
	newmaxK++;

	float *nhp = tpara.nh, *newnh = new float[newmaxK * mydata.L], *newnhp = newnh;
	int const *aszp = mydata.asz;
	for(i = 0; i < newmaxK * mydata.L; i++, newnhp++) *newnhp = 0;
	newnhp = newnh;
	vector<vector<int> >::iterator imaps = maps.begin();
	for(i = 0; i < mydata.L; i++, nhp += tpara.maxK, newnhp += newmaxK, aszp++, imaps++)
	{	vector<PARAUNIT> tparam = tpara.param[i];
		PARAUNIT pu;
		pu.lapseN = pu.popN = 0;
		tpara.param[i].clear(); tpara.param[i].resize(newmaxK);//, pu);
		vector<int>::iterator im = (*imaps).begin();
		vector<PARAUNIT>::iterator itpam = tparam.begin();
		for(j = 0; j < tpara.maxK; j++, im++, itpam++)
			if((*im) >= 0)
			{	tpara.param[i][(*im)] = (*itpam);
				*(newnhp +(*im)) = *(nhp + j);
			}
	}
	delete tpara.nh;
	tpara.nh = newnh;

	for(j = 0; j < newmaxK; j++) tpara.Q[j] = 0;
        bool *rrr = tpara.rec;
        int *ppp = tpara.pop;
        for(k = 0; k < mydata.indN; k++)
        {       for(int l = 0; l < mydata.L; l++, rrr++, ppp++)
                {       if(l==0 || *(rrr))
				tpara.Q[*ppp]+=tpara.idw[k];
                }
        }
	tpara.maxK = newmaxK;
}

void tensorHMMbase::_greedySwitch_new(int indN, int L, int ploidity, TENSORPARA &tpara, vector<bool> const &init)
{
	if(tpara.maxK < 2 || tpara.singlecall || tpara0.popn0 > 0) return;
	int i, j, k, tsz = tpara.maxK * tpara.maxK;
	float *trans = new float[tsz];
	int perm[tpara.maxK];
	bool flag = false;

	for(i = 0; i < tpara.maxK; i++) perm[i] = i;
	float *nhp = tpara.nh + tpara.maxK;
	int *pp[indN], pbak[indN];
	bool *rr[indN];
	for(i = 0; i < indN; i++) 
	{	pp[i] = tpara.pop + i * L + 1;
		rr[i] = tpara.rec + i * L + 1;
		pbak[i] = *(tpara.pop + i * L);
	}
	for(i = 1; i < L; i++, nhp += tpara.maxK)
	{	bool change = false;
		if(tpara.popr[i] > 0.05)
		{	
			int revmap[tpara.maxK], revbak[tpara.maxK], operm[tpara.maxK];
			for(j = 0; j < tpara.maxK; j++) 
			{	revbak[perm[j]] = revmap[perm[j]] = j;
				operm[j] = perm[j];
			}
			int const *pp1[indN], *pp2[indN];
			for(j = 0; j < indN; j++) 
			{	pp1[j] = tpara.pop + j * L + i-1;
				pp2[j] = tpara.pop + j * L + i;
			}
			for(j = 0; j < tsz; j++) trans[j] = 0;
			for(j = 0; j < indN; j++) 
				if((int)init.size() == 0 || !init[j / ploidity]) trans[(*pp1[j]) * tpara.maxK + (*pp2[j])]+=tpara.idw[j];
			_greedyMatch(trans, tpara.maxK, revmap);

			bool check = false;
			for(j = 0; j < tpara.maxK; j++) { if(revbak[j] != revmap[j]) check = true; }
			if(check)
			{	double lp0 = 0, lp1 = 0;
				double r0 = 0, r1 = 0, nnn = 0;
				for(j = 0; j < indN; j++)
					if(!init[j / ploidity])
					{	if(*rr[j]) { r0+=tpara.idw[j]; lp0 += log(tpara.Vh[*pp[j]]/* * tpara.idw[j]*/); }
						if(revmap[*pp[j]] != *(pp[j]-1) || (*pp[j] == pbak[j] && *rr[j])) { r1+=tpara.idw[j]; lp1 += log(tpara.Vh[revmap[*pp[j]]]/**tpara.idw[j]*/); }
						nnn+=tpara.idw[j];
					}
				//lp0=lp1=0;
				lp0 += gsl_sf_lngamma((double)r0 + 1.) + gsl_sf_lngamma((double)nnn - r0 + 1000.);
				lp1 += gsl_sf_lngamma((double)r1 + 1.) + gsl_sf_lngamma((double)nnn - r1 + 1000.);
				double maxp = max(lp0, lp1);
				if(gsl_ran_flat(gammar, 0, 1.) < exp(lp1 - maxp) / (exp(lp1 - maxp) + exp(lp0 - maxp))) 
				{	flag = change = true;
					for(j = 0; j < tpara.maxK; j++) perm[revmap[j]] = j;
				}
			}
		}
		if(flag)
		{	vector<PARAUNIT> tparam = tpara.param[i];
			float tnh[tpara.maxK];
			for(j = 0; j < tpara.maxK; j++) 
			{	tpara.param[i][perm[j]] = tparam[j];	
				tnh[j] = *(nhp + j);
			}
			for(j = 0; j < tpara.maxK; j++) 
				*(nhp + perm[j]) = tnh[j];

			for(j = 0; j < indN; j++)
				if(!init[j / ploidity])
				{	if(*pp[j] != pbak[j]) *rr[j] = false;
					pbak[j] = *pp[j];
					*pp[j] = perm[*pp[j]];
					if(*pp[j] != *(pp[j]-1)) *rr[j] = true;
				}
		}
		else for(j = 0; j < indN; j++) pbak[j] = *pp[j];
	
		for(j = 0; j < indN; j++) 
		{	pp[j] += 1;
			rr[j] += 1;
		}
	}

	delete []trans;

	for(j = 0; j < tpara.maxK; j++) tpara.Q[j] = 0;
        for(j = 0; j < L; j++) tpara.recombN[j] = 0;
	bool *rrr;
	int *ppp;
        for(k = 0; k < indN; k++)
        	if(!init[k / ploidity])
		{       rrr = tpara.rec + k * L;
			ppp = tpara.pop + k * L;
			for(int l = 0; l < L; l++, rrr++, ppp++)
                	{       if(l==0 || *(rrr))
                        	{       tpara.Q[*ppp]+=tpara.idw[k];
                                	tpara.recombN[l]+=tpara.idw[k];
                        	}
                	}
        	}
}

void tensorHMMbase::_greedyMatch(float *trans, int maxK, int *perm)
{
	int j, k, l;
	j = 1;
	bool a_inner = true;
	double alpha = 0.;//checkthis
	while(a_inner)
	{	if(j > 0)
		{	for(k = 0; k < j; k++)
			{	if(trans[j*maxK+perm[j]] + trans[k*maxK+perm[k]] + alpha * ((double)(j==perm[j]) + (double)(k==perm[k]))
					< trans[j*maxK+perm[k]] + trans[k*maxK+perm[j]] + alpha * ((double)(j==perm[k]) + (double)(k==perm[j])))
				{	
					l = perm[j];
					perm[j] = perm[k];	
					perm[k] = l;
					j = k - 1;
					break;
				}
			}
		}
		j++;
		if(j >= maxK) a_inner = false;
	}
}

void tensorHMMbase::_LongtoShort(int *z[], bool const *rec, vector<vector<int> > &zshort, int indN, int L, int st)
{	int i, j;
	zshort.clear();
	zshort.resize(indN);
	for(i = 0; i < indN; i++)
	{	zshort[i].push_back(0 + st);
		zshort[i].push_back(z[i][0]);
		for(j = 1; j < L; j++)
		{	
			if((rec != NULL && rec[j]) || (rec == NULL && z[i][j] != z[i][j - 1]))
			{	zshort[i].push_back(j + st);
				zshort[i].push_back(z[i][j]);
			}
		}
		zshort[i].push_back(L + st);
		zshort[i].push_back(-1);
	}
}

double tensorHMMbase::_logP(MYDATA const &mydata, TENSORPARA const &tpara)
{	double lp = 0, recp;
	int i, j, k;;
	int const *pp, *pp2 = NULL;
	float const *dp, *dp2 = NULL;
	float const *nhp; 
	int const *aszp;
	double const *rp;

//double lp1[20] = {0};
//int nnn[20] = {0};

	for(i = 0; i < mydata.indN; i+=mydata.ploidity)
	{	if((int)mydata.fixState.size() > 0 && (int)mydata.fixState[i].size() > 0) continue;
		pp = tpara.pop + i * mydata.L;
		dp = mydata.data + (mydata.indIndex[(i % mydata.indN)/mydata.ploidity]*mydata.ploidity) * mydata.L;
		
		if(mydata.ploidity == 2) 
		{	pp2 = pp + mydata.L;
			dp2 = dp + mydata.L;
		}
		vector<vector<PARAUNIT> >::const_iterator ipm = tpara.param.begin();
		nhp = tpara.nh;
		aszp = mydata.asz;
		rp = tpara.popr;

		int flipi = 0, flipsz = 0;
		bool curflip = false;
		if((int)tpara.flips.size() > 0) flipsz = (int)tpara.flips[(i%mydata.indN) / mydata.ploidity].size();
		int ed = mydata.L;
		vector<int> tbreaks;
		if((int)mydata.breaks.size() > 0)
		{	tbreaks = mydata.breaks[i / mydata.ploidity];
			ed = (int)tbreaks.size() / 2 - 1;
		}
		for(j = 0; j < ed; j++)
		{	int stpos = j;
			int edpos = j + 1;
			if((int)mydata.breaks.size() > 0)
			{	stpos = mydata.breaks[i / mydata.ploidity][j * 2];
				edpos = mydata.breaks[i / mydata.ploidity][j * 2 + 2];
			}
			while(flipi < flipsz && tpara.flips[(i%mydata.indN)/mydata.ploidity][flipi] < stpos) flipi+=2;
			if(flipsz > 0 && stpos == tpara.flips[(i%mydata.indN)/mydata.ploidity][flipi]) 
				curflip = tpara.flips[(i%mydata.indN)/mydata.ploidity][flipi+1];

			recp =  *rp + 1e-10;
			if(j == 0) lp += log(tpara.Vh[*pp]/* * tpara.idw[i]*/);
			else if(*pp == *(pp - 1)) lp += log(1. - recp + recp * tpara.Vh[*pp]/* * tpara.idw[i]*/);
			else lp += log(recp * tpara.Vh[*pp]/* * tpara.idw[i]*/);
			if(mydata.ploidity == 2)
			{	if(j == 0) lp += log(tpara.Vh[*pp2]/* * tpara.idw[i]*/);
				else if(*pp2 == *(pp2 - 1)) lp += log(1. - recp + recp * tpara.Vh[*pp2]/* * tpara.idw[i]*/);
				else lp += log(recp * tpara.Vh[*pp2]/* * tpara.idw[i]*/);
			}
		//	for(k = 0; k < ided - idst; k++)
			k=0;
			{	double f = _getLp(i, tbreaks, j * (1 + (int)((int)mydata.breaks.size() > 0)), curflip, mydata.ploidity, tpara, dp + mydata.L * mydata.ploidity * k, dp2 + mydata.L * mydata.ploidity * k, pp, pp2, ipm, /*tpara.param.begin() + stpos, */nhp, aszp);
//if(f > -10000. || f < 10000.);else printf("!!!%d:%d,%f\n", i, j, f);
				if(tpara.logscale) lp += f;
				else lp += log(f);
/*if(gID > 0)
{
nnn[(int)dp[j]]++;
lp1[(int)dp[j]]+=f;
}*/

			}

			dp += edpos - stpos;
			pp += edpos - stpos;
			if(mydata.ploidity == 2)
			{	dp2 += edpos - stpos;
				pp2 += edpos - stpos;
			}
			ipm += edpos - stpos;
			nhp += tpara.maxK * (edpos - stpos);
			aszp += (edpos - stpos);
			rp += (edpos - stpos);
		}
	}
//for(i = 0; i < 20; i++) printf("(%d:%d,%f)", i, nnn[i], lp1[i]/((double)nnn[i]+1e-10));

	return lp;
}

double tensorHMMbase::_imputeData(MYDATA const &mydata, TENSORPARA const &tpara, int id, vector<int> const &dbreaks, bool init)
{	
	int j, k;
	int const *pp, *pp2 = NULL;
	float *dp, *dp2 = NULL;
	float const *nhp;
	int const *aszp;
	
	pp = tpara.pop + id * mydata.L;
	dp = mydata.data + mydata.indIndex[id/mydata.ploidity] * mydata.ploidity * mydata.L;
	if(mydata.ploidity == 2) 
	{	pp2 = pp + mydata.L;
		dp2 = dp + mydata.L;
	}
	vector<vector<PARAUNIT> >::const_iterator ipm = tpara.param.begin();
	nhp = tpara.nh;
	aszp = mydata.asz;

	int ed = mydata.L;
	vector<int> tbreaks;
	if((int)mydata.breaks.size() > 0) tbreaks = mydata.breaks[id / mydata.ploidity];
	if((int)dbreaks.size() > 0) ed = (int)dbreaks.size() - 1;
	vector<int>::const_iterator idb = dbreaks.begin();
	
//clock_t cst, ced;
//cst=clock();
	double RT = 0;
	for(j = 0; j < ed; j++, idb++)
	{	int stpos = j;
		int edpos = j + 1;
		int flip = 0;
		if((int)mydata.breaks.size() > 0)
		{	stpos = mydata.breaks[id / mydata.ploidity][*idb];
			edpos = mydata.breaks[id / mydata.ploidity][*(idb+1)];
		}
		if((int)tpara.flips.size() > 0) flip = tpara.flips[id/mydata.ploidity][j*2+1];
		int ttttt=0;
		//for(k = 0; k < ided - idst; k++,ttttt+=tttt)
		k=0;
		{	

			if((int)mydata.fixAllele.size() > 0 && (int)mydata.fixAllele[(id+k)].size() > 0 && !mydata.fixAllele[(id+k)][j]) ;
			else 
			{	RT += _imputeOne(id, tbreaks, *idb, *(idb+1), flip, mydata.ploidity, tpara, dp+ttttt, dp2+ttttt, pp, pp2, ipm, /*tpara.param.begin() + stpos, */nhp, aszp);
			}
		}

		dp += edpos - stpos;
		pp += edpos - stpos;
		if(mydata.ploidity == 2)
		{	dp2 += edpos - stpos;
			pp2 += edpos - stpos;
		}
		ipm += edpos - stpos;
		nhp += tpara.maxK * (edpos - stpos);
		aszp += (edpos - stpos);
	}
//ced=clock();
//printf("impute of %d = %fsec\n", mydata.L, (double)(ced-cst)/1000000.);
	return(RT);
}

void tensorHMMbase::_getDoubleBreaks(vector<int> const &breaks, int L, vector<int> &dbreaks)
{	int i;
	dbreaks.clear();//first of a segment, right to break point, need to contain start position 0
	if((int)breaks.size() > 0)
	{	int bsz = (int)breaks.size();
		for(i = 0; i < bsz; i += 2) 
			if(breaks[i + 1] == 2) dbreaks.push_back(i);
	}
	else 
	{	dbreaks.resize(L + 1);
		for(i = 0; i < L + 1; i++) dbreaks[i] = i;
	}
}

void tensorHMMbase::_readPop(char const *fname, vector<vector<int> > &pop)
{
	pop.clear();
	string str;
	ifstream inFile(fname);
	if(!inFile.is_open()) printf("!!!cannot open file %s\n", fname),exit(0);
	int i, id = 0;
	while(true)
	{	inFile >> str;
		if(inFile.eof()) break;
		if(str[0] == 'i' && str[1] == 'd')
		{	id = atoi(&str[2]);
			if(id < (int)pop.size()) pop[id].clear();
			else pop.resize(id + 1);
		}
		else
		{	for(i = 0; i < str.size(); i++) if(str[i] == ':') break;
			int pos = atoi(&str[0]);
			int state = atoi(&str[i + 1]);
			pop[id].push_back(pos);
			pop[id].push_back(state);
		}
	}
	inFile.close();
}

void tensorHMMbase::_outputPop(char const *fname, int const *pop, int N2, int L, bool const *rec)
{	
	int j, k;
	FILE *f = fopen(fname, "w");
	for(j = 0; j < N2; j++)
       	{       int const *pp = pop + j * L;
		bool const *rr = rec + j * L;
		fprintf(f, "id%d ", j);
                fprintf(f, "0:%d ", *pp);
		pp++;rr++;
               	for(k = 1; k < L; k++,pp++,rr++)
       			//if(*pp != *(pp-1) || (*rr)) fprintf(f, "%d:%d ", k, *pp);
       			if(*pp != *(pp-1)) fprintf(f, "%d:%d ", k, *pp);
                fprintf(f, "\n");
       	}
        fclose(f);
}

void tensorHMMbase::_outputPops(char const *fname, vector<vector<vector<int> > > const &pops)
{	
	int i, j, k;
	FILE *f = fopen(fname, "w");
	for(i = 0; i < (int)pops.size(); i++)
		for(j = 0; j < (int)pops[i].size(); j++)
       		{	fprintf(f, "id%d ", j);
			for(k = 0; k < (int)pops[i][j].size()-2; k += 2)
       				fprintf(f, "%d:%d ", pops[i][j][k], pops[i][j][k+1]);
			fprintf(f, "\n");
       		}
        fclose(f);
}

void tensorHMMbase::_outputPQC(char const *fname, vector<vector<vector<int> > > const &pops, int N2, int L, int mMaxK)
{
	int i, j, k;
	int szx = (int)pops.size();
	int idx[szx * N2];
	vector<vector<int> > cnx(N2, vector<int>(mMaxK, 0));
	FILE *f = fopen(fname, "w");
	for(i = 0; i < L; i++)
	{	int *idxp = idx;
		if(i == 0)
		{	for(j = 0; j < szx; j++)
				for(k = 0; k < N2; k++, idxp++)
				{	cnx[k][pops[j][k][1]]++;
					(*idxp) = 2;
				}
		}
		else
		{	for(j = 0; j < szx; j++)
				for(k = 0; k < N2; k++, idxp++)
				{	if(pops[j][k][(*idxp)] <= i)
					{	cnx[k][pops[j][k][(*idxp)-1]]--;
						cnx[k][pops[j][k][(*idxp)+1]]++;
						(*idxp) += 2;
					}
				}
		}
		for(j = 0; j < N2; j++)
			for(k = 0; k < mMaxK; k++)
				fprintf(f, "%d ", cnx[j][k]);
        fprintf(f, "\n");
	}
	fclose(f);	
}

void tensorHMMbase::summarizePosterior(vector<vector<vector<int> > > &pops, vector<vector<int> > &breaks, int *z, int ploidity, double const *r, int L, int maxK)
{
	int i, j, k, l, t;
	int sz = (int)pops.size();
	if(sz == 0) return;

	int indN = (int)pops[0].size();
	for(i = 0; i < indN; i += ploidity)
	{	k = (int)pow((double)maxK, (double)ploidity);
		vector<vector<int> > cn(L, vector<int>(k, 0));
		for(j = 0; j < sz; j++)
		{	vector<vector<int> > state(ploidity, vector<int>(L, -1));
			for(l = 0; l < ploidity; l++)
			{	for(k = 2; k < (int)pops[j][i + l].size(); k += 2)
				{	for(t = pops[j][i + l][k - 2]; t < pops[j][i + l][k]; t++)
						state[l][t] = pops[j][i + l][k - 1];
				}
			}
			if(ploidity == 1)
			{	for(l = 0; l < L; l++) cn[l][state[0][l]]++;
			}
			else
			{	int d = j * indN / ploidity + i / ploidity;
				if((int)breaks[d].size() == 0)
				{	for(t = 0; t < L; t++)
						cn[t][state[0][t] * maxK + state[1][t]]++;
				}
				else
				{	t = 0;
					for(k = (int)breaks[d].size() - 1; k >= 0; k--)
					{	int s1 = 0, s2 = 0, o = t;
						for(t = o; t <= breaks[d][k]; t++)
						{	s1 += cn[t][state[0][t] * maxK + state[1][t]];
							s2 += cn[t][state[1][t] * maxK + state[0][t]];
						}
						if(s1 >= s2) l = 0;
						else l = 1;
						
						for(t = o; t <= breaks[d][k]; t++)
							cn[t][state[l][t] * maxK + state[1-l][t]]++;
					}
				}
			}
		}
		t = 0;
		int *zp = z + i * L;
		for(j = 0; j < L; j++, zp++)
		{	for(k = 0; k < (int)cn[j].size(); k++) if(cn[j][k] > cn[j][t]) t = k;
			l = t;
			for(k = ploidity - 1; k >= 0; k--)
			{	*(zp + k * L) = l % maxK;
				l /= maxK;
			}
		}	
	}
}

void tensorHMMbase::summarizeHaps(int indN, int L, int ploidity, vector<vector<char> > const &fullsample, vector<vector<int> > const &cumsample, bool singlecall, vector<vector<char> > &haps)
{
	int i, j, k, l;

	haps.clear();
	haps.resize(indN / ploidity * 2, vector<char>(L, 0));
	for(j = 0; j < indN / ploidity; j++)
	{	vector<int> hets(L);
		int hl = 0;
		for(i = 0; i < L; i++)
		{	k = 0;
			int sss = cumsample[j][i * (2 + ploidity)];
			for(l = 1; l < 2 + ploidity; l++)
			{	if(cumsample[j][i * (2 + ploidity) + l] > cumsample[j][i * (2 + ploidity) + k]) k = l;
				sss += cumsample[j][i * (2 + ploidity) + l];
			}
			if(ploidity == 2)
			{	if(cumsample[j][i * 4] >= (cumsample[j][i * 4 + 1] + cumsample[j][i * 4 + 2])/1. && cumsample[j][i * 4] >= cumsample[j][i * 4 + 3]) k = 0;
				else if(cumsample[j][i * 4 + 3] >= (cumsample[j][i * 4 + 1] + cumsample[j][i * 4 + 2])/1. && cumsample[j][i * 4 + 3] >= cumsample[j][i * 4]) k = 3;
				else if(cumsample[j][i * 4 + 1] >= cumsample[j][i * 4 + 2]) k = 1;
				else k = 2;
			}

			haps[j * 2][i] = (char)(k / 2);
			if(ploidity == 1) haps[j * 2 + 1][i] = (char)(k > 0);
			else 
			{	haps[j * 2 + 1][i] = (char)(k % 2);
				if(k > 0 && k < 3) hets[hl++] = i;
			}
		}
		if(ploidity == 2 && hl > 1 && !singlecall)//look for most frequent pairs
		{	int ssz = (int)fullsample.size() / indN;
			int minBound = (int)((double)ssz * 0.5 + 0.9999);
			vector<int> hetcode(ssz, 0);
			vector<int> hcode(ssz, 0);
			int cc = 0, hst = 0;
			for(i = 0; i < hl; i++)
			{	int i1 = hets[i] / 8, i2 = 7 - hets[i] % 8, cn=0;
				int hsz = (int)pow(2., (double)(i-hst+1));
				vector<int> hcn(hsz / 2, 0);
				for(k = 0; k < ssz; k++)
				{	int a1 = (int)((fullsample[k * indN + j * 2][i1] & (1 << i2)) >> i2);
					int a2 = (int)((fullsample[k * indN + j * 2 + 1][i1] & (1 << i2)) >> i2);
					hcode[k] = hcode[k] * 2 + a1;
					hetcode[k] = hetcode[k] * 2 + (int)(a1 != a2);
					if(hetcode[k] == hsz - 1) cn++;
					if(hetcode[k] == hsz - 1) 
						hcn[min(hcode[k], hsz - 1 - hcode[k])]++;
				}
				while(((cn < minBound || hcn[min(cc * 2 + 1,hsz-1-cc*2-1)]+hcn[min(cc * 2,hsz-1-cc*2)] < minBound / 2 + 1) && hst < i - 1) || i - hst > 5)
				{	cn = 0;
					hsz /= 2;
					cc %= (hsz / 2);
					hst++;
					hcn.clear(); hcn.resize(hsz / 2, 0);
					for(k = 0; k < ssz; k++) 
					{	hcode[k] %= hsz;
						hetcode[k] %= hsz;
						if(hetcode[k] == hsz - 1) cn++;
						if(hetcode[k] == hsz - 1) 
							hcn[min(hcode[k], hsz - 1 - hcode[k])]++;
					}
				}

				haps[j * 2][hets[i]] = haps[j * 2 + 1][hets[i]] = 0;
				if(((double)hcn[min(cc * 2 + 1,hsz-1-cc*2-1)] + 0.5) + (double)cumsample[j][hets[i]*4+2]/1000. > ((double)hcn[min(cc * 2,hsz-1-cc*2)] + 0.5) + (double)cumsample[j][hets[i]*4+1]/1000.)
				{	haps[j * 2][hets[i]] = 1;
					cc = cc * 2 + 1;
				}
				else 
				{	haps[j * 2 + 1][hets[i]] = 1;
					cc = cc * 2;
				}
			}
		}	
	}
}

int tensorHMMbase::_sample(double *prob, int sz, char maximum)
{
	int i, j;
	if(maximum != 0)
	{	j = 0;
		for(i = 1; i < sz; i++) if(prob[j] < prob[i]) j = i;
	}
	else
	{   double tprob[sz];
		if(Temp > 1.)
		{	tprob[0] = pow(prob[0], 1./Temp);
			for(i = 1; i < sz; i++) tprob[i] = prob[i] + pow(tprob[i - 1],1./Temp);
		}
		else
		{	tprob[0] = prob[0];
			for(i = 1; i < sz; i++) tprob[i] = prob[i] + tprob[i - 1];
		}
		double un = gsl_ran_flat(gammar, 0, tprob[sz - 1]);
        for(j = 0; j < sz; j++) if(un <= tprob[j]) break;
	}
	return j;
}

void tensorHMMbase::_checkCounts(TENSORPARA const &tpara)
{	int i, j, L = (int)tpara.param.size();
	float k;
	for(i = 0; i < L; i++)
	{	for(j = 0; j < tpara.maxK; j++)
		{	k = 0;
			for(vector<float>::const_iterator ipm = tpara.param[i][j].allele.begin(); ipm < tpara.param[i][j].allele.end(); ipm+=2)
			{	k += (*(ipm+1));
			}
			if(fabs(k - tpara.nh[i * tpara.maxK + j])>NUMPRECISION)
				printf("!!!%d:%d %f vs %f\n", i, j, k, tpara.nh[i*tpara.maxK+j]),exit(0);
		}
	}
}

void tensorHMMbase::initPara(TENSORPARA &tpara)
{
	tpara.rPrior = NULL;
	tpara.priorW = 10.;
	tpara.probAlpha = 1.;
	tpara.defmut = 0.5;
	tpara.maxHapK = 10;
	tpara.A = 1.;
	tpara.burnin = 50;
	tpara.mcmc = 50;
	tpara.shrinkclass = false;
	tpara.singlecall = false;
	tpara.imputedata = true;
	tpara.addsample = true;
	tpara.samplefull = false;
	tpara.changeAlpha = 1;
	tpara.splitstop = 40;
	tpara.refweight = 5;
	tpara.maximum = 0;
	tpara.samplemaximum = false;
	tpara.logscale = true;
	tpara.recrate = 1.;
	tpara.recombcut = 0;
	tpara.heteroVh = 0;

	tpara.pop = NULL;
	tpara.rec = NULL;
	tpara.lapse = NULL;
	tpara.nh = NULL;
	tpara.recombN = NULL;
}

//this function is only designed for ploidity=1
int tensorHMMbase::_matchPopsample(TENSORPARA const &tpara, vector<vector<vector<int> > > &pops, int *z, int maxK, vector<int> &pqc)
{
printf("maxK=%d\n", maxK);fflush(stdout);
	if(maxK < 2) return maxK;
	int i, j, k, l, t;
	int sz = (int)pops.size();
	if(sz == 0) return maxK;
	int indN = (int)pops[0].size();

	vector<int> breaks;
	for(i = 0; i < sz; i++)
		for(j = 0; j < indN; j++)
			for(k = 0; k < (int)pops[i][j].size(); k += 2)
				if(k == 0 || pops[i][j][k + 1] != pops[i][j][k - 1]) breaks.push_back(pops[i][j][k]);
	sort(breaks.begin(), breaks.end());
	for(i = (int)breaks.size() - 1; i > 0; i--)
		if(breaks[i] == breaks[i - 1]) breaks.erase(breaks.begin() + i);
	int bsz = (int)breaks.size();
	int L = breaks[bsz - 1];

	int isz = sz * indN;
	int tsz = maxK * maxK;
	float trans[tsz];
	vector<int> index(isz, 0);
	int *newpop = new int[indN * (bsz - 1)];
	int grevmap[maxK], gperm[maxK];;
	for(k = 0; k < maxK; k++) grevmap[k] = gperm[k] = k;

	pqc.clear();
	pqc.resize(bsz * (1 + indN * maxK), 0);
	for(i = 0; i < bsz - 1; i++)
	{	
//printf("%d,",i);fflush(stdout);
		int states[isz];
		l = 0;
		for(j = 0; j < sz; j++)
		{	for(k = 0; k < indN; k++, l++)
			{	while(index[l] < (int)pops[j][k].size() && pops[j][k][index[l]] < breaks[i]) index[l] += 2;
				if(index[l] >= (int)pops[j][k].size() || pops[j][k][index[l]] > breaks[i]) index[l] -= 2;
				states[l] = pops[j][k][index[l] + 1];
			}
		}

		for(j = sz - 2; j >= 0; j--)
		{	for(k = 0; k < tsz; k++) trans[k] = 0;
			int sj = indN * j, sk = (sz - 1) * indN;
			for(k = sz - 1; k > j; k--, sk -= indN)	
			{	for(l = 0; l < indN; l++)
					trans[states[sk+l] * maxK + states[sj+l]]+=tpara.idw[l];
			} 
			int revmap[maxK], perm[maxK];
			for(k = 0; k < maxK; k++) revmap[k] = k;
			_greedyMatch(trans, maxK, revmap);
			for(k = 0; k < maxK; k++) perm[revmap[k]] = k;
			
			int ts[indN];
			for(k = 0; k < indN; k++) ts[k] = states[sj + k];
			for(k = 0; k < indN; k++) states[sj + k] = perm[ts[k]];
		}
		for(j = 0; j < indN; j++)
		{	vector<double> cn(maxK, 0);
			if(i > 0) cn[gperm[newpop[(i-1)*indN+j]]] = 0.1;
			for(k = 0; k < sz; k++) cn[states[k * indN + j]]++;
			l = 0;
			for(k = 1; k < maxK; k++) if(cn[l] < cn[k]) l = k;
			newpop[i * indN + j] = l;
		}
		
		if(i > 0)
		{	for(k = 0; k < tsz; k++) trans[k] = 0;
			t = i * indN;
			for(l = 0; l < indN; l++)
				trans[newpop[t-indN+l] * maxK + newpop[t+l]]+=tpara.idw[l];
			_greedyMatch(trans, maxK, grevmap);
			for(l = 0; l < maxK; l++) gperm[grevmap[l]] = l;
			int ts[indN];
			for(k = 0; k < indN; k++) ts[k] = newpop[t + k];
			for(k = 0; k < indN; k++) newpop[t + k] = gperm[ts[k]];
		}

		pqc[i * (1 + indN * maxK)] = breaks[i];
		for(j = 0; j < indN; j++)
		{	for(k = 0; k < sz; k++) pqc[i + 1 + (i * indN + j) * maxK + gperm[states[k*indN+j]]]++;
		}
	}
	printf("]\n");fflush(stdout);
	pqc[i * (1 + indN * maxK)] = breaks[i];
	for(j = 0; j < indN * maxK; j++)
	{	pqc[i * (1 + indN * maxK) + 1 + j] = pqc[(i - 1) * (1 + indN * maxK) + 1 + j];
	}

	int rt = 0;
	int *pp = newpop;
	for(i = 0; i < bsz - 1; i++)
	{	//printf("%d;",i);fflush(stdout);
		int *zp = z + breaks[i];
		for(j = 0; j < indN; j++, zp += L, pp++)
		{	int *zpp = zp;
			for(k = breaks[i]; k < breaks[i + 1]; k++, zpp++)
			{	*zpp = *pp;
				rt = max(rt, *pp);
			}
		}
	}
	rt++;

printf("}\n");fflush(stdout);
	delete []newpop;
printf("deleted\n");fflush(stdout);
	return rt;
}
