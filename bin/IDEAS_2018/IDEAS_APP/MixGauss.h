#pragma once
#include "datastructure.h"

class MixGauss
{	public:
		MixGauss(void);
		~MixGauss(void);
	private:
	public:
		double *gausspara0;
		int clustersz0;

		double *gausspara, *gaussprior, *lambda;
		int gausssz, lambdasz, groupn, *groupcode, *coderevmap;
		bool lessoneUpdate;
		const gsl_rng_type *T;
		gsl_rng *gammar;
		int maxxmsz, *xmsz, maxxvsz, maxymsz, *ymsz, maxyvsz, maxyxsz, **xmap, **ymap, *mapspace;
		vector<int> ymap0, xmap0;
		int ymsz0, xmsz0, code0;
		int *neighbor;	
		vector<double> markmean, marksd;
		float *preLP;

	public:
		int clustersz, totalN0;
		double minerr, maxerr, tmpprop[10000];
		bool indc, nb;
		vector<vector<double> > modelparameter; 

	public:
		void computeLP(float *ydata, float *xdata, int id, double priorW, double *lp);
		int* computeLP_subset(float *ydata, float *xdata, int id, int state, double priorW, double *lp);

		void addPara(int g, float *yy, float *xx, int id, float wt);
		void removePara(int g, float *yy, float *xx, int id, float wt);
		void getStateCount(double *cn);
		void getStatePrior(double *&ppp, int &step, double A, double B);

		void updateLambda(double priorW, bool updateprior = true, bool updateNeighbor = true);
		void initializePara(int clusterSZ, float **dataYP, float **dataXP, int totalN, int L, int maxysz, int *ysz, int **ymp, int maxxsz, int *xsz, int **xmp, double priorW, char const *parafile = NULL);
		void clearParameter();
		void rearrangeParameter(int *remap, int newclustersz, float **dataYP, float **dataXP, int totalN, int L, double A);
		void outputParameter(char *fout, double priorW, vector<string> const &fy, vector<string> const &fx);
		void updateParameterPrior(double A);
		
		void simData(int n, int id, int g, float *yp);

		double splitmergeCluster(int type, float **dataYP, float **dataXP, int totalN, int L, float *states, bool *lapse, double priorW, int &mi, int &mj, int &tid, float const *wt);

		void imputeMissingTrack(int id, int st, int ed, float *yp, float *xp, float *yimp, float *ximp, float *states, double diag = 0, bool verbose = false);//double *stateprob);

		void loadGaussPrior0(char const *parafile, vector<string> const &fy, vector<string> const &fx);
		void imputeMissing2(MYDATA const &mydata, TENSORPARA const &tpara, float **dataYP, float **dataXP, double priorW, char const *fname, int gid, float const *wt);
		void preComputeLP(char const *fmixpara, vector<string> const &fy, vector<string> const &fx, float **dataYP, float **dataXP, int totalN, int L);

	private:
		double _dmvnorm(double *lambdap, float *yp, float *xp, int tymsz, int txmsz);
		void _getLambdaOne(double *ll, double *rr, double *tmpspace, double priorW, int i, int id);
		void _cholV(double const *A, int n, double *L);
		void _invL(double const *L, int n, double *iL, bool transpose);
		double _lgammaM(int q, double x);
		void _prepareGroupCode(int totalN);
		void _getMean(double const *pp, float *my, float *mx);
		void _gaussEM(double *pp, double *rr, double priorW, double *rr0);
		void _collapsedPrediction(double const *pp, double const *phi, int id, double *rt, double *tmpspace, int &misyn, int *misy, int &misxn, int *misx, int *map);

		void _XtX(double const *X, int rn, int cn, double *X2);
		void _AtB(double const *A, double const *B, int ca, int cb, int n, double *C);
		void _ABt(double const *A, double const *B, int ra, int rb, int n, double *C);
		void _AB(double const *A, double const *B, int ra, int cb, int n, double *C);

		void _updateVh(double *Q, int K, double A, double *V);
		void _getNeighbor(double priorW);
		void _loadGaussPrior_mix(char const *fmixpara, vector<string> const &fy, vector<string> const &fx, double *&prepara, vector<int> &state);
		void indDist(MYDATA const &mydata, TENSORPARA const &tpara, double *wt);
};
