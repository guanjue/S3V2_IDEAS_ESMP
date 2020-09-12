#pragma once
#include "datastructure.h"

class tensorHMMbase
{
public:
	tensorHMMbase(unsigned int rseed = 0);
	~tensorHMMbase(void);

	const gsl_rng_type *T;
	gsl_rng *gammar;
	clock_t forwardT, backwardT, imputeT;	
	int gID;
	bool outputproportion;
	double Temp;
	TENSORPARA0 tpara0;
	double *tmpM;
///////////////////////////////////////////////////
public: 
	void inferStructure(MYDATA &mydata, TENSORPARA &tpara, char const *outname, bool outqc, vector<bool> const &myinit = vector<bool>());
	void initPara(TENSORPARA &tpara);

protected:
	double _oneRound(MYDATA &mydata, TENSORPARA &tpara, vector<bool> &init, int iter);
	double _updateStructure(MYDATA &mydata, TENSORPARA &tpara, int id, int ploidity, vector<int> const &dbreaks);
	void _forwardSum(MYDATA &mydata, TENSORPARA &tpara, int id, int ploidity, vector<vector<int> > const &ss, vector<vector<int> > const &mlist, vector<int> const &dbreaks, double *Vh, double *Vh2, double *P); 
	double _backwardSample(MYDATA &mydata, TENSORPARA &tpara, int id, int ploidity, vector<vector<int> > const &ss, vector<vector<int> > const &mlist, vector<int> const &dbreaks, double *Vh, double *Vh2, double *P);

	double _logP(MYDATA const &mydata, TENSORPARA const &tpara);

	void _LongtoShort(int *z[], bool const *rec, vector<vector<int> > &zshort, int indN, int L, int st = 0);
	void _readPop(char const *fname, vector<vector<int> > &pop);
	void _outputPop(char const *fname, int const *pop, int N2, int L, bool const *rec);
	void _outputPQC(char const *fname, vector<vector<vector<int> > > const &pops, int N2, int L, int mMaxK);
	void summarizePosterior(vector<vector<vector<int> > > &pops, vector<vector<int> > &breaks, int *z, int ploidity, double const *r, int L, int maxK);
	void summarizeHaps(int indN, int L, int ploidity, vector<vector<char> > const &fullsample, vector<vector<int> > const &cumsample, bool singlecall, vector<vector<char> > &haps);
	int _sample(double *prob, int sz, char maximum);
	void _outputPops(char const *fname, vector<vector<vector<int> > > const &pops);
        void _updateVh(double *Q, int K, double A, double *V, bool verbose = false);
	int _matchPopsample(TENSORPARA const &tpara, vector<vector<vector<int> > > &pops, int *z, int maxK, vector<int> &pqc);

	virtual void _getAsz(MYDATA const &mydata)=0;
	virtual void _updatePara_remove(PARAUNIT &para, float data, float b, int id, int pos)=0;
	virtual void _updatePara_add(PARAUNIT &para, float data, float b, int id, int pos)=0;
	virtual void _sumforwardOne(int id, vector<int> const &breaks, int bst, int bed, TENSORPARA &tpara, int ploidity, float *&dp1, float *&dp2, double const *&rp, float const *&nhp, vector<vector<PARAUNIT> >::const_iterator &ipm, int const *&aszp, double *&lpp, double *&lmm, double *&lss, int SSZ, int MK, vector<vector<int> > const &ss, vector<vector<int> > const &mlist, double *Vh, double *Vh2)=0;
	virtual int _tracebackOne(int id, vector<int> const &breaks, int bst, int bed, TENSORPARA &tpara, int ploidity, int *&pp, int *&pp2, bool *&rr, bool *&rr2, double const *&rp, float const *&nhp, vector<vector<PARAUNIT> >::const_iterator &ipm, double *&lpp, double *&lmm, double *&lss, int SSZ, int MK, vector<vector<int> > const &ss, vector<vector<int> > const &mlist, double *Vh, double *Vh2)=0;
	virtual double _getLp(int id, vector<int> const &breaks, int bid, bool flip, int ploidity, TENSORPARA const &tpara, float const *dp, float const *dp2, int const *pp, int const *pp2, vector<vector<PARAUNIT> >::const_iterator ipm, float const *nhp, int const *aszp)=0;
	virtual double _imputeOne(int id, vector<int> const &breaks, int bst, int bed, int flip, int ploidity, TENSORPARA const &tpara, float *dp, float *dp2, int const *pp, int const *pp2, vector<vector<PARAUNIT> >::const_iterator ipm, float const *nhp, int const *aszp)=0;
	virtual void _updatePrior(vector<bool> const &init, bool oneround = false)=0;
	virtual void outputResult(char const *fname)=0;

private:
	int _shrinkClasses(int indN, int L, float *z);
	void _removeInd(MYDATA const &mydata, TENSORPARA &tpara, int id, float wb = 1., int st = -1, int ed = -1, bool clean = false);
        void _addInd(MYDATA const &mydata, TENSORPARA &tpara, int id, float wb = 1., int st = -1, int ed = -1);
        int _updateRecomb(TENSORPARA &tpara, float trueindN, int L);
        double _updateA(int maxK, double const *Q, double A);
	void _collapseGroups(int indN, int L, int ploidity, TENSORPARA &tpara, vector<bool> const &init);
	void _collapseGroups_new(MYDATA const &mydata, TENSORPARA &tpara);
	void _greedySwitch_new(int indN, int L, int ploidity, TENSORPARA &tpara, vector<bool> const &init);
	void _greedyMatch(float *trans, int maxK, int *perm);
	virtual double _imputeData(MYDATA const &mydata, TENSORPARA const &tpara, int id, vector<int> const &dbreaks, bool init);
	void _getDoubleBreaks(vector<int> const &breaks, int L, vector<int> &dbreaks);
	void _checkCounts(TENSORPARA const &tpara);
};

