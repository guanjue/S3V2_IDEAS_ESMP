#include "datastructure.h"
#include "tensorHMMbase.h"
#include "MixGauss.h"

class genomicTensor : public tensorHMMbase
{
public:

        genomicTensor(unsigned int rs=0);
        ~genomicTensor();

	MYDATA mydata;
	TENSORPARA tpara;
	int posSZ, maxymsz, *ymsz, maxxmsz, *xmsz, **ymap, **xmap;
	float **dataYP, **dataXP;
	int *basePop;
	bool *baseR;
	double *emitMatrix;
	vector<vector<double> > emit0;
	int emitK, maxG, maxPosSZ;
	double *tmpSpace, *TLP, *tmpSpace2, *tpointer, *dpointer;/////////////////////
	int *gClass, gClassSZ, posClassSZ;
	double overallmean, overallsd, overallmeanx, overallsdx, log2;
	bool independentInd, fixstateN, norm;

	MixGauss mylik;
	int nthread;

public:
	void run(char const *finput, char const *fbed, char const *fcov, char const *fcovbed, int burnin, int mcmc, int thread, int maxHapK, int maxGG, int maxPos, double AA, double recr, double heteroVh, bool samplemaximum, char const *foutput, bool sqc, char const *fparam, bool add2, double merr, double mxerr, bool ind, bool indind, int startK, int fixC, vector<vector<int> > const &removelist, bool nb, char const *fmixpara, char const *parafile, char const *statefile, char const *clusterfile, char const *para0file, char const *state0file, char const *cluster0file, char const *profile0file);
	void parseParameter(char const *parafile);
	void removeData(vector<vector<int> > const &list);

protected:
        void _getAsz(MYDATA const &mydata);
        void _updatePara_remove(PARAUNIT &para, float data, float b, int id, int pos);
        void _updatePara_add(PARAUNIT &para, float data, float b, int id, int pos);
        void  _sumforwardOne(int id, vector<int> const &breaks, int bst, int bed, TENSORPARA &tpara, int ploidity, float *&dp1, float *&dp2, double const *&rp, float const *&nhp, vector<vector<PARAUNIT> >::const_iterator &ipm, int const *&aszp, double *&lpp, double *&lmm, double *&lss, int SSZ, int MK, vector<vector<int> > const &ss, vector<vector<int> > const &mlist, double *Vh, double *Vh2);
        int _tracebackOne(int id, vector<int> const &breaks, int bst, int bed, TENSORPARA &tpara, int ploidity, int *&pp, int *&pp2, bool *&rr, bool *&rr2, double const *&rp, float const *&nhp, vector<vector<PARAUNIT> >::const_iterator &ipm, double *&lpp, double *&lmm, double *&lss, int SSZ, int MK, vector<vector<int> > const &ss, vector<vector<int> > const &mlist, double *Vh, double *Vh2);
	double _getLp(int id, vector<int> const &breaks, int bid, bool flip, int ploidity, TENSORPARA const &tpara, float const *dp, float const *dp2, int const *pp, int const *pp2, vector<vector<PARAUNIT> >::const_iterator ipm, float const *nhp, int const *aszp);
	double _imputeOne(int id, vector<int> const &breaks, int bst, int bed, int flip, int ploidity, TENSORPARA const &tpara, float *dp, float *dp2, int const *pp, int const *pp2, vector<vector<PARAUNIT> >::const_iterator ipm, float const *nhp, int const *aszp);
	double _updateLapseData(int id, double *P, int posst, int posed, float *nlapse);

	void readInput(MYDATA &mydata, char const *input, char const *fbed, char const *fcov, char const *fcovbed, bool add2 = false);
	void _loadData(float **&datamatrix, char const *input, vector<string> &sid, vector<string> &fid, vector<SNPINFO> &snpinfo, int &indN, int &L, int &totalN, int &maxmsz, int *&msz, int **&map, int *&indIndex);

	void outputResult(char const *fname);
	void _updatePrior(vector<bool> const &init, bool oneround = false);
	void _rearrangeState(int add, int maxG);
	
	double _lgammaM(int q, double x);
	void _getProp(TENSORPARA const &tpara, int id, int pos, int state, int asz, vector<vector<PARAUNIT> >::const_iterator ipm, float const *nhp, int const *op, double *prop);
	void _getProp_soft(TENSORPARA const &tpara, int id, int pos, int asz, vector<vector<PARAUNIT> >::const_iterator ipm, float const *nhp, int const *op, double *prop);
	void _updateBase();
	void _getMLE_P(int *tpop, int N, int o, double *E, int oK, int K, double *PP, int *tttg, double *priorpp, int priorstep);
	void _getMLE_HyperP(double *logp, double *ak, int sz);

	double _solveDigamma(double x0, double y, double precision = 0.01, int maxitern = 100);
	void _getMLE_HyperP1(double *logp, double *ak, double &M, int sz, double N, double precision = 0.01, int maxitern = 100);

	void _splitmergeCluster(int type);
	void _refreshParameter();

	void groupData(float *data, int index, int msz, double step, vector<int> const &inputlist, vector<vector<double> > &values, vector<vector<int> > &outputlist);

	void _getLpVar(float *ypp, float *xpp, int nst, int ned, int *tymsz, int *txmsz, vector<vector<double> > const &props, double priorW, int asz, vector<double> &rt);
	void _testVariability(vector<int> const &response, double *vprob);
	void _groupComposition(vector<vector<double> > &ggg, vector<double> &nnn, vector<vector<double> > &prop, vector<double> &newnnn, double A, double B);
	void _pairEnrich(vector<vector<double> > const &prop, vector<double> const &nnn, vector<vector<double> > &score);

	void BronKerbosch1(vector<vector<double> > const &score, vector<int> P, vector<int> R, vector<int> X, vector<vector<int> > &clique, double scorecut);
	void _groupComposition1(vector<vector<double> > const &ggg, vector<double> const &nnn, vector<vector<double> > &prop, vector<double> &newnnn);
	void _propErrorCorrect(vector<vector<double> > &ggg, vector<double> &nnn, vector<int> const &response, vector<vector<int> > &glabel, double a);

	void loadPreviousRun(char const *statefile, char const *clusterfile, double posRate);
	void loadOtherResult(char const *state0file, char const *cluster0file, double posRate);
	void _loadData_list(float **&datamatrix, char const *listfile, char const *bedfile, vector<string> &sid, vector<string> &fid, vector<SNPINFO> &snpinfo, int &indN, int &L, int &totalN, int &maxmsz, int *&msz, int **&map, int *&indIndex);
	void dataNorm();

	double _imputeData(MYDATA const &mydata, TENSORPARA const &tpara, int id, vector<int> const &dbreaks, bool init);//reimplemented to allow parallelization
	void _imputeMissingData(char const *fname, int gid);
};

