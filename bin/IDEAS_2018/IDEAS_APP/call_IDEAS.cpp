/*#include "datastructure.h"
#include "genomicTensor.h"
#include "ideas.h"

unsigned int rseed = 0;
bool hpass = false; 
bool splitmerge = true;
bool SA = false;
bool gzip = true;
int bST = -1, bED = -1;
vector<string> imputelist;

///////////////////////////////////////////////////
void IDEAS::call_IDEAS(int argc, char* argv[])
{
	int i, j;
	int minCut = 5, admixN = 30, burnin = 50, mcmc = 50, maxHapK = 0, maxGG = 0, maxPos = 0, thread = 1;
	bool sqc = false, lik = false, samplemaximum = true, outputproportion = false;
	double log2 = -1;
	bool add2 = false, ind = false, indind = false, nb = false, norm = false;
	double error = 0.01, A = 0, recr = 100., heteroVh = 1., minerr = 0.5, maxerr = 100000000.;
    char const *output = argv[1], *fixPop = NULL, *input = NULL, *fbed = NULL, *fcov = NULL, *fcovbed = NULL, *fparam = NULL, *fmixpara = NULL, *parafile = NULL, *statefile = NULL, *clusterfile = NULL, *para0file = NULL, *state0file = NULL, *cluster0file = NULL, *profile0file = NULL;
	int startK=0, fixC=0;
	vector<string> fS, fA;
	vector<vector<int> > remove;

	i = 1;
    input = argv[1];
    if(argc > 2 && argv[2][0] != '-')
    {	fbed = argv[2];
		i = 2;
    }
		
    for(i = i + 1; i < argc; i++)
    {	if(strcmp(argv[i], "-rseed") == 0)
        {	rseed = atoi(argv[i + 1]);
			i++;
		}
        else if(strcmp(argv[i], "-sa") == 0)
		{	SA = true;
		}
        else if(strcmp(argv[i], "-o") == 0)
        {	output = argv[i + 1];
			i++;
		}
        else if(strcmp(argv[i], "-P") == 0)
        {	maxPos = atoi(argv[i + 1]);
			i++;
		}
        else if(strcmp(argv[i], "-K") == 0)
        {	maxHapK = atoi(argv[i + 1]);
			i++;
		}
        else if(strcmp(argv[i], "-G") == 0)
        {	maxGG = atoi(argv[i + 1]);
			startK = max(startK, maxGG);	
			i++;
		}
        else if(strcmp(argv[i], "-c") == 0)
        {	minCut = atoi(argv[i + 1]);
			i++;
		}
        else if(strcmp(argv[i], "-err") == 0)
        {	error = atof(argv[i + 1]);
			i++;
		}
        else if(strcmp(argv[i], "-inv") == 0)
        {	bST = atoi(argv[i + 1]);
            bED = atoi(argv[i + 2]);
			i += 2;
		}
        else if(strcmp(argv[i], "-nogz") == 0)
		{	gzip = false;
		}
        else if(strcmp(argv[i], "-impute") == 0)
        {	int l = (int)strlen(argv[i + 1]);
			int ii, jj;
			jj = 0;
			for(ii = 0; ii < l; ii++) 
            {	if(argv[i + 1][ii] == ',')
                {	argv[i + 1][ii] = 0;
                    imputelist.push_back(&argv[i + 1][jj]);
					jj = ii + 1;
                    argv[i + 1][ii] = ',';
				}
			}
            if(jj < ii) imputelist.push_back(&argv[i + 1][jj]);
			//for(ii = 0; ii < (int)imputelist.size(); ii++)
			//	printf("%s\n", imputelist[ii].c_str());
			i++;
		}
        else if(strcmp(argv[i], "-fixPop") == 0)
        {	fixPop = argv[i + 1];
			i++;
		}
        else if(strcmp(argv[i], "-fixHap") == 0)
		{	if((int)fS.size() != (int)fA.size()) printf("-fixHap option cannot be used together with -fixAllele option.\n"),exit(0);
            fS.push_back(argv[i + 1]);
            fA.push_back(argv[i + 2]);
			i += 2;
		}
        else if(strcmp(argv[i], "-fixAllele") == 0)
		{	if((int)fS.size() > 0) printf("-fixAllele option cannot be used together with -fixHap option.\n"),exit(0);
            fA.push_back(argv[i + 1]);
			i++;
		}
        else if(strcmp(argv[i], "-fixC") == 0)
        {	fixC = atoi(argv[i + 1]);
			i++;
		}
        else if(strcmp(argv[i], "-admixn") == 0)
        {	admixN = atoi(argv[i + 1]);
			i++;
		}
        else if(strcmp(argv[i], "-sqc") == 0)
			sqc = true;
        else if(strcmp(argv[i], "-sample") == 0)
        {	mcmc = atoi(argv[i + 2]);
            burnin = atoi(argv[i + 1]);
			i += 2;
		}
        else if(strcmp(argv[i], "-lik") == 0)
			lik = true;
        else if(strcmp(argv[i], "-A") == 0)
        {	A = atof(argv[i + 1]);
			i ++;
		}
        else if(strcmp(argv[i], "-z") == 0)
        {	fcov = argv[i + 1];
			i++;
            if(i + 1 < argc && argv[i + 1][0] != '-')
            {	fcovbed = argv[i + 1];
				i++;
			}
		}
        else if(strcmp(argv[i], "-max") == 0)
		{	samplemaximum = true;
		}
        else if(strcmp(argv[i], "-outputprop") == 0)
		{	outputproportion = true;
		}
        else if(strcmp(argv[i], "-rec") == 0)
        {	recr = atof(argv[i + 1]) * 1000.;
			i++;
		}
        else if(strcmp(argv[i], "-log2") == 0)
		{	log2 = 1.;	
            if(argc > i + 1 && argv[i + 1][0] >= 48 && argv[i + 1][0] < 58)
            {	log2 = atof(argv[i + 1]);
				i++;
				//printf("log2 constant = %f\n", log2);
			}
		}
        else if(strcmp(argv[i], "-minerr") == 0)
        {	minerr = atof(argv[i + 1]);
			i++;
		}
        else if(strcmp(argv[i], "-maxerr") == 0)
        {	maxerr = atof(argv[i + 1]);
			i++;
		}
        else if(strcmp(argv[i], "-add2") == 0) add2 = true;
        else if(strcmp(argv[i], "-param") == 0)
        {	fparam = argv[i + 1];
			i++;
		}
        else if(strcmp(argv[i], "-mixpara") == 0)
        {	fmixpara = argv[i + 1];
			i++;
		}
        else if(strcmp(argv[i], "-startpara") == 0) //start with specified gauss parameters
        {	parafile = argv[i + 1];
			i++;
		}
        else if(strcmp(argv[i], "-prevrun") == 0) //restore segments from previous run on the same data and then continue to run
        {	statefile = argv[i + 1];
            clusterfile = argv[i + 2];
			i += 2;
		}
        else if(strcmp(argv[i], "-otherpara") == 0) //use other gauss parameters as priors
        {	para0file = argv[i + 1];
			i++;
            if(argc > i + 1 && argv[i + 1][0] != '-')
            {	profile0file = argv[i + 1];
				i++;
			}
		}
        else if(strcmp(argv[i], "-otherstate") == 0) //use other segments as priors
        {	state0file = argv[i + 1];
            cluster0file = argv[i + 2];
			i += 2;
		}
        else if(strcmp(argv[i], "-norm") == 0)
		{	norm = true;
		}
        else if(strcmp(argv[i], "-C") == 0)
        {	startK=atoi(argv[i+1]);
			i++;
		}
        else if(strcmp(argv[i], "-ind") == 0)
		{	ind = true;
		}
        else if(strcmp(argv[i], "-indind") == 0)
		{	indind = true;
		}
        else if(strcmp(argv[i], "-h") == 0)
        {	heteroVh = atof(argv[i + 1]);
			i++;
		}
        else if(strcmp(argv[i], "-hp") == 0)
		{	hpass = true;
		}
        else if(strcmp(argv[i], "-nosm") == 0)
		{	splitmerge = false;
		}
        else if(strcmp(argv[i], "-remove") == 0)
		{	i++;
            int ind = atoi(argv[i]), ll = (int)strlen(argv[i]) - 1;
			vector<int> row;
            for(j = 1; j < ll; j++) if(argv[i][j] == ':') break;
			do {
                row.push_back(atoi(&argv[i][j + 1]));
                for(j = j + 1; j < ll; j++) if(argv[i][j] == ',') break;
			} while(j < ll);
			if((int)remove.size() <= ind) remove.resize(ind + 1, vector<int>());
			remove[ind] = row;
		}
        else if(strcmp(argv[i], "-nb") == 0)
		{	nb = true;
		}
        else if(strcmp(argv[i], "-thread") == 0)
        {	thread = atoi(argv[i + 1]);
			i++;
		}
    }

	time_t tst, ted;
	time(&tst);
	genomicTensor gtm(rseed);
	gtm.outputproportion = outputproportion;
	gtm.log2 = log2;
	gtm.norm = norm;
	if(fixC > 0) startK = fixC;
	gtm.run(input, fbed, fcov, fcovbed, burnin, mcmc, thread, maxHapK, maxGG, maxPos, A, recr, heteroVh, samplemaximum, output, sqc, fparam, add2, minerr, maxerr, ind, indind, startK, fixC, remove, nb, fmixpara, parafile, statefile, clusterfile, para0file, state0file, cluster0file, profile0file); 
	
	time(&ted);
	//printf("total time = %dsec\n", (int)(ted - tst));

    //return 0;
}

*/
