#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <vector>
#include <algorithm>
#include <iterator>
#include <sys/stat.h>

using namespace std;

typedef struct MySortType {
        int st, ed, did, sid, g[10];
	string chr, disease, rsid, orient;
	double r2;
} MySortType;

bool operator<(const MySortType &a, const MySortType &b)
{       return(a.chr < b.chr || (a.chr==b.chr && a.st < b.st) || (a.chr==b.chr && a.st==b.st && a.ed < b.ed));
}

inline bool is_empty_c(ifstream & pFile)
{	return pFile.peek() == ifstream::traits_type::eof();
}
	
time_t rseed;
bool renewmissingfile = true;

char const *tmpfolder = "tmp/";

/*---------------------------------------------------------------*/
void readFiles(char const *fname, vector<vector<string> > &flist);
void prepareBED(vector<string> const &chr, vector<int> const &gsz, int wsz, vector<string> const &exbed, char const *fout);
void readBED(char const *fname, vector<MySortType> &mBED, int center=-1, bool skip1 = false, bool split = false);
void mergeBED(vector<MySortType> &X, vector<MySortType> const &Y, int op);
void prepareMatrix(vector<vector<string> > const &flist, vector<string> &bedfile, char const *gszfile, vector<string> const &chr, int bedcol, char const *fout, bool norm = false, char const *sublist = NULL);

bool DirectoryExists(const char* dirName);
/*---------------------------------------------------------------*/
int main(int argc, char* argv[])
{
	time(&rseed);
	int i, j, k, l; 
	char const *gszfile = NULL;//"extra/hg19.genome";
	char const *bedfile = NULL;
	char const *invfile = NULL;
	char const *sublist = NULL;
	char const *fout = "try.data";
	int wsz = 200;
	int bedcol = 5;//4: sum over covered; 5: mean over all; 6: mean over covered; 7: min; 8: max
	bool norm = false, bychr = false;
	vector<string> chrlist;
	vector<string> exbed;
	for(i = 2; i < argc; i++)
	{	if(strcmp(argv[i], "-gsz") == 0)
		{	gszfile = argv[i + 1];
			i++;
		}
		else if(strcmp(argv[i], "-bed") == 0)
		{	bedfile = argv[i + 1];
			i++;
		}
		else if(strcmp(argv[i], "-inv") == 0)
		{	invfile = argv[i + 1];
			i++;
		}
		else if(strcmp(argv[i], "-wsz") == 0)
		{	wsz = atoi(argv[i + 1]);
			i++;
		}
		else if(strcmp(argv[i], "-exclude") == 0)
		{	exbed.push_back(argv[i + 1]);
			i++;
		}
		else if(strcmp(argv[i], "-bychr") == 0) //output different chromosomes in different files
		{	bychr = true;
		}
		else if(strcmp(argv[i], "-o") == 0)
		{	fout = argv[i + 1];
			i++;
		}
		else if(strcmp(argv[i], "-c") == 0)
		{	bedcol = atoi(argv[i + 1]);
			i++;
		}
		else if(strcmp(argv[i], "-norm") == 0)
		{	norm = true;
		}
		else if(strcmp(argv[i], "-sub") == 0)
		{	sublist = argv[i + 1];
			i++;
		}
		else if(strcmp(argv[i], "-a") == 0)
		{	renewmissingfile = false;
		}
		else if(strcmp(argv[i], "-dir") == 0)
		{	tmpfolder = argv[i + 1];
			i++;
		}
		else if(argv[i][0] == '-' && argv[i][1] == 'c' && argv[i][2] == 'h' && argv[i][3] == 'r')
		{	bool flag = false;
			l = (int)strlen(argv[i]);
			for(j = 4; j < l; j++)
			{	if(argv[i][j] == ',') flag = false;
				else if(argv[i][j] == '-') flag = true;
				else
				{	k = atoi(&argv[i][j]);
					int oj = j;
					for(j = j + 1; j < l; j++) if(argv[i][j] == ',' || argv[i][j] == '-') break;
					if(flag)
					{	int ii, jj = atoi(chrlist[(int)chrlist.size() - 1].c_str());
						for(int ii = jj + 1; ii <= k; ii++)
						{	ostringstream convert;
							convert<<ii;
							chrlist.push_back(convert.str());
						}
					}
					else 
					{	char c;
						if(j < l)
						{	c = argv[i][j];
							argv[i][j] = 0;
						}
						chrlist.push_back(&argv[i][oj]);
						if(j < l) argv[i][j] = c;
					}
					j--;
				}
			}
		}
	}
	
	//printf("data list file: %s\n", argv[1]);
	//if(gszfile != NULL) printf("genome size file: %s\n", gszfile);
	//for(i = 0; i < (int)exbed.size(); i++)
	//	printf("exclude regions in %s\n", exbed[i].c_str());

	if(bedfile == NULL && gszfile == NULL)
	{	printf("Either bed file or genome size file must be provided\n");exit(0);
	}

	if(!DirectoryExists(tmpfolder))
	{	char str[200];
		sprintf(str, "mkdir %s", tmpfolder);
		system(str);
	}

	vector<vector<string> > flist;
	readFiles(argv[1], flist);
	//for(i = 0; i < (int)flist.size(); i++)
	//{	for(j = 0; j < (int)flist[i].size(); j++)
	//	{	printf("<%s> ", flist[i][j].c_str());
	//	}
	//	printf("\n");
	//}

	string tmp;
	vector<string> chr;
	vector<int> gsz;
	if(gszfile != NULL)
	{	ifstream f(gszfile);
		if(f.fail()) printf("Cannot open %s\n", gszfile), exit(0);
		for(i = 0; getline(f, tmp); i++) 
		{	l = (int)tmp.size();
			for(j = 0; j < l; j++) if(tmp[j] == ' ' || tmp[j] == '\t') break;
			chr.push_back(tmp.substr(0,j));
			gsz.push_back(atoi(tmp.substr(j + 1).c_str()));
		}
		f.close();
	}

	if((int)chrlist.size() > 0 && bedfile == NULL)
	{	vector<bool> sel((int)chr.size(), false);
		for(i = 0; i < (int)chrlist.size(); i++)
		{	string tchr = "chr";
			tchr.append(chrlist[i]);
			for(j = 0; j < (int)chr.size(); j++)
				if(chr[j] == tchr) sel[j] = true;
		}
		vector<string> nchr;
		vector<int> ngsz;
		for(i = 0; i < (int)sel.size(); i++)
			if(sel[i])
			{	nchr.push_back(chr[i]);
				if((int)gsz.size() > 0) ngsz.push_back(gsz[i]);
			}
		chr = nchr;
		gsz = ngsz;
	}

	vector<string> mybedfile;
	if(bedfile == NULL)
	{	char tmpbed[1000];
		if(bychr)
		{	for(i = 0; i < (int)chr.size(); i++)
			{	sprintf(tmpbed, "%s%s.%d.bed", tmpfolder, chr[i].c_str(), rseed);
				prepareBED(vector<string>(1,chr[i]), vector<int>(1,gsz[i]), wsz, exbed, tmpbed);
				mybedfile.push_back(tmpbed);
			}
			prepareMatrix(flist, mybedfile, gszfile, chr, bedcol, fout, norm, sublist);
		}
		else
		{	sprintf(tmpbed, "%s%s.%d.bed", tmpfolder, gszfile, rseed);
			prepareBED(chr, gsz, wsz, exbed, tmpbed);
			mybedfile.push_back(tmpbed);
			prepareMatrix(flist, mybedfile, gszfile, vector<string>(), bedcol, fout, norm, sublist);
		}
	}
	else 
	{	mybedfile.push_back(bedfile);
		prepareMatrix(flist, mybedfile, gszfile, vector<string>(), bedcol, fout, norm, sublist);
	}

	return 0;
}

void readFiles(char const *fname, vector<vector<string> > &flist)
{	
	ifstream f(fname);
	int i, j, k, l, oj;
	string tmp;
	if(f.fail()) printf("Cannot open %s\n", fname), exit(0);
//	getline(f, tmp);

	flist.clear();
	for(i = 0; getline(f, tmp); i++)
	{	if(tmp[0] == '#') continue;
		l = tmp.size();
		k = 0;
		vector<string> row;
		oj = 0;
		for(j = 0; j < l; j++)
		{	if(j == l - 1 || tmp[j] == ' ' || tmp[j] == '\t')
			{	row.push_back(tmp.substr(oj, j+(int)(j==l-1)-oj));
				k++;
				if(k >= 3) break;
				oj = j + 1;
			}
		}
		flist.push_back(row);
	}	
	f.close();
}

void prepareBED(vector<string> const &chr, vector<int> const &gsz, int wsz, vector<string> const &exbed, char const *fout)
{	int i, j, k;
	FILE *f = fopen(fout, "w");
	if(f == NULL) printf("Fail to create BED file %s\n", fout), exit(0);
	
	vector<MySortType> mybed;
	for(i = 0; i < (int)chr.size(); i++)
	{	for(j = 0; j < gsz[i]; j += wsz)
		{	MySortType rt;
			rt.chr = chr[i];
			rt.st = j;
			rt.ed = j + wsz;
			mybed.push_back(rt);
		}
	}
	sort(mybed.begin(), mybed.end());
	
	if((int)exbed.size() > 0)
	{	vector<MySortType> rmbed;
		for(i = 0; i < (int)exbed.size(); i++)
		{	vector<MySortType> tbed;
			readBED(exbed[i].c_str(), tbed);
			if(i == 0) rmbed = tbed;
			else mergeBED(rmbed, tbed, 0);
		}
		mergeBED(mybed, rmbed, 1);
	}
	
	for(i = 0; i < (int)mybed.size(); i++)
		fprintf(f, "%s\t%d\t%d\tR%d\n", mybed[i].chr.c_str(), mybed[i].st, mybed[i].ed, i);

	fclose(f);
}

void readBED(char const *fname, vector<MySortType> &mBED, int center, bool skip1, bool split)
{	
	mBED.clear();
	
	FILE *f = fopen(fname, "r");
	char tmp[10000];
	while(fgets(tmp, 10000, f) != NULL)
	{	MySortType rt;
		int i, l = (int)strlen(tmp) - 1, j;
		if(tmp[0] != 'c' && (tmp[0] < 48 || tmp[0] >= 58)) { continue; }
		i = 0;
		for(j = i; j < l; j++) { if(tmp[j] == ' ' || tmp[j] == '\t') break; };
		if(skip1) { i = j; i++; for(j = i; j < l; j++) { if(tmp[j] == ' ' || tmp[j] == '\t') break; } }
		if(tmp[i] != 'c') continue;
		tmp[j] = 0;
		rt.chr = &tmp[i];
		tmp[j] = ' ';

		for(i = i + 1; i < l; i++) if(tmp[i] == ' ' || tmp[i] == '\t' || tmp[i] == ':') break;
		rt.st = atoi(&tmp[i + 1]);
		for(i = i + 1; i < l; i++) if(tmp[i] == ' ' || tmp[i] == '\t' || tmp[i] == '-') break;
		rt.ed = atoi(&tmp[i + 1]);
		for(i = i + 1; i < l; i++) if(tmp[i] == ' ' || tmp[i] == '\t') break;
		for(i = i + 1; i < l; i++) if(tmp[i] == ' ' || tmp[i] == '\t') break;
		for(j = i + 1; j < l; j++) if(tmp[j] == '\t' || tmp[j] == ' ') break;
		tmp[j] = 0;
		rt.rsid = &tmp[i + 1];
		tmp[j] = ' ';
		i = j;
		for(j = i + 1; j < l; j++) if(tmp[j] == '\t' || tmp[j] == ' ') break;
		tmp[j] = 0;
		rt.orient = &tmp[i + 1];
		
		if(center > 0)
		{	rt.st = (int)((double)rt.st / (double)center + 0.5) * center;
			rt.ed =	(int)((double)rt.ed / (double)center + 0.5) * center;
		}
		if(rt.ed > rt.st) 
		{	if(!split) mBED.push_back(rt);
			else
			{	for(j = 0; j < (rt.ed - rt.st) / center; j++)
				{	MySortType trt = rt;
					trt.st = rt.st + j * center; trt.ed = trt.st + center;
					mBED.push_back(trt);
				}
			}
		}
		if(rt.ed<rt.st) printf("#!!!%d\t%s %d %d | %s\n", (int)mBED.size(), rt.chr.c_str(), rt.st, rt.ed, tmp);
	}
	fclose(f);

	sort(mBED.begin(), mBED.end());
}

//assume the positions have been sorted
//op= 0: add, 1: minus, 2: intersection
void mergeBED(vector<MySortType> &X, vector<MySortType> const &Y, int op)
{
	int i, j, k, istx = 0, iedx, isty = 0, iedy;
	vector<string> uniquechr;
	k = -1;
	for(i = 0; i < (int)X.size(); i++)
	{	if(k >= 0 && X[i].chr == uniquechr[k]) continue;
		for(j = 0; j < (int)uniquechr.size(); j++)
			if(uniquechr[j] == X[i].chr) break;
		k = j;
		if(j >= (int)uniquechr.size()) uniquechr.push_back(X[i].chr);
	}
	k = -1;
	for(i = 0; i < (int)Y.size(); i++)
	{	if(k >= 0 && Y[i].chr == uniquechr[k]) continue;
		for(j = 0; j < (int)uniquechr.size(); j++)
			if(uniquechr[j] == Y[i].chr) break;
		k = j;
		if(j >= (int)uniquechr.size()) uniquechr.push_back(Y[i].chr);
	}

	vector<MySortType> Z;
	for(i = 0; i < (int)uniquechr.size(); i++)
	{	vector<int> pos;
		for(j = istx; j < (int)X.size(); j++)
			if(X[j].chr == uniquechr[i]) 
			{	pos.push_back(X[j].st);
				pos.push_back(X[j].ed);
			}
			else break;
		iedx = j;
		for(j = isty; j < (int)Y.size(); j++)
			if(Y[j].chr == uniquechr[i])
			{	pos.push_back(Y[j].st);
				pos.push_back(Y[j].ed);
			}
			else break;
		iedy = j;
		sort(pos.begin(), pos.end());
		k = 1;
		for(j = 1; j < (int)pos.size(); j++)
			if(pos[j] != pos[j - 1]) pos[k++] = pos[j];
		pos.resize(k);

		vector<int> selX((int)pos.size(), -1), selY = selX;
		int a, b = 0;
		for(j = istx; j < iedx; j++)
		{	for(b; b < (int)pos.size(); b++)
			{	if(pos[b] == X[j].st) break;
			}
			for(a = b; a < (int)pos.size(); a++)
			{	if(pos[a] < X[j].ed) selX[a] = j;
				else break;
			}
		}
		b = 0;
		for(j = isty; j < iedy; j++)
		{	for(b; b < (int)pos.size(); b++)
			{	if(pos[b] == Y[j].st) break;
			}
			for(a = b; a < (int)pos.size(); a++)
			{	if(pos[a] < Y[j].ed) selY[a] = j;
				else break;
			}
		}
		
		for(j = 0; j < (int)pos.size(); j++)
		{	for(k = j + 1; k < (int)pos.size(); k++)
				if(selX[k] != selX[k - 1] || selY[k] != selY[k - 1]) break;
			MySortType nx;
			if(selX[j] >= 0) nx = X[selX[j]];
			else if(selY[j] >= 0) nx = Y[selY[j]];
			else {	j = k - 1; continue;	}

			if(k >= (int)pos.size()) k--;
			nx.st = pos[j];
			nx.ed = pos[k];
			if(op == 0 || (op == 1 && selY[j] < 0) || (op == 2 && selX[j] >= 0 && selY[j] >= 0))
			{	Z.push_back(nx);
			}
			
			j = k - (int)(k < (int)pos.size() - 1);
		}

		istx = iedx;
		isty = iedy;
	}
	
	X = Z;
}

void prepareMatrix(vector<vector<string> > const &flist, vector<string> &bedfile, char const *gszfile, vector<string> const &chr, int bedcol, char const *fout, bool norm, char const *sublist)
{
	int i, j, k;
	char tmpstr[100000];
	vector<int> subindex;
	if(sublist != NULL)
	{	FILE *f = fopen(sublist, "r");
		char tmp[1000];
		while(fgets(tmp, 1000, f) != NULL)
			subindex.push_back(atoi(&tmp[0]) - 1);
		fclose(f);
	}
	bool missingflag = true;
	FILE *f;
	if(renewmissingfile)
	{	char str[200];
		sprintf(str, "%smissingdata.txt", tmpfolder);
		f = fopen(str, "w");
		fclose(f);
	}

	vector<string> obedfile = bedfile;
	vector<int> bedsz;
	for(j = 0; j < (int)bedfile.size(); j++)
	{		FILE *tf = fopen(bedfile[j].c_str(), "r");
			char ttmp[1000];
			sprintf(ttmp, "%stmpbed%d.%d", tmpfolder, j, rseed);
			bedfile[j] = ttmp;
			FILE *nf = fopen(ttmp, "w");
			int iii;
			for(iii = 0; fgets(ttmp, 1000, tf) != NULL; iii++)
			{	int k, l = (int)strlen(ttmp) - 1;
				if(l < 4) continue;
				for(k = 0; k < l; k++) if(ttmp[k]==' ' || ttmp[k] == '\t') break;
				for(k = k + 1; k < l; k++) if(ttmp[k]==' ' || ttmp[k] == '\t') break;
				for(k = k + 1; k < l; k++) if(ttmp[k]==' ' || ttmp[k] == '\t') break;
				ttmp[k] = 0;
				fprintf(nf, "%s\t%d\n", ttmp, iii);
			}
			bedsz.push_back(iii);
			fclose(tf);
			fclose(nf);
	}

	vector<string> bednames;
	for(i = 0; i < (int)flist.size(); i++)
	{	printf("Fetching %s ... ", flist[i][2].c_str());fflush(stdout);
		
		bool url = false;
		sprintf(tmpstr, "%s", flist[i][2].c_str());
		if(tmpstr[0]=='h' && tmpstr[1]=='t' && tmpstr[2]=='t' && tmpstr[3]=='p') url=true;
		if(tmpstr[0]=='f' && tmpstr[1]=='t' && tmpstr[2]=='p') url=true;
		
		size_t found = flist[i][2].find_last_of("/");
		string name = flist[i][2].substr(found+1), name0 = name;
		int l = (int)name.size();
		while(l > 0 && (name[l - 1] == ' ' || name[l - 1] == '\t')) l--;
		name[l] = 0;
		if(l > 3 && name[l - 3] == '.' && name[l - 2] == 'g' && name[l - 1] == 'z') 
		{	name0 = name = name.substr(0, l - 3); 
			ostringstream convert1;
			convert1 << name << rseed;
			name = convert1.str();
			if(url) sprintf(tmpstr, "wget -q --no-check-certificate O - %s | gunzip > %s%s", flist[i][2].c_str(), tmpfolder, name.c_str());
			else sprintf(tmpstr, "gunzip -c %s > %s%s", flist[i][2].c_str(), tmpfolder, name.c_str());
		}
		else 
		{	ostringstream convert1;
			convert1 << name << rseed;
			name = convert1.str();
			if(url) sprintf(tmpstr, "wget -q --no-check-certificate -O %s%s %s", tmpfolder, name.c_str(), flist[i][2].c_str());
			else sprintf(tmpstr, "cp %s %s%s", flist[i][2].c_str(), tmpfolder, name.c_str());
		}
		system(tmpstr);

		char ttt[1000];
		sprintf(ttt, "%s%s", tmpfolder, name.c_str());
		ifstream ifile(ttt);
        if(is_empty_c(ifile))
		{	char str[200];
			sprintf(str, "%smissingdata.txt", tmpfolder);
			f = fopen(str, "a");
			fprintf(f, "%s\n", flist[i][2].c_str());
			fclose(f);
			if(ifile) 
			{	sprintf(tmpstr, "rm %s%s", tmpfolder, name.c_str());
				system(tmpstr);
				ifile.close();
				printf("Failed\n");
			}
			missingflag = false;
			continue;
		}

		int l0 = (int)name0.size();
		if(l0 > 4 && name0[l0 - 4] == '.' && name0[l0 - 3] == 'b' && name0[l0 - 2] == 'a' && name0[l0 - 1] == 'm')
		{	if(gszfile == NULL) printf("BAM file cannot be processed because genome size file is not found.\n"),exit(0);
			
			sprintf(tmpstr, "bedtools bamtobed -i %s%s > %s%s.tb", tmpfolder, name.c_str(), tmpfolder, name.c_str());
			//printf("%s\n", tmpstr);fflush(stdout);
			system(tmpstr);
			sprintf(tmpstr, "rm %s%s", tmpfolder, name.c_str());
			system(tmpstr);

			sprintf(tmpstr, "genomeCoverageBed -i %s%s.tb -bg -g %s > %s%s.tc", tmpfolder, name.c_str(), gszfile, tmpfolder, name.c_str());
			//printf("%s\n", tmpstr);fflush(stdout);
			system(tmpstr);
			sprintf(tmpstr, "rm %s%s.tb", tmpfolder, name.c_str());
			system(tmpstr);

			sprintf(tmpstr, "bedSort %s%s.tc %s%s.tb", tmpfolder, name.c_str(), tmpfolder, name.c_str());
			//printf("%s\n", tmpstr);fflush(stdout);
			system(tmpstr);
			sprintf(tmpstr, "rm %s%s.tc", tmpfolder, name.c_str());
			system(tmpstr);

			sprintf(tmpstr, "bedGraphToBigWig %s%s.tb %s %s%s", tmpfolder, name.c_str(), gszfile, tmpfolder, name.c_str());
			//printf("%s\n", tmpstr);fflush(stdout);
			system(tmpstr);
			sprintf(tmpstr, "rm %s%s.tb", tmpfolder, name.c_str());
			system(tmpstr);
		}
		
		for(j = 0; j < (int)bedfile.size(); j++)
		{	ostringstream convert;
			if((int)chr.size() == 0) convert<<tmpfolder<<flist[i][0]<<"."<<flist[i][1]<<".bed";
			else convert<<tmpfolder<<flist[i][0]<<"."<<flist[i][1]<<chr[j].c_str()<<".bed";
		
			string bedname = convert.str();
			//sprintf(tmpstr, "bigWigAverageOverBed -minMax tmp/%s %s stdout | cut -f%d > %s", name.c_str(), bedfile[j].c_str(), bedcol, bedname.c_str());
			sprintf(tmpstr, "bigWigAverageOverBed -minMax %s%s %s %s", tmpfolder, name.c_str(), bedfile[j].c_str(), bedname.c_str());
			system(tmpstr);
			bednames.push_back(bedname);
	
			ifstream tf(bedname.c_str());
			vector<double> score(bedsz[j], 0);
			string ttmp;
			double sum = 0, ts;
			for(int iii = 0; getline(tf, ttmp); iii++)
			{	istringstream buf(ttmp);
		                istream_iterator<string> beg(buf), end;
		                vector<string> tokens(beg, end);
				ts = atof(tokens[bedcol - 1].c_str());
			//	if(bedcol == 5 || bedcol==6) ts *= 200.;
				score[atoi(tokens[0].c_str())] = ts;
				sum += ts;
			}
			tf.close();
			if(norm)
			{	for(int iii = 0; iii < bedsz[j]; iii++)
				{	score[iii] = score[iii] / sum * 100000000.;
				}
			}
			f = fopen(bedname.c_str(), "w");
			if((int)subindex.size() > 0)
			{	for(int jjj = 0; jjj < (int)subindex.size(); jjj++)
				{	fprintf(f, "%3.2f\n", score[subindex[jjj]]);
				}
			}
			else 
			{	for(int jjj = 0; jjj < bedsz[j]; jjj++)
				{	fprintf(f, "%3.2f\n", score[jjj]);
				}
			}
			fclose(f);

			sprintf(tmpstr, "gzip -f %s", bedname.c_str());
			system(tmpstr);
		}
		sprintf(tmpstr, "rm %s%s", tmpfolder, name.c_str());
		system(tmpstr);
	}
	for(i = 0; i < (int)bedfile.size(); i++)
	{	sprintf(tmpstr, "rm -f %s", bedfile[i].c_str());
		system(tmpstr);
	}
	bedfile = obedfile;
return;
	
	ostringstream legend;
	legend << "#chr pos_st pos_ed ";
	for(i = 0; i < (int)flist.size(); i++)
	{	legend << flist[i][0] << "." << flist[i][1] << " ";
	}
	for(j = 0; j < (int)bedfile.size(); j++)
	{	printf("merge files %d ... ", j);fflush(stdout);
		sprintf(tmpstr, "cut -f1-3 %s | paste -", bedfile[j].c_str());
		for(i = j; i < (int)bednames.size(); i+=(int)bedfile.size())
		{	sprintf(tmpstr, "%s %s", tmpstr, bednames[i].c_str());
		}
		char myfout[1000];
		if((int)chr.size() == 0) sprintf(myfout, "%s", fout);
		else sprintf(myfout, "%s.%s", fout, chr[j].c_str());
		char tfout[1000];
		sprintf(tfout, "%st%s", tmpfolder, myfout);
		sprintf(tmpstr, "%s > %s", tmpstr, tfout);
		system(tmpstr);

		sprintf(tmpstr, "echo \"%s\" | cat - %s > %s%s", legend.str().c_str(), tfout, tmpfolder, myfout);
		system(tmpstr);
		sprintf(tmpstr, "rm %s", tfout);
		system(tmpstr);
	
		for(i = j; i < (int)bednames.size(); i+=(int)bedfile.size())
		{	sprintf(tmpstr, "rm %s", bednames[i].c_str());
			system(tmpstr);
		}
		printf("done\n");fflush(stdout);
	}
}

bool DirectoryExists(const char* dirName)
{	//DWORD attribs = ::GetFileAttributesA(dirName);
	//if(attribs == INVALID_FILE_ATTRIBUTES) return false;
	//else return (attribs & FILE_ATTRIBUTE_DIRECTORY);
	struct stat st;
	if(stat(dirName, &st)==0 && S_ISDIR(st.st_mode)) return true;
	else return false;
}
