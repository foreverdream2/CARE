#include "util.h"

#include <iostream>
#include <sstream>
#include <cstdlib>
#include <limits>
#include <algorithm>
using namespace std;

double convert_to_double(const string n, const bool NAflag)
{
	istringstream iss(n);
	double v;

	if(NAflag && n == "NA"){
		v = numeric_limits<double>::quiet_NaN();
	}else{
		iss >> v;

		if(iss.fail())
		{
			cerr << "Cannot convert double for " << n << endl;
			exit(1);
		}
	}

	return v;
}

void split(string line, vector<string> &vec, const char sep)
{
	istringstream iss(line);
	for(vec.clear(); getline(iss,line,sep); vec.push_back(line));
}

void intersect(const vector<string> &v1, const vector<string> &v2, vector<string> &result)
{
	result.clear();

	vector<string>::const_iterator iter1 = v1.begin(), iter2 = v2.begin();

	while (iter1 != v1.end() && iter2 != v2.end())
	{
		if (*iter1 < *iter2)
			iter1++;
		else if (*iter2 < *iter1)
			iter2++;
		else {
			result.push_back(*iter1);
			iter1++;
			iter2++;
		}
	}
}


class rank_node
{
public:
	rank_node(double v, size_t i):v(v),i(i){}
	bool operator < (const rank_node &r) const {return v < r.v;}

	double v;
	size_t i;
};

void Benjamini_Hochberg(double FDR[], double pvalue[], const size_t N)
{
	size_t i;
	double qvalue = 1;
	vector<rank_node> sortvec;

	for(i=0;i<N;i++) sortvec.push_back(rank_node(pvalue[i],i));

	sort(sortvec.begin(), sortvec.end());

	for(i=0;i<N;i++)
	{
		qvalue = min<double>(qvalue, sortvec[N-1-i].v*N/(N-i));
		FDR[sortvec[N-1-i].i] = qvalue;
	}
}
