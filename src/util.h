#ifndef UTIL_H_
#define UTIL_H_

#include <string>
#include <vector>
using namespace std;

#define check(flag,msg) if(flag){cerr << msg << endl; exit(1);}

// split string into vector by separator
void split(string line, vector<string> &vec, const char sep);

// convert string to double and monitor NA
double convert_to_double(const string n, const bool NAflag);

// get overlap between sorted vectors (Please sort vectors before)
void intersect(const vector<string> &v1, const vector<string> &v2, vector<string> &result);

// FDR correction for p values
void Benjamini_Hochberg(double FDR[], double pvalue[], const size_t N);

#endif
