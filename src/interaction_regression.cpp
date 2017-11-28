extern "C" {
#include "lm.h"
}

#include "util.h"
#include "Matrix_Map.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <set>
#include <limits>
#include <cstdlib>
#include <cstring>
#include <cmath>

using namespace std;


// general format of output result matrix
void output_matrix(ofstream &fout, const string &title, const vector<string> &rownames, const double mat[], const size_t N)
{
	size_t i;
	vector<string>::const_iterator iter = rownames.begin();

	for(; iter != rownames.end(); iter++, mat += N)
	{
		fout << title << '\t' << *iter;
		for(i=0;i<N;i++) fout << '\t' << mat[i];
		fout << '\n';
	}
}


// data: column major order
void normalize_features(double data[], const size_t n, const vector<string> &features)
{
	size_t i, j, p = features.size();

	// flag for non-category variables (not mutation)
	double aver, sigma;

	for(i=0; i<p; i++, data += n)
	{
		// don't normalize mutation categorical columns
		if(features[i].find("Mutation") != string::npos) continue;

		for(aver=sigma=0, j=0;j<n;j++)
		{
			aver += data[j];
			sigma += SQR(data[j]);
		}

		aver/= n;
		sigma = sigma/(n-1) - SQR(aver) * n/(n-1);

		if(fabs(sigma) < EPS) continue;

		sigma = sqrt(sigma);

		for(j=0;j<n;j++) data[j] = (data[j] - aver)/sigma;
	}
}



int main(int argc, char *argv[])
{
	string type, line, input_X, input_Y, input_T, input_B, output;
	size_t i, j, k, parseCnt = (argc-1)/2, N_data, N_target, N_dim, N_response, N_feature, N_maxtarget = 1, N_background = 0;

	// show all error messages
	bool verbose = false, partial = false;

	double r1 = 1, r2 = 1;

	int flag, *flag_NA;

	ifstream fin;
	ofstream fout;

	istringstream iss;

	// column and row names
	vector<string> features, targets, responses, common_rownames;
	vector<string>::iterator viter;

	// matrix map is used to load in data first
	Matrix_Map *response_map, *feature_map, *background_map = NULL;

	// feature name to index
	map<string, size_t> inxmap_feature;
	map<string, size_t>::iterator miter;

	// target indices
	set<size_t> targets_inx;
	set<size_t>::iterator siter;

	// response to targets map
	map<string, set<size_t> > target_map;
	map<string, set<size_t> >::iterator titer;

	// cache matrix for features
	double *data_feature, *data_response, *data_background = NULL, *data_ptr, *data_pivot, *data_target,
		*tvalue, *pvalue, *X, *Y, *I, *pred, *result_t, *result_p;

	if (argc < 7)
	{
		if(argc == 2 && string(argv[1]) == "-help")
		{
			cout << endl;
			cout << "interaction regression\n" << endl;

			cout << "Usage: interaction_regression -x Feature -y Response -o Output\n" << endl;

			cout << "\t-t target annotation. Default: empty (response name as target)" << endl;
			cout << "\t-b background matrix. Default: empty" << endl;
			cout << "\t-p partial regression (0 or 1). Default: " << (partial?1:0) << endl;
			cout << "\t-r1 combination ratio 1. Default: " << r1 << endl;
			cout << "\t-r2 combination ratio 2. Default: " << r2 << endl;
			cout << "\t-v verbose (0 or 1). Default: " << (verbose?1:0) << endl;

			cout << "\nReport bugs to Peng Jiang (peng.jiang.software@gmail.com)\n" <<endl;
			exit(0);
		}else{
			cerr << "Insufficient number of arguments, do \"interaction_regression -help\" for help."<<endl;
			exit(1);
		}
	}

	// read in all parameters
	for(i=0;i<parseCnt;i++)
	{
		type = argv[2*i+1];
		line = argv[2*i+2];

		if(type == "-x"){
			input_X = line;

		}else if (type == "-y"){
			input_Y = line;

		}else if (type == "-t"){
			input_T = line;

		}else if (type == "-b"){
			input_B = line;

		}else if (type == "-o"){
			output = line;

		}else if (type == "-v"){
			check(line != "0" && line != "1", "Please input 1 or 0 for verbose flag.")
			verbose = (line != "0");

		}else if (type == "-p"){
			check(line != "0" && line != "1", "Please input 1 or 0 for partial regression flag.")
			partial = (line != "0");

		}else if (type == "-r1"){
			r1 = convert_to_double(line, false);

		}else if (type == "-r2"){
			r2 = convert_to_double(line, false);

		}else if (type == "-help"){
			cerr << "Please don't use \"-help\" as parameter input." << endl;
			exit(1);

		}else{
			cerr << "Cannot recognize parameter \""<< type << "\"." << endl;
			exit(1);
		}
	}

	check(input_X.empty(), "Cannot find Feature")
	check(input_Y.empty(), "Cannot find Response")
	check(output.empty(), "Cannot find output name")

	// load in all data matrix first
	response_map = new Matrix_Map(input_Y);
	feature_map = new Matrix_Map(input_X);

	intersect(response_map->rownames, feature_map->rownames, common_rownames);

	// has background
	if(!input_B.empty())
	{
		background_map = new Matrix_Map(input_B);
		intersect(common_rownames, background_map->rownames, features);
		common_rownames = features;

		// background map to column major matrix
		N_background = background_map->header.size();
		data_background = background_map->convert_to_array(common_rownames);
		delete background_map;
	}

	// response map to column major matrix
	responses = response_map->header;
	N_response = responses.size();
	data_response = response_map->convert_to_array(common_rownames);
	delete response_map;

	// features map to column major matrix
	features = feature_map->header;
	N_feature = features.size();
	data_feature = feature_map->convert_to_array(common_rownames);
	delete feature_map;

	// row dimension of all matrices
	N_data = common_rownames.size();

	// normalize features to interpret low order coefficients
	normalize_features(data_feature, N_data, features);


	// build feature index map
	for(i=0, viter = features.begin(); viter != features.end(); i++, viter++)
	{
		check(inxmap_feature.find(*viter) != inxmap_feature.end(), "Duplicated feature " << *viter)
		inxmap_feature[*viter] = i;
	}

	// load in target annotations
	if(!input_T.empty())
	{
		fin.open(input_T.c_str());
		check(fin.fail(), "Cannot open target file " << input_T)

		while(getline(fin, line, '\n'))
		{
			iss.str(line);
			getline(iss, line, '\t');
			getline(iss, input_T, '\t');
			iss.clear();

			if(input_T.empty()) continue;

			// convert feature name to index
			split(input_T, targets, ',');

			for(viter=targets.begin(); viter!= targets.end(); viter++)
			{
				miter = inxmap_feature.find(*viter);

				if (miter!=inxmap_feature.end())
					targets_inx.insert(miter->second);
				else if (verbose)
					cerr << "Cannot find " << *viter << " for " << line << endl;
			}

			if(!targets_inx.empty())
			{
				target_map[line] = targets_inx;
				N_maxtarget = max<size_t>(N_maxtarget, targets_inx.size());
				targets_inx.clear();
			}
		}

		fin.close();
		fin.clear();

	}else{
		// no target annotation, use response name as target
		for(viter=responses.begin(); viter!= responses.end(); viter++)
		{
			miter = inxmap_feature.find(*viter);

			if (miter!=inxmap_feature.end())
			{
				targets_inx.insert(miter->second);
				target_map[*viter] = targets_inx;
				targets_inx.clear();

			}else if(verbose)
				cerr << "Cannot find " << *viter << endl;
		}
	}

	// X: intercept, background, targets, pivot, targets * pivot
	N_dim = 1 + N_background + N_maxtarget + 1 + N_maxtarget;

	// data points should be larger than dimension
	check(N_data < N_dim + 1, "Only "<< N_data << " samples. At least " << N_dim + 1 << " required.")


	// allocate space for linear regression
	X = new double[N_data * N_dim];
	for(i=0;i<N_data;i++) X[i] = 1;	// intercept

	I = new double[N_dim*N_dim];
	pred = new double[N_data];
	tvalue = new double[N_dim];
	pvalue = new double[N_dim];
	flag_NA = new int[N_data];	// marker of NA value in Y

	// allocate space for result output: possible targets
	result_t = new double[N_feature * (N_maxtarget + 1)];
	result_p = new double[N_feature * (N_maxtarget + 1)];

	// process each response one by one
	fout.open(output.c_str());
	check(fout.fail(), "Cannot open output " << output)

	// output header
	fout << "Response\tType\tTarget";
	for(viter = features.begin(); viter!=features.end(); viter++) fout << '\t' << *viter;
	fout << '\n';

	// progress step
	parseCnt = max<size_t>(N_response/100,1);
	gsl_set_error_handler_off();

	for(i=0, Y=data_response; i<responses.size(); i++, Y+=N_data)
	{
		// progress bar
		if(i%parseCnt==0) cout << round(100.0*float(i+1)/N_response) << '%' << endl;

		// no target annotation
		titer = target_map.find(responses[i]);
		if(titer == target_map.end()) continue;
		set<size_t> &s = titer->second;

		N_target = s.size();
		check(N_target==0, "Impossible zero targets")

		// clear up result cache: targets + self covariate
		N_dim = N_feature * (N_target+1);

		for(j=0;j<N_dim;j++)
		{
			result_t[j] = 0;
			result_p[j] = 1;
		}

		// dimension of current regression
		if(partial){
			// intercept, background, target, feature
			N_dim = 1 + N_background + N_target + 1;
		}else{
			// intercept, background, target, feature, target * feature
			N_dim = 1 + N_background + N_target + 1 + N_target;
		}

		// count NA values in Y
		for(j=0; j<N_data; j++)
		{
			if(isnan(Y[j])){
				flag_NA[j] = 1;
				Y[j] = 0;
			}else
				flag_NA[j] = 0;
		}

		// prepare X: intercept, background, targets, pivot, pivot*targets

		// background
		if(N_background > 0)
			memcpy(X + N_data, data_background, N_data * N_background * sizeof(double));

		// targets
		data_target = X + (1+N_background) * N_data;

		targets.clear();

		if(partial) targets.push_back("Self");

		for(siter=s.begin(); siter != s.end(); siter++, data_target += N_data)
		{
			j = *siter;
			memcpy(data_target, data_feature + j*N_data, N_data*sizeof(double));

			if(!partial) targets.push_back(features[j]);
		}

		// current location of data target is for data pivot
		data_pivot = data_target;

		// each pivot as feature
		for(j=0;j<N_feature;j++)
		{
			// jump all targets in features
			if(s.find(j) != s.end()) continue;

			// pivot
			memcpy(data_pivot, data_feature + j*N_data, N_data*sizeof(double));

			// start of target region
			data_target = X + (1+N_background)*N_data;

			// targets*pivot
			for(k=0; k<N_target; k++, data_target += N_data)
			{
				data_ptr = data_target + (1+N_target) *N_data;
				mul_vec(data_target, data_pivot, data_ptr, N_data);
			}

			// linear regression
			if(partial)
				flag = linear_regression(X, Y, N_data, N_dim, flag_NA, tvalue, pvalue, I, pred);
			else
				flag = linear_regression_combo(X, Y, N_data, N_dim, N_target, r1, r2, flag_NA, tvalue, pvalue, I, pred);

			if(flag == OLS_SUCCESS)
			{
				if(partial)
				{
					result_t[j] = tvalue[N_dim-1];
					result_p[j] = pvalue[N_dim-1];

				}else{
					for(k=0; k < N_target; k++)
					{
						result_t[k*N_feature + j] = tvalue[N_dim - N_target + k];
						result_p[k*N_feature + j] = pvalue[N_dim - N_target + k];
					}
				}

			}else if(verbose)
				cerr << "Jump failure " << responses[i] << '\t' << features[j] << endl;
		}

		// output results
		output_matrix(fout, responses[i] + "\tt", targets, result_t, N_feature);
		output_matrix(fout, responses[i] + "\tp", targets, result_p, N_feature);

		// FDR correction
		for(j=0,data_ptr=result_p; j < targets.size(); j++, data_ptr+=N_feature) Benjamini_Hochberg(data_ptr, data_ptr, N_feature);

		output_matrix(fout, responses[i] + "\tfdr", targets, result_p, N_feature);
	}

	fout.close();

	if(N_background > 0) delete[] data_background;

	delete[] data_response;
	delete[] data_feature;
	delete[] X;
	delete[] I;
	delete[] pred;
	delete[] tvalue;
	delete[] pvalue;
	delete[] flag_NA;
	delete[] result_t;
	delete[] result_p;

	return 0;
}
