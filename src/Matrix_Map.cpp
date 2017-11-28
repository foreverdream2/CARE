#include "Matrix_Map.h"
#include "util.h"

#include <iostream>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <cstdlib>
using namespace std;

Matrix_Map::Matrix_Map(const string file_name)
{
	size_t i, N;
	string line;
	istringstream iss;

	double *data;

	ifstream fin(file_name.c_str());
	check(fin.fail(), cerr <<  "Cannot open file " << file_name)

	getline(fin, line, '\n');
	split(line, header, '\t');
	N = header.size();

	while(getline(fin, line, '\n'))
	{
		iss.str(line);
		getline(iss, line, '\t');

		check(data_map.find(line) != data_map.end(), "Duplicated response ID " << line)

		rownames.push_back(line);
		data = data_map[line] = new double[N];

		for(i=0;i<N;i++)
		{
			getline(iss, line, '\t');
			data[i] = convert_to_double(line, true);
		}

		iss.clear();
	}

	fin.close();

	// order row names for intersection operation
	sort(rownames.begin(), rownames.end());
}


Matrix_Map::~Matrix_Map()
{
	for (map<string, double*>::iterator iter=data_map.begin(); iter!=data_map.end(); iter++)
		delete[] iter->second;
}


double *Matrix_Map::convert_to_array(const vector<string> &row_order)
{
	size_t i, j, N_col = header.size(), N_row = row_order.size();
	vector<string>::const_iterator iter;
	map<string, double*>::const_iterator miter;

	double *data = new double[N_row * N_col];

	for (iter = row_order.begin(), i=0; iter != row_order.end(); iter++, i++)
	{
		miter = data_map.find(*iter);
		check(miter == data_map.end(), "Cannot find " << *iter)

		for(j=0;j<N_col;j++) data[j*N_row + i] = miter->second[j];
	}

	return data;
}
