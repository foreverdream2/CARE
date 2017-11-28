#ifndef MATRIX_MAP_H_
#define MATRIX_MAP_H_

#include <string>
#include <vector>
#include <map>
using namespace std;


// row major order matrix map
class Matrix_Map {
public:
	Matrix_Map(const string file_name);
	~Matrix_Map();

	vector<string> header, rownames;
	map<string, double*> data_map;

	// convert data map to column major order matrix
	double *convert_to_array(const vector<string> &row_order);

};

#endif
