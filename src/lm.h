#ifndef LM_H_
#define LM_H_

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_linalg.h>

#define EPS 1e-10
#define OLS_SUCCESS 1
#define OLS_FAIL 0

#define SQR(x) ((x)*(x))

void mul_vec(const double A[], const double B[], double result[], const size_t N);
double dot_vec(const double A[], const double B[], const size_t N);

int linear_regression(
	// input
	const double X[], const double Y[], const size_t N, const size_t p,

	// position of NA values in Y
	const int flag_NA[],

	// output
	double tvalue[], double pvalue[],

	// workspace
	double I[], double pred[]
	);


// specialized function for CARE analysis
int linear_regression_combo(
	// input
	const double X[], const double Y[], const size_t N, const size_t p, const size_t m,

	// combination ratios
	const double r1, const double r2,

	// position of NA values in Y
	const int flag_NA[],

	// output
	double tvalue[], double pvalue[],

	// workspace
	double I[], double pred[]
	);



void print_matrix(const gsl_matrix *X, FILE *fp);
void print_vector(const double A[], const size_t N);

#endif
