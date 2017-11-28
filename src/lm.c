#include "lm.h"
#include <math.h>
#include <string.h>

inline void mul_vec(const double A[], const double B[], double result[], const size_t N)
{
	size_t i;
	for(i=0;i<N;i++) result[i] = A[i]*B[i];
}

inline double dot_vec(const double A[], const double B[], const size_t N)
{
	size_t i;
	double result=0;

	for(i=0;i<N;i++) result += A[i]*B[i];
	return result;
}


int linear_regression(const double X[], const double Y[], const size_t N, const size_t p, const int flag_NA[], double tvalue[], double pvalue[], double I[], double pred[])
{
	size_t i,j, degree = N-p, N_NA=0;
	const double *arr;
	double sigma = 0, err;

	for(i=0;i<N;i++) N_NA += flag_NA[i];
	degree -= N_NA;

	for(i=0, arr=X; i<p; i++, arr+=N)
	{
		memcpy(pred, arr, N*sizeof(double));

		if(N_NA > 0)
		{
			for(j=0;j<N;j++)
			{
				if(flag_NA[j] == 1) pred[j] = 0;
			}
		}

		// X'Y
		pvalue[i] = dot_vec(pred, Y, N);

		// X'X to lower triangle
		for(j=0;j<=i;j++) I[i*p + j] = dot_vec(pred, X + j*N, N);
	}

	gsl_matrix_view Iv = gsl_matrix_view_array(I, p, p);

	// I = (X'X)^-1
	if(gsl_linalg_cholesky_decomp(&Iv.matrix) == GSL_EDOM) return OLS_FAIL;

	gsl_linalg_cholesky_invert(&Iv.matrix);

	// beta = I * X'Y
	for(i=0,arr=I; i<p; i++,arr+=p) tvalue[i] = dot_vec(arr, pvalue, p);

	// Y pred = X * beta
	memset(pred, 0, sizeof(double)*N);

	for(i=0, arr=X; i<p; i++,arr+=N)
	{
		err = tvalue[i];
		for(j=0;j<N;j++) pred[j] += arr[j]*err;
	}

	// calculate standard error
	for(i=0; i<N; i++)
	{
		if(flag_NA[i]==0){
			err = Y[i] - pred[i];
			sigma += err*err;
		}
	}

	sigma /= degree;

	for(i=0;i<p;i++)
	{
		err = sigma* I[i*(p+1)];

		if(fabs(err) < EPS)
			err = EPS;
		else
			err = sqrt(err);

		tvalue[i] /= err;
		pvalue[i] = 2*gsl_cdf_tdist_Q(fabs(tvalue[i]), degree);
	}

	return OLS_SUCCESS;
}





int linear_regression_combo(
	const double X[], const double Y[], const size_t N, const size_t p, const size_t m,
	const double r1, const double r2,
	const int flag_NA[], double tvalue[], double pvalue[], double I[], double pred[])
{
	size_t i,j, degree = N-p, N_NA=0;
	const double *arr;
	double sigma = 0, err;

	for(i=0;i<N;i++) N_NA += flag_NA[i];
	degree -= N_NA;

	for(i=0, arr=X; i<p; i++, arr+=N)
	{
		memcpy(pred, arr, N*sizeof(double));

		if(N_NA > 0)
		{
			for(j=0;j<N;j++)
			{
				if(flag_NA[j] == 1) pred[j] = 0;
			}
		}

		// X'Y
		pvalue[i] = dot_vec(pred, Y, N);

		// X'X to lower triangle
		for(j=0;j<=i;j++) I[i*p + j] = dot_vec(pred, X + j*N, N);
	}

	gsl_matrix_view Iv = gsl_matrix_view_array(I, p, p);

	// I = (X'X)^-1
	if(gsl_linalg_cholesky_decomp(&Iv.matrix) == GSL_EDOM) return OLS_FAIL;

	gsl_linalg_cholesky_invert(&Iv.matrix);

	// beta = I * X'Y
	for(i=0,arr=I; i<p; i++,arr+=p) tvalue[i] = dot_vec(arr, pvalue, p);

	// Y pred = X * beta
	memset(pred, 0, sizeof(double)*N);

	for(i=0, arr=X; i<p; i++,arr+=N)
	{
		err = tvalue[i];
		for(j=0;j<N;j++) pred[j] += arr[j]*err;
	}

	// calculate standard error
	for(i=0; i<N; i++)
	{
		if(flag_NA[i]==0){
			err = Y[i] - pred[i];
			sigma += err*err;
		}
	}

	sigma /= degree;

	j = p - m - 1;

	for(i=p-m; i<p; i++)
	{
		// var(j) + var(i) + cov(j, i)
		err = SQR(r1) * I[j*(p+1)] + SQR(r2) * I[i*(p+1)] + 2*r1*r2*I[j*p + i];
		err *= sigma;

		if(fabs(err) < EPS)
			err = EPS;
		else
			err = sqrt(err);

		tvalue[i] = (r1*tvalue[j] + r2*tvalue[i])/err;
		pvalue[i] = 2*gsl_cdf_tdist_Q(fabs(tvalue[i]), degree);
	}

	return OLS_SUCCESS;
}



void print_matrix(const gsl_matrix *X, FILE *fp)
{
	size_t i,j;

	for(i=0;i<X->size1;i++)
	{
		for(j=0;j<X->size2;j++) fprintf(fp,"%f\t,\t",gsl_matrix_get(X,i,j));

		fprintf(fp,"\n");
	}
}

void print_vector(const double A[], const size_t N)
{
	size_t i;

	for(i=0;i<N;i++) fprintf(stdout,"%f\t,\t",A[i]);
	fprintf(stdout, "\n");
}
