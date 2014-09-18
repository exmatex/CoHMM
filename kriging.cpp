#include <cstdio>
#include <vector>
#include <cstdlib>
#include <cmath>

//MKL
#ifdef HAVE_MKL
#include <mkl.h>
#else
#include <lapacke.h>
#include <cblas.h>
#endif

enum variogramApprox_t
{
	SPHERE,
	EXPON,
	GAUSS
};

/*
 * Inputs
 *   int n: the number of points we care about
 *   double [n * strideX] row: the row
 *   int strideX: the stride between points in the row
 *   double [n * strideY] col: the col
 *   int strideY: the stride between points in the col
 *   variogramApprox_t approx approximationMode
 * Output
 *   double: the variogram'd value
 */
double variogram(int n, double *row, int strideX, double * col, int strideY, variogramApprox_t approx)
{
	const double c0 = 0.0;
	const double c1 = 1.0;
	const double a = 0.25;
	
	///Calculate the distance between sample datas
	double l = 0.0;
	for(int i = 0; i < n; i++)
	{
		double x = row[i*strideX];
		double y = col[i*strideY];

		l += pow(x-y, 2.0);
	}
	double h = sqrt(l);

	double V = 0.0;

	///Approximation mode select
	switch(approx)
	{
	case SPHERE:
		if(h < a)
		{
			V = c0 + c1 * ( (1.5 * (h/a)) - (0.5 * pow(h/a, 3)));
		}
		else
		{
			V = c0 + c1;
		}
		break;
	case EXPON:
		V = c0 + c1 * (1 - exp(-3*h/a));
		break;
	case GAUSS:
		V = c0 + c1 * (1 - exp(-pow(3*h,2)/(a*a)));
		break;
	default:
		break;
	}
	return V;
}

/*
 * Inputs
 *   double [pointSize]: w
 *   int pointSize: the number of dimensions (for us, that is 7)
 *   vector<double[pointSize]>: oldWs
 *   vector<double>: outputs corresponding to oldWs
 * Output
 *   double[2]: Estimated value for W and error
 */
int kriging(double * w, int pointDims, std::vector<double *> oldWs, std::vector<double> oldVals, double *retVal)
{
	if(oldWs.size() == 0)
	{
		retVal[1] = 9.9999999;
		return 0;
	}
	int n = oldVals.size() + 1;
    ///Z(x_i) known values of surrounding points + 1 
	//Copy oldVals into a new array, add a 0 to the end, and transpose
#ifndef HAVE_MKL
	double *Z = (double *)malloc(sizeof(double) * n);
	//double *bufferZ = (double *)malloc(sizeof(double) * n + 7);
	//double *Z = (void*)( ((uintptr_t)bufferZ+63) & ~ (uintptr_t)0x0F );
#else
	double * Z = (double *)mkl_malloc(sizeof(double) * n, 64);
#endif
	
	for(int i = 0; i < n-1; i++)
	{
		Z[i] = oldVals[i];
	}
	Z[n-1] = 0.0;

	///Prepare matrices with Variogram model
#ifndef HAVE_MKL
	double *K = (double *)calloc(n*n, sizeof(double));
	//double *bufferK = (double *)calloc(n*n + 7, sizeof(double));
	//double *K = (void*)( ((uintptr_t)bufferK+63) & ~ (uintptr_t)0x0F );
#else
	double * K = (double *)mkl_calloc(n*n, sizeof(double), 64);
#endif

	//Variogram the inside
	for(int i = 0; i < n-1; i++)
	{
		double * col = oldWs[i];
		for(int j = 0; j < n-1; j++)
		{
			double * row = oldWs[j];
			K[i + j*(n)] = variogram(pointDims, row, 1, col, 1, SPHERE);
		}
	}
	//Set the outsides to 1s and 0s
	for(int j = 0; j < n; j++)
	{
		K[n-1 + j*(n)] = 1.0;
	}
	for(int i = 0; i < n; i++)
	{
		K[i + (n-1)*(n)] = 1.0;
	}
	K[n-1 + (n-1)*(n)] = 0.0;

#ifndef HAVE_MKL
	double * M = (double *)malloc(sizeof(double) * n);
	//double * bufferM = (double *)malloc(sizeof(double) * n + 7);
	//double * M = (void*)( ((uintptr_t)bufferM+63) & ~ (uintptr_t)0x0F );
#else
	double * M = (double *)mkl_malloc(sizeof(double) * n, 64);
#endif
#ifndef HAVE_MKL
	double * Mo = (double *)malloc(sizeof(double) * n);
	//double * bufferMo = (double *)malloc(sizeof(double) * n+7);
	//double * Mo = (void*)( ((uintptr_t)bufferMo+64) & ~ (uintptr_t)0x0F );
#else
	double * Mo = (double *)mkl_malloc(sizeof(double) * n, 64);
#endif
	for(int i = 0; i < n-1; i++)
	{
		M[i] = variogram(pointDims, oldWs[i], 1, w, 1, SPHERE);
		Mo[i] = M[i];
	}
	M[n-1] = 1.0;
	Mo[n-1] = 1.0;

	///Proesses Kriging system KW = M
	//W = np.linalg.solve(K,M)
	//Need to LU fatorize K
	int * ipiv = (int *)malloc(sizeof(int) * n);
	int info;
	info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, K, n, ipiv);
	if(info != 0)
	{
		fprintf(stderr, "LAPACK dgetrf failed on %d \n", info);

		return 0;
	}
	//KW=M, solve for W
	info = LAPACKE_dgetrs(LAPACK_ROW_MAJOR, 'N', n, 1, K, n, ipiv, M, 1);

	//info = LAPACKE_dgetrs(LAPACK_ROW_MAJOR, 'N', n, 1, K, n, M, n, 1, ipiv);
	if(info != 0)
	{
		fprintf(stderr, "LAPACK dgetrs failed on %d \n", info);
		return 0;
	}
	double * W = M;
	



	M = Mo;
	///Process results
/*
void cblas_dgemm
(
	const CBLAS_ORDER Order,
	const CBLAS_TRANSPOSE TransA,
	const CBLAS_TRANSPOSE TransB,
	const MKL_INT M,
	const MKL_INT N,
	const MKL_INT K,
	const double alpha,
	const double *A,
	const MKL_INT lda,
	const double *B,
	const MKL_INT ldb,
	const double beta,
	double *C,
	const MKL_INT ldc
);
*/

	//Value = W.T*Z
	cblas_dgemm
	(
		CblasRowMajor,
		CblasTrans,
		CblasNoTrans,
		1, //Rows of A and C
		1, //Cols of B and C
		n, //Cols of A and rows of B
		1.0, //alpha
		W, //A
		1, //lda before transpose?
		Z, //B
		1, //ldb before transpose?
		0.0, //beta
		&retVal[0],	//C
		1	//ldc
	);
	//Error = W.T*M
	cblas_dgemm
	(
		CblasRowMajor,
		CblasTrans,
		CblasNoTrans,
		1, //Rows of A and C
		1, //Cols of B and C
		n, //Cols of A and rows of B
		1.0, //alpha
		W, //A
		1, //lda before transpose?
		M, //B
		1, //ldb before transpose?
		0.0, //beta
		&retVal[1],	//C
		1	//ldc
	);

	//Free stuff
#ifndef HAVE_MKL
	free(Z);
	free(K);
	free(M);
	free(W);
	//free(bufferZ);
	//free(bufferK);
	//free(bufferM);
	//free(bufferMo);
#else
	mkl_free(Z);
	mkl_free(K);
	mkl_free(M);
	mkl_free(W);
#endif

	//Return stuff what is important
	return 1;
}



