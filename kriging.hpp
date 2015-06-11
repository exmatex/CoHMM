/** header file for kriging.cpp
 * **/
#ifndef KRIGING_HPP
#define KRIGING_HPP

#include <vector>

int kriging(double * w, int pointDims, std::vector<double *> oldWs, std::vector<double> oldVals, double* retVal);

//#ifndef HAVE_MKL
extern "C"
{ 
void dgetrs_(char* trans, int* n, int* nrhs, double* A, int* lda, int* ipiv, double* B, int* ldb, int* info);   
void dgetrf_(int* m, int* n, double* A, int* lda, int* ipiv, int* info);   
void dgemm_(char *transa, char *transb, int *m, int *n, int *k, double *alpha, double *a, int *lda, double *b, int *ldb, double *beta, double *c, int *ldc);
void* MKL_calloc(size_t num, size_t size, int align);
#define mkl_calloc MKL_calloc
void* MKL_malloc(size_t size, int align);
#define mkl_malloc MKL_malloc
void MKL_free(void *ptr);
#define mkl_free MKL_free
}
//#endif

#endif
