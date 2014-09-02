/** header file for kriging.cpp
 * **/
#ifndef KRIGING_HPP
#define KRIGING_HPP

#include <vector>

int kriging(double * w, int pointDims, std::vector<double *> oldWs, std::vector<double> oldVals, double* retVal);

#ifndef HAVE_MKL
extern "C"
{ 
void dgetrs_(char* trans, int* n, int* nrhs, double* A, int* lda, int* ipiv, double* B, int* ldb, int* info);   
void dgetrf_(int* m, int* n, double* A, int* lda, int* ipiv, int* info);   
}
#endif

#endif
