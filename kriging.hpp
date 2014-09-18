#ifndef KRIGING_HPP
#define KRIGING_HPP

#include <vector>

int kriging(double * w, int pointDims, std::vector<double *> oldWs, std::vector<double> oldVals, double* retVal);

#endif

