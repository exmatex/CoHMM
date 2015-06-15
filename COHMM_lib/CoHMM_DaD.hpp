#ifndef COHMM_SWIFTT_HPP
#define COHMM_SWIFTT_HPP

#include <hiredis.h>

#include "2DKriging.hpp"

//Can't use default arguments because of swig issues
bool initEverything(bool doKriging, bool doCoMD, int * dims, double * dt, double * delta, double * gamma);
bool initEverything(bool doKriging, bool doCoMD, int * dims, double * dt, double * delta, double * gamma, const char * redis_host);
int prepFirstFlux(bool doKriging, bool doCoMD, int * dims, double * dt, double * delta, double * gamma, int curStep);
int prepFirstFlux(bool doKriging, bool doCoMD, int * dims, double * dt, double * delta, double * gamma, int curStep, const char * redis_host);
int prepSecondFlux(bool doKriging, bool doCoMD, int * dims, double * dt, double * delta, double * gamma, int curStep);
int prepSecondFlux(bool doKriging, bool doCoMD, int * dims, double * dt, double * delta, double * gamma, int curStep, const char * redis_host);
int prepThirdFlux(bool doKriging, bool doCoMD, int * dims, double * dt, double * delta, double * gamma, int curStep);
int prepThirdFlux(bool doKriging, bool doCoMD, int * dims, double * dt, double * delta, double * gamma, int curStep, const char * redis_host);
int prepLastFlux(bool doKriging, bool doCoMD, int * dims, double * dt, double * delta, double * gamma, int curStep);
int prepLastFlux(bool doKriging, bool doCoMD, int * dims, double * dt, double * delta, double * gamma, int curStep, const char * redis_host);
int finishStep(bool doKriging, bool doCoMD, int * dims, double * dt, double * delta, double * gamma, int curStep);
int finishStep(bool doKriging, bool doCoMD, int * dims, double * dt, double * delta, double * gamma, int curStep, const char * redis_host);
bool cloudFlux(bool doKriging, bool doCoMD, int curStep, int phase, int taskID);
bool cloudFlux(bool doKriging, bool doCoMD, int curStep, int phase, int taskID, const char * redis_host);
bool outputVTK(bool doKriging, bool doCoMD, int * dims, double * dt, double * delta, double * gamma, int curStep, const char * redis_host);
bool outputVTK(bool doKriging, bool doCoMD, int * dims, double * dt, double * delta, double * gamma, int curStep);

FluxOut fluxFn(bool doKriging, bool doCoMD, FluxIn * input, redisContext * headRedis);

#endif

