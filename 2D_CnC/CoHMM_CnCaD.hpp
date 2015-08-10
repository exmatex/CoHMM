#ifndef COHMM_CNCAD_HPP
#define COHMM_CNCAD_HPP

#include <hiredis.h>

#include "2DKriging.hpp"

#include "CoHMM_CnC.hpp"

bool initEverything(bool doKriging, bool doCoMD, int * dims, double * dt, double * delta, double * gamma, const char * redis_host, CnCaDContext &ctxt);
int prepFirstFlux(bool doKriging, bool doCoMD, int * dims, double * dt, double * delta, double * gamma, int curStep, const char * redis_host, CnCaDContext &ctxt);
int prepSecondFlux(bool doKriging, bool doCoMD, int * dims, double * dt, double * delta, double * gamma, int curStep, const char * redis_host, CnCaDContext &ctxt);
int prepThirdFlux(bool doKriging, bool doCoMD, int * dims, double * dt, double * delta, double * gamma, int curStep, const char * redis_host, CnCaDContext &ctxt);
int prepLastFlux(bool doKriging, bool doCoMD, int * dims, double * dt, double * delta, double * gamma, int curStep, const char * redis_host, CnCaDContext &ctxt);
int finishStep(bool doKriging, bool doCoMD, int * dims, double * dt, double * delta, double * gamma, int curStep, const char * redis_host, CnCaDContext &ctxt);
bool cloudFlux(bool doKriging, bool doCoMD, int curStep, int phase, int taskID, const char * redis_host, CnCaDContext &ctxt);
bool outputVTK(bool doKriging, bool doCoMD, int * dims, double * dt, double * delta, double * gamma, int curStep, const char * redis_host, CnCaDContext &ctxt);


FluxOut fluxFn(bool doKriging, bool doCoMD, FluxIn * input, redisContext * headRedis, CnCaDContext &ctxt);
FluxOut randomCoMDImbalance(FluxIn * input);

bool checkRedisHost(const char * inHost);
char * getRedisHost(const char * filePath);

#endif
