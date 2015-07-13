%module cohmm_swiftt
%{
extern bool initEverything(bool doKriging, bool doCoMD, int * dims, double * dt, double * delta, double * gamma);

extern int prepFirstFlux(bool doKriging, bool doCoMD, int * dims, double * dt, double * delta, double * gamma, int curStep);
extern int prepSecondFlux(bool doKriging, bool doCoMD, int * dims, double * dt, double * delta, double * gamma, int curStep);
extern int prepThirdFlux(bool doKriging, bool doCoMD, int * dims, double * dt, double * delta, double * gamma, int curStep);
extern int prepLastFlux(bool doKriging, bool doCoMD, int * dims, double * dt, double * delta, double * gamma, int curStep);
extern int finishStep(bool doKriging, bool doCoMD, int * dims, double * dt, double * delta, double * gamma, int curStep);

extern bool cloudFlux(bool doKriging, bool doCoMD, int curStep, int phase, int taskID);

extern bool outputVTK(bool doKriging, bool doCoMD, int * dims, double * dt, double * delta, double * gamma, int curStep);

extern bool tryShortCircuit(int * dims, int curStep);

extern int checkStepForFaults(int * dims, int curPhase, int curStep);
extern bool retryCloudFlux(bool doKriging, bool doCoMD, int curStep, int phase, int taskID, int round);
%}

extern bool initEverything(bool doKriging, bool doCoMD, int * dims, double * dt, double * delta, double * gamma);

extern int prepFirstFlux(bool doKriging, bool doCoMD, int * dims, double * dt, double * delta, double * gamma, int curStep);
extern int prepSecondFlux(bool doKriging, bool doCoMD, int * dims, double * dt, double * delta, double * gamma, int curStep);
extern int prepThirdFlux(bool doKriging, bool doCoMD, int * dims, double * dt, double * delta, double * gamma, int curStep);
extern int prepLastFlux(bool doKriging, bool doCoMD, int * dims, double * dt, double * delta, double * gamma, int curStep);
extern int finishStep(bool doKriging, bool doCoMD, int * dims, double * dt, double * delta, double * gamma, int curStep);

extern bool cloudFlux(bool doKriging, bool doCoMD, int curStep, int phase, int taskID);

extern bool outputVTK(bool doKriging, bool doCoMD, int * dims, double * dt, double * delta, double * gamma, int curStep);

extern bool tryShortCircuit(int * dims, int curStep);

extern int checkStepForFaults(int * dims, int curPhase, int curStep);
extern bool retryCloudFlux(bool doKriging, bool doCoMD, int curStep, int phase, int taskID, int round);
