#ifndef COHMM_CNC_HPP
#define COHMM_CNC_HPP

#include "CoHMM.hpp"
#include "CoHMM_Context.hpp"
#include "RedisWrapper.hpp"

class CoHMM_CnC : public CoHMM
{
public:
    CoHMM_CnC(unsigned int dimX, unsigned int dimY, unsigned int nSteps, bool doKriging, bool doCoMD, const char * redisHost = "localhost");
    CoHMM_CnC(unsigned int dims[2], double delta[2], double dt[2], double gamma[3], unsigned int nSteps, bool doKriging, bool doCoMD, const char * redisHost = "localhost");
    virtual ~CoHMM_CnC();

    virtual void processFluxes();

protected:
    virtual bool queueTask(FluxIn &task, unsigned int taskID);
    #ifdef CNC_DIST
	   CnC::dist_cnc_init<CoHMMContext> dinit;
    #endif
	CoHMMContext ctxt;
    char redis_host[RedisWrapper::MAX_HOST_LENGTH];


};

#endif
