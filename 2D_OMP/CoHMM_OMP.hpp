#ifndef COHMM_OMP_HPP
#define COHMM_OMP_HPP

#include "CoHMM.hpp"

class CoHMM_OMP : public CoHMM
{
public:
    CoHMM_OMP(unsigned int dimX, unsigned int dimY, unsigned int nSteps, bool doKriging, bool doCoMD, const char * redisHost = "localhost");
    CoHMM_OMP(unsigned int dims[2], double delta[2], double dt[2], double gamma[3], unsigned int nSteps, bool doKriging, bool doCoMD, const char * redisHost = "localhost");
    virtual ~CoHMM_OMP();

    virtual void processFluxes();

protected:
    virtual bool queueTask(FluxIn &task, unsigned int taskID);

    std::vector<FluxIn> taskQueue;
    std::vector<FluxOut> taskResults;
};

#endif
