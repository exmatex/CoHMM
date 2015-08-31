#ifndef COHMM_CHT_HPP
#define COHMM_CHT_HPP

#include "chunks_and_tasks.h"

#include "CoHMM.hpp"
#include "RedisWrapper.hpp"

class CoHMM_CHT : public CoHMM
{
public:
    CoHMM_CHT(unsigned int dimX, unsigned int dimY, unsigned int nSteps, bool doKriging, bool doCoMD, const char * redisHost = "localhost");
    CoHMM_CHT(unsigned int dims[2], double delta[2], double dt[2], double gamma[3], unsigned int nSteps, bool doKriging, bool doCoMD, const char * redisHost = "localhost");
    virtual ~CoHMM_CHT();

    virtual void processFluxes();

protected:
    virtual bool queueTask(FluxIn &task, unsigned int taskID);

    std::vector<cht::ChunkID> taskHandles;
	std::vector<cht::ChunkID> chunkHandles;
    std::vector<FluxOut> taskResults;
    char redis_host[RedisWrapper::MAX_HOST_LENGTH];
};

#endif
