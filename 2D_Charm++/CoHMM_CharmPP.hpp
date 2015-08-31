#ifndef COHMM_CHARMPP_HPP
#define COHMM_CHARMPP_HPP

#include "2DKriging.hpp"
#include "CoHMM.hpp"
#include "FluxPupWrappers.hpp"

#include "CoHMM_TBA_Charm.decl.h"

//ReadOnly vars
extern char gRedis_host[48];
extern bool gDoKriging;
extern bool gDoCoMD;
extern CProxy_CoHMM_CharmPP gMainProxy;

class CoHMM_CharmPP : public CoHMM, public CBase_CoHMM_CharmPP
{
    CoHMM_CharmPP_SDAG_CODE
public:
    CoHMM_CharmPP();
    CoHMM_CharmPP(CkArgMsg * m);
    CoHMM_CharmPP(unsigned int dimX, unsigned int dimY, unsigned int nSteps, bool doKriging, bool doCoMD, const char * redisHost = "localhost");
    CoHMM_CharmPP(unsigned int dims[2], double delta[2], double dt[2], double gamma[3], unsigned int nSteps, bool doKriging, bool doCoMD, const char * redisHost = "localhost");
    virtual ~CoHMM_CharmPP();

    virtual void processFluxes();

    void fluxCallBack(FluxOut result, int tid);

    void actualCheckAndIncOrSet(int potSetVal);

protected:
    virtual bool queueTask(FluxIn &task, unsigned int taskID);

    void cohmmDriver(unsigned int phaseOffset = 0);

    void processFutures();

    std::vector<FluxOut> results;
    unsigned int doneTasks;
    int expectedTasks;

    static const size_t TASK_QUEUE_SIZE = 40;
};

///TODO: Separate hpp?
class Flux: public CBase_Flux
{
public:
    Flux(FluxIn task, int tid);
};

#include "CoHMM_TBA_Charm.def.h"

#endif
