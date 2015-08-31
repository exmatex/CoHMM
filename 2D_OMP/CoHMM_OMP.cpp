#include "CoHMM_OMP.hpp"

#include <cstring>

#include <omp.h>

CoHMM_OMP::CoHMM_OMP(unsigned int dimX, unsigned int dimY, unsigned int nSteps, bool doKriging, bool doCoMD, const char * redisHost)
:
CoHMM(dimX, dimY, nSteps, doKriging, doCoMD, redisHost)
{

}

CoHMM_OMP::CoHMM_OMP(unsigned int dims[2], double delta[2], double dt[2], double gamma[3], unsigned int nSteps, bool doKriging, bool doCoMD, const char * redisHost)
:
CoHMM(dims, delta, dt, gamma, nSteps, doKriging, doCoMD, redisHost)
{

}

CoHMM_OMP::~CoHMM_OMP()
{

}

void CoHMM_OMP::processFluxes()
{
    //Get the number of tasks
    unsigned int nTasks = this->taskQueue.size();
    //Make room for the results
    this->taskResults.resize(nTasks);
    //Compute the results
    #pragma omp parallel for
    for(unsigned int i = 0; i < nTasks; i++)
    {
        this->taskResults[i] = fluxFn(this->doKriging, this->doCoMD, &this->taskQueue[i]);
    }

    //Use futures and phase to write to appropriate field
    unsigned int destField;
    if(this->curPhase == 0 || this->curPhase == 3)
    {
        //Write to fields[0]
        destField = 0;
    }
    else
    {
        //Write to fields[1]
        destField = 1;
    }
    for(unsigned int i = 0; i < this->futures.size(); i++)
    {
        //Do we need to grab the results from the taskResults?
        if(futures[i].alreadyComputed = false)
        {
            //Yes
            memcpy(this->fields[destField][i].f.f, this->taskResults[futures[i].taskID].f, sizeof(double)*7);
            memcpy(this->fields[destField][i].g.f, this->taskResults[futures[i].taskID].g, sizeof(double)*7);
        }
        else
        {
            //No, so get it from the future
            memcpy(this->fields[destField][i].f.f, this->futures[i].f.f, sizeof(double)*7);
            memcpy(this->fields[destField][i].g.f, this->futures[i].g.f, sizeof(double)*7);
        }
    }

    //Clear taskMap, taskQueue, and taskResults
    this->taskQueue.clear();
    this->taskResults.clear();
    this->taskMap.clear();
    //And done
    return;
}


bool CoHMM_OMP::queueTask(FluxIn &task, unsigned int taskID)
{
    this->taskQueue.push_back(task);
    return true;
}
