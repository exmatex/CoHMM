#include "CoHMM_CnC.hpp"

#include <cstring>



CoHMM_CnC::CoHMM_CnC(unsigned int dimX, unsigned int dimY, unsigned int nSteps, bool doKriging, bool doCoMD, const char * redisHost)
:
CoHMM(dimX, dimY, nSteps, doKriging, doCoMD, redisHost)
{
    strcpy(this->redis_host, redisHost);
}

CoHMM_CnC::CoHMM_CnC(unsigned int dims[2], double delta[2], double dt[2], double gamma[3], unsigned int nSteps, bool doKriging, bool doCoMD, const char * redisHost)
:
CoHMM(dims, delta, dt, gamma, nSteps, doKriging, doCoMD, redisHost)
{
    strcpy(this->redis_host, redisHost);
}

CoHMM_CnC::~CoHMM_CnC()
{

}

void CoHMM_CnC::processFluxes()
{
    //All tasks already pushed
    //Wait for tasks to complete
	this->ctxt.wait();

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
    Flux_Tag tag;
    // < <step, phase>, task>
    tag.first.first = this->curStep;
    tag.first.second = this->curPhase;

    for(unsigned int i = 0; i < this->futures.size(); i++)
    {
        //Do we need to grab the results from the CnC run
        if(futures[i].alreadyComputed = false)
        {
            //Yes, so get it from the context
            tag.second = futures[i].taskID;
        	Flux_Result fluxResult;
            this->ctxt.fluxResults.get(tag, fluxResult);
            //Now copy the data
            memcpy(this->fields[destField][i].f.f, fluxResult.result.f, sizeof(double)*7);
            memcpy(this->fields[destField][i].f.f, fluxResult.result.g, sizeof(double)*7);
        }
        else
        {
            //No, so get it from the future
            memcpy(this->fields[destField][i].f.f, this->futures[i].f.f, sizeof(double)*7);
            memcpy(this->fields[destField][i].g.f, this->futures[i].g.f, sizeof(double)*7);
        }
    }

    //Clean up CnC: Nothing to do

    //Clear taskMap
    this->taskMap.clear();
    //And done
    return;
}


bool CoHMM_CnC::queueTask(FluxIn &task, unsigned int taskID)
{
    //this->taskQueue.push_back(task);
    //Set up tag
    Flux_Tag tag;
    // < <step, phase>, task>
    tag.first.first = this->curStep;
    tag.first.second = this->curPhase;
    tag.second = taskID;
    //Set up item
    Flux_Item item;
    item.task = task;
    strcpy(item.redis_host, this->redis_host);
    item.doKriging = this->doKriging;
    item.doCoMD = this->doCoMD;
    //Put to item collection
    this->ctxt.fluxItems.put(tag, item);
    //Put to step collection
    this->ctxt.fluxTags.put(tag);

    return true;
}
