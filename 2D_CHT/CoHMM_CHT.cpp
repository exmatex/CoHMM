#include "CoHMM_CHT.hpp"

#include <cstring>

#include "FluxInChunk.hpp"
#include "FluxOutChunk.hpp"
#include "FluxTask.hpp"


CoHMM_CHT::CoHMM_CHT(unsigned int dimX, unsigned int dimY, unsigned int nSteps, bool doKriging, bool doCoMD, const char * redisHost)
:
CoHMM(dimX, dimY, nSteps, doKriging, doCoMD, redisHost)
{
    strcpy(this->redis_host, redisHost);
}

CoHMM_CHT::CoHMM_CHT(unsigned int dims[2], double delta[2], double dt[2], double gamma[3], unsigned int nSteps, bool doKriging, bool doCoMD, const char * redisHost)
:
CoHMM(dims, delta, dt, gamma, nSteps, doKriging, doCoMD, redisHost)
{
    strcpy(this->redis_host, redisHost);
}

CoHMM_CHT::~CoHMM_CHT()
{
    for(unsigned int i = 0; i < this->chunkHandles.size(); i++)
    {
        cht::deleteChunk(chunkHandles[i]);
    }
    for(unsigned int i = 0; i < this->taskHandles.size(); i++)
    {
		cht::deleteChunk(taskHandles[i]);
    }
}

void CoHMM_CHT::processFluxes()
{
    //For each task, block, copy to result buffer, and free
    this->taskResults.resize(this->taskHandles.size());
    for(unsigned int i = 0; i < this->taskHandles.size(); i++)
	{
		cht::shared_ptr<FluxOutChunk const> resChunk;
		cht::getChunk(this->taskHandles[i], resChunk);
        this->taskResults[i] = resChunk->output;
		cht::deleteChunk(chunkHandles[i]);
		cht::deleteChunk(taskHandles[i]);
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
    this->taskHandles.clear();
    this->taskResults.clear();
    this->chunkHandles.clear();
    //And done
    return;
}


bool CoHMM_CHT::queueTask(FluxIn &task, unsigned int taskID)
{
    //this->taskQueue.push_back(task);
    //Make FluxInObject
    FluxInChunk inObject;
    inObject.input = task;
    inObject.taskID = taskID;
    inObject.doKriging = this->doKriging;
    inObject.doCoMD = this->doCoMD;
    strcpy(inObject.redisHost, this->redis_host);

    //Make dummy IDs: This may be less memory efficient than a constant reservation
    cht::ChunkID dummyId;
    this->chunkHandles.push_back(dummyId);
    this->taskHandles.push_back(dummyId);
    //Push chunk to CHT
    this->chunkHandles.back() = cht::registerChunk(new FluxInChunk(inObject));
    //Push task to CHT
    this->taskHandles.back() = cht::executeMotherTask<FluxTask>(chunkHandles.back());

    return true;
}
