#include "CoHMM_CharmPP.hpp"

#include <cstring>

#include <omp.h>

#include "RedisWrapper.hpp"

//ReadOnly vars
char gRedis_host[48];
bool gDoKriging;
bool gDoCoMD;
CProxy_CoHMM_CharmPP gMainProxy;

CoHMM_CharmPP::CoHMM_CharmPP(unsigned int dimX, unsigned int dimY, unsigned int nSteps, bool doKriging, bool doCoMD, const char * redisHost)
:
CoHMM(dimX, dimY, nSteps, doKriging, doCoMD, redisHost)
{
    this->results.resize(dimX*dimY);
    this->doneTasks = 0;
    this->expectedTasks = -1;
}

CoHMM_CharmPP::CoHMM_CharmPP(unsigned int dims[2], double delta[2], double dt[2], double gamma[3], unsigned int nSteps, bool doKriging, bool doCoMD, const char * redisHost)
:
CoHMM(dims, delta, dt, gamma, nSteps, doKriging, doCoMD, redisHost)
{
    this->results.resize(dims[0]*dims[1]);
    this->doneTasks = 0;
    this->expectedTasks = -1;
}

CoHMM_CharmPP::CoHMM_CharmPP()
:
CoHMM()
{
    //Do nothing here because a separeate method does what we want
}

//Actually the main
CoHMM_CharmPP::CoHMM_CharmPP(CkArgMsg * m)
{
    if( m->argc != 5)
    {
        ckerr <<  "./2D_DaDTest <dim_x> <dim_y> <nsteps> <redis_server>" << endl;
        CkExit();
    }
    //Init stuff
    //Set up parameters
	const bool fineGrainFT = false;
	const bool doKriging = true;
    gDoKriging = doKriging;
	const bool doCoMD = false;
    gDoCoMD = doCoMD;
	unsigned int dims[2] = {(unsigned int) atoi(m->argv[1]), (unsigned int) atoi(m->argv[2])};
	double dt[2] = {0.1, 0.1};
	double delta[2] = {1.0, 1.0};
	double gamma[3];
	gamma[0] = 0; //mom_gamma
	gamma[1] = gamma[0]; //strain_gamma
	gamma[2] = 0.1 * gamma[1];//en_gamma
    strcpy(gRedis_host, m->argv[4]);

	unsigned int numSteps = atoi(m->argv[3]);
    //Call lateInit
    this->lateInit(dims, delta, dt, gamma, numSteps, doKriging, doCoMD, m->argv[4]);
    //Actually Init
    {
        double startTime = CkWallTimer();
        FILE * stFile = fopen("stime.dat", "w");
        fprintf(stFile, "%f", startTime);
        fclose(stFile);
    }
	ckout << "Initializing " << dims[0] << " by " << dims[1] << " grid" << endl;
    this->initializeConservedFields(InitialConditions_e::CENTRALIZED);
	ckout << "Initialized" << endl;
	ckout << "Running for " << numSteps << " iterations" << endl;
    //Start Driver
    this->cohmmDriver();
}

CoHMM_CharmPP::~CoHMM_CharmPP()
{
    //Nothing, everything is handled through STL
}

void CoHMM_CharmPP::cohmmDriver(unsigned int phaseOffset)
{
    //Are we done?
    if(this->curStep == this->nSteps)
    {
        //One last print
		ckout << this->nSteps << ": Vising to Verifying" << endl;
        this->outputVTK();
        //We are done
        ckout << "Ran for " << nSteps << " iterations" << endl;
        {
            double endTime = CkWallTimer();
            FILE * eFile = fopen("etime.dat", "w");
            fprintf(eFile, "%f", endTime);
            fclose(eFile);
        }
        CkExit();
    }
    //Nope, so check what phase we are in
    else
    {
        //At this point, any fluxes have been processed and we are ready to proceed
        unsigned int nTasks;
        unsigned int switchPhase = this->curPhase + phaseOffset;
        switch(switchPhase)
        {
            case 0:
                ckout << this->curStep << ": Vising to Verifying" << endl;
                this->outputVTK();
                ckout << this->curStep <<  ": First Flux" << endl;
                this->prepFirstFlux();
                nTasks = this->getNumberOfTasks();
                ckout << this->curStep <<  ": Doing " << nTasks << " fluxes" << endl;
                this->processFluxes();
            break;
            case 1:
                ckout << this->curStep <<  ": Second Flux" << endl;
                this->prepSecondFlux();
                nTasks = this->getNumberOfTasks();
                ckout << this->curStep <<  ": Doing " << nTasks << " fluxes" << endl;
                this->processFluxes();
            break;
            case 2:
                ckout << this->curStep <<  ": Third Flux" << endl;
                this->prepThirdFlux();
                nTasks = this->getNumberOfTasks();
                ckout << this->curStep <<  ": Doing " << nTasks << " fluxes" << endl;
                this->processFluxes();
            break;
            case 3:
                ckout << this->curStep <<  ": Last Flux" << endl;
                this->prepLastFlux();
                nTasks = this->getNumberOfTasks();
                ckout << this->curStep <<  ": Doing " << nTasks << " fluxes" << endl;
                this->processFluxes();
            break;
            case 4:
                ckout << this->curStep <<  ": Finish Step, no Fluxes" << endl;
                this->finishStep();
                //Recursivevly call again: finishStep handled phase and step
                this->cohmmDriver();
            break;
        }
    }
}

void CoHMM_CharmPP::processFluxes()
{
    //Get the number of tasks
    unsigned int nTasks = this->getNumberOfTasks();
    //Don't need to queue or do anything as futures are already done
    //So just set the number of expected tasks
    gMainProxy.checkAndIncOrSet(nTasks);
    //And that is it. We will recover control later
    return;
}

void CoHMM_CharmPP::processFutures()
{
    //Based on phase, write to a specific destination
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
    //Now just get the results
    for(unsigned int i = 0; i < this->futures.size(); i++)
    {
        //Do we need to grab the results from the Charm++ run
        if(futures[i].alreadyComputed = false)
        {
            //Yes, so get it from results
            //Now copy the data
            memcpy(this->fields[destField][i].f.f, this->results[futures[i].taskID].f, sizeof(double)*7);
            memcpy(this->fields[destField][i].f.f, this->results[futures[i].taskID].g, sizeof(double)*7);
        }
        else
        {
            //No, so get it from the future
            memcpy(this->fields[destField][i].f.f, this->futures[i].f.f, sizeof(double)*7);
            memcpy(this->fields[destField][i].g.f, this->futures[i].g.f, sizeof(double)*7);
        }
    }

    //Clear taskMap, but not futures or results as we will need them later
    this->taskMap.clear();
    //Now return to the solver with a phaseOffst of 1 so we can continue without breaking everything else
    this->cohmmDriver(1);
    return;
}


bool CoHMM_CharmPP::queueTask(FluxIn &task, unsigned int taskID)
{
    //Resize if needed: Try to avoid doing this because it can mess things up horribly
    ///WARNING: May need to serialize/atomic
    if(taskID >= this->results.size())
    {
        this->results.resize(this->results.size() + CoHMM_CharmPP::TASK_QUEUE_SIZE);
    }
    //Actually queue task
    CProxy_Flux::ckNew(task, taskID);
    return true;
}


void CoHMM_CharmPP::fluxCallBack(FluxOut result, int tid)
{
    //Set result value: CAN BE VERY BAD DEPENDING ON HOW REALLOCATION WORKS!!!
    this->results[tid] = result;

    //Daisy chain method to checkAndIncOrSet with -1 to increment
    gMainProxy.checkAndIncOrSet(-1);
}

void CoHMM_CharmPP::actualCheckAndIncOrSet(int potSetVal)
{
    //Charm++ allegedly guarantees this will all be atomic if we call via checkAndIncOrSet
    //Did we pass in a non-negative?
    if(potSetVal >= 0)
    {
        //Yup, so this is the number of tasks
        this->expectedTasks = potSetVal;
    }
    else
    {
        //Nope, so a task just finished: Increment
        this->doneTasks++;
    }

    //See if we are done
    if(this->expectedTasks == this->doneTasks)
    {
       //We are, so reset the variables
       this->doneTasks = 0;
       this->expectedTasks = -1;
       //Then fire the callBack that returns control to the solver
       this->processFutures();

    }
    //We aren't, so wait for the next call
}

Flux::Flux(FluxIn task, int tid)
{
    //Build redis
    RedisWrapper::getContext(gRedis_host);
    //Compute flux
    FluxOut res= fluxFn(gDoKriging, gDoCoMD, &task);
    //Do callback to write result
    gMainProxy.fluxCallBack(res, tid);
}
