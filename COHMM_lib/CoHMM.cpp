#include <cstring>
#include <iostream>

#include "RedisWrapper.hpp"

#include "CoHMM.hpp"
#include "2DKriging.hpp"
#include "kriging.hpp"
#include "redisBuckets.hpp"

CoHMM::CoHMM(unsigned int dimX, unsigned int dimY, unsigned int nSteps, bool doKriging, bool doCoMD, const char * redisHost)
{
    this->dims[0] = dimX;
    this->dims[1] = dimY;
    this->nSteps = nSteps;
    this->doKriging = doKriging;
    this->doCoMD = doCoMD;
    this->delta[0] = 1.0;
    this->delta[1] = 1.0;
    this->dt[0] = 0.1;
    this->dt[1] = 0.1;
    this->gamma[0] = 0;
    this->gamma[1] = this->gamma[0];
    this->gamma[2] = 0.1 * this->gamma[1];

    this->futures.resize(this->dims[0]*this->dims[1]);
    this->fields[0].resize(this->dims[0]*this->dims[1]);
    this->fields[1].resize(this->dims[0]*this->dims[1]);


    this->curStep = 0;
    this->curPhase = 0;

    //Init redis
    RedisWrapper::getContext(redisHost);
}

CoHMM::CoHMM(unsigned int dims[2], double delta[2], double dt[2], double gamma[3], unsigned int nSteps, bool doKriging, bool doCoMD, const char * redisHost)
{
    this->dims[0] = dims[0];
    this->dims[1] = dims[1];
    this->nSteps = nSteps;
    this->doKriging = doKriging;
    this->doCoMD = doCoMD;
    this->delta[0] = delta[0];
    this->delta[1] = delta[1];
    this->dt[0] = dt[0];
    this->dt[1] = dt[1];
    this->gamma[0] = gamma[0];
    this->gamma[1] = gamma[1];
    this->gamma[2] = gamma[2];

    this->futures.resize(this->dims[0]*this->dims[1]);
    this->fields[0].resize(this->dims[0]*this->dims[1]);
    this->fields[1].resize(this->dims[0]*this->dims[1]);

    this->curStep = 0;
    this->curPhase = 0;

    //Connect to redis
    RedisWrapper::getContext(redisHost);
}

CoHMM::CoHMM()
{

}

CoHMM::~CoHMM()
{

}

void CoHMM::lateInit(unsigned int dimX, unsigned int dimY, unsigned int nSteps, bool doKriging, bool doCoMD, const char * redisHost)
{
    this->dims[0] = dimX;
    this->dims[1] = dimY;
    this->nSteps = nSteps;
    this->doKriging = doKriging;
    this->doCoMD = doCoMD;
    this->delta[0] = 1.0;
    this->delta[1] = 1.0;
    this->dt[0] = 0.1;
    this->dt[1] = 0.1;
    this->gamma[0] = 0;
    this->gamma[1] = this->gamma[0];
    this->gamma[2] = 0.1 * this->gamma[1];

    this->futures.resize(this->dims[0]*this->dims[1]);
    this->fields[0].resize(this->dims[0]*this->dims[1]);
    this->fields[1].resize(this->dims[0]*this->dims[1]);

    this->curStep = 0;
    this->curPhase = 0;

    //Init redis
    RedisWrapper::getContext(redisHost);
}

void CoHMM::lateInit(unsigned int dims[2], double delta[2], double dt[2], double gamma[3], unsigned int nSteps, bool doKriging, bool doCoMD, const char * redisHost)
{
    this->dims[0] = dims[0];
    this->dims[1] = dims[1];
    this->nSteps = nSteps;
    this->doKriging = doKriging;
    this->doCoMD = doCoMD;
    this->delta[0] = delta[0];
    this->delta[1] = delta[1];
    this->dt[0] = dt[0];
    this->dt[1] = dt[1];
    this->gamma[0] = gamma[0];
    this->gamma[1] = gamma[1];
    this->gamma[2] = gamma[2];

    this->futures.resize(this->dims[0]*this->dims[1]);
    this->fields[0].resize(this->dims[0]*this->dims[1]);
    this->fields[1].resize(this->dims[0]*this->dims[1]);

    this->curStep = 0;
    this->curPhase = 0;

    //Connect to redis
    RedisWrapper::getContext(redisHost);
}

void CoHMM::initializeConservedFields(InitialConditions_e initCondition)
{
    init_conserved_fields(this->fields[0], this->dims, this->dims[0]*this->dims[1], initCondition);
}

void CoHMM::outputVTK()
{
    //Build lattice for sanity's sake
	Lattice l;
	l.dim_x = this->dims[0];
	l.dim_y = this->dims[1];
	l.dx = this->delta[0];
	l.dy = this->delta[1];
	l.dt_x = this->dt[0];
	l.dt_y = this->dt[1];
	//VTK output it
	printf_fields_vtk(this->curStep, this->fields[0], l, this->dims[0]*this->dims[1]);
}

void CoHMM::prepFirstFlux()
{
    //Phase 0
    //Generate flux tasks of conserved fields in 0
    int numTasks =  prepTasks(this->fields[0]);
    return;
}

void CoHMM::prepSecondFlux()
{
    //Phase 1
    this->curPhase++;
    //Generates phase 1's w
    wNSqrt(this->fields[0], this->fields[1], dims, dt, delta);
    //and the phase 1 tasks that are associated with them
    int numTasks =  prepTasks(this->fields[1]);
    return;
}

void CoHMM::prepThirdFlux()
{
    //Phase 2
    this->curPhase++;
    //Now do the jiang tambor stuff
    wSummation(this->fields[0], this->fields[1], dims, dt, delta);
    //Now we ge to start the next half-step
    int numTasks =  prepTasks(this->fields[1]);
    return;
}

void CoHMM::prepLastFlux()
{
    //Phase 3
    this->curPhase++;
    //Generates phase 1's w
    wNSqrt(this->fields[1], this->fields[0], dims, dt, delta);
    //and the phase 1 tasks that are associated with them
    int numTasks =  prepTasks(this->fields[0]);
    return;
}

void CoHMM::finishStep()
{
    //Phase 4 (and 0)
    this->curPhase++;
    //Now do the jiang tambor stuff
    wSummation(this->fields[1], this->fields[0], dims, dt, delta);
    //Now do a shift_back (in place)
    inPlace_shift_back(this->fields[0].data(), this->dims);
    //Reset phase and increment step
    this->curPhase = 0;
    this->curStep++;
    return;
}

int CoHMM::prepTasks(std::vector<Node> &field)
{
	//Task counter
	unsigned int taskCounter = 0;

	//Iterate over points to see who is getting krig'd and who is getting comd'd
	for(int y = 0; y < dims[1]; y++)
	{
		for(int x = 0; x < dims[0]; x++)
		{
			//Grab comd database
			std::vector<double *> wVec;
			std::vector<double *> fVec;
			std::vector<double *> gVec;
			getSortedSubBucketNearZero(field[x + dims[0]*y].w.w, (char *)"comd", comdDigits, 2, &wVec, &fVec, &gVec, zeroThresh);

			//Check for exact value
			bool useDB = ifConservedFieldsMatch(field[x+dims[0]*y].w.w, &wVec, dbT);
			//We actually found it
			if(useDB == true)
			{
				//Write the result to the appropriate future
				memcpy(&futures[x+dims[0]*y].f, fVec[0], sizeof(double)*7);
				memcpy(&futures[x+dims[0]*y].g, gVec[0], sizeof(double)*7);
				//Signal that we don't have to do a look up later
				futures[x+dims[0]*y].alreadyComputed = true;
			}
			else
			{
				//We did not
				//Check gradient
				bool smallGradient;
				if(doKriging == true)
				{
					smallGradient = checkGradient(x, y, field, dims, delta);
				}
				else
				{
					smallGradient = false;
				}
				//Can we krig?
				if(smallGradient == true)
				{
					//Grab krig database
					std::vector<double *> wVecK;
					std::vector<double *> fVecK;
					std::vector<double *> gVecK;
					getSortedSubBucketNearZero(field[x + dims[0]*y].w.w, (char *)"krig", krigDigits, 1, &wVecK, &fVecK, &gVecK, zeroThresh);
					bool useDB = ifConservedFieldsMatch(field[x+dims[0]*y].w.w, &wVecK, 0.0);
					//Did we already krig this?
					if(useDB == true)
					{
						//Write the result to the appropriate future
						memcpy(&futures[x+dims[0]*y].f, fVecK[0], sizeof(double)*7);
						memcpy(&futures[x+dims[0]*y].g, gVecK[0], sizeof(double)*7);
						//Signal that we don't have to do a look up later
						futures[x+dims[0]*y].alreadyComputed = true;
					}
					else
					{
						//Do we have at least two w's?
						if (wVec.size() >= 2)
						{
							//We do, so kriging task
							//Is there already a task?
							if(taskMap.find(field[x+dims[0]*y].w) == taskMap.end())
							{
								//Not previously added to map
								//So add it
								taskMap[field[x+dims[0]*y].w] = taskCounter;
								//Now build a FluxIn
								FluxIn actualTask;
								memcpy(&actualTask.fields, &field[x+dims[0]*y].w, sizeof(Conserved));
								actualTask.tryKriging = true;
								//Now enqueue it
                                this->queueTask(actualTask, taskCounter);
								//Increment task counter
								taskCounter++;
							}
							//Set future to task ID and indicate we need to do a look up
							futures[x+dims[0]*y].alreadyComputed = false;
							futures[x+dims[0]*y].taskID = taskMap[field[x+dims[0]*y].w];
						}
						//Nope, so CoMD
						else
						{
							//Is there already a task?
							if(taskMap.find(field[x+dims[0]*y].w) == taskMap.end())
							{
								//Not previously added to map
								//So add it
								taskMap[field[x+dims[0]*y].w] = taskCounter;
								//Now build a FluxIn
								FluxIn actualTask;
								memcpy(&actualTask.fields, &field[x+dims[0]*y].w, sizeof(Conserved));
								actualTask.tryKriging = false;
								//Now enqueue it
                                this->queueTask(actualTask, taskCounter);
								//Increment task counter
								taskCounter++;
							}
							//Set future to task ID and indicate we need to do a look up
							futures[x+dims[0]*y].alreadyComputed = false;
							futures[x+dims[0]*y].taskID = taskMap[field[x+dims[0]*y].w];
						}
					}
					freeClear(wVecK);
					freeClear(fVecK);
					freeClear(gVecK);
				}
				//No krig, so comd
				else
				{
					//Is there already a task?
					if(taskMap.find(field[x+dims[0]*y].w) == taskMap.end())
					{
						//Not previously added to map
						//So add it
						taskMap[field[x+dims[0]*y].w] = taskCounter;
						//Now build a FluxIn
						FluxIn actualTask;
						memcpy(&actualTask.fields, &field[x+dims[0]*y].w, sizeof(Conserved));
						actualTask.tryKriging = false;
						//Now enqueue it
                        this->queueTask(actualTask, taskCounter);
						//Increment task counter
						taskCounter++;
					}
					//Set future to task ID and indicate we need to do a look up
					futures[x+dims[0]*y].alreadyComputed = false;
					futures[x+dims[0]*y].taskID = taskMap[field[x+dims[0]*y].w];
				}
			}
			freeClear(wVec);
			freeClear(fVec);
			freeClear(gVec);
		}
	}
	return taskCounter;
}

int CoHMM::getNumberOfTasks()
{
    return this->taskMap.size();
}
