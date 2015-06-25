#include <map>
#include <vector>
#include <string>
#include <cassert>
#include <cstring>
#include <fstream>

#include <hiredis.h>

///TODO: Find a more OS neutral way to do this
#include <unistd.h>

extern "C"
{
#include <CoMD_lib.h>
}

#include "redisShmem.hpp"
#include "2DKriging.hpp"
#include "redisBuckets.hpp"
#include "kriging.hpp"
#include "CoHMM_DaD.hpp"


bool backToTheFuture(Node * fields, FluxFuture * futures, int * dims, int curStep, int curPhase, redisContext * headRedis)
{
	std::map<unsigned int, FluxOut> retMap;
	for(int i = 0; i < dims[0]*dims[1]; i++)
	{
		//Did we compute the results
		if(futures[i].alreadyComputed == false)
		{
			//We did
			unsigned int taskID = futures[i].taskID;
			//Check if we already fetched this output
			if(retMap.find(taskID) == retMap.end())
			{
				//We did not, so fetch it
				FluxOut res;
				getSingle<FluxOut>(&res, curStep, curPhase, taskID, headRedis, "RESULT");
				///TODO: Verify this is deep copy, but it should be
				retMap[taskID] = res;
			}
			//Copy result to fields
			memcpy(fields[i].f.f, retMap[taskID].f, sizeof(double)*7);
			memcpy(fields[i].g.f, retMap[taskID].g, sizeof(double)*7);
		}
		else
		{
			//We did not, so it is already in the future
			memcpy(fields[i].f.f, futures[i].f.f, sizeof(double)*7);
			memcpy(fields[i].g.f, futures[i].g.f, sizeof(double)*7);
		}
	}
	return true;
}

int prepTasks(Node * fields, FluxFuture * futures, bool doKriging, int * dims, double * dt, double * delta,  int curStep, int curPhase, redisContext * headRedis)
{
	//Task map
	std::map<Conserved, unsigned int> taskMap;
	unsigned int taskCounter = 0;

	//Iterate over points to see who is getting krig'd and who is getting comd'd
	for(int y = 0; y < dims[1]; y++)
	{
		for(int x = 0; x < dims[0]; x++)
		{
			///TODO: Refactor to avoid repeated code and to use individual methods for swift/t
			//Grab comd database
			std::vector<double *> wVec;
			std::vector<double *> fVec;
			std::vector<double *> gVec;
			getSortedSubBucketNearZero(fields[x + dims[0]*y].w.w, (char *)"comd", headRedis, comdDigits, 2, &wVec, &fVec, &gVec, zeroThresh);

			//Check for exact value
			bool useDB = ifConservedFieldsMatch(fields[x+dims[0]*y].w.w, &wVec, dbT);
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
					smallGradient = checkGradient(x, y, fields, dims, delta);
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
					getSortedSubBucketNearZero(fields[x + dims[0]*y].w.w, (char *)"krig", headRedis, krigDigits, 1, &wVecK, &fVecK, &gVecK, zeroThresh);
					bool useDB = ifConservedFieldsMatch(fields[x+dims[0]*y].w.w, &wVecK, 0.0);
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
							if(taskMap.find(fields[x+dims[0]*y].w) == taskMap.end())
							{
								//Not previously added to map
								//So add it
								taskMap[fields[x+dims[0]*y].w] = taskCounter;
								//Now build a FluxIn
								FluxIn actualTask;
								memcpy(&actualTask.fields, &fields[x+dims[0]*y].w, sizeof(Conserved));
								actualTask.tryKriging = true;
								//Now enqueue it
								putSingle<FluxIn>(&actualTask, curStep, curPhase, taskCounter, headRedis, "TASK");
								//Increment task counter
								taskCounter++;
							}
							//Set future to task ID and indicate we need to do a look up
							futures[x+dims[0]*y].alreadyComputed = false;
							futures[x+dims[0]*y].taskID = taskMap[fields[x+dims[0]*y].w];
						}
						//Nope, so CoMD
						else
						{
							//Is there already a task?
							if(taskMap.find(fields[x+dims[0]*y].w) == taskMap.end())
							{
								//Not previously added to map
								//So add it
								taskMap[fields[x+dims[0]*y].w] = taskCounter;
								//Now build a FluxIn
								FluxIn actualTask;
								memcpy(&actualTask.fields, &fields[x+dims[0]*y].w, sizeof(Conserved));
								actualTask.tryKriging = false;
								//Now enqueue it
								putSingle<FluxIn>(&actualTask, curStep, curPhase, taskCounter, headRedis, "TASK");
								//Increment task counter
								taskCounter++;
							}
							//Set future to task ID and indicate we need to do a look up
							futures[x+dims[0]*y].alreadyComputed = false;
							futures[x+dims[0]*y].taskID = taskMap[fields[x+dims[0]*y].w];
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
					if(taskMap.find(fields[x+dims[0]*y].w) == taskMap.end())
					{
						//Not previously added to map
						//So add it
						taskMap[fields[x+dims[0]*y].w] = taskCounter;
						//Now build a FluxIn
						FluxIn actualTask;
						memcpy(&actualTask.fields, &fields[x+dims[0]*y].w, sizeof(Conserved));
						actualTask.tryKriging = false;
						//Now enqueue it
						putSingle<FluxIn>(&actualTask, curStep, curPhase, taskCounter, headRedis, "TASK");
						//Increment task counter
						taskCounter++;
					}
					//Set future to task ID and indicate we need to do a look up
					futures[x+dims[0]*y].alreadyComputed = false;
					futures[x+dims[0]*y].taskID = taskMap[fields[x+dims[0]*y].w];
				}
			}
			freeClear(wVec);
			freeClear(fVec);
			freeClear(gVec);
		}
	}
	return taskCounter;
}

//Initialize all fields and store in database at KEY_0_0
bool initEverything(bool doKriging, bool doCoMD, int * dims, double * dt, double * delta, double * gamma, const char * redis_host)
{
	const char * redisHostName;
	//Check redis host
	bool needFree = checkRedisHost(redis_host);
	if(needFree == true)
	{
		redisHostName = getRedisHost(redis_host);
	}
	else
	{
		redisHostName = redis_host;
	}
	//Connect to redis
	redisContext * headRedis = redisConnect(redisHostName, 6379);
	if(headRedis == NULL || headRedis->err)
	{
		printf("Redis error: %s\n", headRedis->errstr);
		return false;
	}
	//Allocate fields buffer
	Node * field = new Node[dims[0]*dims[1]]();
	//Initialize fields
	init_conserved_fields(field, dims, dims[0]*dims[1]);
	//Put to DB
	putBlocks<Node>(field, dims[0], dims[1], 0, 0, headRedis, "FIELD");
	//cleanup redis
	redisFree(headRedis);
	//Free memory
	delete [] field;
	if(needFree == true)
	{
		delete [] redisHostName;
	}
	//And we're done
	return true;
}

//Essentially run through to the first flux of the first half-step
int prepFirstFlux(bool doKriging, bool doCoMD, int * dims, double * dt, double * delta, double * gamma, int curStep, const char * redis_host)
{
	const char * redisHostName;
	//Check redis host
	bool needFree = checkRedisHost(redis_host);
	if(needFree == true)
	{
		redisHostName = getRedisHost(redis_host);
	}
	else
	{
		redisHostName = redis_host;
	}
	//Connect to redis
	redisContext * headRedis = redisConnect(redisHostName, 6379);
	if(headRedis == NULL || headRedis->err)
	{
		printf("Redis error: %s\n", headRedis->errstr);
		return false;
	}
	//Allocate fields buffer
	Node * field = new Node[dims[0]*dims[1]]();
	//Get field data from previous step
	getBlocks<Node>(field, dims[0], dims[1], curStep, 0, headRedis, "FIELD");
	//Prep futures
	FluxFuture * futures = new FluxFuture[dims[0]*dims[1]]();
	int numTasks =  prepTasks(field, futures, doKriging , dims, dt, delta, curStep, 0, headRedis);
	//Write futures
	putBlocks<FluxFuture>(futures, dims[0], dims[1], curStep, 0, headRedis, "FUTS");
	//cleanup redis
	redisFree(headRedis);
	//Don't save fields, we didn't do anything with it
	delete [] field;
	delete [] futures;
	if(needFree == true)
	{
		delete [] redisHostName;
	}
	//Return the number of tasks
	return numTasks;
}

int prepSecondFlux(bool doKriging, bool doCoMD, int * dims, double * dt, double * delta, double * gamma, int curStep, const char * redis_host)
{
	const char * redisHostName;
	//Check redis host
	bool needFree = checkRedisHost(redis_host);
	if(needFree == true)
	{
		redisHostName = getRedisHost(redis_host);
	}
	else
	{
		redisHostName = redis_host;
	}
	//Connect to redis
	redisContext * headRedis = redisConnect(redisHostName, 6379);
	if(headRedis == NULL || headRedis->err)
	{
		printf("Redis error: %s\n", headRedis->errstr);
		return false;
	}
	//Dependent on phase 0's w, and the output of phase 0's fluxes
	//Allocate fields buffer
	Node * field = new Node[dims[0]*dims[1]]();
	//Get field data from previous step
	getBlocks<Node>(field, dims[0], dims[1], curStep, 0, headRedis, "FIELD");
	//Get futures from previous step
	FluxFuture * futures = new FluxFuture[dims[0]*dims[1]]();
	getBlocks<FluxFuture>(futures, dims[0], dims[1], curStep, 0, headRedis, "FUTS");
	backToTheFuture(field, futures, dims, curStep, 0, headRedis);
	//Generates phase 1's w
	wNSqrt(field, dims, dt, delta);
	//and the phase 1 tasks that are associated with them
	int numTasks =  prepTasks(field, futures, doKriging , dims, dt, delta, curStep, 1, headRedis);
	//Write futures for phase 1
	putBlocks<FluxFuture>(futures, dims[0], dims[1], curStep, 1, headRedis, "FUTS");
	//cleanup redis
	redisFree(headRedis);
	//Don't save fields, we only need the original w's and the final f's and g's
	delete [] field;
	delete [] futures;
	if(needFree == true)
	{
		delete [] redisHostName;
	}
	//Return the number of tasks
	return numTasks;
}

int prepThirdFlux(bool doKriging, bool doCoMD, int * dims, double * dt, double * delta, double * gamma, int curStep, const char * redis_host)
{
	//Need phase 0's w's and phase 1's f's and g's
	const char * redisHostName;
	//Check redis host
	bool needFree = checkRedisHost(redis_host);
	if(needFree == true)
	{
		redisHostName = getRedisHost(redis_host);
	}
	else
	{
		redisHostName = redis_host;
	}
	//Connect to redis
	redisContext * headRedis = redisConnect(redisHostName, 6379);
	if(headRedis == NULL || headRedis->err)
	{
		printf("Redis error: %s\n", headRedis->errstr);
		return false;
	}
	//Dependent on phase 0's w, and the output of phase 0's fluxes
	//Allocate fields buffer
	Node * aField = new Node[dims[0]*dims[1]]();
	Node * bField = new Node[dims[0]*dims[1]]();
	//Get field data from phase 0
	getBlocks<Node>(aField, dims[0], dims[1], curStep, 0, headRedis, "FIELD");
	//Get futures from previous phase
	FluxFuture * futures = new FluxFuture[dims[0]*dims[1]]();
	getBlocks<FluxFuture>(futures, dims[0], dims[1], curStep, 1, headRedis, "FUTS");
	backToTheFuture(bField, futures, dims, curStep, 1, headRedis);
	//Now do the jiang tambor stuff
	wSummation(aField, bField, dims, dt, delta);
	//Now we ge to start the next half-step
	int numTasks =  prepTasks(bField, futures, doKriging , dims, dt, delta, curStep, 2, headRedis);
	//Write futures
	putBlocks<FluxFuture>(futures, dims[0], dims[1], curStep, 2, headRedis, "FUTS");
	//Write field for later use
	putBlocks<Node>(bField, dims[0], dims[1], curStep, 2, headRedis, "FIELD");
	//cleanup redis
	redisFree(headRedis);
	//Free Memory
	delete[] aField;
	delete[] bField;
	delete[] futures;
	if(needFree == true)
	{
		delete [] redisHostName;
	}
	//Return the number of tasks
	return numTasks;
}

int prepLastFlux(bool doKriging, bool doCoMD, int * dims, double * dt, double * delta, double * gamma, int curStep, const char * redis_host)
{
	//Essentially the same as secondFlux
	//Really should refactor
	//Connect to redis
	const char * redisHostName;
	//Check redis host
	bool needFree = checkRedisHost(redis_host);
	if(needFree == true)
	{
		redisHostName = getRedisHost(redis_host);
	}
	else
	{
		redisHostName = redis_host;
	}
	redisContext * headRedis = redisConnect(redisHostName, 6379);
	if(headRedis == NULL || headRedis->err)
	{
		printf("Redis error: %s\n", headRedis->errstr);
		return false;
	}
	//Dependent on phase 2's w, and the output of phase 2's fluxes
	//Allocate fields buffer
	Node * field = new Node[dims[0]*dims[1]]();
	//Get field data from previous step
	getBlocks<Node>(field, dims[0], dims[1], curStep, 2, headRedis, "FIELD");
	//Get futures from previous step
	FluxFuture * futures = new FluxFuture[dims[0]*dims[1]]();
	getBlocks<FluxFuture>(futures, dims[0], dims[1], curStep, 2, headRedis, "FUTS");
	backToTheFuture(field, futures, dims, curStep, 2, headRedis);
	//Generates phase 3's w
	wNSqrt(field, dims, dt, delta);
	//and the phase 3 tasks that are associated with them
	int numTasks =  prepTasks(field, futures, doKriging , dims, dt, delta, curStep, 3, headRedis);
	//Write futures for phase 3
	putBlocks<FluxFuture>(futures, dims[0], dims[1], curStep, 3, headRedis, "FUTS");
	//cleanup redis
	redisFree(headRedis);
	//Don't save fields, we only need the original w's and the final f's and g's
	delete [] field;
	delete [] futures;
	if(needFree == true)
	{
		delete [] redisHostName;
	}
	//Return the number of tasks
	return numTasks;
}

int finishStep(bool doKriging, bool doCoMD, int * dims, double * dt, double * delta, double * gamma, int curStep, const char * redis_host)
{
	//Basically phase 3, but replace the last flux call with a shiftback
	//Need phase 2's w's and phase 3's f's and g's
	//Connect to redis
	const char * redisHostName;
	//Check redis host
	bool needFree = checkRedisHost(redis_host);
	if(needFree == true)
	{
		redisHostName = getRedisHost(redis_host);
	}
	else
	{
		redisHostName = redis_host;
	}
	redisContext * headRedis = redisConnect(redisHostName, 6379);
	if(headRedis == NULL || headRedis->err)
	{
		printf("Redis error: %s\n", headRedis->errstr);
		return false;
	}
	//Dependent on phase 2's w, and the output of phase 3's fluxes
	//Allocate fields buffer
	Node * fieldA = new Node[dims[0]*dims[1]]();
	Node * fieldB = new Node[dims[0]*dims[1]]();
	//Get field data from phase 2
	getBlocks<Node>(fieldA, dims[0], dims[1], curStep, 2, headRedis, "FIELD");
	//Get futures from previous phase
	FluxFuture * futures = new FluxFuture[dims[0]*dims[1]]();
	getBlocks<FluxFuture>(futures, dims[0], dims[1], curStep, 3, headRedis, "FUTS");
	backToTheFuture(fieldB, futures, dims, curStep, 3, headRedis);
	//Now do the jiang tambor stuff
	wSummation(fieldA, fieldB, dims, dt, delta);
	//Now we do a shift back
	//COPY RIGHT TO LEFT
	shift_back(fieldA, dims[0]*dims[1], dims, fieldB);
	//Write field for later use
	putBlocks<Node>(fieldA, dims[0], dims[1], curStep+1, 0, headRedis, "FIELD");
	//cleanup redis
	redisFree(headRedis);
	//Free Memory
	delete [] fieldA;
	delete [] fieldB;
	delete [] futures;
	if(needFree == true)
	{
		delete [] redisHostName;
	}
	//Return the number of tasks (0)
	return 0;
}

FluxOut fluxFn(bool doKriging, bool doCoMD, FluxIn * input, redisContext * headRedis)
{
	FluxOut output;

	//See if we are kriging or solving
	bool doSolver = (not doKriging) or (not input->tryKriging);
	bool tryKriging = not doSolver;
	//Do we bother kriging?
	if(tryKriging == true)
	{
		//We do
		//Get data for kriging
		std::vector<double *> oldWs;
		std::vector<double *> oldFs;
		std::vector<double *> oldGs;
		getSortedSubBucketNearZero(input->fields.w, (char *)"comd", headRedis, comdDigits, 10, &oldWs, &oldFs, &oldGs, zeroThresh);
		//Call Kriging on each point
		double resF[2];
		double resG[2];
		double error;
		for(int i = 0; i < 7; i++)
		{
			int info;
			///TODO: Consider refactoring this to kriging.cpp
			std::vector<double> oldF;
			std::vector<double> oldG;
			for(int j = 0; j < oldFs.size(); j++)
			{
				oldF.push_back(oldFs[j][i]);
				oldG.push_back(oldGs[j][i]);
			}
			info = kriging(input->fields.w, 7, oldWs, oldF, resF);
			info = kriging(input->fields.w, 7, oldWs, oldG, resG);
			//Write result
			output.f[i] = resF[0];
			output.g[i] = resG[0];
			//Set error
			if(resF[1] > error)
			{
				error = resF[1];
			}
			if(resG[1] > error)
			{
				error = resG[1];
			}
		}
		freeClear(oldWs);
		freeClear(oldFs);
		freeClear(oldGs);
		//Was error too high?
		if(error > errorThresh)
		{
			//It was, so fall back to comd
			doSolver = true;
		}
		else
		{
			//It was not, so put it to the kriging db
			putData(input->fields.w, output.f, output.g, (char *)"krig", headRedis, krigDigits);
		}
	}
	//Either kriging failed or we never tried
	if(doSolver == true)
	{
		//Call CoMD
		//4 strains
		double strain_xx = input->fields.w[0];
		double strain_xy = input->fields.w[1];
		double strain_yx = input->fields.w[2];
		double strain_yy = input->fields.w[3];
		//2 momentum_fluxes
		double momentum_x = input->fields.w[4];
		double momentum_y = input->fields.w[5];
		//enery_flux
		double energy = input->fields.w[6];
		//CoMD or approximation?
		if(doCoMD == true)
		{
			CoMD_input theInput;
			///FIXME: Fix potdir
			strcpy(theInput.potDir,"../pots");
			strcpy(theInput.potName,"Cu01.eam.alloy");
			strcpy(theInput.potType,"setfl");
			theInput.doeam = 1;
			theInput.nx = 6;
			theInput.ny = 6;
			theInput.nz = 6;
			theInput.nSteps = 1000;
			theInput.printRate = 1;
			//MUST SPECIFY THE FOLLOWING
			theInput.dt = 10.0;
			theInput.lat = 3.6186;
			theInput.temperature = 0;
			theInput.initialDelta = 0.0;
			theInput.defGrad[0] = strain_xx;
			theInput.defGrad[1] = strain_xy;
			theInput.defGrad[2] = strain_yx;
			theInput.defGrad[3] = strain_yy;
			theInput.enDens = energy;
			theInput.momDens[0] = momentum_x;
			theInput.momDens[1] = momentum_y;
			theInput.momDens[2] = 0.0;
			//theInput.rank = rank;
			//theInput.calls = calls;
			theInput.rank = 0;
			theInput.calls = 0;
			//Call it
			CoMD_return theRet = CoMD_lib(&theInput);
			//4 stresses
			double rho = 1.0;
			output.f[0] = momentum_x/rho;
			output.f[1] = momentum_y/rho;
			output.f[2] = 0.0;
			output.f[3] = 0.0;
			output.f[4] = theRet.stressXX;
			output.f[5] = theRet.stressXY;
			output.g[0] = 0.0;
			output.g[1] = 0.0;
			output.g[2] = momentum_x/rho;
			output.g[3] = momentum_y/rho;
			output.g[4] = theRet.stressYX;
			output.g[5] = theRet.stressYY;
			output.f[6] = -theRet.energyDensX;
			output.g[6] = -theRet.energyDensY;
		}
		else
		{
			//4 stresses
			double rho = 1.0;
			//out.f[0] = theRet.momX;
			output.f[0] = momentum_x/rho;
			output.f[1] = momentum_y/rho;
			output.f[2] = 0.0;
			//o->t.f[2] = theRet.momY;
			output.f[3] = 0.0;
			output.f[4] = strain_xx-1 + 0.75*(strain_yy-1);
			output.f[5] = 1.9*strain_xy;
			output.g[0] = 0.0;
			//o->t.g[1] = theRet.momX;
			output.g[1] = 0.0;
			output.g[2] = momentum_x/rho;
			//o->t.g[3] = theRet.momY;
			output.g[3] = momentum_y/rho;
			output.g[4] = 1.9*strain_yx;
			output.g[5] = strain_yy-1 + 0.75*(strain_xx-1);
			//f->6] = g[6] = -0.295;
			output.f[6] = -output.f[0]*sqrt(output.f[4]*output.f[4] + output.f[5]*output.f[5]);
			output.g[6] = -output.g[3]*sqrt(output.g[4]*output.g[4] + output.g[5]*output.g[5]);
		}
		//Put result to DB for future use: Warning, flush if we switch to comd as this is horrible
		putData(input->fields.w, output.f, output.g, (char *)"comd", headRedis, comdDigits);
	}
	return output;
}

//Searches for FluxInput at KEY_curStep_phase_taskID
bool cloudFlux(bool doKriging, bool doCoMD, int curStep, int phase, int taskID, const char * redis_host)
{
	FluxIn input;
	FluxOut output;
	//Connect to redis
	const char * redisHostName;
	//Check redis host
	bool needFree = checkRedisHost(redis_host);
	if(needFree == true)
	{
		redisHostName = getRedisHost(redis_host);
	}
	else
	{
		redisHostName = redis_host;
	}
	redisContext * headRedis = redisConnect(redisHostName, 6379);
	if(headRedis == NULL || headRedis->err)
	{
		printf("Redis error: %s\n", headRedis->errstr);
	}
	//Grab task
	getSingle<FluxIn>(&input, curStep, phase, taskID, headRedis, "TASK");
	//Call fluxFn with input
	output = fluxFn(doKriging, doCoMD, &input, headRedis);

	//Write result to DB
	putSingle<FluxOut>(&output, curStep, phase, taskID, headRedis, "RESULT");
	//cleanup redis
	redisFree(headRedis);
	if(needFree == true)
	{
		delete [] redisHostName;
	}
	return true;
}

bool outputVTK(bool doKriging, bool doCoMD, int * dims, double * dt, double * delta, double * gamma, int curStep, const char * redis_host)
{
	//Connect to redis
	const char * redisHostName;
	//Check redis host
	bool needFree = checkRedisHost(redis_host);
	if(needFree == true)
	{
		redisHostName = getRedisHost(redis_host);
	}
	else
	{
		redisHostName = redis_host;
	}
	redisContext * headRedis = redisConnect(redisHostName, 6379);
	if(headRedis == NULL || headRedis->err)
	{
		printf("Redis error: %s\n", headRedis->errstr);
		return false;
	}
	//Allocate fields buffer
	Node * field = new Node[dims[0]*dims[1]]();
	//Get field data step
	getBlocks<Node>(field, dims[0], dims[1], curStep, 0, headRedis, "FIELD");
	//Build lattice for sanity's sake
	Lattice l;
	l.dim_x = dims[0];
	l.dim_y = dims[1];
	l.dx = delta[0];
	l.dy = delta[1];
	l.dt_x = dt[0];
	l.dt_y = dt[1];
	//VTK output it
	printf_fields_vtk(curStep, field, l, dims[0]*dims[1]);
	//cleanup redis
	redisFree(headRedis);
	//Don't save fields, we didn't do anything with it
	delete [] field;
	if(needFree == true)
	{
		delete [] redisHostName;
	}
	return true;
}

bool tryShortCircuit(int * dims, int curStep, const char * redis_host)
{
	//For short circuiting purposes, we do it on a per-step basis
	//Just check for all tiles  of phase 0 of the next timestep
	///WARNING: This could fail if the job crashes partway through a puts
	//Connect to redis
	const char * redisHostName;
	//Check redis host
	bool needFree = checkRedisHost(redis_host);
	if(needFree == true)
	{
		redisHostName = getRedisHost(redis_host);
	}
	else
	{
		redisHostName = redis_host;
	}
	redisContext * headRedis = redisConnect(redisHostName, 6379);
	if(headRedis == NULL || headRedis->err)
	{
		printf("Redis error: %s\n", headRedis->errstr);
		return false;
	}
	//Check each tile
	bool retBool = true;
	unsigned int nTiles = getNumBlocks(dims[0], dims[1]);
	for(unsigned int i = 0; i < nTiles; i++)
	{
		//Build Key
		char keyBuffer[maxKeyLength];
		buildBlockKey(keyBuffer, curStep+1, 0, i, dims[0], dims[1], "FIELD");
		//Use Key to check for existence
		redisReply * reply;
		reply = (redisReply *) redisCommand(headRedis, "EXISTS %s", keyBuffer);
		assert(reply->type == REDIS_REPLY_INTEGER);
		if(reply->integer == 0)
		{
			retBool = false;
			i = nTiles + 1;
		}
		freeReplyObject(reply);
	}
	//Kill redis
	redisFree(headRedis);
	if(needFree == true)
	{
		delete [] redisHostName;
	}
	return retBool;
}

inline bool checkRedisHost(const char * inHost)
{
	//Check if we passed a path
	//For santiy's sake
	//	we assume all paths are either relative and start with a "."
	//	or absolute and start with a "/"
	if( (inHost[0] == '.') or (inHost[0] == '/') )
	{
		return true;
	}
	else
	{
		return false;
	}
}

char * getRedisHost(const char * filePath)
{
	char * retHost = nullptr;
	//Open file we passed in
	std::ifstream redisFile;
	redisFile.open(filePath);
	if(redisFile.is_open())
	{
		//Get hostname of self
		char cHost[64];
		gethostname(cHost, 64);
		std::string strHost(cHost);
		bool found = false;
		std::string line;
		//Iterate through file we passed in
		while( std::getline(redisFile, line) && found == true )
		{
			//Check if this is the line that corresponds to us
			if(line.compare(0, strHost.length() - 1, strHost) == 0)
			{
				//This is it, use the host after the tab
				std::string ourRedis = line.substr(line.find('\t')+1);
				strcpy(cHost, ourRedis.c_str());
				//End the loop
				found = true;
			}
		}
		//We either found it or we didn't, close the file
		redisFile.close();
	}
	//Hopefully return the hostname, not nullptr
	return retHost;
}

bool initEverything(bool doKriging, bool doCoMD, int * dims, double * dt, double * delta, double * gamma)
{
	return initEverything(doKriging, doCoMD, dims, dt, delta, gamma, "localhost");
}
int prepFirstFlux(bool doKriging, bool doCoMD, int * dims, double * dt, double * delta, double * gamma, int curStep)
{
	return prepFirstFlux(doKriging, doCoMD, dims, dt, delta, gamma, curStep, "localhost");
}
int prepSecondFlux(bool doKriging, bool doCoMD, int * dims, double * dt, double * delta, double * gamma, int curStep)
{
	return prepSecondFlux(doKriging, doCoMD, dims, dt, delta, gamma, curStep, "localhost");
}
int prepThirdFlux(bool doKriging, bool doCoMD, int * dims, double * dt, double * delta, double * gamma, int curStep)
{
	return prepThirdFlux(doKriging, doCoMD, dims, dt, delta, gamma, curStep, "localhost");
}
int prepLastFlux(bool doKriging, bool doCoMD, int * dims, double * dt, double * delta, double * gamma, int curStep)
{
	return prepLastFlux(doKriging, doCoMD, dims, dt, delta, gamma, curStep, "localhost");
}
int finishStep(bool doKriging, bool doCoMD, int * dims, double * dt, double * delta, double * gamma, int curStep)
{
	return finishStep(doKriging, doCoMD, dims, dt, delta, gamma, curStep, "localhost");
}
bool cloudFlux(bool doKriging, bool doCoMD, int curStep, int phase, int taskID)
{
	return cloudFlux(doKriging, doCoMD, curStep, phase, taskID, "localhost");
}
bool outputVTK(bool doKriging, bool doCoMD, int * dims, double * dt, double * delta, double * gamma, int curStep)
{
	return outputVTK(doKriging, doCoMD, dims, dt, delta, gamma, curStep, "localhost");
}
bool tryShortCircuit(int * dims, int curStep)
{
	return tryShortCircuit(dims, curStep, "localhost");
}
