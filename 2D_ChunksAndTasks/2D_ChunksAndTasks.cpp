#include <iostream>
#include <cstdlib>
#include <cstring>
#include <vector>

#include "chunks_and_tasks.h"

#include "CoHMM_DaD.hpp"

#include "2D_ChunksAndTasks.hpp"
#include "FluxTask.hpp"
#include "FluxChunk.hpp"
#include "RetryTask.hpp"
#include "RetryChunk.hpp"
#include "CppBoolChunk.hpp"

bool iterativeRetry(int * dims, bool doKriging, bool doCoMD, unsigned int step, unsigned int phase, const char * redis_host)
{
	//Make base RetryChunk
	RetryChunk baseTask;
	memcpy(baseTask.redisHost, redis_host, sizeof(char) * MAX_HOST_LENGTH);
	baseTask.doKriging = doKriging;
	baseTask.doCoMD = doCoMD;
	baseTask.step = step;
	baseTask.phase = phase;

	int curRound = 0;
	int nTasks = checkStepForFaults(dims, step, phase, curRound, redis_host);
	//Ierate until it is not failed
	while(nTasks != 0)
	{
		std::cout << step << ": Redoing " << nTasks << " Tasks" << std::endl;
		//Make vectpr of chunk ids of length nTasks
		std::vector<cht::ChunkID> taskHandles(nTasks);
		std::vector<cht::ChunkID> chunkHandles(nTasks);
		baseTask.round = curRound;
		//For each task, enqueue it
		for(unsigned int i = 0; i < nTasks; i++)
		{
			//Set chunk's ID
			baseTask.taskID = i;
			//Push chunk
			chunkHandles[i] = cht::registerChunk(new RetryChunk(baseTask));
			//Push task
			taskHandles[i] = cht::executeMotherTask<RetryTask>(chunkHandles[i]);
		}
		//For each task, block and free
		for(unsigned int i = 0; i < nTasks; i++)
		{
			cht::shared_ptr<CppBoolChunk const> resBool;
			cht::getChunk(taskHandles[i], resBool);
			if(resBool->retVal != true)
			{
				std::cerr << "Error: CHT Task Failed Somehow" << std::endl;
			}
			cht::deleteChunk(chunkHandles[i]);
			cht::deleteChunk(taskHandles[i]);
		}
		//See if we are done
		curRound++;
		nTasks = checkStepForFaults(dims, step, phase, curRound, redis_host);
	}
}

bool fluxParallelFor(bool doKriging, bool doCoMD, unsigned int step, unsigned int phase, unsigned int nTasks, const char * redisHost)
{
	//Make base FluxChunk
	FluxChunk baseTask;
	memcpy(baseTask.redisHost, redisHost, sizeof(char) * MAX_HOST_LENGTH);
	baseTask.doKriging = doKriging;
	baseTask.doCoMD = doCoMD;
	baseTask.step = step;
	baseTask.phase = phase;
	//Make vectpr of chunk ids of length nTasks
	std::vector<cht::ChunkID> taskHandles(nTasks);
	std::vector<cht::ChunkID> chunkHandles(nTasks);
	//For each task, enqueue it
	for(unsigned int i = 0; i < nTasks; i++)
	{
		//Set chunk's ID
		baseTask.taskID = i;
		//Push chunk
		chunkHandles[i] = cht::registerChunk(new FluxChunk(baseTask));
		//Push task
		taskHandles[i] = cht::executeMotherTask<FluxTask>(chunkHandles[i]);
	}
	//For each task, block and free
	for(unsigned int i = 0; i < nTasks; i++)
	{
		cht::shared_ptr<CppBoolChunk const> resBool;
		cht::getChunk(taskHandles[i], resBool);
		if(resBool->retVal != true)
		{
			std::cerr << "Error: CHT Task Failed Somehow" << std::endl;
		}
		cht::deleteChunk(chunkHandles[i]);
		cht::deleteChunk(taskHandles[i]);
	}
	return true;
}

int main(int argc, char ** argv)
{
	//<dim_x> <dim_y> <nsteps> <redis_server> <database error threshold> <Kriging error threshold> <Gaussian noise strength>
	//dimX dimY nSteps redis_server
	if( argc != 6)
	{
		std::cerr <<  "mpirun -np 1 ./2D_ChunksAndTasks <dim_x> <dim_y> <nsteps> <redis_server> <nWorkers>" << std::endl;
		return 1;
	}
	//Start CHT
	std::cout << "Trying to start CHT" << std::endl;
	cht::start();
	cht::extras::setNWorkers(atoi(argv[5]));
	std::cout << "Succeeded to start CHT" << std::endl;
	//Set up parameters
	const bool doKriging = true;
	const bool doCoMD = false;
	const bool fineGrainFT = false;
	int dims[2] = {atoi(argv[1]), atoi(argv[2])};
	double dt[2] = {0.1, 0.1};
	double delta[2] = {1.0, 1.0};
	double gamma[3];
	gamma[0] = 0; //mom_gamma
	gamma[1] = gamma[0]; //strain_gamma
	gamma[2] = 0.1 * gamma[1];//en_gamma

	unsigned int numSteps = atoi(argv[3]);

	//Initialize
	std::cout << "Initializing " << dims[0] << " by " << dims[1] << " grid" << std::endl;
	initEverything(doKriging, doCoMD, dims, dt, delta, gamma);
	std::cout << "Initialized" << std::endl;
	//Loop
	std::cout << "Running for " << numSteps << " iterations" << std::endl;
	for(unsigned int t = 0; t < numSteps; t++)
	{
		int nTasks;
		std::cout << t << ": Vising to Verifying" << std::endl;
		outputVTK(doKriging, doCoMD, dims, dt, delta, gamma, t, argv[4]);
		//Do a short circuit test
		if(tryShortCircuit(dims, t, argv[4]))
		{
			//Short circuit succeeded
			std::cout << t << ": Short Circuit Successful, on to the next step!" << std::endl;
		}
		else
		{
			std::cout << t << ": First Flux" << std::endl;
			nTasks = prepFirstFlux(doKriging, doCoMD, dims, dt, delta, gamma, t, argv[4]);
			std::cout << t << ": Doing " << nTasks << " fluxes" << std::endl;
			fluxParallelFor(doKriging, doCoMD, t, 0, nTasks, argv[4]);
			if(fineGrainFT == true)
			{
				std::cout << t << ": Checking First Flux" << std::endl;
				iterativeRetry(dims, doKriging, doCoMD, t, 0, argv[4]);
			}
			std::cout << t << ": Second Flux" << std::endl;
			nTasks = prepSecondFlux(doKriging, doCoMD, dims, dt, delta, gamma, t, argv[4]);
			std::cout << t << ": Doing " << nTasks << " fluxes" << std::endl;
			fluxParallelFor(doKriging, doCoMD, t, 1, nTasks, argv[4]);
			if(fineGrainFT == true)
			{
				std::cout << t << ": Checking Second Flux" << std::endl;
				iterativeRetry(dims, doKriging, doCoMD, t, 1, argv[4]);
			}
			std::cout << t << ": Third Flux" << std::endl;
			nTasks = prepThirdFlux(doKriging, doCoMD, dims, dt, delta, gamma, t, argv[4]);
			std::cout << t << ": Doing " << nTasks << " fluxes" << std::endl;
			fluxParallelFor(doKriging, doCoMD, t, 2, nTasks, argv[4]);
			if(fineGrainFT == true)
			{
				std::cout << t << ": Checking Third Flux" << std::endl;
				iterativeRetry(dims, doKriging, doCoMD, t, 2, argv[4]);
			}
			std::cout << t << ": Last Flux" << std::endl;
			nTasks = prepLastFlux(doKriging, doCoMD, dims, dt, delta, gamma, t, argv[4]);
			std::cout << t << ": Doing " << nTasks << " fluxes" << std::endl;
			fluxParallelFor(doKriging, doCoMD, t, 3, nTasks, argv[4]);
			if(fineGrainFT == true)
			{
				std::cout << t << ": Checking Last Flux" << std::endl;
				iterativeRetry(dims, doKriging, doCoMD, t, 3, argv[4]);
			}
			std::cout << t << ": Finish Step, no Fluxes" << std::endl;
			finishStep(doKriging, doCoMD, dims, dt, delta, gamma, t, argv[4]);
		}
	}
	//Final vis
	std::cout << numSteps << ": Vising to Verifying" << std::endl;
	outputVTK(doKriging, doCoMD, dims, dt, delta, gamma, numSteps, argv[4]);
	std::cout << "Ran for " << numSteps << " iterations" << std::endl;
	//End CHT
	cht::stop();
	return 0;
}
