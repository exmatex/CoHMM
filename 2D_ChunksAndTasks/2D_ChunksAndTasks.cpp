#include <iostream>
#include <cstdlib>
#include <cstring>
#include <vector>

#include "chunks_and_tasks.h"

#include "CoHMM_DaD.hpp"

#include "2D_ChunksAndTasks.hpp"
#include "FluxTask.hpp"
#include "FluxChunk.hpp"
#include "CppBoolChunk.hpp"


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
	if( argc != 5)
	{
		std::cerr <<  "./2D_DaDTest <dim_x> <dim_y> <nsteps> <redis_server>" << std::endl;
		return 1;
	}
	//Set up parameters
	bool doKriging = true;
	bool doCoMD = false;
	int dims[2] = {atoi(argv[1]), atoi(argv[2])};
	double dt[2] = {0.1, 0.1};
	double delta[2] = {1.0, 1.0};
	double gamma[3];
	gamma[0] = 0; //mom_gamma
	gamma[1] = gamma[0]; //strain_gamma
	gamma[2] = 0.1 * gamma[1];//en_gamma

	unsigned int numSteps = 10;

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
		std::cout << t << ": First Flux" << std::endl;
		nTasks = prepFirstFlux(doKriging, doCoMD, dims, dt, delta, gamma, t, argv[4]);
		std::cout << t << ": Doing " << nTasks << " fluxes" << std::endl;
		#pragma omp parallel for
		for(unsigned int i = 0; i < nTasks; i++)
		{
			cloudFlux(doKriging, doCoMD, t, 0, i, argv[4]);
		}
		std::cout << t << ": Second Flux" << std::endl;
		nTasks = prepSecondFlux(doKriging, doCoMD, dims, dt, delta, gamma, t, argv[4]);
		std::cout << t << ": Doing " << nTasks << " fluxes" << std::endl;
		#pragma omp parallel for
		for(unsigned int i = 0; i < nTasks; i++)
		{
			cloudFlux(doKriging, doCoMD, t, 1, i, argv[4]);
		}
		std::cout << t << ": Third Flux" << std::endl;
		nTasks = prepThirdFlux(doKriging, doCoMD, dims, dt, delta, gamma, t, argv[4]);
		std::cout << t << ": Doing " << nTasks << " fluxes" << std::endl;
		#pragma omp parallel for
		for(unsigned int i = 0; i < nTasks; i++)
		{
			cloudFlux(doKriging, doCoMD, t, 2, i, argv[4]);
		}
		std::cout << t << ": Last Flux" << std::endl;
		nTasks = prepLastFlux(doKriging, doCoMD, dims, dt, delta, gamma, t, argv[4]);
		std::cout << t << ": Doing " << nTasks << " fluxes" << std::endl;
		#pragma omp parallel for
		for(unsigned int i = 0; i < nTasks; i++)
		{
			cloudFlux(doKriging, doCoMD, t, 3, i, argv[4]);
		}
		std::cout << t << ": Finish Step, no Fluxes" << std::endl;
		finishStep(doKriging, doCoMD, dims, dt, delta, gamma, t, argv[4]);
	}
	//Final vis
	std::cout << numSteps << ": Vising to Verifying" << std::endl;
	outputVTK(doKriging, doCoMD, dims, dt, delta, gamma, numSteps, argv[4]);
	std::cout << "Ran for " << numSteps << " iterations" << std::endl;

	return 0;
}

