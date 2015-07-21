#include <iostream>
#include <cstdlib>

#include "CoHMM_DaD.hpp"

int main(int argc, char ** argv)
{
	//<dim_x> <dim_y> <nsteps> <redis_server> <database error threshold> <Kriging error threshold> <Gaussian noise strength>
	//dimX dimY nSteps redis_server
	if( argc != 5)
	{
		std::cerr <<  "./2D_OMP <dim_x> <dim_y> <nsteps> <redis_server>" << std::endl;
		return 1;
	}
	//Set up parameters
	const bool fineGrainFT = false;
	const bool doKriging = true;
	const bool doCoMD = false;
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
		int curRound;
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
			#pragma omp parallel for
			for(unsigned int i = 0; i < nTasks; i++)
			{
				cloudFlux(doKriging, doCoMD, t, 0, i, argv[4]);
			}
			if(fineGrainFT == true)
			{
				std::cout << t << ": Checking First Flux" << std::endl;
				curRound = 0;
				nTasks = checkStepForFaults(dims, t, 0, curRound, argv[4]);
				while(nTasks != 0)
				{
					std::cout << t << ": Redoing " << nTasks << " Tasks" << std::endl;
					#pragma omp parallel
					for (unsigned int i = 0; i < nTasks; i++)
			 		{
						retryCloudFlux(doKriging, doCoMD, t, 0, i, curRound, argv[4]);
					}
					//See if we are done
					curRound++;
					nTasks = checkStepForFaults(dims, t, 0, curRound, argv[4]);
				}
			}
			std::cout << t << ": Second Flux" << std::endl;
			nTasks = prepSecondFlux(doKriging, doCoMD, dims, dt, delta, gamma, t, argv[4]);
			std::cout << t << ": Doing " << nTasks << " fluxes" << std::endl;
			#pragma omp parallel for
			for(unsigned int i = 0; i < nTasks; i++)
			{
				cloudFlux(doKriging, doCoMD, t, 1, i, argv[4]);
			}
			if(fineGrainFT == true)
			{
				std::cout << t << ": Checking Second Flux" << std::endl;
				curRound = 0;
				nTasks = checkStepForFaults(dims, t, 1, curRound, argv[4]);
				while(nTasks != 0)
				{
					std::cout << t << ": Redoing " << nTasks << " Tasks" << std::endl;
					#pragma omp parallel
					for (unsigned int i = 0; i < nTasks; i++)
			 		{
						retryCloudFlux(doKriging, doCoMD, t, 1, i, curRound, argv[4]);
					}
					//See if we are done
					curRound++;
					nTasks = checkStepForFaults(dims, t, 1, curRound, argv[4]);
				}
			}
			std::cout << t << ": Third Flux" << std::endl;
			nTasks = prepThirdFlux(doKriging, doCoMD, dims, dt, delta, gamma, t, argv[4]);
			std::cout << t << ": Doing " << nTasks << " fluxes" << std::endl;
			#pragma omp parallel for
			for(unsigned int i = 0; i < nTasks; i++)
			{
				cloudFlux(doKriging, doCoMD, t, 2, i, argv[4]);
			}
			if(fineGrainFT == true)
			{
				std::cout << t << ": Checking Third Flux" << std::endl;
				curRound = 0;
				nTasks = checkStepForFaults(dims, t, 2, curRound, argv[4]);
				while(nTasks != 0)
				{
					std::cout << t << ": Redoing " << nTasks << " Tasks" << std::endl;
					#pragma omp parallel
					for (unsigned int i = 0; i < nTasks; i++)
			 		{
						retryCloudFlux(doKriging, doCoMD, t, 2, i, curRound, argv[4]);
					}
					//See if we are done
					curRound++;
					nTasks = checkStepForFaults(dims, t, 2, curRound, argv[4]);
				}
			}
			std::cout << t << ": Last Flux" << std::endl;
			nTasks = prepLastFlux(doKriging, doCoMD, dims, dt, delta, gamma, t, argv[4]);
			std::cout << t << ": Doing " << nTasks << " fluxes" << std::endl;
			#pragma omp parallel for
			for(unsigned int i = 0; i < nTasks; i++)
			{
				cloudFlux(doKriging, doCoMD, t, 3, i, argv[4]);
			}
			if(fineGrainFT == true)
			{
				std::cout << t << ": Checking Last Flux" << std::endl;
				curRound = 0;
				nTasks = checkStepForFaults(dims, t, 3, curRound, argv[4]);
				while(nTasks != 0)
				{
					std::cout << t << ": Redoing " << nTasks << " Tasks" << std::endl;
					#pragma omp parallel
					for (unsigned int i = 0; i < nTasks; i++)
			 		{
						retryCloudFlux(doKriging, doCoMD, t, 3, i, curRound, argv[4]);
					}
					//See if we are done
					curRound++;
					nTasks = checkStepForFaults(dims, t, 3, curRound, argv[4]);
				}
			}
			std::cout << t << ": Finish Step, no Fluxes" << std::endl;
			finishStep(doKriging, doCoMD, dims, dt, delta, gamma, t, argv[4]);
		}
	}
	//Final vis
	std::cout << numSteps << ": Vising to Verifying" << std::endl;
	outputVTK(doKriging, doCoMD, dims, dt, delta, gamma, numSteps, argv[4]);
	std::cout << "Ran for " << numSteps << " iterations" << std::endl;

	return 0;
}
