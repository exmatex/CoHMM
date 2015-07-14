#include <iostream>
#include <cstdlib>

#include <mpi.h>

#include "CoHMM_DaD.hpp"

int main(int argc, char ** argv)
{
	//Start up MPI
	MPI_Init(&argc, &argv);
	int pid, nProcs;
	MPI_Comm_rank(MPI_COMM_WORLD, &pid);
	MPI_Comm_size(MPI_COMM_WORLD, &nProcs);

	//dimX dimY nSteps redis_server
	if( argc != 5)
	{
		std::cerr <<  "mpirun -np <nProcs> ./2D_MPI <dim_x> <dim_y> <nsteps> <redis_server>" << std::endl;
		return 1;
	}
	//Set up parameters
	const bool fineGrainFT = true;
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
	if(pid == 0)
	{
		std::cout << "Initializing " << dims[0] << " by " << dims[1] << " grid" << std::endl;
		initEverything(doKriging, doCoMD, dims, dt, delta, gamma);
		std::cout << "Initialized" << std::endl;
		//Loop
		std::cout << "Running for " << numSteps << " iterations" << std::endl;
	}
	for(unsigned int t = 0; t < numSteps; t++)
	{
		int nTasks;
		bool shortCircuit;
		int curRound;
		if(pid == 0)
		{
			std::cout << t << ": Vising to Verifying" << std::endl;
			outputVTK(doKriging, doCoMD, dims, dt, delta, gamma, t, argv[4]);
			//Do a short circuit test
			shortCircuit = tryShortCircuit(dims, t, argv[4]);
		}
		MPI_Bcast(&shortCircuit, 1, MPI::BOOL, 0, MPI_COMM_WORLD);
		if(shortCircuit)
		{
			//Short circuit succeeded
			std::cout << t << ": Short Circuit Successful, on to the next step!" << std::endl;
		}
		else
		{
			if(pid == 0)
			{
				std::cout << t << ": First Flux" << std::endl;
				nTasks = prepFirstFlux(doKriging, doCoMD, dims, dt, delta, gamma, t, argv[4]);
				std::cout << t << ": Doing " << nTasks << " fluxes" << std::endl;
			}
			MPI_Bcast(&nTasks, 1, MPI::INT, 0, MPI_COMM_WORLD);
			for(unsigned int i = 0; i < nTasks; i++)
			{
				if(i % nProcs == pid)
				{
					cloudFlux(doKriging, doCoMD, t, 0, i, argv[4]);
				}
			}
			if(fineGrainFT == true)
			{
				curRound = 0;
				if(pid == 0)
				{
					std::cout << t << ": Checking First Flux" << std::endl;
					nTasks = checkStepForFaults(dims, t, 0, curRound, argv[4]);
				}
				MPI_Bcast(&nTasks, 1, MPI::INT, 0, MPI_COMM_WORLD);
				while(nTasks != 0)
				{
					if(pid == 0)
					{
						std::cout << t << ": Redoing " << nTasks << " Tasks" << std::endl;
					}
					for (unsigned int i = 0; i < nTasks; i++)
			 		{
						if(i % nProcs == pid)
						{
							retryCloudFlux(doKriging, doCoMD, t, 0, i, curRound, argv[4]);
						}
					}
					//See if we are done
					curRound++;
					if(pid == 0)
					{
						nTasks = checkStepForFaults(dims, t, 0, curRound, argv[4]);
					}
					MPI_Bcast(&nTasks, 1, MPI::INT, 0, MPI_COMM_WORLD);
				}
			}
			if(pid == 0)
			{
				std::cout << t << ": Second Flux" << std::endl;
				nTasks = prepSecondFlux(doKriging, doCoMD, dims, dt, delta, gamma, t, argv[4]);
				std::cout << t << ": Doing " << nTasks << " fluxes" << std::endl;
			}
			MPI_Bcast(&nTasks, 1, MPI::INT, 0, MPI_COMM_WORLD);
			for(unsigned int i = 0; i < nTasks; i++)
			{
				if(i % nProcs == pid)
				{
					cloudFlux(doKriging, doCoMD, t, 1, i, argv[4]);
				}
			}
			if(fineGrainFT == true)
			{
				curRound = 0;
				if(pid == 0)
				{
					std::cout << t << ": Checking Second Flux" << std::endl;
					nTasks = checkStepForFaults(dims, t, 1, curRound, argv[4]);
				}
				MPI_Bcast(&nTasks, 1, MPI::INT, 0, MPI_COMM_WORLD);
				while(nTasks != 0)
				{
					if(pid == 0)
					{
						std::cout << t << ": Redoing " << nTasks << " Tasks" << std::endl;
					}
					for (unsigned int i = 0; i < nTasks; i++)
			 		{
						if(i % nProcs == pid)
						{
							retryCloudFlux(doKriging, doCoMD, t, 1, i, curRound, argv[4]);
						}
					}
					//See if we are done
					curRound++;
					if(pid == 0)
					{
						nTasks = checkStepForFaults(dims, t, 1, curRound, argv[4]);
					}
					MPI_Bcast(&nTasks, 1, MPI::INT, 0, MPI_COMM_WORLD);
				}
			}
			if(pid == 0)
			{
				std::cout << t << ": Third Flux" << std::endl;
				nTasks = prepThirdFlux(doKriging, doCoMD, dims, dt, delta, gamma, t, argv[4]);
				std::cout << t << ": Doing " << nTasks << " fluxes" << std::endl;
			}
			MPI_Bcast(&nTasks, 1, MPI::INT, 0, MPI_COMM_WORLD);
			for(unsigned int i = 0; i < nTasks; i++)
			{
				if(i % nProcs == pid)
				{
					cloudFlux(doKriging, doCoMD, t, 2, i, argv[4]);
				}
			}
			if(fineGrainFT == true)
			{
				curRound = 0;
				if(pid == 0)
				{
					std::cout << t << ": Checking Third Flux" << std::endl;
					nTasks = checkStepForFaults(dims, t, 2, curRound, argv[4]);
				}
				MPI_Bcast(&nTasks, 1, MPI::INT, 0, MPI_COMM_WORLD);
				while(nTasks != 0)
				{
					if(pid == 0)
					{
						std::cout << t << ": Redoing " << nTasks << " Tasks" << std::endl;
					}
					for (unsigned int i = 0; i < nTasks; i++)
			 		{
						if(i % nProcs == pid)
						{
							retryCloudFlux(doKriging, doCoMD, t, 2, i, curRound, argv[4]);
						}
					}
					//See if we are done
					curRound++;
					if(pid == 0)
					{
						nTasks = checkStepForFaults(dims, t, 2, curRound, argv[4]);
					}
					MPI_Bcast(&nTasks, 1, MPI::INT, 0, MPI_COMM_WORLD);
				}
			}
			if(pid == 0)
			{
				std::cout << t << ": Last Flux" << std::endl;
				nTasks = prepLastFlux(doKriging, doCoMD, dims, dt, delta, gamma, t, argv[4]);
				std::cout << t << ": Doing " << nTasks << " fluxes" << std::endl;
			}
			MPI_Bcast(&nTasks, 1, MPI::INT, 0, MPI_COMM_WORLD);
			for(unsigned int i = 0; i < nTasks; i++)
			{
				if(i % nProcs == pid)
				{
					cloudFlux(doKriging, doCoMD, t, 3, i, argv[4]);
				}
			}
			if(fineGrainFT == true)
			{
				curRound = 0;
				if(pid == 0)
				{
					std::cout << t << ": Checking Last Flux" << std::endl;
					nTasks = checkStepForFaults(dims, t, 3, curRound, argv[4]);
				}
				MPI_Bcast(&nTasks, 1, MPI::INT, 0, MPI_COMM_WORLD);
				while(nTasks != 0)
				{
					if(pid == 0)
					{
						std::cout << t << ": Redoing " << nTasks << " Tasks" << std::endl;
					}
					for (unsigned int i = 0; i < nTasks; i++)
			 		{
						if(i % nProcs == pid)
						{
							retryCloudFlux(doKriging, doCoMD, t, 3, i, curRound, argv[4]);
						}
					}
					//See if we are done
					curRound++;
					if(pid == 0)
					{
						nTasks = checkStepForFaults(dims, t, 3, curRound, argv[4]);
					}
					MPI_Bcast(&nTasks, 1, MPI::INT, 0, MPI_COMM_WORLD);
				}
			}
			if(pid == 0)
			{
				std::cout << t << ": Finish Step, no Fluxes" << std::endl;
				finishStep(doKriging, doCoMD, dims, dt, delta, gamma, t, argv[4]);
			}
		}
	}
	if(pid == 0)
	{
		//Final vis
		std::cout << numSteps << ": Vising to Verifying" << std::endl;
		outputVTK(doKriging, doCoMD, dims, dt, delta, gamma, numSteps, argv[4]);
		std::cout << "Ran for " << numSteps << " iterations" << std::endl;
	}
	//End MPI
	MPI_Finalize();

	return 0;
}
