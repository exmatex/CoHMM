#include <iostream>
#include <cstdlib>

#include "chunks_and_tasks.h"

#include "CoHMM_CHT.hpp"

int main(int argc, char ** argv)
{
	//<dim_x> <dim_y> <nsteps> <redis_server> <database error threshold> <Kriging error threshold> <Gaussian noise strength>
	//dimX dimY nSteps redis_server
	if( argc != 6)
	{
		std::cerr <<  "mpirun -np 1 ./2D_CHT <dim_x> <dim_y> <nsteps> <redis_server> <nWorkers>" << std::endl;
		return 1;
	}
	//Start CHT
	std::cout << "Trying to start CHT with " << argv[5] << " workers" << std::endl;
	cht::extras::setNWorkers(atoi(argv[5]));
	cht::start();
	std::cout << "Succeeded to start CHT" << std::endl;

	//Set up parameters
	const bool fineGrainFT = false;
	const bool doKriging = true;
	const bool doCoMD = false;
	unsigned int dims[2] = {(unsigned int) atoi(argv[1]), (unsigned int) atoi(argv[2])};
	double dt[2] = {0.1, 0.1};
	double delta[2] = {1.0, 1.0};
	double gamma[3];
	gamma[0] = 0; //mom_gamma
	gamma[1] = gamma[0]; //strain_gamma
	gamma[2] = 0.1 * gamma[1];//en_gamma

	unsigned int numSteps = atoi(argv[3]);

	CoHMM_CHT solver(dims, delta, dt, gamma, numSteps, doKriging, doCoMD, argv[4]);

	//Initialize
	std::cout << "Initializing " << dims[0] << " by " << dims[1] << " grid" << std::endl;
	solver.initializeConservedFields();
	std::cout << "Initialized" << std::endl;
	//Loop
	std::cout << "Running for " << numSteps << " iterations" << std::endl;
	for(unsigned int t = 0; t < numSteps; t++)
	{
		int nTasks;
		std::cout << t << ": Vising to Verifying" << std::endl;
		solver.outputVTK();
		std::cout << t << ": First Flux" << std::endl;
		solver.prepFirstFlux();
		nTasks = solver.getNumberOfTasks();
		std::cout << t << ": Doing " << nTasks << " fluxes" << std::endl;
		solver.processFluxes();
		std::cout << t << ": Second Flux" << std::endl;
		solver.prepSecondFlux();
		nTasks = solver.getNumberOfTasks();
		std::cout << t << ": Doing " << nTasks << " fluxes" << std::endl;
		solver.processFluxes();
		std::cout << t << ": Third Flux" << std::endl;
		solver.prepThirdFlux();
		nTasks = solver.getNumberOfTasks();
		std::cout << t << ": Doing " << nTasks << " fluxes" << std::endl;
		solver.processFluxes();
		std::cout << t << ": Last Flux" << std::endl;
		solver.prepLastFlux();
		nTasks = solver.getNumberOfTasks();
		std::cout << t << ": Doing " << nTasks << " fluxes" << std::endl;
		solver.processFluxes();
		std::cout << t << ": Finish Step, no Fluxes" << std::endl;
		solver.finishStep();
	}
	//Final vis
	std::cout << numSteps << ": Vising to Verifying" << std::endl;
	solver.outputVTK();
	std::cout << "Ran for " << numSteps << " iterations" << std::endl;

	//End CHT
	cht::stop();
	return 0;
}
