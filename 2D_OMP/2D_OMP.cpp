#include <iostream>
#include <cstdlib>

#include "CoHMM_OMP.hpp"

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
	unsigned int dims[2] = {(unsigned int) atoi(argv[1]), (unsigned int) atoi(argv[2])};
	double dt[2] = {0.1, 0.1};
	double delta[2] = {1.0, 1.0};
	double gamma[3];
	gamma[0] = 0; //mom_gamma
	gamma[1] = gamma[0]; //strain_gamma
	gamma[2] = 0.1 * gamma[1];//en_gamma

	unsigned int numSteps = atoi(argv[3]);

	CoHMM_OMP solver(dims, delta, dt, gamma, numSteps, doKriging, doCoMD, argv[4]);

	//Initialize
	std::cout << "Initializing " << dims[0] << " by " << dims[1] << " grid" << std::endl;
	solver.initializeConservedFields(InitialConditions_e::X);
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
	return 0;
}