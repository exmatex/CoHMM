#include <iostream>
#include <cstdlib>

#include "CoHMM_DaD.hpp"

#include "2D_CnC.hpp"

DaDContext::DaDContext()
	:
	CnC::context<DaDContext>(),
	fluxTags(*this, "fluxTags"),
	fluxTask(*this, "fluxTasks"),
	globalItem(*this, "globalItem")
{
	//Indicate which tags correspond to which steps
	fluxTags.prescribes(fluxTask, *this);

	//Indicate which step consumes which item (ha ha)
	fluxTask.consumes(globalItem);
}

int Flux_Task::execute(const Flux_Tag &tag, DaDContext &c) const
{
	//Get The Item
	Flux_Item runConfig;
	c.globalItem.get(0, runConfig);
	//Get params
	unsigned int step = std::get<0>(tag);
	unsigned int phase = std::get<1>(tag);
	unsigned int task = std::get<2>(tag);
	//Call
	cloudFlux(runConfig.doKriging, runConfig.doCoMD, step, phase, task, runConfig.redis_host);
	//Return
	return CnC::CNC_Success;
}

int main(int argc, char ** argv)
{
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
	}
	//Final vis
	std::cout << numSteps << ": Vising to Verifying" << std::endl;
	outputVTK(doKriging, doCoMD, dims, dt, delta, gamma, numSteps, argv[4]);
	std::cout << "Ran for " << numSteps << " iterations" << std::endl;

	return 0;
}

