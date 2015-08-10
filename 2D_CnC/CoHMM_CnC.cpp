#include <iostream>
#include <cstdlib>

#include "CoHMM_CnC.hpp"
#include "CoHMM_CnCaD.hpp"

CnCaDContext::CnCaDContext()
	:
	CnC::context<CnCaDContext>(),
	fluxTags(*this, "fluxTags"),
	fluxTask(*this, "fluxTasks"),
	globalItem(*this, "globalItem"),
	blockItems(*this, "blockItems"),
	taskItems(*this, "taskItems"),
	futureItems(*this, "futureItems")
{
	//Indicate which tags correspond to which steps
	fluxTags.prescribes(fluxTask, *this);

	//Indicate which step consumes which item (ha ha)
	fluxTask.consumes(globalItem);

	//Nothing else because we are doing this dirty
}

int Flux_Task::execute(const Flux_Tag &tag, CnCaDContext &c) const
{
	//Get The Item
	Global_Singleton_Item runConfig;
	c.globalItem.get(0, runConfig);
	//Get params
	unsigned int step = tag.first.first;
	unsigned int phase = tag.first.second;
	unsigned int task = tag.second;
	//Call
	cloudFlux(runConfig.doKriging, runConfig.doCoMD, step, phase, task, runConfig.redis_host);
	//Return
	return CnC::CNC_Success;
}

void parallelFor(unsigned int step, unsigned int phase, unsigned int nTasks, CnCaDContext &ctxt)
{
	for(unsigned int i = 0; i < nTasks; i++)
	{
		//Make task tag
		Flux_Tag tag;
		tag.first.first = step;
		tag.first.second = phase;
		tag.second = i;
		//Put task tag
		ctxt.fluxTags.put(tag);
	}
	ctxt.wait();
	return;
}

int main(int argc, char ** argv)
{
	//dimX dimY nSteps redis_server
	if( argc != 5)
	{
		std::cerr <<  "./2D_DaDTest <dim_x> <dim_y> <nsteps> <redis_server>" << std::endl;
		return 1;
	}
	//Initialize CnC
#ifdef CNC_DIST
	CnC::dist_cnc_init<CnCaDContext> dinit;
#endif
	CnCaDContext ctxt;

	//Set up parameters
	const bool doKriging = true;
	const bool doCoMD = false;
	int dims[2] = {atoi(argv[1]), atoi(argv[2])};
	double dt[2] = {0.1, 0.1};
	double delta[2] = {1.0, 1.0};
	double gamma[3];
	gamma[0] = 0; //mom_gamma
	gamma[1] = gamma[0]; //strain_gamma
	gamma[2] = 0.1 * gamma[1];//en_gamma

	//Put the one item
	Global_Singleton_Item runConfig;
	strcpy(runConfig.redis_host, argv[4]);
	runConfig.doKriging = doKriging;
	runConfig.doCoMD = doCoMD;
	ctxt.globalItem.put(0, runConfig);

	unsigned int numSteps = atoi(argv[3]);

	//Initialize
	std::cout << "Initializing " << dims[0] << " by " << dims[1] << " grid" << std::endl;
	initEverything(doKriging, doCoMD, dims, dt, delta, gamma, argv[4]);
	std::cout << "Initialized" << std::endl;
	//Loop
	std::cout << "Running for " << numSteps << " iterations" << std::endl;
	for(unsigned int t = 0; t < numSteps; t++)
	{
		int nTasks;
		std::cout << t << ": Vising to Verifying" << std::endl;
		outputVTK(doKriging, doCoMD, dims, dt, delta, gamma, t, argv[4]);
		//Do a short circuit test
		{
			std::cout << t << ": First Flux" << std::endl;
			nTasks = prepFirstFlux(doKriging, doCoMD, dims, dt, delta, gamma, t, argv[4]);
			std::cout << t << ": Doing " << nTasks << " fluxes" << std::endl;
			parallelFor(t, 0, nTasks, ctxt);
			std::cout << t << ": Second Flux" << std::endl;
			nTasks = prepSecondFlux(doKriging, doCoMD, dims, dt, delta, gamma, t, argv[4]);
			std::cout << t << ": Doing " << nTasks << " fluxes" << std::endl;
			parallelFor(t, 1, nTasks, ctxt);
			std::cout << t << ": Third Flux" << std::endl;
			nTasks = prepThirdFlux(doKriging, doCoMD, dims, dt, delta, gamma, t, argv[4]);
			std::cout << t << ": Doing " << nTasks << " fluxes" << std::endl;
			parallelFor(t, 2, nTasks, ctxt);
			std::cout << t << ": Last Flux" << std::endl;
			nTasks = prepLastFlux(doKriging, doCoMD, dims, dt, delta, gamma, t, argv[4]);
			std::cout << t << ": Doing " << nTasks << " fluxes" << std::endl;
			parallelFor(t, 3, nTasks, ctxt);
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
