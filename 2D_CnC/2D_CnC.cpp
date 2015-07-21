#include <iostream>
#include <cstdlib>

#include "CoHMM_DaD.hpp"

#include "2D_CnC.hpp"

Flux_Tag::Flux_Tag(int step, int phase, int task)
{
	this->step = step;
	this->phase = phase;
	this->task = task;
}

Flux_Tag::Flux_Tag()
{
	this->step = -1;
	this->phase = -1;
	this->task = -1;
}

Retry_Tag::Retry_Tag(int step, int phase, int task, int round)
{
	this->step = step;
	this->phase = phase;
	this->task = task;
	this->round = round;
}

Retry_Tag::Retry_Tag()
{
	this->step = -1;
	this->phase = -1;
	this->task = -1;
	this->round = -1;
}

DaDContext::DaDContext()
	:
	CnC::context<DaDContext>(),
	fluxTags(*this, "fluxTags"),
	retryTags(*this, "retryTags"),
	fluxTask(*this, "fluxTasks"),
	retryTask(*this, "retryTasks"),
	globalItem(*this, "globalItem")
{
	//Indicate which tags correspond to which steps
	fluxTags.prescribes(fluxTask, *this);
	retryTags.prescribes(retryTask, *this);

	//Indicate which step consumes which item (ha ha)
	fluxTask.consumes(globalItem);
	retryTask.consumes(globalItem);
}

int Flux_Task::execute(const Flux_Tag &tag, DaDContext &c) const
{
	//Get The Item
	Flux_Item runConfig;
	c.globalItem.get(0, runConfig);
	//Get params
	unsigned int step = tag.step;
	unsigned int phase = tag.phase;
	unsigned int task = tag.task;
	//Call
	cloudFlux(runConfig.doKriging, runConfig.doCoMD, step, phase, task, runConfig.redis_host);
	//Return
	return CnC::CNC_Success;
}

int Retry_Task::execute(const Retry_Tag &tag, DaDContext &c) const
{
	//Get The Item
	Flux_Item runConfig;
	c.globalItem.get(0, runConfig);
	//Get params
	unsigned int step = tag.step;
	unsigned int phase = tag.phase;
	unsigned int task = tag.task;
	unsigned int curRound = tag.round;
	//Call
	retryCloudFlux(runConfig.doKriging, runConfig.doCoMD, step, phase, task, curRound, runConfig.redis_host);
	//Return
	return CnC::CNC_Success;
}

void parallelFor(unsigned int step, unsigned int phase, unsigned int nTasks, DaDContext &ctxt)
{
	for(unsigned int i = 0; i < nTasks; i++)
	{
		//Make task tag
		Flux_Tag tag(step, phase, i);
		//Put task tag
		ctxt.fluxTags.put(tag);
	}
	ctxt.wait();
	return;
}

void iterativeRetry(int * dims, unsigned int step, unsigned int phase, char * redis_host, DaDContext &ctxt)
{
	int curRound = 0;
	int nTasks = checkStepForFaults(dims, step, phase, curRound, redis_host);
	while(nTasks != 0)
	{
		std::cout << step << ": Redoing " << nTasks << " Tasks" << std::endl;
		for(unsigned int i = 0; i < nTasks; i++)
		{
			//Make tag
			Retry_Tag tag(step, phase, i, curRound);
			//Put tag
			ctxt.retryTags.put(tag);
		}
		ctxt.wait();
		//See if we are done
		curRound++;
		nTasks = checkStepForFaults(dims, step, phase, curRound, redis_host);
	}
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
	CnC::dist_cnc_init<DaDContext> dinit;
#endif
	DaDContext ctxt;
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

	//Put the one item
	Flux_Item runConfig;
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
			parallelFor(t, 0, nTasks, ctxt);
			std::cout << t << ": Checking First Flux" << std::endl;
			if(fineGrainFT == true)
			{
				iterativeRetry(dims, t, 0, argv[4], ctxt);
			}
			std::cout << t << ": Second Flux" << std::endl;
			nTasks = prepSecondFlux(doKriging, doCoMD, dims, dt, delta, gamma, t, argv[4]);
			std::cout << t << ": Doing " << nTasks << " fluxes" << std::endl;
			parallelFor(t, 1, nTasks, ctxt);
			std::cout << t << ": Checking Second Flux" << std::endl;
			if(fineGrainFT == true)
			{
				iterativeRetry(dims, t, 1, argv[4], ctxt);
			}
			std::cout << t << ": Third Flux" << std::endl;
			nTasks = prepThirdFlux(doKriging, doCoMD, dims, dt, delta, gamma, t, argv[4]);
			std::cout << t << ": Doing " << nTasks << " fluxes" << std::endl;
			parallelFor(t, 2, nTasks, ctxt);
			std::cout << t << ": Checking Third Flux" << std::endl;
			if(fineGrainFT == true)
			{
				iterativeRetry(dims, t, 2, argv[4], ctxt);
			}
			std::cout << t << ": Last Flux" << std::endl;
			nTasks = prepLastFlux(doKriging, doCoMD, dims, dt, delta, gamma, t, argv[4]);
			std::cout << t << ": Doing " << nTasks << " fluxes" << std::endl;
			parallelFor(t, 3, nTasks, ctxt);
			std::cout << t << ": Checking Last Flux" << std::endl;
			if(fineGrainFT == true)
			{
				iterativeRetry(dims, t, 3, argv[4], ctxt);
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
