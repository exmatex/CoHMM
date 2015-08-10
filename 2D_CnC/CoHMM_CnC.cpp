#include <iostream>
#include <cstdlib>

#include "CoHMM_CnC.hpp"


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
	unsigned int step = tag.step;
	unsigned int phase = tag.phase;
	unsigned int task = tag.task;
	//Call
	cloudFlux(runConfig.doKriging, runConfig.doCoMD, step, phase, task, runConfig.redis_host);
	//Return
	return CnC::CNC_Success;


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

	return 0;
}
