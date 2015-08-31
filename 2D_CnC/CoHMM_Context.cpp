#include "CoHMM_Context.hpp"

#include "RedisWrapper.hpp"

CoHMMContext::CoHMMContext()
	:
	CnC::context<CoHMMContext>(),
	fluxTags(*this, "fluxTags"),
	fluxTasks(*this, "fluxTasks"),
	fluxItems(*this, "fluxItems"),
    fluxResults(*this, "fluxResults")
{
	//Indicate which tags correspond to which steps
	fluxTags.prescribes(fluxTasks, *this);

	//Indicate which step consumes which item (ha ha)
	fluxTasks.consumes(fluxItems);

    //Indicate which step produces which item
    fluxTasks.produces(fluxResults);
}

int Flux_Task::execute(const Flux_Tag &tag, CoHMMContext &c) const
{
	//Get The Item
	Flux_Item runConfig;
	c.fluxItems.get(tag, runConfig);
	//Get params
	unsigned int step = tag.first.first;
	unsigned int phase = tag.first.second;
	unsigned int task = tag.second;
	//Prep output
	Flux_Result res;
    //Potentially init RedisWrapper
    RedisWrapper::getContext(runConfig.redis_host);
	//Call
    res.result = fluxFn(runConfig.doKriging, runConfig.doCoMD, &runConfig.task);
    //Put item to result
	c.fluxResults.put(tag, res);
	//Return
	return CnC::CNC_Success;
}
