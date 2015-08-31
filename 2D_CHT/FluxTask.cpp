#include "FluxTask.hpp"

CHT_TASK_TYPE_IMPLEMENTATION((FluxTask));
cht::ID FluxTask::execute(FluxInChunk const &task)
{
	//Prep return value
	FluxOutChunk resChunk;
	//Use Input to call flux
    resChunk.output = fluxFn(task.doKriging, task.doCoMD, (FluxIn *)&task.input);
	//Return output
	cht::ChunkID cid_result = registerChunk(new FluxOutChunk(resChunk), cht::persistent);
	return cid_result;
}
