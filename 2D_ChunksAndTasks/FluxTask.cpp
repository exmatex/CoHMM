#include "FluxTask.hpp"

#include "CoHMM_DaD.hpp"

CHT_TASK_TYPE_IMPLEMENTATION((FluxTask));
cht::ID FluxTask::execute(FluxChunk const &task)
{
	//Use task to call cloudFlux... the only part of this that matters
	cloudFlux(task.doKriging, task.doCoMD, task.step, task.phase, task.taskID, task.redisHost); 
	//Return dummy chunk val
	CppBoolChunk resChunk;
	resChunk.retVal = true;
	cht::ChunkID cid_result = registerChunk(new CppBoolChunk(resChunk), cht::persistent);
	return cid_result;
}

