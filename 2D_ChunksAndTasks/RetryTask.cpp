#include "RetryTask.hpp"

#include "CoHMM_DaD.hpp"

CHT_TASK_TYPE_IMPLEMENTATION((RetryTask));
cht::ID RetryTask::execute(RetryChunk const &task)
{
	//Use task to call cloudRetry... the only part of this that matters
	retryCloudFlux(task.doKriging, task.doCoMD, task.step, task.phase, task.taskID, task.round, task.redisHost);
	//Return dummy chunk val
	CppBoolChunk resChunk;
	resChunk.retVal = true;
	cht::ChunkID cid_result = registerChunk(new CppBoolChunk(resChunk), cht::persistent);
	return cid_result;
}
