#ifndef RETRYTASK_HPP
#define RETRYTASK_HPP

#include "chunks_and_tasks.h"

#include "RetryChunk.hpp"
#include "CppBoolChunk.hpp"


struct RetryTask: public cht::Task
{
	cht::ID execute(RetryChunk const &);
	CHT_TASK_INPUT((RetryChunk));
	//Dummy output so we can return these
	CHT_TASK_OUTPUT((CppBoolChunk));
	CHT_TASK_TYPE_DECLARATION;
};


#endif
