#ifndef FLUXTASK_HPP
#define FLUXTASK_HPP

#include "chunks_and_tasks.h"

#include "FluxChunk.hpp"
#include "CppBoolChunk.hpp"


struct FluxTask: public cht::Task
{
	cht::ID execute(FluxChunk const &);
	CHT_TASK_INPUT((FluxChunk));
	//Dummy output so we can return these
	CHT_TASK_OUTPUT((CppBoolChunk));
	CHT_TASK_TYPE_DECLARATION;
};


#endif

