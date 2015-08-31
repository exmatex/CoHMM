#ifndef FLUXTASK_HPP
#define FLUXTASK_HPP

#include "chunks_and_tasks.h"

#include "FluxInChunk.hpp"
#include "FluxOutChunk.hpp"

struct FluxTask: public cht::Task
{
	cht::ID execute(FluxInChunk const &);
	CHT_TASK_INPUT((FluxInChunk));
	CHT_TASK_OUTPUT((FluxOutChunk));
	CHT_TASK_TYPE_DECLARATION;
};


#endif
