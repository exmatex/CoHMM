#ifndef RETRYCHUNK_HPP
#define RETRYCHUNK_HPP

#include "chunks_and_tasks.h"

#include "2D_ChunksAndTasks.hpp"

struct RetryChunk: public cht::Chunk
{
	public:
		//Attributes
		char redisHost[MAX_HOST_LENGTH];
		unsigned int taskID;
		unsigned int step;
		unsigned int phase;
		unsigned int round;
		bool doKriging;
		bool doCoMD;
		//Members
		void writeToBuffer(char * dataBuffer, size_t const bufferSize) const;
		void assignFromBuffer(char const * dataBuffer, size_t const bufferSize);
		size_t getSize() const;
		size_t memoryUsage() const;
	private:
		CHT_CHUNK_TYPE_DECLARATION;
};

#endif
