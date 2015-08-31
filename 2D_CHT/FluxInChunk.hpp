#ifndef FLUXINCHUNK_HPP
#define FLUXINCHUNK_HPP

#include "chunks_and_tasks.h"

#include "RedisWrapper.hpp"
#include "2DKriging.hpp"

struct FluxInChunk: public cht::Chunk
{
	public:
		//Attributes
		FluxIn input;
		char redisHost[RedisWrapper::MAX_HOST_LENGTH];
		unsigned int taskID;
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
