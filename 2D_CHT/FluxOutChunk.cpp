#include "FluxOutChunk.hpp"

#include <cstring>

#include "chunks_and_tasks.h"

CHT_CHUNK_TYPE_IMPLEMENTATION((FluxOutChunk));
void FluxOutChunk::writeToBuffer(char * dataBuffer, size_t const bufferSize) const
{
	if(bufferSize != this->getSize())
	{
		throw std::runtime_error("Wrong buffer size to FluxOutChunk::writeToBuffer.");
	}
	else
	{
		//Copy stuff, one at a time
		size_t offset = 0;
		memcpy(&dataBuffer[offset], &this->output, sizeof(FluxOut));
	}
}

void FluxOutChunk::assignFromBuffer(char const * dataBuffer, size_t const bufferSize)
{
	if (bufferSize != this->getSize())
	{
		throw std::runtime_error("Wrong buffer size to FluxOutChunk::assign_from_buffer.");
	}
	else
	{
		//Copy stuff, one at a time
		size_t offset = 0;
		memcpy((void *)&this->output, &dataBuffer[offset], sizeof(FluxOut));
	}

}

size_t FluxOutChunk::getSize() const
{
	size_t retSize = 0;
	retSize += sizeof(FluxOut);
	return retSize;
}

size_t FluxOutChunk::memoryUsage() const
{
	return this->getSize();
}
