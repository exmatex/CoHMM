#include "FluxInChunk.hpp"

#include <cstring>

#include "chunks_and_tasks.h"

CHT_CHUNK_TYPE_IMPLEMENTATION((FluxInChunk));
void FluxInChunk::writeToBuffer(char * dataBuffer, size_t const bufferSize) const
{
	if(bufferSize != this->getSize())
	{
		throw std::runtime_error("Wrong buffer size to FluxInChunk::writeToBuffer.");
	}
	else
	{
		//Copy stuff, one at a time
		size_t offset = 0;
		memcpy(&dataBuffer[offset], &this->input, sizeof(FluxIn) );
		offset += sizeof(FluxIn);
		memcpy(&dataBuffer[offset], this->redisHost, sizeof(char)*RedisWrapper::MAX_HOST_LENGTH);
		offset += sizeof(char) * RedisWrapper::MAX_HOST_LENGTH;
		memcpy(&dataBuffer[offset], &this->taskID, sizeof(unsigned int));
		offset += sizeof(unsigned int);
		memcpy(&dataBuffer[offset], &this->doKriging, sizeof(bool));
		offset += sizeof(bool);
		memcpy(&dataBuffer[offset], &this->doCoMD, sizeof(bool));
	}
}

void FluxInChunk::assignFromBuffer(char const * dataBuffer, size_t const bufferSize)
{
	if (bufferSize != this->getSize())
	{
		throw std::runtime_error("Wrong buffer size to FluxInChunk::assign_from_buffer.");
	}
	else
	{
		//Copy stuff, one at a time
		size_t offset = 0;
		memcpy((void *)&this->input, &dataBuffer[offset], sizeof(FluxIn));
		offset += sizeof(FluxIn);
		memcpy((void *)this->redisHost, &dataBuffer[offset],  sizeof(char)*RedisWrapper::MAX_HOST_LENGTH);
		offset += sizeof(char) * RedisWrapper::MAX_HOST_LENGTH;
		memcpy((void *)&this->taskID, &dataBuffer[offset], sizeof(unsigned int));
		offset += sizeof(unsigned int);
		memcpy((void *)&this->doKriging, &dataBuffer[offset],  sizeof(bool));
		offset += sizeof(bool);
		memcpy((void *)&this->doCoMD, &dataBuffer[offset], sizeof(bool));
	}

}

size_t FluxInChunk::getSize() const
{
	size_t retSize = 0;
	retSize += sizeof(FluxIn);
	retSize += sizeof(char) * RedisWrapper::MAX_HOST_LENGTH;
	retSize += sizeof(unsigned int);
	retSize += sizeof(bool);
	retSize += sizeof(bool);
	return retSize;
}

size_t FluxInChunk::memoryUsage() const
{
	return this->getSize();
}
