#include "FluxChunk.hpp"

#include <cstring>

#include "chunks_and_tasks.h"

CHT_CHUNK_TYPE_IMPLEMENTATION((FluxChunk));
void FluxChunk::writeToBuffer(char * dataBuffer, size_t const bufferSize) const
{
	if(bufferSize != this->getSize())
	{
		throw std::runtime_error("Wrong buffer size to FluxChunk::writeToBuffer.");
	}
	else
	{
		//Copy stuff, one at a time
		size_t offset = 0;
		memcpy(&dataBuffer[offset], this->redisHost, sizeof(char)*MAX_HOST_LENGTH);
		offset += sizeof(char) * MAX_HOST_LENGTH;
		memcpy(&dataBuffer[offset], &this->taskID, sizeof(unsigned int));
		offset += sizeof(unsigned int);
		memcpy(&dataBuffer[offset], &this->step, sizeof(unsigned int));
		offset += sizeof(unsigned int);
		memcpy(&dataBuffer[offset], &this->phase, sizeof(unsigned int));
		offset += sizeof(unsigned int);
		memcpy(&dataBuffer[offset], &this->doKriging, sizeof(bool));
		offset += sizeof(bool);
		memcpy(&dataBuffer[offset], &this->doCoMD, sizeof(bool));
	}
}

void FluxChunk::assignFromBuffer(char const * dataBuffer, size_t const bufferSize)
{
	if (bufferSize != this->getSize())
	{
		throw std::runtime_error("Wrong buffer size to FluxChunk::assign_from_buffer.");
	}
	else
	{
		//Copy stuff, one at a time
		size_t offset = 0;
		memcpy((void *)this->redisHost, &dataBuffer[offset],  sizeof(char)*MAX_HOST_LENGTH);
		offset += sizeof(char) * MAX_HOST_LENGTH;
		memcpy((void *)&this->taskID, &dataBuffer[offset], sizeof(unsigned int));
		offset += sizeof(unsigned int);
		memcpy((void *)&this->step, &dataBuffer[offset],  sizeof(unsigned int));
		offset += sizeof(unsigned int);
		memcpy((void *)&this->phase, &dataBuffer[offset], sizeof(unsigned int));
		offset += sizeof(unsigned int);
		memcpy((void *)&this->doKriging, &dataBuffer[offset],  sizeof(bool));
		offset += sizeof(bool);
		memcpy((void *)&this->doCoMD, &dataBuffer[offset], sizeof(bool));

	}

}

size_t FluxChunk::getSize() const
{
	size_t retSize = 0;
	retSize += sizeof(char) * MAX_HOST_LENGTH;
	retSize += sizeof(unsigned int);
	retSize += sizeof(unsigned int);
	retSize += sizeof(unsigned int);
	retSize += sizeof(bool);
	retSize += sizeof(bool);
	return retSize;
}

size_t FluxChunk::memoryUsage() const
{
	return this->getSize();
}


