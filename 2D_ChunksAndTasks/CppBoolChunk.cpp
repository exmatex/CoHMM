#include "CppBoolChunk.hpp"

#include <cstring>

CHT_CHUNK_TYPE_IMPLEMENTATION((CppBoolChunk));
void CppBoolChunk::writeToBuffer(char * dataBuffer, size_t const bufferSize) const
{
	if(bufferSize != this->getSize())
	{
		throw std::runtime_error("Wrong buffer size to CppBoolChunk::writeToBuffer.");
	}
	else
	{
		memcpy(dataBuffer, &this->retVal, sizeof(bool)); 
	}
}

void CppBoolChunk::assignFromBuffer(char const * dataBuffer, size_t const bufferSize)
{
	if (bufferSize != this->getSize())
	{
		throw std::runtime_error("Wrong buffer size to CppBoolChunk::assign_from_buffer.");
	}
	else
	{
		memcpy(&this->retVal, dataBuffer, sizeof(bool));
	}
}

size_t CppBoolChunk::getSize() const
{
	return sizeof(bool);
}

size_t CppBoolChunk::memoryUsage() const
{
	return this->getSize();
}

