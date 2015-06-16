#ifndef CPPBOOLCHUNK_HPP
#define CPPBOOLCHUNK_HPP

#include "chunks_and_tasks.h"

struct CppBoolChunk: public cht::Chunk
{
	public:
		//Attributes
		bool retVal;
		//Members
		void writeToBuffer(char * dataBuffer, size_t const bufferSize) const;
		void assignFromBuffer(char const * dataBuffer, size_t const bufferSize);
		size_t getSize() const;
		size_t memoryUsage() const;
	private:
		CHT_CHUNK_TYPE_DECLARATION;

};

#endif

