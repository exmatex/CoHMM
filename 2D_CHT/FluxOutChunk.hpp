#ifndef FLUXOUTCHUNK_HPP
#define FLUXOUTCHUNK_HPP

#include "chunks_and_tasks.h"

#include "2DKriging.hpp"

struct FluxOutChunk: public cht::Chunk
{
	public:
		//Attributes
		FluxOut output;
		//Members
		void writeToBuffer(char * dataBuffer, size_t const bufferSize) const;
		void assignFromBuffer(char const * dataBuffer, size_t const bufferSize);
		size_t getSize() const;
		size_t memoryUsage() const;
	private:
		CHT_CHUNK_TYPE_DECLARATION;
};

#endif
