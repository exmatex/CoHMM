#include <cassert>
#include <cstdio>
#include <cstring>

#include "redisShmem.hpp"

void buildKey(char * key, int curStep, int curPhase, int ID, const char * tag)
{
	sprintf(key, "%s:%s:%04d:%04d:%04d", tag, rshmem_tag, curStep, curPhase, ID);
}

//Determine how many blocks the field has been chunked in to
unsigned int getNumBlocks(int dimX, int dimY)
{
	unsigned int numBlocks = (dimX * dimY) / fieldBlockSize;
	if(fieldBlockSize % (dimX*dimY) != 0)
	{
		numBlocks++;
	}

    return numBlocks;
}

