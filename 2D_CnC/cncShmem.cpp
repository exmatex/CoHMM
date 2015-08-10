#include "cncShmem.hpp"

std::string buildBlockKey(int curStep, int curPhase, int ID, int dimX, int dimY, const char * tag)
{
	char key[maxKeyLength];
	sprintf(key, "%s:%s:%04d:%04d:%04d:%04d:%04d", tag, rshmem_tag, dimX, dimY, curStep, curPhase, ID);
	std::string retString(key);
	return retString;
}

std::string buildSingleKey(int curStep, int curPhase, int ID, const char * tag)
{
	char key[maxKeyLength];
	sprintf(key, "%s:%s:%04d:%04d:%04d", tag, rshmem_tag, curStep, curPhase, ID);
	std::string retString(key);
	return retString;
}

//Determine how many blocks the field has been chunked in to
unsigned int getNumBlocks(int dimX, int dimY)
{
	unsigned int numBlocks = (dimX * dimY) / fieldBlockSize;
	if((dimX*dimY) % fieldBlockSize != 0)
	{
		numBlocks++;
	}
    return numBlocks;
}

bool putNodes(Node * field, int dimX, int dimY, int curStep, int curPhase, CnCaDContext &c)
{
	unsigned int nBlocks = getNumBlocks(dimX, dimY);
	//Push the first nBlocks-1 blocks as they are easy
	for(unsigned int i = 0; i < (nBlocks - 1); i++)
	{
		//Build key
		std::string key = buildBlockKey(curStep, curPhase, i, dimX, dimY, "FIELD");
		//Build data
		NodeBlock block;
		memcpy(block.block, &field[i*fieldBlockSize], sizeof(Node)*fieldBlockSize);
		//Put to cnc
		c.blockItems.put(key, block);
	}
	//See if last block is full
	unsigned int i = nBlocks - 1;
	std::string key = buildBlockKey(curStep, curPhase, i, dimX, dimY, "FIELD");
	unsigned int lastBlock = (dimX*dimY) % fieldBlockSize;
	if( lastBlock == 0)
	{
		//It was, so same as above
		//Build data
		NodeBlock block;
		memcpy(block.block, &field[i*fieldBlockSize], sizeof(Node)*fieldBlockSize);
		//Put to cnc
		c.blockItems.put(key, block);
	}
	else
	{
		//It was not, so only copy what we need
		unsigned int lastBlock = (dimX*dimY) % fieldBlockSize;
		//Build data
		NodeBlock block;
		memcpy(block.block, &field[i*fieldBlockSize], sizeof(Node)*lastBlock);
		//Put to cnc
		c.blockItems.put(key, block);
	}
	//Success
	return true;
}

bool getNodes(Node * field, int dimX, int dimY, int curStep, int curPhase, CnCaDContext &c)
{
	unsigned int nBlocks = getNumBlocks(dimX, dimY);
	size_t blockSize = sizeof(Node) * fieldBlockSize;

	//Push the first nBlocks-1 blocks as they are easy
	for(unsigned int i = 0; i < (nBlocks - 1); i++)
	{
		//Build key
		std::string key = buildBlockKey(curStep, curPhase, i, dimX, dimY, "FIELD");
		//Get the data
		NodeBlock block;
		//Doing this super unsafe-like
		bool gotIt = false;
		while(gotIt == false)
		{
			gotIt = c.blockItems.unsafe_get(key, block);
		}
		memcpy(&field[i*fieldBlockSize], block.block, blockSize);
	}
	//Last one is slightly easier with GET
	//Do a redis pull
	int i = nBlocks - 1;
	std::string key = buildBlockKey(curStep, curPhase, i, dimX, dimY, "FIELD");
	NodeBlock block;
	bool gotIt = false;
	while(gotIt == false)
	{
		gotIt = c.blockItems.unsafe_get(key, block);
	}
	if((dimX*dimY) % fieldBlockSize != 0)
	{
		blockSize = ((dimX*dimY) % fieldBlockSize) * sizeof(Node);
	}
	memcpy(&field[i*fieldBlockSize], block.block, blockSize);
	//Success
	return true;
}

bool putFutures(FluxFuture * field, int dimX, int dimY, int curStep, int curPhase, CnCaDContext &c)
{
	unsigned int nBlocks = getNumBlocks(dimX, dimY);
	//Push the first nBlocks-1 blocks as they are easy
	for(unsigned int i = 0; i < (nBlocks - 1); i++)
	{
		//Build key
		std::string key = buildBlockKey(curStep, curPhase, i, dimX, dimY, "FUTS");
		//Build data
		FutureBlock block;
		memcpy(block.block, &field[i*fieldBlockSize], sizeof(FluxFuture)*fieldBlockSize);
		//Put to cnc
		c.futureItems.put(key, block);
	}
	//See if last block is full
	unsigned int i = nBlocks - 1;
	std::string key = buildBlockKey(curStep, curPhase, i, dimX, dimY, "FUTS");
	FutureBlock block;
	unsigned int lastBlock = (dimX*dimY) % fieldBlockSize;
	if( lastBlock == 0)
	{
		//It was, so same as above
		//Build data
		FutureBlock block;
		memcpy(block.block, &field[i*fieldBlockSize], sizeof(FluxFuture)*fieldBlockSize);
		//Put to cnc
		c.futureItems.put(key, block);
	}
	else
	{
		//It was not, so only copy what we need
		unsigned int lastBlock = (dimX*dimY) % fieldBlockSize;
		//Build data
		FutureBlock block;
		memcpy(block.block, &field[i*fieldBlockSize], sizeof(FluxFuture)*lastBlock);
		//Put to cnc
		c.futureItems.put(key, block);
	}
	//Success
	return true;
}

bool getFutures(FluxFuture * field, int dimX, int dimY, int curStep, int curPhase, CnCaDContext &c)
{
	unsigned int nBlocks = getNumBlocks(dimX, dimY);
	size_t blockSize = sizeof(FluxFuture) * fieldBlockSize;

	//Push the first nBlocks-1 blocks as they are easy
	for(unsigned int i = 0; i < (nBlocks - 1); i++)
	{
		//Build key
		std::string key = buildBlockKey(curStep, curPhase, i, dimX, dimY, "FUTS");
		//Get the data
		FutureBlock block;
		//Doing this super unsafe-like
		bool gotIt = false;
		while(gotIt == false)
		{
			gotIt = c.futureItems.unsafe_get(key, block);
		}
		memcpy(&field[i*fieldBlockSize], block.block, blockSize);
	}
	//Last one is slightly easier with GET
	//Do a redis pull
	int i = nBlocks - 1;
	std::string key = buildBlockKey(curStep, curPhase, i, dimX, dimY, "FUTS");
	FutureBlock block;
	bool gotIt = false;
	while(gotIt == false)
	{
		gotIt = c.futureItems.unsafe_get(key, block);
	}
	if((dimX*dimY) % fieldBlockSize != 0)
	{
		blockSize = ((dimX*dimY) % fieldBlockSize) * sizeof(FluxFuture);
	}
	memcpy(&field[i*fieldBlockSize], block.block, blockSize);
	//Success
	return true;
}




bool putTask(FluxIn * item, int curStep, int curPhase, int ID,  CnCaDContext &c)
{
	//Get Key
	std::string key = buildSingleKey(curStep, curPhase, ID,  "TASK");
	//Put data
	c.taskItems.put(key, *item);
	//Success
	return true;
}

bool putResult(FluxOut * item, int curStep, int curPhase, int ID,  CnCaDContext &c)
{
	//Get Key
	std::string key = buildSingleKey(curStep, curPhase, ID,  "RESULT");
	//Put data
	c.resultItems.put(key, *item);
	//Success
	return true;
}

bool getTask(FluxIn * item, int curStep, int curPhase, int ID,   CnCaDContext &c)
{
	//Get Key
	std::string key = buildSingleKey(curStep, curPhase, ID,  "TASK");
	FluxIn task;
	bool gotIt = false;
	while(gotIt == false)
	{
		gotIt = c.taskItems.unsafe_get(key, task);
	}
	*item = task;
	//Success
	return true;
}

bool getResult(FluxOut * item, int curStep, int curPhase, int ID,   CnCaDContext &c)
{
	//Get Key
	std::string key = buildSingleKey(curStep, curPhase, ID,  "RESULT");
	FluxOut task;
	bool gotIt = false;
	while(gotIt == false)
	{
		gotIt = c.resultItems.unsafe_get(key, task);
	}
	*item = task;
	//Success
	return true;
}
