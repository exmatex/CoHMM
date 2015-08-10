#ifndef CNCSHMEM_HPP
#define CNCSHMEM_HPP

#include <string>
#include <cstring>

#include "CoHMM_CnC.hpp"

const unsigned int maxKeyLength = 48;
const char rshmem_tag[] = "SWIFTT_TEST";

std::string buildSingleKey(int curStep, int curPhase, int ID, const char * tag);
std::string buildBlockKey(int curStep, int curPhase, int ID, int dimX, int dimY, const char * tag);
bool putNodes(Node * field, int dimX, int dimY, int curStep, int curPhase, CnCaDContext &c);
bool getNodes(Node * field, int dimX, int dimY, int curStep, int curPhase, CnCaDContext &c);
bool putFutures(FluxFuture * field, int dimX, int dimY, int curStep, int curPhase, CnCaDContext &c);
bool getFutures(FluxFuture * field, int dimX, int dimY, int curStep, int curPhase, CnCaDContext &c);

bool putTask(FluxIn * item, int curStep, int curPhase, int ID,  CnCaDContext &c);
bool putResult(FluxOut * item, int curStep, int curPhase, int ID,  CnCaDContext &c);

bool getTask(FluxIn * item, int curStep, int curPhase, int ID,   CnCaDContext &c);
bool getResult(FluxOut * item, int curStep, int curPhase, int ID,   CnCaDContext &c);




//Determine how many blocks the field has been chunked in to
unsigned int getNumBlocks(int dimX, int dimY);


#endif
