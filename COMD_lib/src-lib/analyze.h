
#ifndef __ANALYZE_H_
#define __ANALYZE_H_

#include "CoMDTypes.h"
#include "CoMD_lib.h"

CoMD_return printTensor(int step, real_t* mat9);
CoMD_return printAverageTensor(int step, real_t* mat9);
CoMD_return printAverageFlux(int step, real_t* avFlux);
int calcStress(SimFlat* s);
void printThingsToFile(FILE* file, SimFlat* s, int iStep, double elapsedTime, real_t avStress, int counter);
void calcAverageTensor(real_t* avStress, real_t* mat9);
void calcAverageEnergyFlux(real_t* avFlux, real_t* flux);
void calcAveragePosition(SimFlat* s, real3 *avPos, int* counter);
void writeVtk(char* fileName, SimFlat* s);
void writeAverageVtk(char* fileName, SimFlat* s, real3* avPos, int counter);
void stressAutoCorrelation(SimFlat* s, real_t* mat9, real_t* g, real_t** stress, int* nTimes);
void initStressAutoCorr(real_t** g, real_t*** stress, int length);
void printStressAutoCorr(char* fileName, real_t* g, int* nTimes, int printRate, real_t dt);
#endif
