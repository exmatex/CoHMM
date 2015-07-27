#ifndef __DEFORMATION_H_
#define __DEFORMATION_H_

#include "CoMDTypes.h"
#include "CoMD_lib.h"

FILE* stressOut;
FILE* averageStressOut;

#if 0
CoMD_return printTensor(int step, real_t* mat9);
CoMD_return printAverageTensor(int step, real_t* mat9);
void calcAverageTensor(real_t* avStress, real_t* mat9);
#endif
void matVec3(real_t *mat, real_t *vec);
void matInv3x3 (real_t *in, real_t *out);
void forwardDeformation(SimFlat *s);
void reverseDeformation(SimFlat *s);

Deformation* initDeformation(SimFlat *sim, real_t* defGrad);
void destroyDeformation(Deformation* defInfo);

#endif
