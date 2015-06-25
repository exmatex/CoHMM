/* file includes analyze functions to generate output of physical quantities */
/* for the lib version do NOT introduce global variables */
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <memory.h>
#include <string.h>
#include <unistd.h>

#include "CoMDTypes.h"
#include "memUtils.h"
#include "constants.h"
#include "eam.h"
#include <assert.h>
#include <math.h>

#include "CoMD_lib.h"
#include "ljForce.h"

#define POT_SHIFT 1.0

CoMD_return printTensor(int step, real_t* mat9)
{
	int index = step;
	CoMD_return retVal;

	retVal.simIter = step;
	retVal.stressXX = mat9[0];
	retVal.stressXY = mat9[1];
	retVal.stressXZ = mat9[2];
	retVal.stressYX = mat9[3];
	retVal.stressYY = mat9[4];
	retVal.stressYZ = mat9[5];
	retVal.stressZX = mat9[6];
	retVal.stressZY = mat9[7];
	retVal.stressZZ = mat9[8];
	
/*
   //fprintf(stressOut, "[ %g %g %g] \n [%g %g %g] \n [%g %g %g] \n",
   fprintf(stdout, "%6d % 12.6g % 12.6g % 12.6g % 12.6g % 12.6g % 12.6g % 12.6g % 12.6g % 12.6g \n",
         step,
         mat9[0], mat9[1], mat9[2],
         mat9[3], mat9[4], mat9[5],
         mat9[6], mat9[7], mat9[8]);
*/
	return retVal;
}
//averaged stress tensor -> communicated to the macrosolver 
CoMD_return printAverageTensor(int step, real_t* mat9)
{

	CoMD_return retVal;
	//int index = step;

	retVal.simIter = step;
	retVal.stressXX = mat9[0]/step;
	retVal.stressXY = mat9[1]/step;
	retVal.stressXZ = mat9[2]/step;
	retVal.stressYX = mat9[3]/step;
	retVal.stressYY = mat9[4]/step;
	retVal.stressYZ = mat9[5]/step;
	retVal.stressZX = mat9[6]/step;
	retVal.stressZY = mat9[7]/step;
	retVal.stressZZ = mat9[8]/step;

	return retVal;
}
//averaged stress tensor -> communicated to the macrosolver 
CoMD_return printAverageFlux(int step, real_t* avFlux)
{

	CoMD_return retVal;
	//int index = step;

	retVal.energyDensX = avFlux[0]/step;
	retVal.energyDensY = avFlux[1]/step;
	retVal.energyDensZ = avFlux[2]/step;

	return retVal;
}
// averaging of the stress tensor to reduce fluct
//in fact just a sum of the stresses, dividing by the number of measurement steps done in printAverageTensor 
void calcAverageTensor(real_t* avStress, real_t* mat9)
{
    for(int i=0; i<9; ++i){
        avStress[i] += mat9[i];
    }
}
//averaging of the energy flux to reduce fluct
//in fact just a sum of the stresses, dividing by the number of measurement steps done in printAverageTensor 
void calcAverageEnergyFlux(real_t* avFlux, real_t* flux)
{
    for(int i=0; i<3; ++i){
        avFlux[i] += flux[i];
    }
}
//averaged stress for LJ pot only!!!
int calcStressLJ(SimFlat* s)
{
   LjPotential* pot = (LjPotential *) s->pot;
   real_t sigma = pot->sigma;
   real_t epsilon = pot->epsilon;
   real_t rCut = pot->cutoff;
   real_t rCut2 = rCut*rCut;

   // zero forces and energy
   int fSize = s->boxes->nTotalBoxes*MAXATOMS;

   real_t s6 = sigma*sigma*sigma*sigma*sigma*sigma;

   real_t rCut6 = s6 / (rCut2*rCut2*rCut2);
   real_t eShift = POT_SHIFT * rCut6 * (rCut6 - 1.0);

   //init virial stress here
   real_t localStress[9];
   for (int m=0; m<9; m++) 
   {
      localStress[m] = 0.0;
   }

   int nbrBoxes[27];
   // loop over local boxes
   for (int iBox=0; iBox<s->boxes->nLocalBoxes; iBox++)
   {
      int nIBox = s->boxes->nAtoms[iBox];
      if ( nIBox == 0 ) continue;
      int nNbrBoxes = getNeighborBoxes(s->boxes, iBox, nbrBoxes);
      // loop over neighbors of iBox
      for (int jTmp=0; jTmp<nNbrBoxes; jTmp++)
      {
         int jBox = nbrBoxes[jTmp];

         assert(jBox>=0);

         int nJBox = s->boxes->nAtoms[jBox];
         if ( nJBox == 0 ) continue;

         // calculate energy contribution based on whether
         // the neighbor box is local or remote
         real_t localScale;
         if (jBox < s->boxes->nLocalBoxes)
            localScale = 1.0;
         else
            localScale = 0.5;

         // loop over atoms in iBox
         for (int iOff=iBox*MAXATOMS,ii=0; ii<nIBox; ii++,iOff++)
         {
            int iId = s->atoms->gid[iOff];
            // kinetic energy contribution to the virial stress
            int iSpecies = s->atoms->iSpecies[iOff];
            real_t invMass = 1.0/s->species[iSpecies].mass;
#if 1
            for (int ii=0; ii<3; ii++) 
            {
               for (int jj=0; jj<3; jj++) 
               {
                  int m = 3*ii + jj;
                  localStress[m] -= s->atoms->p[iOff][ii]*s->atoms->p[iOff][jj]*invMass;
               }
            }
#endif
            // loop over atoms in jBox
            for (int jOff=MAXATOMS*jBox,ij=0; ij<nJBox; ij++,jOff++)
            {
               real_t dr[3];
               int jId = s->atoms->gid[jOff];  
               if (jBox < s->boxes->nLocalBoxes && jId <= iId )
                  continue; // don't double count local-local pairs.
               real_t r2 = 0.0;
               for (int m=0; m<3; m++)
               {
                  dr[m] = s->atoms->r[iOff][m]-s->atoms->r[jOff][m];
                  r2+=dr[m]*dr[m];
               }

               if ( r2 > rCut2) continue;

               // Important note:
               // from this point on r actually refers to 1.0/r
               r2 = 1.0/r2;
               real_t r6 = s6 * (r2*r2*r2);

               // different formulation to avoid sqrt computation
               real_t fr = - 4.0*epsilon*r6*r2*(12.0*r6 - 6.0);

               // stress computation
               for (int ii=0; ii<3; ii++) 
               {
                  for (int jj=0; jj<3; jj++) 
                  {
                     int m = 3*ii + jj;
                     localStress[m] -= localScale*fr*dr[ii]*dr[jj];
                  }
               }
            } // loop over atoms in jBox
         } // loop over atoms in iBox
      } // loop over neighbor boxes
   } // loop over local boxes in system

   // renormalize stress
   for (int m=0; m<9; m++) 
   {
      // divide by the undeformed volume
      localStress[m] = localStress[m]/s->defInfo->globalVolume;
      // normalize the volume accounting for the strain
      // write to the global array
      s->defInfo->stress[m] = localStress[m]/s->defInfo->strain[9];
   }

   return 0;
}
/* so far the stress is calcuted within the eam force calculation function (see eam.c) but it may be moved over here  */

void printThingsToFile(FILE* file, SimFlat* s, int iStep, double elapsedTime, real_t avStress, int counter)
{
   // keep track previous value of iStep so we can calculate number of steps.
   static int iStepPrev = -1; 

   int nEval = iStep - iStepPrev; // gives nEval = 1 for zeroth step.
   iStepPrev = iStep;
   
   //if (! printRank() )
   //   return;
   //   FIXME print info only at first print out
#if 0
   if (counter == 0)
   {
      fprintf(file, 
       "#  Loop   Time(fs)       Total Energy   Potential Energy     Kinetic Energy  	Temperature  	Stress      Average Stress        # Atoms\n");
   }
#endif
   real_t time = iStep*s->dt;
   real_t eTotal = (s->ePotential+s->eKinetic) / s->atoms->nGlobal;
   //real_t eTotal = (s->ePotential+s->eKinetic+s->eKineticCom) / s->atoms->nGlobal;
   real_t eK = s->eKinetic / s->atoms->nGlobal;
   real_t eU = s->ePotential / s->atoms->nGlobal;
   real_t Temp = (s->eKinetic / s->atoms->nGlobal) / (kB_eV * 1.5);
   real_t stressXX = s->defInfo->stress[0];

   double timePerAtom = 1.0e6*elapsedTime/(double)(nEval*s->atoms->nLocal);

   if(counter>0){
   fprintf(file, " %6d %10.2f %18.12f %18.12f %18.12f %18.12f %18.12f %18.12f %12d\n",
           iStep, time, eTotal, eU, eK, Temp, stressXX, avStress/counter, s->atoms->nGlobal);
   }else{
   fprintf(file, " %6d %10.2f %18.12f %18.12f %18.12f %18.12f %18.12f %18.12f %12d\n",
           iStep, time, eTotal, eU, eK, Temp, stressXX, avStress, s->atoms->nGlobal);

   }
}

void writeVtk(char* fileName, SimFlat* s){

   //Temp output file for physical data
   //char fileName[255];
   //sprintf(fileName,"energy%i_%i.dat",cmd.calls, cmd.rank);
   real3 pos[s->atoms->nGlobal];
   
   FILE *partOut = fopen(fileName, "w");
   if (partOut == NULL) {
       printf("Error writing datafile < %s >!\n", fileName);
       exit(0);
   }
   // This is the CoMD main loop
    fprintf(partOut, "# vtk DataFile Version 2.0\nparticles\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS %i floats\n", s->atoms->nGlobal);

   for (int iBox=0; iBox<s->boxes->nLocalBoxes; ++iBox)
   {
      for (int iOff=MAXATOMS*iBox, ii=0; ii<s->boxes->nAtoms[iBox]; ++ii, ++iOff)
      {
        pos[s->atoms->gid[iOff]][0] = s->atoms->r[iOff][0];
        pos[s->atoms->gid[iOff]][1] = s->atoms->r[iOff][1];
        pos[s->atoms->gid[iOff]][2] = s->atoms->r[iOff][2];
      }
   }
   for (int iBox=0; iBox<s->boxes->nLocalBoxes; ++iBox)
   {
      for (int iOff=MAXATOMS*iBox, ii=0; ii<s->boxes->nAtoms[iBox]; ++ii, ++iOff)
      {
        fprintf(partOut, "%g %g %g\n",pos[s->atoms->gid[iOff]][0], pos[s->atoms->gid[iOff]][1], pos[s->atoms->gid[iOff]][2]);
      }
   }
   fclose(partOut);
}

void writeAverageVtk(char* fileName, SimFlat* s, real3* avPos, int counter){

   //Temp output file for physical data
   //char fileName[255];
   //sprintf(fileName,"energy%i_%i.dat",cmd.calls, cmd.rank);
   printf("counter%u\n", counter);
   FILE *partOut = fopen(fileName, "w");
   if (partOut == NULL) {
       printf("Error writing datafile < %s >!\n", fileName);
       exit(0);
   }
   // This is the CoMD main loop
    fprintf(partOut, "# vtk DataFile Version 2.0\nparticles\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS %i floats\n", s->atoms->nGlobal);

   for (int iBox=0; iBox<s->boxes->nLocalBoxes; ++iBox)
   {
      for (int iOff=MAXATOMS*iBox, ii=0; ii<s->boxes->nAtoms[iBox]; ++ii, ++iOff)
      {
        fprintf(partOut, "%g %g %g\n",avPos[s->atoms->gid[iOff]][0]/counter, avPos[s->atoms->gid[iOff]][1]/counter, avPos[s->atoms->gid[iOff]][2]/counter);
      }
   }
   fclose(partOut);
}
//simple averaging of the postions to reduce fluct for vtk visu
void calcAveragePosition(SimFlat* s, real3* avPos, int* counter)
{
   for (int iBox=0; iBox<s->boxes->nLocalBoxes; ++iBox)
   {
      for (int iOff=MAXATOMS*iBox, ii=0; ii<s->boxes->nAtoms[iBox]; ++ii, ++iOff)
      {
        avPos[s->atoms->gid[iOff]][0] += s->atoms->r[iOff][0];
        avPos[s->atoms->gid[iOff]][1] += s->atoms->r[iOff][1];
        avPos[s->atoms->gid[iOff]][2] += s->atoms->r[iOff][2];
      }
   }
   *counter+=1;
}

void initStressAutoCorr(real_t** g, real_t*** stress, int length){
  g[0] = (real_t*)comdMalloc(length*sizeof(real_t));
  stress[0] = (real_t**)comdMalloc(9*sizeof(real_t*));
  for (int i=0; i<9; ++i){
    stress[0][i] = (real_t*)comdMalloc(length*sizeof(real_t));
    for (int j=0; j<length; ++j){
        //roberts bootstraps
      stress[0][i][j] = 0.0;
    }
  }
}

void stressAutoCorrelation(SimFlat* s, real_t* mat9, real_t* g, real_t** stress, int* nTimes){

  *nTimes += 1;
  int steps = *nTimes;
  real_t temperature = (s->eKinetic/s->atoms->nGlobal)/kB_eV/1.5;
  //
  //	stressXX = mat9[0];
  //	stressXY = mat9[1];
  //	stressXZ = mat9[2];
  //	stressYX = mat9[3];
  //	stressYY = mat9[4];
  //	stressYZ = mat9[5];
  //	stressZX = mat9[6];
  //	stressZY = mat9[7];
  //	stressZZ = mat9[8];
  stress[0][steps-1] = mat9[0];
  stress[1][steps-1] = mat9[1];
  stress[2][steps-1] = mat9[2];
  stress[3][steps-1] = mat9[3];
  stress[4][steps-1] = mat9[4];
  stress[5][steps-1] = mat9[5];
  stress[6][steps-1] = mat9[6];
  stress[7][steps-1] = mat9[7];
  stress[8][steps-1] = mat9[8];

//calc summs over time
  real_t stressAvXY = 0;
  real_t stressAvYZ = 0;
  real_t stressAvZX = 0;
  real_t meanstressAvXY = 0;
  real_t meanstressAvYZ = 0;
  real_t meanstressAvZX = 0;
  real_t nAvXY = 0;
  real_t nAvXZ = 0;
  real_t nAvYZ = 0;
  real_t meannAvXY = 0;
  real_t meannAvXZ = 0;
  real_t meannAvYZ = 0;


  for (int time=0; time<steps; ++time) {
    meanstressAvXY += stress[1][time];
    meanstressAvYZ += stress[5][time];
    meanstressAvZX += stress[6][time];
    meannAvXY += (stress[0][time] - stress[4][time]);
    meannAvXZ += (stress[0][time] - stress[8][time]);
    meannAvYZ += (stress[4][time] - stress[8][time]);
  }
  meanstressAvXY /= steps;
  meanstressAvYZ /= steps;
  meanstressAvZX /= steps;
  meannAvXY /= steps; 
  meannAvXZ /= steps; 
  meannAvYZ /= steps; 
  for (int time=0; time<steps; ++time) {
    stressAvXY += (stress[1][time] - meanstressAvXY)*(stress[1][0] - meanstressAvXY);
    stressAvYZ += (stress[5][time] - meanstressAvYZ)*(stress[5][0] - meanstressAvYZ);
    stressAvZX += (stress[6][time] - meanstressAvZX)*(stress[6][0] - meanstressAvZX);
    nAvXY += ((stress[0][time] - stress[4][time]) - meannAvXY)*((stress[0][0] - stress[4][0]) - meannAvXY);
    nAvXZ += ((stress[0][time] - stress[8][time]) - meannAvXZ)*((stress[0][0] - stress[8][0]) - meannAvXZ);
    nAvYZ += ((stress[4][time] - stress[8][time]) - meannAvYZ)*((stress[4][0] - stress[8][0]) - meannAvYZ);
  }
//normalize by number of measurements
  stressAvXY /= steps;
  stressAvYZ /= steps;
  stressAvZX /= steps;
  nAvXY /= steps; 
  nAvXZ /= steps; 
  nAvYZ /= steps; 
  g[steps-1] = s->defInfo->globalVolume/(5*kB_eV*temperature)*(stressAvXY + stressAvYZ + stressAvZX);
  //g[steps-1] = s->defInfo->globalVolume/(5*kB_eV*temperature)*(stressAvXY - (meanstressAvXY*meanstressAvXY) + stressAvYZ - (meanstressAvYZ*meanstressAvYZ) + stressAvZX -(meanstressAvZX*meanstressAvZX));
  g[steps-1] += s->defInfo->globalVolume/(30*kB_eV*temperature)*(nAvXY + nAvXZ + nAvYZ);
  //g[steps-1] += s->defInfo->globalVolume/(30*kB_eV*temperature)*(nAvXY - (meannAvXY*meannAvXY) + nAvXZ - (meannAvXZ*meannAvXZ) + nAvYZ - (meannAvYZ*meannAvYZ));
  //printf("autocorr %g step %i\n", g[steps-1], steps);

}

void printStressAutoCorr(char* fileName, real_t* g, int* nTimes, int printRate, double dt){

  int steps = *nTimes;
  FILE *corrOut = fopen(fileName, "w");
  if (corrOut == NULL) {
      printf("Error writing datafile < %s >!\n", fileName);
      exit(0);
  }
  for (int time=0; time<steps; ++time) {
    fprintf(corrOut, "%g %g\n", time*printRate*dt, g[time]); 
  }
  fclose(corrOut);
}
