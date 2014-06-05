
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <math.h>
#include <string.h>



#define GNUPLOT "/usr/bin/gnuplot -persist"

//Simulation parameters

int xp;

int microResponse, nSteps,hmmPrintRate, dimX, CoMDnx, CoMDny, CoMDnz, CoMDnSteps,interpolationType,plot, genuine;
double rho0, dx, dt, initialDeformation, initialEnergy, energyPerturbation, deformationStart, CoMDdt, deformationEnd, latticeParameter,refinementThreshold, coarseningThreshold;



void plotError(FILE *gp, double stressErrorMax[nSteps]) {
  int i,j;
  fprintf(gp, "set xrange [0:%i]\n",nSteps);
  fprintf(gp, "set xlabel 'Continuum time step \n");
  fprintf(gp, "set ylabel 'Average relative stress error'\n");
  fprintf(gp, "plot '-' u 1:2 w l notitle\n\n");
  for (j = 0; j < 1; j++) { // gnuplot forces us to repeat ourselves                                            
    for (i = 0; i < nSteps; i++) {
      fprintf(gp, "%i %.16lf\n",i,stressErrorMax[i]);
    }
    fprintf(gp, "e\n");
  }
  fflush(gp);
}

void plotErrorPos(FILE *gp2, double stressErrorPos[dimX]) {
  int i,j;
  
  fprintf(gp2, "set xrange [0:%i]\n",dimX);
  fprintf(gp2, "set xlabel 'Continuum grid' \n");
  fprintf(gp2, "set yrange [0.0:0.000000003]\n");
  fprintf(gp2, "set ylabel 'MSress error'\n");
  fprintf(gp2, "plot '-' u 1:2 w lp notitle\n\n");
  for (j = 0; j < 1; j++) { // gnuplot forces us to repeat ourselves                                                         
    for (i = 0; i < dimX; i++) {
      fprintf(gp2, "%i %f\n",i,stressErrorPos[i]);
    }
    fprintf(gp2, "e\n");
  }
  fflush(gp2);
}

void readInput(int xp){

  char inputname[1000];
  sprintf(inputname,"src/input%i.txt",xp);
  FILE *input = fopen(inputname, "r");
  if (input==NULL) {
      printf("Could not find an input file! Closing :-( \n", GNUPLOT);
      exit(0);
  }
  
  rewind (input);
  fscanf (input, "CoMD 1 Analytical 0 %d\n",&microResponse);
  //fscanf (input, "Wave speed for the analytical response %lf\n",&c);
  fscanf (input, "Reference material's density %lf\n",&rho0);
  fscanf (input, "Number of HMM time steps %d\n", &nSteps);
  fscanf (input, "Number of HMM grid points %d\n", &dimX);
  fscanf (input, "HMM print rate %d\n",&hmmPrintRate);
  fscanf (input, "HMM spatial resolution dx in A %lf\n", &dx);
  fscanf (input, "HMM time resolution dt in fs %lf\n", &dt);
  fscanf (input, "Initial deformation %lf\n", &initialDeformation);
  fscanf (input, "Initial energy density where undeformed in eV per A cube %lf\n", &initialEnergy);
  fscanf (input, "Energy added with deformation in eV per A cube %lf\n", &energyPerturbation);
  fscanf (input, "Deformation start and end in percentage of HMM mesh %lf %lf\n", &deformationStart, &deformationEnd);
  fscanf (input, "CoMD parameters:\nLattice parameter in A %lf\nNumber of lattices nx ny nz %d %d %d\nTime resolution dt in fs %lf\nNumber of MD time steps %d\n", &latticeParameter, &CoMDnx,&CoMDny,&CoMDnz,&CoMDdt,&CoMDnSteps);
  fscanf (input, "Spatial adaptive sampling parameters:\nRefinement threshold %lf\nCoarsening threshold %lf\nInterpolation: linear 0 Akima splines 1 %d\nPlot the CoHMM trajectory 0 no 1 yes %d\n", &refinementThreshold, &coarseningThreshold,&interpolationType,&plot);
  fscanf (input, "Avoid using the adaptive sampling 1 yes 0 no %d\n", &genuine);
  

  coarseningThreshold=(initialDeformation)*coarseningThreshold/dx;

  fclose(input);
}


int main(int argc, char **argv) {
  int i,j,k,l;


  if(argc != 2)
    {
      fprintf(stderr, "./error.o (input number)\n");
      return 1;
    }
  

  xp=  coarseningThreshold =atoi(argv[1]);

  readInput(xp);

  FILE *gp;
  gp = popen(GNUPLOT,"w");

  fprintf(gp, "set terminal postscript \n");
  fprintf(gp, "set output \'data/simulation%i/error/errorTrajectory_CoMD%i_dimX%i_nSteps%i_epsR%g_epsC%g.ps\' \n",xp,microResponse,dimX,nSteps,refinementThreshold,coarseningThreshold);
  //FILE *gp2;
  //gp2 = popen(GNUPLOT,"w");

  //  fprintf(gp2, "set terminal gif animate delay 10 \n");
  //fprintf(gp2, "set output \'data/error/trajectory_CoMD%i_dimX%i_nSteps%i_epsR%g_epsC%g.gif\' \n",microResponse,dimX,nSteps,refinementThreshold,coarseningThreshold);
  FILE *fp_approx;
  char filename[1000];
  char filename_exact[1000];
  //memset(filename_exact, 0, sizeof(filename_exact));
  sprintf(filename_exact,"data/simulation%i/trajectory/trajectory_CoMD%i_genuine.dat",xp,microResponse);
  FILE *fp_exact;
   fp_exact=fopen(filename_exact,"r");
  if(fp_exact==NULL){
    printf("Could not open genuine .dat file, closing...\n");
    exit(0);
  }

  FILE *fp_error;
  char filename_error[1000];
  //memset(filename_error, 0, sizeof(filename_error));
  sprintf(filename_error,"data/simulation%i/error/error_CoMD%i_dimX%i_nSteps%i_epsR%g_epsC%g.dat",xp,microResponse,dimX,nSteps,refinementThreshold,coarseningThreshold);
  
  fp_error=fopen(filename_error,"a");
  if(fp_error==NULL){
    printf("Could not open error .dat file, closing...\n");
    exit(0);
  }

  double stressError=0.0;
  double stressErrorMax[100];
  double stressErrorAvg[nSteps];

  double avgErrorStress=0.0;
  double stdDeviationErrorStress=0.0;


  double stressErrorPos[dimX];

  double strainError=0.0;
  double strainErrorMax[100];
  double strainErrorAvg[nSteps];


  double stressExact;
  double stressApprox;

  double strainExact;
  double strainApprox;

  double position;

  double stressMax=0.0;
  double strainMax=0.0;

  rewind(fp_exact);
  for(i=0;i<nSteps;i++){
    for(j=0;j<dimX;j++){
      fscanf(fp_exact,"%lf%lf%lf\n",&position,&strainExact,&stressExact);
      if(strainExact>strainMax) strainMax=strainExact;
      if(stressExact>stressMax) stressMax=stressExact;
    }
  } 

  printf("\n maximum strain = %g maximum stress = %g\n",strainMax, stressMax);
  //results.adaptiveInterpolation.1000.dat
 


  for(k=1;k>=1;k--){ 
    for(l=1;l<2;l++){
      
      stressErrorMax[k]=0.0;
      strainErrorMax[k]=0.0;
      avgErrorStress=0.0;
      stdDeviationErrorStress=0.0;

    
      // memset(filename, 0, sizeof(filename));
      sprintf(filename,"data/simulation%i/trajectory/trajectory_CoMD%i_adaptiveSampling_dimX%i_nSteps%i_epsR%g_epsC%g.dat",xp,microResponse,dimX,nSteps,refinementThreshold,coarseningThreshold);
      fp_approx=fopen(filename,"r");
 
 if(fp_approx==NULL){
   printf("Could not open adaptive .dat file %s, closing...\n", filename);
    exit(0);
  }

      printf("************** Opening %s **************\n\n",filename);
      rewind(fp_exact);
      for(i=0;i<nSteps;i++){
	stressErrorAvg[i]=0;

	for(j=0;j<dimX;j++){
	  
	  fscanf(fp_exact,"%lf%lf%lf\n",&position,&strainExact,&stressExact);
	  fscanf(fp_approx,"%lf%lf%lf\n",&position,&strainApprox,&stressApprox);
	  	  
	  stressError=fabs((stressExact-stressApprox))/stressMax;
	  avgErrorStress+=stressError/(nSteps*dimX);
	  if (stressError>stressErrorMax[k]){ 
	    stressErrorMax[k]=stressError;
	  }
	  stressErrorPos[j]=stressError;
	  stressErrorAvg[i]+=stressError/dimX;
	  
	  strainError=fabs((strainExact-strainApprox))/strainMax;
	  if (strainError>strainErrorMax[k]) {
	    strainErrorMax[k]=strainError;
	  }
	  strainErrorAvg[i]+=strainError/dimX;
	  
	  if(stressError>1) printf("Stress error of %f percent at time step %i and position %i\n\n",stressError,i,j);
	}
	
      }
    
      
      printf("Maximum relative strain error: %.8g \n",strainErrorMax[k]);
      printf("Maximum relative stress error :%.8g \n\n",stressErrorMax[k]);
      
      printf("Average relative stress error :%.8g \n\n",avgErrorStress);
    

      fclose(fp_approx);
      
      fprintf(fp_error,"epsR %g, epsC  %g max stress error %g average stress error %g\n",refinementThreshold, coarseningThreshold, stressErrorMax[k], avgErrorStress);
    }
  }
  plotError(gp,stressErrorAvg);
  
  pclose(gp);  
 
  fclose(fp_error);
  fclose(fp_exact);
}

