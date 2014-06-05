#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <CoMD_lib.h>
//#include <gsl/gsl_spline.h>
//#include <gsl/gsl_errno.h>


#define dimXMax 1000000
#define L 1000


FILE *fp1;
FILE *fp2;
FILE *fp3;
FILE *fperror;
FILE *input;
#ifdef GNUPLOT
FILE *gp;
#endif // GNUPLOT
FILE *sp;

float timer;

int dimX;
int nSteps;

int hmmPrintRate;

double stressTensor[dimXMax];
double stressTensorPast1[dimXMax];
double stressTensorPast2[dimXMax];

int numCoMDPast1;
int numCoMDPast2;
int numCoMD;

int microResponse; //1 to use atomistic response, 0 to use analytical response

int interpolationType;

int necessaryCoMD1[dimXMax];
int necessaryCoMD2[dimXMax];

int CoMD[dimXMax]; //position of CoMD calls, 1 for CoMD call, 0 otherwise
int printCoMD[dimXMax];
double fluxesTensor[dimXMax];

double dx;
double dt;

double CoMDdt;
int CoMDnSteps;
int CoMDnx;
int CoMDny;
int CoMDnz;

int checkPointRate=10;

struct stat st = {0};

//Lots of global variables that should be put as args of the main and passed to the various functions :)


double c; //1.0// wave speed for linear stress-strain relationship
double rho0;//221.815;//1.33 // mass / undeformed MD volume

double rho;

double latticeParameter;//3.609;//zero temp equilibrium lattice//3.6186845 is the equilibrium lattice for an energy density of -0.295
double CoMDTemp;
double t    = 0.0; // time (evolving)


 // Notation: A = A_11 -- deformation gradient along the x axis
//           p        -- momentum density
//           e        -- energy density

double * A0; // initial condition for A
double * e0; // initial conditions for e


double initialDeformation;//0.14 max with CoMD. Linear regime up to 0.05
double initialEnergy;//-0.2963; //-0.295 for ~300K initial temperature

double energyPerturbation;//0.005;//0.0;//0.01 to compare with full MD
double refinementThreshold;//error threshold of the refinement step, relative to the initial stress on the shock front //0.000001
double coarseningThreshold; //sensitivity parameter of the coarsening step//0.01

int linearInterpolation;
//int splineInterpolation;

int plot; // set to 1 to save the full trajectory. Can LOTS of data for big simulations!

int genuine; //set to 0 tu use the adaptive interpolation, to 1 to use the genuine CoHMM

int perturbationBeginning;
int perturbationEnd;
double deformationStart, deformationEnd;

//int meshCoarsening=5;

double initialPerturbation;


double * A;
double * p;
double * e;

double * temp1, * temp2, * temp3;
double * temp4, * temp5, * temp6;



//double kA= 100000; //10000 for shock and spike, 10 for gaussian and sinus // parameters affecting the probability of evaluating the stress in each grid point
//double kE= 100000; //10000 for shock ans spike, 10 for gaussian and sinus// the higher the more likely        

int CoMDcalls=0;
int refinementCalls=0;
int coarseningCalls=0;
//tolerance for interpolation error 
double interpolationThreshold;

double coarseningThreshold;
typedef struct task_s
{
	int startInc;
	int endExc;
} task_t;


char * updatedStress;
task_t * taskQueue;

task_t * taskQueuePast1;
task_t * taskQueuePast2;

// The elastic energy corresponding to the deformation gradient A at zero temperature
// (for the linear stress-strain relationship)
double zeroTempEnergyDensity(double A) {
    return (A-1)*(A-1)/2;

}

double gaussianRand(double m, double s){	/* normal random variate generator */
				        /* mean m, standard deviation s */
	double x1, x2, w, y1;
	//static double y2;
	//static int use_last = 0;

	//	if (use_last){		        /* use value from previous call */
	//  y1 = y2;
	//  use_last = 0;
	//}
	//	else{
	do {
	  x1 = 2.0 * (double)rand() / RAND_MAX - 1.0;
	  x2 = 2.0 * (double)rand() / RAND_MAX - 1.0;
	  w = x1 * x1 + x2 * x2;
	  
	} while ( w >= 1.0 );
	
	w = sqrt( (-2.0 * log( w ) ) / w );
	y1 = x1 * w;
	//y2 = x2 * w;
	//use_last = 1;
	//	}
	return( m + y1 * s );
}

/*
double splineInterpolate(double * grid, double * stress, int N, double x){
  double y;
  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  //change the following line to have whatever interpolation function from the gnu scientific library    
  const gsl_interp_type *t = gsl_interp_akima;//_periodic;
  gsl_spline *spline = gsl_spline_alloc (t, N);

 
  gsl_spline_init (spline, grid, stress, N);

  y=gsl_spline_eval (spline, x, acc);

  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);

  return y;
}
*/

double linearInterpolate(double ya, double yb, int xa, int xb, int x){
  //arguments:
  //index of the beginning and ending interpolation points, index of the interpolated point
  //former = 1 to interpolate on the stress from the previous step
  //former = 0 to interpolate the current stress

  double y;

  y=ya+1.0*(yb-ya)*(x-xa)/(xb-xa);
  
  return y;
}

char checkInterpolation(int xa, int xb, int x, int number){
  //function that checks how the interpolation is affected at x when removing the CoMD call at x
  double * gridPoints;
  double * stressPoints;
  
  int check=1;
  int i,j,k;
  //arguments:                                          
  //xa = index of the beginning, xb = index of the last interpolated points, x index of the interpolation result            

  //if checkInterpolation returns 1 then the interpolation is not good enough and CoMD needs to be called at x
  //if checkInterpolation returns 0 then the interpolation is good enough and CoMD doesn't need to be called at x
  double *w[3] = {A, p, e};

  if(number==1){
    gridPoints=malloc(sizeof(double) * numCoMDPast1);
    stressPoints=malloc(sizeof(double) * numCoMDPast1);
    

    double y=stressTensorPast1[x];
      // FIXME
    if(linearInterpolation==1){
      if(fabs(linearInterpolate(stressTensorPast1[xa],stressTensorPast1[xb],xa,xb,x)-y)>interpolationThreshold){
	check=1;
      }else{
	check=0;
      }
    }
/*
    if(splineInterpolation==1){
      i=0;
      for(j=0;j<numCoMDPast1;j++){
	if(taskQueuePast1[j].startInc!=x){
	  gridPoints[i]=taskQueuePast1[j].startInc*dx;
	  stressPoints[i]=stressTensorPast1[taskQueuePast1[j].startInc];       
	  i++;
	}
      }
      
      if(fabs(splineInterpolate(gridPoints,stressPoints,i,x*dx)-y)>interpolationThreshold){//x*dx
	check=1;
      }else{
      check=0;
      }      
    }
*/
  }

  if(number==2){
    gridPoints=malloc(sizeof(double) * numCoMDPast2);
    stressPoints=malloc(sizeof(double) * numCoMDPast2);
    double y=stressTensorPast2[x];
    if(linearInterpolation==1){
      if(fabs(linearInterpolate(stressTensorPast2[xa],stressTensorPast2[xb],xa,xb,x)-y)>interpolationThreshold){
        check=1;
      }else{
        check=0;
      }
    }

    /*if(splineInterpolation==1){
      i=0;
      for(j=0;j<numCoMDPast2;j++){
	if(taskQueuePast2[j].startInc!=x){
	  gridPoints[i]=taskQueuePast2[j].startInc*dx;
	  stressPoints[i]=stressTensorPast2[taskQueuePast2[j].startInc];
	  i++;
        }       
      }
      
      if(fabs(splineInterpolate(gridPoints,stressPoints,i,x*dx)-y)>interpolationThreshold){//x*dx
        check=1;
      }else{
	check=0;
      }
    }*/
  }
  free(gridPoints);
  free(stressPoints);
  return check;
}

char Coarsening(double *strainArr, double *energyArr, double * momentumArr, int index, int step)
{

  //  if((index == 0))
  //	{
  //		updatedStress[index] = 1;
  //		return 1;
  //	}

  //calculate the variations of the continuum fields
  double dALeft = (strainArr[index-1]-strainArr[index])/dx;
  double dARight = (strainArr[index+1]-strainArr[index])/dx;
  double dELeft = (energyArr[index-1]-energyArr[index])/dx;
  double dERight = (energyArr[index+1]-energyArr[index])/dx;

  //periodic boundary conditions should be put here on the variations.
  if(index == 0)
    {
      if(step == 0){
	updatedStress[index] = 1;                                                                                      
	return 1;
      }
      dALeft=0;
      dELeft=0;
    }
  if(index==dimX-1){    
    if(step == 0){
      updatedStress[index] = 1;
      return 1;
    }
    dARight=0;
    dERight=0;
  }

  //evaluate change in conditions 
  //Convenient mapping of the change in the continuum fields to determine where CoMDs are obviously necessary
  // double condChange = (step==0) ? 1-exp(-1000*kA*(fabs(dALeft)+fabs(dARight))-1000*kE*(fabs(dELeft)+fabs(dERight))) : 1-exp(-kA*(fabs(dALeft)+fabs(dARight))-kE*(fabs(dELeft)+fabs(dERight)));
  //This mapping is not nessecary and thresholds can be defined instead: above a certain variation in the continuum fields, a finer-scale  simulation is obviously necessary, with no need of a fine knowledge as regards the relation between the continuum fields and the finer\cale response.
  
  double condChange = (step==0) ? 1000*(fabs(dALeft)+1000*fabs(dARight))+2000*(fabs(dELeft)+2000*fabs(dERight)) : (fabs(dALeft)+fabs(dARight)) + 2*(fabs(dELeft)+2*fabs(dERight));
  if(condChange > coarseningThreshold) //>0.5 for the mapping
    {
      updatedStress[index] = 1;
      return 1;
    }
  else
    {
      updatedStress[index] = 0;
      return 0;
    }
}


double stressFn(double deformation, double energy, double momentum, int step, int gridPoint) {
  if(microResponse==1){
	CoMD_input theInput;
	
	//Set the various inputs of CoMD

    strcpy(theInput.potDir,"CoMDLib/potentials");
	strcpy(theInput.potName,"Cu01.eam.alloy");
	//strcpy(theInput.potName,"Cu_u6.eam");
	strcpy(theInput.potType,"setfl");
	//strcpy(theInput.potType,"funcfl");
	theInput.doeam = 1; 
	theInput.nx = CoMDnx;
	theInput.ny = CoMDny;
	theInput.nz = CoMDnz;
	theInput.nSteps = CoMDnSteps;//50 for full MD with no temp, 1000 for equilibrated stress
	theInput.printRate = 1; 
	//MUST SPECIFY THE FOLLOWING               
	theInput.dt = CoMDdt;
	theInput.lat = latticeParameter;
	theInput.energy = energy;
	theInput.initialDelta = 0.0;
	theInput.defGrad = deformation;
	theInput.temperature=CoMDTemp;//1200;//1000.0;


	//printf("Library[%d]: Calling stress/strain: %g -> \n", index, strain);
	//fflush(stdout);
	int retSteps;

	//	CoMD_return * theRet = CoMD_lib(&theInput, &retSteps);
	CoMD_return theRet = CoMD_lib(&theInput);

	printf("Library: Calculated strain -> stress: %g -> %g energy -> temperature: %g -> %g\n",deformation-1,(theRet.stressXX),energy,(theRet.temp));

	//printf("Library: Calculated strain -> stress: %g -> %g\n",deformation-1,(theRet.stressXX)); //Uncomment to print CoMD results to terminal



	stressTensor[gridPoint]=(theRet.stressXX);

	return (theRet.stressXX);


}else{

    double stressXX = 1.0364*rho*c*c*(deformation-1)-0*(energy-initialEnergy)+gaussianRand(0,0.000001);
  //printf("Calculated strain stress: %g -> %g\n", deformation-1, stressXX);
  stressTensor[gridPoint]=stressXX;
  return stressXX; // For this toy problem, just the derivative of the zero temp energy wrt A
  
 }
}

int mod(int x, int y) {
    int rem = x % y;
    return rem < 0 ? rem + y : rem;
}

// Exact solution for A(t), assuming linear stress/strain relation
double exactA(int i) {
    double del_i = (int) ((t * c) / dx);
    int i1 = mod(i+del_i, dimX);
    int i2 = mod(i-del_i, dimX);
    return 0.5*(A0[i1] + A0[i2]);
}


void refinement(int xa,int xb, int number, int step){//Algorithm 3 from the paper
  //Defines where the finer-scale simulation have to be called during the current time step, based on their importance during the previous one
  int j;
  int necessaryLeft=0;
  int necessaryRight=0;
  if(number==1){
    for(j=0;j<numCoMDPast1;j++){
      if(necessaryCoMD1[j]==1){ 
	if((taskQueuePast1[j].startInc>=xa-1)&&(taskQueuePast1[j].startInc<=xa+(xb-xa)/2)){
	  necessaryLeft=1;
	}
	if((taskQueuePast1[j].startInc>=xa+(xb-xa)/2)&&(taskQueuePast1[j].startInc<=xb+1)){
	  necessaryRight=1;
	}
      }
    }
  }
  if(number==2){
    for(j=0;j<numCoMDPast2;j++){
      if(necessaryCoMD2[j]==1){
        if((taskQueuePast2[j].startInc>=xa-1)&&(taskQueuePast2[j].startInc<=xa+(xb-xa)/2)){
          necessaryLeft=1;
        }
        if((taskQueuePast2[j].startInc>=xa+(xb-xa)/2)&&(taskQueuePast2[j].startInc<=xb+1)){
          necessaryRight=1;
        }
      }
    }
  }

  if(necessaryLeft==1){
    CoMD[xa+(xb-xa)/2]=1;    
    //if((linearInterpolation==1)){
    CoMD[xa+(xb-xa)/2-1]=1;
    CoMD[xa+(xb-xa)/2+1]=1;
    //}
    if(abs((xb-xa)/2)>1) refinement(xa,xa+(xb-xa)/2,number, step); 
  }
  if(necessaryRight==1){
    CoMD[xa+(xb-xa)/2]=1;
    //if((linearInterpolation==1)){
    CoMD[xa+(xb-xa)/2-1]=1;
    CoMD[xa+(xb-xa)/2+1]=1;
    //}
    if(abs(xb-xa-(xb-xa)/2)>1) refinement(xa+(xb-xa)/2,xb,number,step);
  }
} 

void findPoints(int number, int step){//Algorithm 2 of the paper
  //function that identifies where CoMD was necessary to maintain the interpolation errors under a certain threshold
  int x,xa,xb,j;
  if(number==1){
    necessaryCoMD1[0]=1;
    necessaryCoMD1[numCoMDPast1-1]=1;
    for(j=1;j<numCoMDPast1-1;j++){
      xa=taskQueuePast1[j-1].startInc;
      xb=taskQueuePast1[j+1].startInc;
      x=taskQueuePast1[j].startInc;
      if(checkInterpolation(xa,xb,x,1)==1){
	necessaryCoMD1[j]=1;
      }else{
      necessaryCoMD1[j]=0;
      }
    }
  }

  if(number==2){
    necessaryCoMD2[0]=1;
    necessaryCoMD2[numCoMDPast2-1]=1;
    for(j=1;j<numCoMDPast2-1;j++){
      xa=taskQueuePast2[j-1].startInc;
      xb=taskQueuePast2[j+1].startInc;
      x=taskQueuePast2[j].startInc;
      if(checkInterpolation(xa,xb,x,2)==1){
        necessaryCoMD2[j]=1;
      }else{
	necessaryCoMD2[j]=0;
      }
    }
  }
  
  //recursive function that refines the grid until the interpolation makes errors under the defined threshold 
  refinement(0,dimX-1,number, step);

}

void findPointsTest(){
  //crash test of HMM:
  //function that tests the sensivity to interpolation errors, putting a static location of microscopic calls

  int i;
  for(i=0;i<dimX;i++){
    CoMD[i]=(i%10);
  }
} 

// inputs: conserved fields A, p, e
// outputs: fluxes f_A, f_p, f_e

void fluxes(double *l_A, double *l_p, double *l_e, double *f_A, double *f_p, double *f_e, int step, int number) {
  int i,j;
  int xa=0;
  int xb=dimX-1;
  int x=dimX/2;
  char ret;
  

  //  printf("\n\n\n %i I \n\n\n",step);

  numCoMD = 0;
  //Identify which points require a finer-scale simulation based on the double adaptive interpolation scheme
  for(i = 0; i < dimX; i++){
    ret=0;

    if((CoMD[i]==1)||(Coarsening(l_A, l_e, l_p, i,step)==1)) ret=1;  
    if(Coarsening(l_A, l_e, l_p, i,step)==1){
      //coarsening step
      coarseningCalls++;
    }else{
      if(CoMD[i]==1){
	//refinement step
	refinementCalls++;
      }
    }

    if(genuine==1){
      ret=1;
      refinementCalls=0;
      coarseningCalls=0;
    }

    if(ret == 1){
      printCoMD[i]=1;
      updatedStress[i]=1;
      //Add a task
      taskQueue[numCoMD].startInc = i;
      taskQueue[numCoMD].endExc = i+1;
      numCoMD++;
      CoMDcalls++;
    }
  }
  
  //Distribute the CoMD calls on the cores
  //omp_set_dynamic(0);     // Explicitly disable dynamic teams
  //omp_set_num_threads(1);
  #pragma omp parallel for private(i)
  for(i = 0; i < numCoMD; i++)
    {      
     
      double stress = stressFn(l_A[taskQueue[i].startInc],l_e[taskQueue[i].startInc], l_p[taskQueue[i].startInc], step, taskQueue[i].startInc);
           
      double v = l_p[taskQueue[i].startInc] / rho;
      
      f_A[taskQueue[i].startInc] = -v;
      f_p[taskQueue[i].startInc] = -stress;
      f_e[taskQueue[i].startInc] = -stress*v;
    }

  int x1,x2;  
  double y1,y2;
  
  double * gridPoints;
  double * stressPoints;
  double * momentumPoints;

  gridPoints=malloc(sizeof(double) * numCoMD);
  stressPoints=malloc(sizeof(double) * numCoMD);
  momentumPoints=malloc(sizeof(double) * numCoMD);
  if(genuine==0){
    for(j=0;j<numCoMD;j++){
      //update the microscopic fields with the CoMD results
      gridPoints[j]=taskQueue[j].startInc*dx;
      stressPoints[j]=stressTensor[taskQueue[j].startInc];
      momentumPoints[j]=l_p[taskQueue[j].startInc];
    }
  }
  /*
  gsl_interp_accel *accStress;
  gsl_interp_accel *accMomentum;

  //Choose any interpolation function from the gnu scientific library here                                                              

  const gsl_interp_type *intType = gsl_interp_akima;//_periodic;
  gsl_spline *splineStress;
  gsl_spline *splineMomentum;

  if((splineInterpolation==1)&&(genuine==0)){
    accStress = gsl_interp_accel_alloc ();
    accMomentum = gsl_interp_accel_alloc ();    
        
    splineStress = gsl_spline_alloc (intType, numCoMD);
    splineMomentum = gsl_spline_alloc (intType, numCoMD);
    gsl_spline_init (splineStress, gridPoints, stressPoints, numCoMD);
    gsl_spline_init (splineMomentum, gridPoints, momentumPoints, numCoMD);
  }
*/

  if(genuine==0){
    //for(i = 0; i < dimX; i++){
    for(i = 1; i < dimX-1; i++){
      if(updatedStress[i] == 0){     
	if((linearInterpolation==1)||(step==0)){
	  for(j=0;j<numCoMD;j++){
	    if(taskQueue[j].startInc<i){
	      x1=taskQueue[j].startInc;
	  }
	    if(taskQueue[j].startInc>i){
	      x2=taskQueue[j].startInc;
	      j=numCoMD;
	    }	  
	  }       	
	  y1=-f_p[x1];
	  y2=-f_p[x2];
	  //linear interpolation of the CoMD results
	  double stress = linearInterpolate(y1,y2,x1,x2,i);	
	  stressTensor[i]=stress;
	  double v = linearInterpolate(l_p[x1],l_p[x2],x1,x2,i) / rho;	
	  f_A[i] = -v;
	  f_p[i] = -stress;
	  f_e[i] = -stress*v;
	}
	/*if((splineInterpolation==1)&&(step>0)){
	  //Akima spline interpolation of the CoMD results
	  double stress=gsl_spline_eval (splineStress, i*dx, accStress);
	  
	  stressTensor[i]=stress;	
	  double v=gsl_spline_eval (splineMomentum, i*dx, accMomentum)/rho;
	  
	  f_A[i] = -v;
	  f_p[i] = -stress;
	  f_e[i] = -stress*v;
	  
	}*/
      }
    }
  }
  
  //Copy the former CoMD results to later check the importance thay had in for the quality of the interpolation
  if(number==1){
    numCoMDPast1=numCoMD;
    for(j=0;j<dimX;j++){
      stressTensorPast1[j]=stressTensor[j];
    }

    for(j=0;j<numCoMDPast1;j++){ 
      taskQueuePast1[j].startInc = taskQueue[j].startInc;
      taskQueuePast1[j].endExc = taskQueue[j].endExc;
    }
  }

  if(number==2){
    numCoMDPast2=numCoMD;
    for(j=0;j<dimX;j++){
      stressTensorPast2[j]=stressTensor[j];
    }

    for(j=0;j<numCoMDPast2;j++){
      taskQueuePast2[j].startInc = taskQueue[j].startInc;
      taskQueuePast2[j].endExc = taskQueue[j].endExc;
    }
  }
  free(gridPoints);
  free(stressPoints);
  free(momentumPoints);
  /*if((splineInterpolation==1)&&(genuine==0)){
    gsl_spline_free (splineStress);
    gsl_spline_free (splineMomentum);
    gsl_interp_accel_free (accStress);
    gsl_interp_accel_free (accMomentum);
  }*/
}


void initializedConservedFields() {
  int i;
  for (i = 0; i < dimX; i++) {

    //various initial conditions
    
    //A0[i] = A[i] = ((i < perturbationEnd)) ? 1.0+initialDeformation : 1.0; //shock       
	    
     A0[i] = A[i] = ((i < perturbationEnd)&&(i >= perturbationBeginning)) ? 1.0+initialDeformation : 1.0; //spike
      
    //    A0[i] = A[i] = 1.0 + initialDeformation*exp(-0.5*pow(1.0*(i-dimX/2)/(dimX/10),2)); //gaussian
      
    // A0[i] = A[i] =1.0 + initialDeformation*sin(i*2*3.1415926/dimX); //sinus

      // small initial step in deformation gradient
      //e0[i] = e[i] = initialEnergy;
      p[i] = 0;

      //if(microResponse==0){
      //initialEnergy=0;
      //energyPerturbation=zeroTempEnergyDensity(1+initialDeformation);
      //}
      
      e0[i] = e[i] = ((i < perturbationEnd)&&(i >= perturbationBeginning)) ? initialEnergy+energyPerturbation : initialEnergy; // spike
      // e0[i] = e[i] = initialEnergy+energyPerturbation*exp(-0.5*pow(1.0*(i-dimX/2)/(dimX/10),2)); //gaussian

  }
}

void initializedConservedFields_old() {
  int i;
  for (i = 0; i < dimX; i++) {
    A0[i] = A[i] = (i < dimX/4) ? 1.01 : 1.0; // small initial step in deformation gradient                                 
    p[i] = 0;
    e[i] = zeroTempEnergyDensity(A[i]);
  }
}

double netEnergy() {
    double acc = 0.0;
    int i;
    for (i = 0; i < dimX; i++) {
        acc += e[i];
    }
    return acc;
}

double minmod(double x, double y) {
  return (x > 0 == y > 0) ? 0.5*(fabs(x)+fabs(y)) : 0;
}


double TadmorMinmod(double x, double y){
  double result;
  if((x > 0) != (y > 0)){
    result=0;
  }else{
    if(x>0){
      result=fmin(x,y);
    }else{
      result=fmax(x,y);
    }
  }
  return result;
}

void modulatedDerivative(double *w, double *ret) {
    int i;
    for (i = 0; i < dimX; i++) {
        int im = (i-1+dimX)%dimX;
        int ip = (i+1)%dimX;
        //ret[i] = minmod(w[ip]-w[i], w[i]-w[im]) / dx;
	ret[i] = TadmorMinmod(w[ip]-w[i], w[i]-w[im]) / dx;                                                                      
    }
}

void halfStep2ndOrder(int step) {
  int i,j;
  double *ws[3] = {A, p, e};
  double *fs[3] = {temp1, temp2, temp3};
  double *wps[3] = {temp4, temp5, temp6};
 

  //adaptive interpolation algorithm called to find the points where CoMD is needed for the first half step
  if((step>0)&&(genuine==0)) findPoints(1,step);
  // calculate f^(n+1/2)                                                                                                                 
  fluxes(ws[0], ws[1], ws[2], fs[0], fs[1], fs[2], step, 1);
  // numCoMDPast=numCoMD;  

  
  // calculate w^(n+1/2)                                                                                        
  for (j = 0; j < 3; j++) {
    double *w  = ws[j];
    double *wp = wps[j];
    double *f  = fs[j];
    modulatedDerivative(f, wp);
    for (i = 0; i < dimX; i++)
      wp[i] = w[i] - (0.5*dt/(2*dx)) * wp[i];
  }

 
  //adaptive interpolation algorithm called to find the points where CoMD is needed for the second half step
  if((step>0)&&(genuine==0)) findPoints(2,step); 

  // calculate f^(n+1/2) 
  fluxes(wps[0], wps[1], wps[2], fs[0], fs[1], fs[2], step,2);


  // calculate w^(n+1)                                                                                          
  for (j = 0; j < 3; j++) {
    double *w  = ws[j];
    double *f = fs[j];
    double *dw = temp4;
    modulatedDerivative(w, dw);
    double w0 = w[0];
    for (i = 0; i < dimX; i++) {
      int ip = (i+1)%dimX;
      double wn = (i == dimX-1) ? w0 : w[i+1];
      w[i] = 0.5*(w[i] + wn) + (dx/8) * (dw[i] - dw[ip]) - (0.5*dt/dx) * (f[ip] - f[i]);
    }
  }
}


void fullStep(int step) {
    int i,j;
    // two half steps on a staggered grid


    // if(step>=meshCoarsening) findPoints();

    halfStep2ndOrder(step);
   
    halfStep2ndOrder(step);
    
    // rotate back to original index <-> coordinate map
    double *ws[3] = {A, p, e};
    for (j = 0; j < 3; j++) {
        double *w = ws[j];
        double w_last = w[dimX-1];
	for (i = dimX-1; i > 0; i--){
            w[i] = w[i-1];	 
	}
        w[0] = w_last;
    }
    
    double stressTensorLast=stressTensor[dimX-1];
    for (i = dimX-1; i > 0; i--){            
      stressTensor[i]=stressTensor[i-1];
    }    
    stressTensor[0]=stressTensorLast;
    t += dt;
}

#ifdef GNUPLOT
void plotFields(FILE *gp) {
  int i,j,print;
  print=1;
  fprintf(gp, "set yrange [-0.001:0.101]\n");  
  fprintf(gp, "set xrange [0:%i]\n",dimX);
  fprintf(gp, "set xlabel 'Continuum grid'\n");
  fprintf(gp, "set key right top\n");    
  fprintf(gp, "plot '-' u 1:2 ti 'Strain (delta l/L)' w l lw 2 linecolor rgb \"red\", '-' u 1:3 ti 'Microscopic simulation' w points pt 6 ps 1 lw 1 linecolor rgb \"green\", '-' u 1:4 ti 'Energy density' w l lw 2 linecolor rgb \"blue\"\n\n");  
    for (j = 0; j < 3; j++) { // gnuplot forces us to repeat ourselves
        for (i = 0; i < dimX; i++) {
	  if(i%1==0){
	    if(i>0) print=printCoMD[i-1];
	    if(print==1){
	      fprintf(gp, "%i %f %f %f\n", i, A[i]-1, A[i]-1, 10*(e[i]-initialEnergy));
	    }else{
	      //fprintf(gp, "%i %f -\n", i, A[i]-1); //!!!!!!!!!!!!!!
	    }
	  }
	}
        fprintf(gp, "e\n");
    }
    fflush(gp);
}
#endif // GNUPLOT

void storeFields(int step, int xp) {
  int i,j,print;
  char storename[255];
  FILE *sp;
  if(genuine==0){
    sprintf(storename,"output/simulation%i/storage/storage%i.adaptiveInterpolation.%i.%g.dat",xp,step,dimX,refinementThreshold);
  }else{
    sprintf(storename,"output/simulation%i/storage/storage%i.dat",xp,step);
  }
  sp = fopen(storename,"w");

 
  print=1;
  for (j = 0; j < 1; j++) {
    for (i = 0; i < dimX; i++) {
      if(i%1==0){
	if(i>0) print=printCoMD[i-1];
	if(print==1){
	  fprintf(sp, "%i %f %f\n", i, stressTensor[i], stressTensor[i]);
	}else{
	  fprintf(sp, "%i %f -\n", i, stressTensor[i]);
	}
      }
    }
  }
  fclose(sp);
}

void checkPoint(int step, int xp){
  int i,j,print;
  char checkpointname[255];
  FILE *cp;
  
    sprintf(checkpointname,"output/simulation%i/storage/checkpoint%i.dat",xp,step);
  
  
  
  cp = fopen(checkpointname,"w");

 
  print=1;
  for (j = 0; j < 1; j++) {
    for (i = 0; i < dimX; i++) {
      fprintf(cp, "%i %f %f %f %f\n", i, A[i],stressTensor[i], e[i], p[i]);
    }
  }
  fclose(cp);
}


static inline void loadBar(int x, int n, int r, int w)
{
  // Only update r times.
  if ( x % (n/r) != 0 ) return;
 
  // Calculuate the ratio of complete-to-incomplete.
  float ratio = x/(float)n;
  int   c     = ratio * w;
 
  //Custom 2013 summer school loading bar, essential to the code

  if((int)(ratio*100)==2)   printf("\n\n:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n");
  if((int)(ratio*100)==4)   printf(":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n");
  if((int)(ratio*100)==8)   printf(":::::::::::::::::::::::::::::::::::::::::::::-'    `-::::::::::::::::::\n");
  if((int)(ratio*100)==12)  printf("::::::::::::::::::::::::::::::::::::::::::-'          `::::::::::::::::\n");
  if((int)(ratio*100)==16)  printf(":::::::::::::::::::::::::::::::::::::::-  '   /(_M_)\\  `:::::::::::::::\n");
  if((int)(ratio*100)==20)  printf(":::::::::::::::::::::::::::::::::::-'        |       |  :::::::::::::::\n");
  if((int)(ratio*100)==24)  printf("::::::::::::::::::::::::::::::::-         .   \\/~V~\\/  ,:::::::::::::::\n");
  if((int)(ratio*100)==28)  printf("::::::::::::::::::::::::::::-'             .          ,::::::::::::::::\n");
  if((int)(ratio*100)==32)  printf(":::::::::::::::::::::::::-'                 `-.    .-::::::::::::::::::\n");
  if((int)(ratio*100)==36)  printf(":::::::::::::::::::::-'                  _,,-::::::::::::::::::::::::::\n");
  if((int)(ratio*100)==40)  printf("::::::::::::::::::-'                _,--:::::::::::::::::::::::::::::::\n");
  if((int)(ratio*100)==44)  printf("::::::::::::::-'               _.--:::::::::::::::::::::::#####::::::::\n");
  if((int)(ratio*100)==48)  printf("::::::::-'    ##     ###.-:::::###::::::::::::::::::::::::#####::::####\n");
  if((int)(ratio*100)==52)  printf("::::-'       ###_.::######:::::###::::::::::::::#####:##########:::####\n");
  if((int)(ratio*100)==56)  printf(":'         .-###::########:::::###::::::::::::::#####:##########:::####\n");
  if((int)(ratio*100)==60)  printf("      ..--:::###::########:::::###:::::######:::#####:##########:::####\n");
  if((int)(ratio*100)==64)  printf(" _.-':::##:::###:#########:::::###:::::######:::#####:#################\n");
  if((int)(ratio*100)==68)  printf("'#########:::###:#########::#########::######:::#####:#################\n");
  if((int)(ratio*100)==72)  printf(":#########:::#############::#########::######:::#######################\n");
  if((int)(ratio*100)==76)  printf("##########:::########################::################################\n");
  if((int)(ratio*100)==80)  printf("##########:::##########################################################\n");
  if((int)(ratio*100)==84)  printf("##########:::##########################################################\n");
  if((int)(ratio*100)==88)  printf("#######################################################################\n");
  if((int)(ratio*100)==92)  printf("#######################################################################\n");
  if((int)(ratio*100)==96)  printf("##################################################### Co Design 2013 ##\n");
  if((int)(ratio*100)==98)  printf("#######################################################################\n\n\n");
}

void storeErrors(FILE *fperror) {
  //prints outputs to compare adaptive mesh and database scheme against the original macrosolver                                           
  int i,j;

  for (i = 0; i < dimX; i++) {
    fprintf(fperror, "%.16lf %.16lf %.16lf\n", i*dx, A[i]-1,stressTensor[i]);
  }


  fflush(fperror);
}

void readInput(int xp){
  char inputname[1000];
  sprintf(inputname,"input/input%i.txt",xp);
  input = fopen(inputname, "r");
  if (input==NULL) {
      printf("Could not find an input file! Closing :-( \n");
      exit(0);
  }
  
  rewind (input);
  fscanf (input, "CoMD 1 Analytical 0 %d\n",&microResponse);
  //fscanf (input, "Wave speed for the analytical response %lf\n",&c);
  // fscanf (input, "Reference material's density %lf\n",&rho0);
  fscanf (input, "Number of HMM time steps %d\n", &nSteps);
  fscanf (input, "Number of HMM grid points %d\n", &dimX);
  fscanf (input, "HMM print rate %d\n",&hmmPrintRate);
  fscanf (input, "HMM spatial resolution dx in A %lf\n", &dx);
  fscanf (input, "HMM time resolution dt in fs %lf\n", &dt);
  fscanf (input, "Initial deformation %lf\n", &initialDeformation);
  fscanf (input, "Initial energy density where undeformed in eV per A cube %lf\n", &initialEnergy);
  fscanf (input, "Energy added with deformation in eV per A cube %lf\n", &energyPerturbation);
  fscanf (input, "Deformation start and end in percentage of HMM mesh %lf %lf\n", &deformationStart, &deformationEnd);
  fscanf (input, "CoMD parameters:\nTemperature %lf\nLattice parameter in A %lf\nNumber of lattices nx ny nz %d %d %d\nTime resolution dt in fs %lf\nNumber of MD time steps %d\n", &CoMDTemp,&latticeParameter, &CoMDnx,&CoMDny,&CoMDnz,&CoMDdt,&CoMDnSteps);
  fscanf (input, "Spatial adaptive sampling parameters:\nRefinement threshold %lf\nCoarsening threshold %lf\nInterpolation: linear 0 Akima splines 1 %d\nPlot the CoHMM trajectory 0 no 1 yes %d\n", &refinementThreshold, &coarseningThreshold,&interpolationType,&plot);
  fscanf (input, "Avoid using the adaptive sampling 1 yes 0 no %d\n", &genuine);
  
  // KMB: enforce linear interpolation to remove dependency on GSL
  assert("Only linear interpolation is supported!" && interpolationType == 0);
  /*if(interpolationType==1){
    splineInterpolation=1;
    linearInterpolation=0;
  }else{
    splineInterpolation=0;
    linearInterpolation=1;
  }
  */
  linearInterpolation=1;

   if(genuine==1){
     linearInterpolation=0;
     //splineInterpolation=0;
   }
   //rho=rho0;
   //rho=1.33;
   //c=1;
   rho=4*6592.9/pow(latticeParameter,3); //copper density in CoMD units   
   c=sqrt(1/rho);
   
   printf("\n\n\n\n\n\n\n\n#############################################################\n Launch of HMM simulation with the following parameters:\n#############################################################\n\n");
   
   printf( "CoMD 1 Analytical 0 %d\n",microResponse);
   printf ( "Wave speed for the analytical response %lf\n",c);
   printf ( "Reference material's density %lf\n",rho);
   printf ( "Number of HMM time steps %d\n", nSteps);
   printf ( "Number of HMM grid points %d\n", dimX);
   printf ( "HMM print rate %d\n",hmmPrintRate);
   printf ( "HMM spatial resolution dx in A %lf\n", dx);
   printf ( "HMM time resolution dt in fs %lf\n", dt);
   printf ( "Initial deformation %lf\n", initialDeformation);
   printf ( "Initial energy density where undeformed in eV per A cube %lf\n", initialEnergy);
   printf ( "Energy added with deformation in eV per A cube %lf\n", energyPerturbation);
   printf ( "Deformation start and end in percentage of HMM mesh %lf %lf\n", deformationStart, deformationEnd);
   printf ("CoMD parameters:\nTemperature %lf\nLattice parameter in A %lf\nNumber of lattices nx ny nz %d %d %d\nTime resolution dt in fs %lf\nNumber of MD time steps %d\n", CoMDTemp,latticeParameter, CoMDnx,CoMDny,CoMDnz,CoMDdt,CoMDnSteps);
   printf ( "Spatial adaptive sampling parameters:\nRefinement threshold %lf\nCoarsening threshold %lf\nInterpolation: linear 0 Akima splines 1 %d\nPlot the CoHMM trajectory 0 no 1 yes %d\n", refinementThreshold, coarseningThreshold,interpolationType,plot);
   printf ( "Avoid using the adaptive sampling 1 yes 0 no %d\n", genuine);
   
   perturbationBeginning=(int)(deformationStart*dimX/100); //full MD initial conditions
   perturbationEnd=(int)(deformationEnd*dimX/100);
   
   //perturbationBeginning=dimX/2-dimX/40; //speedup conditions
   //perturbationEnd=dimX/2+dimX/40;
  
   interpolationThreshold=(1.03*rho*c*c*initialDeformation)*refinementThreshold;//(initialDeformation)*refinementThreshold;
   coarseningThreshold=(1.03*rho*c*c*initialDeformation)*coarseningThreshold/dx;//(initialDeformation)*coarseningThreshold/dx;
   
}

void createFiles(int xp){

  char dirname[1000];
  char dirname1[1000];
  char dirname2[1000];
  char dirname3[1000];

  sprintf(dirname,"output/simulation%i",xp);
  sprintf(dirname1,"output/simulation%i/storage",xp);
  sprintf(dirname2,"output/simulation%i/adaptive",xp);
  sprintf(dirname3,"output/simulation%i/trajectory",xp);
  
  if (stat(dirname, &st) == -1) {
    mkdir(dirname, 0700);
    mkdir(dirname1, 0700);
    mkdir(dirname2, 0700);
    mkdir(dirname3, 0700);
  }

  char filename1[1000];
  char filename2[1000];
  char filename3[1000];
  char errorname[1000];
  char gnuplotname[1000];
  
  // open pipe to gnuplot instance
#ifdef GNUPLOT
  gp = popen(GNUPLOT,"w");
  if (gp==NULL) {
    printf("Error opening pipe to GNU plot < %s > t! \n", GNUPLOT);
    exit(0);
  }
#endif // GNUPLOT
    
  if(genuine==0){
    sprintf(errorname,"output/simulation%i/trajectory/trajectory_CoMD%i_adaptiveSampling_dimX%i_nSteps%i_epsR%g_epsC%g.dat",xp,microResponse,dimX,nSteps,refinementThreshold,coarseningThreshold);	
  }else{	  
    sprintf(errorname,"output/simulation%i/trajectory/trajectory_CoMD%i_genuine.dat",xp,microResponse);	  
  }
  fperror = fopen(errorname,"w");
  
  /*if(splineInterpolation==1){
    sprintf(filename1,"output/simulation%i/adaptive/s_totalCoMDcalls_dimX%i_nSteps%i_epsR%g_espC%g.dat",xp,dimX,nSteps,refinementThreshold, coarseningThreshold );
    sprintf(filename2,"output/simulation%i/adaptive/s_coMDcallsPerTimeStep_dimX%i_nSteps%i_epsR%g_espC%g.dat",xp,dimX,nSteps,refinementThreshold, coarseningThreshold );     
    sprintf(filename3,"output/simulation%i/adaptive/s_percentageOfCoMDcallsPerTimeStep_dimX%i_nSteps%i_epsR%g_espC%g.dat",xp,dimX,nSteps,refinementThreshold, coarseningThreshold );
    sprintf(gnuplotname, "set output \'output/simulation%i/s_cohmm_dimX%i_nSteps%i_epsR%g_espC%g.ps\' \n",xp,dimX,nSteps,refinementThreshold, coarseningThreshold);
  }*/
  if(linearInterpolation==1){
    sprintf(filename1,"output/simulation%i/adaptive/l_totalCoMDcalls_dimX%i_nSteps%i_epsR%g_espC%g.dat",xp,dimX,nSteps,refinementThreshold, coarseningThreshold );
    sprintf(filename2,"output/simulation%i/adaptive/l_coMDcallsPerTimeStep_dimX%i_nSteps%i_epsR%g_espC%g.dat",xp,dimX,nSteps,refinementThreshold, coarseningThreshold );
    sprintf(filename3,"output/simulation%i/adaptive/l_percentageOfCoMDcallsPerTimeStep_dimX%i_nSteps%i_epsR%g_espC%g.dat",xp,dimX,nSteps,refinementThreshold, coarseningThreshold );
    sprintf(gnuplotname, "set output \'output/simulation%i/l_cohmm_dimX%i_nSteps%i_epsR%g_espC%g.ps\' \n",xp,dimX,nSteps,refinementThreshold, coarseningThreshold);
  }
  
  if(genuine==1){
    sprintf(filename1,"output/simulation%i/adaptive/g_totalCoMDcalls_dimX%i_nSteps%i.dat",xp,dimX,nSteps);
    sprintf(filename2,"output/simulation%i/adaptive/g_coMDcallsPerTimeStep_dimX%i_nSteps%i.dat",xp,dimX,nSteps);
    sprintf(filename3,"output/simulation%i/adaptive/g_percentageOfCoMDcallsPerTimeStep_dimX%i_nSteps%i.dat",xp,dimX,nSteps);
    sprintf(gnuplotname, "set output \'output/simulation%i/g_cohmm_dimX%i_nSteps%i.ps\' \n",xp,dimX,nSteps);
  }
  
  
#ifdef GNUPLOT
  //fprintf(gp, "set terminal gif animate delay 10 \n");
  fprintf(gp, "set terminal postscript \n");
  
  fprintf(gp, "%s\n", gnuplotname);
  //fprintf(gp, "set title \'Multi-Scale Modeling\'\n");
#endif // GNUPLOT
  
  fp1=fopen(filename1,"w");
  if (fp1==NULL) {
    printf("Error opening file/n");
    exit(0);
  }
  fp2=fopen(filename2,"w");
  if (fp2==NULL) {
    printf("Error opening file/n");
    exit(0);
  }
  fp3=fopen(filename3,"w");
  if (fp3==NULL) {
    printf("Error opening file/n");
    exit(0);
  } 

}


void cleanMemory(){
    fclose(fp1);
    fclose(fp2);
    fclose(fp3);
    fclose(fperror);
    fclose(input);
#ifdef GNUPLOT
    pclose(gp);
#endif // GNUPLOT
    
    free(A0);
    free(A);
    free(p);
    free(e);
    free(e0);
    free(temp1);
    free(temp2);
    free(temp3);
    free(temp4);
    free(temp5);
    free(temp6);
    free(updatedStress);
    free(taskQueue);
    free(taskQueuePast1);
    free(taskQueuePast2);

}

void epilogue(int xp){

  FILE *summary;
  
  char summaryname[1000];
 
  sprintf(summaryname,"output/simulation%i/summary.dat",xp); 
  
  summary=fopen(summaryname,"a");
  

  fprintf(summary, "CoMD 1 Analytical 0 %d\n",microResponse);
  fprintf (summary, "Wave speed for the analytical response %lf\n",c);
  fprintf (summary, "Reference material's density %lf\n",rho0);
  fprintf (summary, "Number of HMM time steps %d\n", nSteps);
  fprintf (summary, "Number of HMM grid points %d\n", dimX);
  fprintf (summary, "HMM spatial resolution dx in A %lf\n", dx);
  fprintf (summary, "HMM time resolution dt in fs %lf\n", dt);
  fprintf (summary, "Initial deformation %lf\n", initialDeformation);
  fprintf (summary, "Initial energy density where undeformed in eV per A cube %lf\n", initialEnergy);
  fprintf (summary, "Energy added with deformation in eV per A cube %lf\n", energyPerturbation);
  fprintf (summary, "Deformation start and end in percentage of HMM mesh %lf %lf\n", deformationStart, deformationEnd);
  fprintf (summary,"CoMD parameters:\nLattice parameter in A %lf\nNumber of lattices nx ny nz %d %d %d\nTime resolution dt in fs %lf\nNumber of MD time steps %d\n", latticeParameter, CoMDnx,CoMDny,CoMDnz,CoMDdt,CoMDnSteps);
  fprintf ( summary,"Spatial adaptive sampling parameters:\nRefinement threshold %lf\nCoarsening threshold %lf\nInterpolation: linear 0 Akima splines 1 %d\nPlot the CoHMM trajectory 0 no 1 yes %d\n", refinementThreshold, coarseningThreshold,interpolationType,plot);
  fprintf (summary, "Avoid using the adaptive sampling 1 yes 0 no %d\n\n", genuine);




    fprintf(summary,"\nNumber of CoMD calls                                  %i \n",CoMDcalls);
    fprintf(summary,"Number of CoMD calls avoided with the adaptive mesh   %i \n\n\n",(nSteps*dimX*4-CoMDcalls));

    fprintf(summary,"Number of CoMD calls from the mesh coarsening scheme  %i \n",coarseningCalls);
    fprintf(summary,"Number of CoMD calls from the mesh refinement scheme  %i \n",refinementCalls);




  fprintf(summary,"\n==> user time of %g seconds\n\n",timer);
  fprintf(summary,"\n###################################################################\n\n\n");
  
  fclose(summary);
  
}

int main(int argc, char **argv) {
    int i,j;
    if(argc != 2){
      fprintf(stderr, "./main (input file number)");
      return 1;
    }
    int xp=atoi(argv[1]);

    time_t t1, t2;
    t1=time(NULL);
 
    readInput(xp);
    
    createFiles(xp);
    
    A0 = malloc(sizeof(double) * dimX);
    A = malloc(sizeof(double) * dimX);
    p = malloc(sizeof(double) * dimX);
    e = malloc(sizeof(double) * dimX);
    e0 = malloc(sizeof(double) * dimX);
    temp1 = malloc(sizeof(double) * dimX);
    temp2 = malloc(sizeof(double) * dimX);
    temp3 = malloc(sizeof(double) * dimX);
    temp4 = malloc(sizeof(double) * dimX);
    temp5 = malloc(sizeof(double) * dimX);
    temp6 = malloc(sizeof(double) * dimX);
    
    updatedStress = malloc(sizeof(char) * dimX);
    
    taskQueue = malloc(sizeof(task_t) * dimX);
    taskQueuePast1 = malloc(sizeof(task_t) * dimX);
    taskQueuePast2 = malloc(sizeof(task_t) * dimX);
    
    
    initializedConservedFields();
    int CoMDcallsPast,CoMDcallsPerTimeStep;
    int coarseningCallsPast,coarseningCallsPerTimeStep;
    int refinementCallsPast,refinementCallsPerTimeStep;
        
    
    for (i = 0; i < nSteps; i++) {
      CoMDcallsPast=CoMDcalls;
      coarseningCallsPast=coarseningCalls;
      refinementCallsPast=refinementCalls;
      //if(microResponse==1){ 
      //	printf("\n\n################ Macrosolver step %i ################\n",i+1);
      // }else{
      //	if(i>0) loadBar(i,nSteps,50,100);
      // }
      
      for(j=0;j<dimX;j++){ 	
	CoMD[j]=0;
	printCoMD[j]=0;
      }      
      fullStep(i);     
      
      CoMDcallsPerTimeStep=CoMDcalls-CoMDcallsPast;
      coarseningCallsPerTimeStep=coarseningCalls-coarseningCallsPast;
      refinementCallsPerTimeStep=refinementCalls-refinementCallsPast;
      if(i%hmmPrintRate==0){
	//printf("Adaptive plot !\n");
	fprintf(fp1,"%i %i %i %i %i\n",i+1,CoMDcalls, coarseningCalls, refinementCalls, (i+1)*dimX*4);
	//fprintf(fp2,"%i %i %i %i %i\n",i+1,CoMDcallsPerTimeStep, coarseningCallsPerTimeStep, refinementCallsPerTimeStep, dimX*4);
	fprintf(fp3,"%i %f %f %f\n",i+1,100.0*CoMDcallsPerTimeStep/(dimX*4), 100.0*coarseningCallsPerTimeStep/(dimX*4), 100.0*refinementCallsPerTimeStep/(dimX*4));
	if(microResponse==1){
	  printf("\nCumulative: total calls, gradient calls, refinement calls, max calls\n");
	  printf("%i %i %i %i %i\n",i+1,CoMDcalls, coarseningCalls, refinementCalls, (i+1)*dimX*4);
	  printf("\nCurrent step: total calls, gradient calls, refinement calls, max calls\n");
	  printf("%i %i %i %i %i\n",i+1,CoMDcallsPerTimeStep, coarseningCallsPerTimeStep, refinementCallsPerTimeStep, dimX*4);
	
	}else{
	  if(i>0) loadBar(i,nSteps,50,100);
	}
	//printf("%i %f %f %f\n",i+1,100.0*CoMDcallsPerTimeStep/(dimX*4), 100.0*coarseningCallsPerTimeStep/(dimX*4), 100.0*refinementCallsPerTimeStep/(dimX*4));


      }
      if (plot==1) storeErrors(fperror); //<====

#ifdef GNUPLOT
      if((i>0)&&(i%hmmPrintRate==0)) plotFields(gp);
#endif // GNUPLOT

      if(plot==1){
	storeFields(i,xp);
      }
      
      if(i%checkPointRate==0) checkPoint(i,xp);

 
    }
    
    printf("Number of CoMD calls                                  %i \n",CoMDcalls);
    printf("Number of CoMD calls avoided with the adaptive mesh   %i \n\n\n",(nSteps*dimX*4-CoMDcalls));

    printf("Number of CoMD calls from the mesh coarsening scheme  %i \n",coarseningCalls);
    printf("Number of CoMD calls from the mesh refinement scheme  %i \n",refinementCalls);

    cleanMemory();


    t2=time(NULL);
    timer = (float)(t2-t1);


    epilogue(xp);

    
    printf("Wall clock execution time %g\n",timer);
    
    return 0;
}

