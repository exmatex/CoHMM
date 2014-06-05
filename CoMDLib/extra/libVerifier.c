#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "../CoMD_lib.h"

inline double GetCurrentTime( void )
{
	struct timeval tp;
	double time;

	gettimeofday( &tp, 0 );
	time = tp.tv_usec;
	time *= 1.0e-6;
	time += tp.tv_sec;

	return( time );
}

int main(int argc, char ** argv)
{

	CoMD_input theInput;

	theInput.potDir[0] = 0;
	theInput.potName[0] = 0;
	theInput.potType[0] = 0;
	theInput.doeam = -1;
	theInput.nx = 4;
	theInput.ny = 4;
	theInput.nz = 4;
	theInput.nSteps = 100;
	theInput.printRate = 10;
	//MUST SPECIFY THE FOLLOWING
	theInput.dt = 1.0;
	theInput.lat = -1.0;
	theInput.temperature = 300;
	theInput.initialDelta = 0.0;
	theInput.defGrad = 1.1;

	if(argc == 3)
	{
		//We assume input of the form ./libDriver $(nSteps) $(cubicRoot of nx*ny*nz)
		theInput.nSteps = atoi(argv[1]);
		theInput.nx = atoi(argv[2]);
		theInput.ny = atoi(argv[2]);
		theInput.nz = atoi(argv[2]);
	}
	else
	{
		if(argc == 2)
		{
			//We assume input of the form ,/libDriver $(nSteps)
			theInput.nSteps = atoi(argv[1]);
		}
		else
		{
			if(argc == 4)
			{
				theInput.nSteps = atoi(argv[1]);
				theInput.nx = atoi(argv[2]);
				theInput.ny = atoi(argv[2]);
				theInput.nz = atoi(argv[2]);
				theInput.defGrad = atof(argv[3]);
			}
			else
			{
				if(argc != 1)
				{
					fprintf(stderr, "Error: ./libDriver, ./libDriver $(nSteps), or ./libDriver $(nSteps) $(cubicRoot of nx*ny*nz), or ./libVerifier $(nSteps) $(cubicRoot of nx*ny*nz) $(defGrad) \n");
					return 1;
				}
			}
		}
	}

	
	int i;
	int nSteps = theInput.nSteps;
	printf("Step\teTot\tePot\teKin\ttemp\tStress\n");
	for(i = 0; i < nSteps; i+= theInput.printRate)
	{
		theInput.nSteps = i;
		CoMD_return  theRet = CoMD_lib(&theInput);
		printf("%d\t%f\t%f\t%f\t%f\t%f\n", theRet.simIter, theRet.eTot, theRet.ePot, theRet.eKin, theRet.temp, theRet.stressXX);
	}


	//printf("%d: %f, %f, %f, %f, %f, %f, %f, %f, %f \n", theRet.simIter, theRet.stressXX, theRet.stressXY, theRet.stressXZ, theRet.stressYX, theRet.stressYY, theRet.stressYZ, theRet.stressZX, theRet.stressZY, theRet.stressZZ);
	

	return 0;
}

