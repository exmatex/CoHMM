#include <iostream>
#include <cstdlib>
#include <cstring>

#include "cohmm_dad.decl.h"

#include "CoHMM_DaD.hpp"

//ReadOnly vars
char gRedis_host[48];
bool gDoKriging;
bool gDoCoMD;
int gDims[2];
double gDt[2];
double gGamma[3];
double gDelta[2];
CProxy_Main mainProxy;

class Flux : public CBase_Flux
{
	public:
	Flux(unsigned int step, unsigned int phase, unsigned int tid)
	{
		cloudFlux(gDoKriging, gDoCoMD, step, phase, tid, gRedis_host); 
		mainProxy.countCallBacks();
	}
};

class Main : public CBase_Main
{
	private:
		unsigned int nTasks;
		unsigned int curPhase;
		unsigned int curStep;
		unsigned int nSteps;

		void spawnTasks()
		{
			ckout << curStep << ": Doing " << nTasks << " fluxes" << endl;
			for(unsigned int i = 0; i < nTasks; i++)
			{
				CProxy_Flux::ckNew(curStep, curPhase, i);
			}
		}

	public:
	Main(CkArgMsg * m)
	{
		//<dim_x> <dim_y> <nsteps> <redis_server> <database error threshold> <Kriging error threshold> <Gaussian noise strength>
		//dimX dimY nSteps redis_server 
		if( m->argc != 5)
		{
			std::cerr <<  "./2D_DaDTest <dim_x> <dim_y> <nsteps> <redis_server>" << endl;
			CkExit();
		}
		//Set up parameters, mostly as globals
		mainProxy = thisProxy;
		gDoKriging = false;
		gDoCoMD = false;
		gDims[0] = atoi(m->argv[1]);
		gDims[1] = atoi(m->argv[2]);
		gDt[0] = 0.1;
		gDt[1] = 0.1;
		gDelta[0] = 1.0;
		gDelta[1] = 1.0;
		gGamma[0] = 0; //mom_gamma
		gGamma[1] = gGamma[0]; //strain_gamma
		gGamma[2] = 0.1 * gGamma[1];//en_gamma
		strcpy(gRedis_host, m->argv[4]);
		nSteps = atoi(m->argv[3]);

		//Initialize
		ckout << "Initializing " << gDims[0] << " by " << gDims[1] << " grid" << endl;
		initEverything(gDoKriging, gDoCoMD, gDims, gDt, gDelta, gGamma, gRedis_host);
		ckout << "Initialized" << endl;
		ckout << "Running for " << nSteps << " iterations" << endl;
	
		//Call the countCallBacks method
		countCallBacks();
	};

	void countCallBacks()
	{
		//Decrrement if needed
		if(nTasks != 0)
		{
			nTasks--;
		}
		//See if we move on
		if(nTasks == 0)
		{
			//Big switch
			switch(curPhase)
			{
				case 0:
					//Visualize Step
					ckout << nSteps << ": Vising to Verifying" << endl;
					outputVTK(gDoKriging, gDoCoMD, gDims, gDt, gDelta, gGamma, curStep, gRedis_host);
					//Check timestep
					if(curStep == nSteps)
					{
						//We are done
						ckout << "Ran for " << nSteps << " iterations" << endl;
						CkExit();
					}
					else
					{
						//Try short circuit
						if(tryShortCircuit(gDims, curStep, gRedis_host))
						{
							//Success
							ckout << curStep << ": Short Circuit Successful, on to the next step!" << endl;
							//Increment step
							curStep++;
							//Recursively calll
							///WARNING: A sufficiently long run that is recovering from a fault could be very bad due to too much recursion
							countCallBacks();
						}
						else
						{
							//Spawn phase 0
							ckout << curStep << ": First Flux" << endl;
							nTasks = prepFirstFlux(gDoKriging, gDoCoMD, gDims, gDt, gDelta, gGamma, curStep, gRedis_host);
							spawnTasks();
							curPhase = 1;
						}
					}
					break;
				case 1:
					//Phase 1
					ckout << curStep << ": Second Flux" << endl;
					nTasks = prepSecondFlux(gDoKriging, gDoCoMD, gDims, gDt, gDelta, gGamma, curStep, gRedis_host);
					spawnTasks();
					curPhase = 2;
					break;
				case 2:
					//Phase 2
					ckout << curStep << ": Third Flux" << endl;
					nTasks = prepThirdFlux(gDoKriging, gDoCoMD, gDims, gDt, gDelta, gGamma, curStep, gRedis_host);
					spawnTasks();
					curPhase = 3;
					break;
				case 3:
					//Phase 3
					ckout << curStep << ": Last Flux" << endl;
					nTasks = prepLastFlux(gDoKriging, gDoCoMD, gDims, gDt, gDelta, gGamma, curStep, gRedis_host);
					spawnTasks();
					curPhase = 4;
					break;
				case 4:
					ckout << curStep << ": Finish Step, no Fluxes" << endl;
					finishStep(gDoKriging, gDoCoMD, gDims, gDt, gDelta, gGamma, curStep, gRedis_host);
					//Increment step
					curStep++;
					//Reset phase
					curPhase = 0;
					//Recursively call this again
					countCallBacks();
					break;
				default:
					break;
			}
		}
	}
};

#include "cohmm_dad.def.h"

