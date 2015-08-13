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
CProxy_Main gMainProxy;

class Flux : public CBase_Flux
{
	public:
	Flux(unsigned int step, unsigned int phase, unsigned int tid)
	{
		cloudFlux(gDoKriging, gDoCoMD, step, phase, tid, gRedis_host);
		gMainProxy.fluxCallBack();
	}
};

class Retry : public CBase_Retry
{
	public:
	Retry(unsigned int step, unsigned int phase, unsigned int round, unsigned int tid)
	{
		retryCloudFlux(gDoKriging, gDoCoMD, step, phase, tid, round, gRedis_host);
		gMainProxy.fluxCallBack();
	}
};

class Main : public CBase_Main
{
	Main_SDAG_CODE
	private:
		unsigned int nTasks;
		unsigned int curPhase;
		unsigned int curStep;
		int curRound;
		unsigned int nSteps;
		#ifdef FINERFT
			const bool fineGrainFT = true;
		#else
			const bool fineGrainFT = false;
		#endif
		int expectedTasks;
		unsigned int doneTasks;

		void actualCheckAndIncOrSet(int potSetVal)
		{
		    //Charm++ allegedly guarantees this will all be atomic if we call via checkAndIncOrSet
		    //Did we pass in a non-negative?
		    if(potSetVal >= 0)
		    {
		        //Yup, so this is the number of tasks
		        this->expectedTasks = potSetVal;
		    }
		    else
		    {
		        //Nope, so a task just finished: Increment
		        this->doneTasks++;
		    }

		    //See if we are done
		    if(this->expectedTasks == this->doneTasks)
		    {
		       //We are, so reset the variables
		       this->doneTasks = 0;
		       this->expectedTasks = -1;
		       //Then fire the callBack that returns control to the solver
			   //Call driver with phase offset by 0
		       this->cohmmDriver(0);

		    }
		    //We aren't, so wait for the next call
		}


		void spawnFluxes()
		{
			//ckout << curStep << ": Doing " << nTasks << " fluxes" << endl;
			for(unsigned int i = 0; i < nTasks; i++)
			{
				CProxy_Flux::ckNew(curStep, curPhase, i);
			}
			gMainProxy.checkAndIncOrSet(nTasks);
		}

		void spawnRetries()
		{
			//ckout << curStep << ": Redoing " << nTasks << " Tasks" << endl;
			for(unsigned int i = 0; i < nTasks; i++)
			{
				CProxy_Retry::ckNew(curStep, curPhase - 1, curRound, i);
			}
			gMainProxy.checkAndIncOrSet(nTasks);
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
		gMainProxy = thisProxy;
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

		//Init vars
		curPhase = 0;
		nTasks = 0;
		curStep = 0;
		curRound = -1;
		this->doneTasks = 0;
		this->expectedTasks = -1;
		//Call the countCallBacks method
		this->cohmmDriver(0);
	};

	void fluxCallBack()
	{
		gMainProxy.checkAndIncOrSet(-1);
	}

	void cohmmDriver(unsigned int phaseOffset)
	{
		//Are we done?
	    if(this->curStep == this->nSteps)
	    {
	        //One last print
			ckout << this->nSteps << ": Vising to Verifying" << endl;
	        outputVTK(gDoKriging, gDoCoMD, gDims, gDt, gDelta, gGamma, curStep, gRedis_host);
	        ckout << "Ran for " << this->nSteps << " iterations" << endl;
	        //We are done
	        {
	            double endTime = CkWallTimer();
	            FILE * eFile = fopen("etime.dat", "w");
	            fprintf(eFile, "%f", endTime);
	            fclose(eFile);
	        }
	        CkExit();
		}
		//Nope, so check what phase we are in
		else
		{
			//At this point, any fluxes have been processed and we are ready to proceed
			unsigned int switchPhase = this->curPhase + phaseOffset;
			switch(switchPhase)
			{
				case 0:
					ckout << this->curStep << ": Vising to Verifying" << endl;
					outputVTK(gDoKriging, gDoCoMD, gDims, gDt, gDelta, gGamma, curStep, gRedis_host);
					ckout << this->curStep <<  ": First Flux" << endl;
					nTasks = prepFirstFlux(gDoKriging, gDoCoMD, gDims, gDt, gDelta, gGamma, curStep, gRedis_host);
					ckout << this->curStep <<  ": Doing " << nTasks << " fluxes" << endl;
					spawnFluxes();
					curPhase = 1;
				break;
				case 1:
					ckout << this->curStep <<  ": Second Flux" << endl;
					nTasks = prepSecondFlux(gDoKriging, gDoCoMD, gDims, gDt, gDelta, gGamma, curStep, gRedis_host);
					ckout << this->curStep <<  ": Doing " << nTasks << " fluxes" << endl;
					spawnFluxes();
					curPhase = 2;
				break;
				case 2:
					ckout << this->curStep <<  ": Third Flux" << endl;
					nTasks = prepThirdFlux(gDoKriging, gDoCoMD, gDims, gDt, gDelta, gGamma, curStep, gRedis_host);
					ckout << this->curStep <<  ": Doing " << nTasks << " fluxes" << endl;
					spawnFluxes();
					curPhase = 3;
				break;
				case 3:
					ckout << this->curStep <<  ": Last Flux" << endl;
					nTasks = prepLastFlux(gDoKriging, gDoCoMD, gDims, gDt, gDelta, gGamma, curStep, gRedis_host);
					ckout << this->curStep <<  ": Doing " << nTasks << " fluxes" << endl;
					spawnFluxes();
					curPhase = 4;
				break;
				case 4:
					ckout << this->curStep <<  ": Finish Step, no Fluxes" << endl;
					finishStep(gDoKriging, gDoCoMD, gDims, gDt, gDelta, gGamma, curStep, gRedis_host);
					curStep++;
					curPhase = 0;
					//Recursivevly call again: finishStep handled phase and step
					cohmmDriver(0);
				break;
			}
		}
	}
/*
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
					ckout << curStep << ": Vising to Verifying" << endl;
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
						}
						else
						{
							//Spawn phase 0
							ckout << curStep << ": First Flux" << endl;
							nTasks = prepFirstFlux(gDoKriging, gDoCoMD, gDims, gDt, gDelta, gGamma, curStep, gRedis_host);
							spawnFluxes();
							curPhase = 1;
						}
					}
					break;
				case 1:
					//Phase 1
					curRound++;
					if(fineGrainFT == true)
					{
						nTasks = checkStepForFaults(gDims, curStep, curPhase - 1, curRound, gRedis_host);
					}
					if(nTasks != 0)
					{
						ckout << curStep << ": Redoing " << nTasks << " Tasks" << endl;
						spawnRetries();
					}
					else
					{
						curRound = -1;
						ckout << curStep << ": Second Flux" << endl;
						nTasks = prepSecondFlux(gDoKriging, gDoCoMD, gDims, gDt, gDelta, gGamma, curStep, gRedis_host);
						spawnFluxes();
						curPhase = 2;
					}
					break;
				case 2:
					//Phase 2
					curRound++;
					if(fineGrainFT == true)
					{
						nTasks = checkStepForFaults(gDims, curStep, curPhase - 1, curRound, gRedis_host);
					}
					if(nTasks != 0)
					{
						ckout << curStep << ": Redoing " << nTasks << " Tasks" << endl;
						spawnRetries();
					}
					else
					{
						curRound = -1;
						ckout << curStep << ": Third Flux" << endl;
						nTasks = prepThirdFlux(gDoKriging, gDoCoMD, gDims, gDt, gDelta, gGamma, curStep, gRedis_host);
						spawnFluxes();
						curPhase = 3;
					}
					break;
				case 3:
					//Phase 3
					curRound++;
					if(fineGrainFT == true)
					{
						nTasks = checkStepForFaults(gDims, curStep, curPhase - 1, curRound, gRedis_host);
					}
					if(nTasks != 0)
					{
						ckout << curStep << ": Redoing " << nTasks << " Tasks" << endl;
						spawnRetries();
					}
					else
					{
						curRound = -1;
						ckout << curStep << ": Last Flux" << endl;
						nTasks = prepLastFlux(gDoKriging, gDoCoMD, gDims, gDt, gDelta, gGamma, curStep, gRedis_host);
						spawnFluxes();
						curPhase = 4;
					}
					break;
				case 4:
					curRound++;
					if(fineGrainFT == true)
					{
						nTasks = checkStepForFaults(gDims, curStep, curPhase - 1, curRound, gRedis_host);
					}
					if(nTasks != 0)
					{
						ckout << curStep << ": Redoing " << nTasks << " Tasks" << endl;
						spawnRetries();
					}
					else
					{
						curRound = -1;
						ckout << curStep << ": Finish Step, no Fluxes" << endl;
						finishStep(gDoKriging, gDoCoMD, gDims, gDt, gDelta, gGamma, curStep, gRedis_host);
						//Increment step
						curStep++;
						//Reset phase
						curPhase = 0;
					}
					break;
				default:
					break;
			}
			//Are we waiting for any tasks?
			if(nTasks == 0)
			{
				//Nope, so recursively call
				///WARNING: THis could be bad for the stack if we are recovering from a fault after a lot of timesteps
				countCallBacks();
			}
		}
	}
	*/
};

#include "cohmm_dad.def.h"
