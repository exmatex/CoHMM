#include <iostream>
#include <cstdlib>

#include "cohmm_dad.decl.h"

#include "CoHMM_DaD.hpp"

class Main : public CBase_Main
{
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
		//Set up parameters
		bool doKriging = true;
		bool doCoMD = false;
		int dims[2] = {atoi(m->argv[1]), atoi(m->argv[2])};
		double dt[2] = {0.1, 0.1};
		double delta[2] = {1.0, 1.0};
		double gamma[3];
		gamma[0] = 0; //mom_gamma
		gamma[1] = gamma[0]; //strain_gamma
		gamma[2] = 0.1 * gamma[1];//en_gamma

		unsigned int numSteps = atoi(m->argv[3]);

		//Initialize
		ckout << "Initializing " << dims[0] << " by " << dims[1] << " grid" << endl;
		initEverything(doKriging, doCoMD, dims, dt, delta, gamma);
		ckout << "Initialized" << endl;
		//Loop
		ckout << "Running for " << numSteps << " iterations" << endl;
		for(unsigned int t = 0; t < numSteps; t++)
		{
			int nTasks;
			ckout << t << ": Vising to Verifying" << endl;
			outputVTK(doKriging, doCoMD, dims, dt, delta, gamma, t, m->argv[4]);
			//Do a short circuit test
			if(tryShortCircuit(dims, t, m->argv[4]))
			{
				//Short circuit succeeded
				ckout << t << ": Short Circuit Successful, on to the next step!" << endl;
			}
			else
			{
				ckout << t << ": First Flux" << endl;
				nTasks = prepFirstFlux(doKriging, doCoMD, dims, dt, delta, gamma, t, m->argv[4]);
				ckout << t << ": Doing " << nTasks << " fluxes" << endl;
				#pragma omp parallel for
				for(unsigned int i = 0; i < nTasks; i++)
				{
					cloudFlux(doKriging, doCoMD, t, 0, i, m->argv[4]);
				}
				ckout << t << ": Second Flux" << endl;
				nTasks = prepSecondFlux(doKriging, doCoMD, dims, dt, delta, gamma, t, m->argv[4]);
				ckout << t << ": Doing " << nTasks << " fluxes" << endl;
				#pragma omp parallel for
				for(unsigned int i = 0; i < nTasks; i++)
				{
					cloudFlux(doKriging, doCoMD, t, 1, i, m->argv[4]);
				}
				ckout << t << ": Third Flux" << endl;
				nTasks = prepThirdFlux(doKriging, doCoMD, dims, dt, delta, gamma, t, m->argv[4]);
				ckout << t << ": Doing " << nTasks << " fluxes" << endl;
				#pragma omp parallel for
				for(unsigned int i = 0; i < nTasks; i++)
				{
					cloudFlux(doKriging, doCoMD, t, 2, i, m->argv[4]);
				}
				ckout << t << ": Last Flux" << endl;
				nTasks = prepLastFlux(doKriging, doCoMD, dims, dt, delta, gamma, t, m->argv[4]);
				ckout << t << ": Doing " << nTasks << " fluxes" << endl;
				#pragma omp parallel for
				for(unsigned int i = 0; i < nTasks; i++)
				{
					cloudFlux(doKriging, doCoMD, t, 3, i, m->argv[4]);
				}
				ckout << t << ": Finish Step, no Fluxes" << endl;
				finishStep(doKriging, doCoMD, dims, dt, delta, gamma, t, m->argv[4]);
			}
		}
		//Final vis
		ckout << numSteps << ": Vising to Verifying" << endl;
		outputVTK(doKriging, doCoMD, dims, dt, delta, gamma, numSteps, m->argv[4]);
		ckout << "Ran for " << numSteps << " iterations" << endl;

		CkExit();
	};
};

#include "cohmm_dad.def.h"

