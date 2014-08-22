/** flux file to call CoMD or calc of analytic response
 * **/
#ifdef CHARM
#include "main.decl.h"
#endif
#include <ctime> 
//#include <random>
#include "2DKriging.hpp"
#include "redisBuckets.hpp"
#include "kriging.hpp"
extern "C"
{
#include <CoMD_lib.h>
#include <hiredis.h>
#include <stdlib.h>
}
#include "types.h"
#define DB
//#define CoMD
//#define C_RAND
//FIXME get rid of that:
char  potDir[] = "../pots";


/** flux computation executed in parallel 
 * @param *in         flux input contains information about field on specific node (input)
 * @return *out       returns flux output information about computed fluxes on specific node (output)
 * @param *dbCache    database cache pointer (input)
 * **/
#ifdef CHARM
void fluxFn(fluxInput *in, fluxOutput *out, Input inp)
#elif CNC
int fluxFn::execute(const int & id, flux_context & fluxText) const
#else
void fluxFn(fluxInput *in, fluxOutput *out, std::map<std::string, std::vector<char *> > *dbCache, double* startKr, double* stopKr, double* startCo, double* stopCo, Input inp)
#endif
{
#ifdef CNC
	fluxInput inVal;
	fluxText.fluxInp.get(id, inVal);
	fluxInput * in = &inVal;
	fluxOutput outVal;
	fluxOutput* out = &outVal;
    Input inp;
    inp.head_node = in->headNode;
    //printf("id %d cnc input %s\n", id, inp.head_node.c_str());
    //printf("id %d cnc input %lf\n", id, in->w.w[0]);
#endif//CNC
//#ifdef CHARM
#if 1
	//Prep outputs
	out->error = 0.0;
    //printf("input %s\n", in->headNode);
    //printf("input %lf\n", in->w.w[0]);

#ifdef DB
	//Make redis context
	redisContext * rTask;
	rTask = redisConnect(inp.head_node.c_str(), 6379);
	//CkPrintf("Redis headnode: %s\n", headNode.c_str());
	if(rTask == NULL || rTask->err)
	{
		//printf("Redis error: %s\n", rTask->errstr);
		//printf("Redis error: %s\n", in->headNode);
		printf("Redis connection delay\n");
        //sleep(20);
	    rTask = redisConnect(inp.head_node.c_str(), 6379);
		//return NULL;
	}
    //if(diff1 > 1.0) printf("CONNECTION timing %lf\n", diff1);
	if(rTask == NULL || rTask->err)
	{
		printf("Redis error: %s\n", rTask->errstr);
		printf("Redis error: %s\n", inp.head_node.c_str());
#ifdef CHARM
        CkExit();
#else
        exit(0);
#endif
		//return NULL;
	}
#endif
	//Is this CoMD?
	if(in->callCoMD == true)
	{
        double startCo = getUnixTime();
		//Call CoMD
		//4 strains
		double strain_xx = in->w.w[0];
		double strain_xy = in->w.w[1];
		double strain_yx = in->w.w[2];
		double strain_yy = in->w.w[3];
		//2 momentum_fluxes
		double momentum_x = in->w.w[4];
		double momentum_y = in->w.w[5];
		//enery_flux
		double energy = in->w.w[6];

#ifdef CoMD
		CoMD_input theInput;

		strcpy(theInput.potDir,potDir);
		strcpy(theInput.potName,"Cu01.eam.alloy");
		strcpy(theInput.potType,"setfl");
		theInput.doeam = 1; 
		theInput.nx = 6;
		theInput.ny = 6;
		theInput.nz = 6;
		theInput.nSteps = 1000;
		theInput.printRate = 1; 
		//MUST SPECIFY THE FOLLOWING               
		theInput.dt = 10.0;
		theInput.lat = 3.6186;
		theInput.temperature = 0;
		theInput.initialDelta = 0.0;
		theInput.defGrad[0] = strain_xx;
		theInput.defGrad[1] = strain_xy;
		theInput.defGrad[2] = strain_yx;
		theInput.defGrad[3] = strain_yy;
		theInput.enDens = energy;
		theInput.momDens[0] = momentum_x;
		theInput.momDens[1] = momentum_y;
		theInput.momDens[2] = 0.0;
		//theInput.rank = rank;
		//theInput.calls = calls;
		theInput.rank = 0;
		theInput.calls = 0;

		//Call it
		CoMD_return theRet = CoMD_lib(&theInput);

		//4 stresses
		double rho = 1.0;
		out->f[0] = momentum_x/rho;
		out->f[1] = momentum_y/rho;
		out->f[2] = 0.0;
		out->f[3] = 0.0;
		out->f[4] = theRet.stressXX;
		out->f[5] = theRet.stressXY;
		out->g[0] = 0.0;
		out->g[1] = 0.0;
		out->g[2] = momentum_x/rho;
		out->g[3] = momentum_y/rho;
		out->g[4] = theRet.stressYX;
		out->g[5] = theRet.stressYY;
		out->f[6] = -theRet.energyDensX;
		out->g[6] = -theRet.energyDensY;
#elif defined(C_RAND)
        //gaussian noise, gamma = strength
        unsigned int tid = omp_get_thread_num();
        unsigned s_seed = 1103515245 * tid + floor(1000*(*startCo));
        std::srand(s_seed);
        int seed = rand()%10000; 
        //fake comd values
	    //4 stresses
	    double rho = 1.0;
	    //out.f[0] = theRet.momX;
        double *num;
        num = new double[2];
        gaussian_random(&seed, num);
        //printf("num %lf %lf\n", num[0], num[1]);
        out->f[0] = momentum_x/rho + (inp.noise*num[0]);
        out->f[1] = momentum_y/rho + (inp.noise*num[1]);
        out->f[2] = 0.0;
        //o->t.f[2] = theRet.momY;
        out->f[3] = 0.0;
        gaussian_random(&seed, num);
        out->f[4] = strain_xx-1 + 0.75*(strain_yy-1) + (strain_gamma*num[0]);
        out->f[5] = 1.9*strain_xy + (strain_gamma*num[1]);
        out->g[0] = 0.0;
        //o->t.g[1] = theRet.momX;
        out->g[1] = 0.0;
        gaussian_random(&seed, num);
	    out->g[2] = momentum_x/rho + (inp.noise*num[0]);
	    //o->t.g[3] = theRet.momY;
	    out->g[3] = momentum_y/rho + (inp.noise*num[1]);
        gaussian_random(&seed, num);
	    out->g[4] = 1.9*strain_yx + (strain_gamma*num[0]);
	    out->g[5] = strain_yy-1 + 0.75*(strain_xx-1) + (strain_gamma*num[1]);
	    //f->6] = g[6] = -0.295;
	    out->f[6] = -out->f[0]*sqrt(out->f[4]*out->f[4] + out->f[5]*out->f[5]);
	    out->g[6] = -out->g[3]*sqrt(out->g[4]*out->g[4] + out->g[5]*out->g[5]);
        delete[] num;
#else
    //fake comd values
		//4 stresses
		double rho = 1.0;
		//out.f[0] = theRet.momX;
		out->f[0] = momentum_x/rho;
		out->f[1] = momentum_y/rho;
		out->f[2] = 0.0;
		//o->t.f[2] = theRet.momY;
		out->f[3] = 0.0;
		out->f[4] = strain_xx-1 + 0.75*(strain_yy-1);
		out->f[5] = 1.9*strain_xy;
		out->g[0] = 0.0;
		//o->t.g[1] = theRet.momX;
		out->g[1] = 0.0;
		out->g[2] = momentum_x/rho;
		//o->t.g[3] = theRet.momY;
		out->g[3] = momentum_y/rho;
		out->g[4] = 1.9*strain_yx;
		out->g[5] = strain_yy-1 + 0.75*(strain_xx-1);
		//f->6] = g[6] = -0.295;
		out->f[6] = -out->f[0]*sqrt(out->f[4]*out->f[4] + out->f[5]*out->f[5]);
		out->g[6] = -out->g[3]*sqrt(out->g[4]*out->g[4] + out->g[5]*out->g[5]);
#endif//COMD
#ifdef DB
		//Write result to database
		putData(in->w.w, out->f, out->g, (char *)"comd", rTask, comdDigits);		
        double stopCo = getUnixTime();
        //if(diff1 > 1.0) printf("put timing %lf\n", diff1);
        //out->diffCo = stopCo - startCo;
#endif
	}
	//It is kriging time!
	else
	{
        double startKr = getUnixTime();
        //printf("krigingtime\n");
		//Get data for kriging
		std::vector<double *> oldWs;
		std::vector<double *> oldFs;
		std::vector<double *> oldGs;
		//getCachedSortedSubBucketNearZero(in->w.w, (char *)"comd", rTask, comdDigits, 10, &oldWs, &oldFs, &oldGs, zeroThresh, dbCache); 
#ifdef DB
		getSortedSubBucketNearZero(in->w.w, (char*)"comd", rTask, comdDigits, 10, &oldWs, &oldFs, &oldGs, zeroThresh); 
#endif

		//Call Kriging on each point
        double *resF = new double[2];
        double *resG = new double[2];
        int info;
		for(int i = 0; i < 7; i++)
		{
			///TODO: Consider refactoring this in to kriging.cpp
			std::vector<double> oldF;
			std::vector<double> oldG;
			for(int j = 0; j < oldFs.size(); j++)
			{
				oldF.push_back(oldFs[j][i]);
				oldG.push_back(oldGs[j][i]);
			}
			info = kriging(in->w.w, 7, oldWs, oldF, resF);	
			info = kriging(in->w.w, 7, oldWs, oldG, resG);	
			//Write result
			out->f[i] = resF[0];
			out->g[i] = resG[0];
			//Set error
			if(resF[1] > out->error)
			{
				out->error = resF[1];
			}
			if(resG[1] > out->error)
			{
				out->error = resG[1];
			}
		}
        freeClear(oldWs);
        freeClear(oldFs);
        freeClear(oldGs);
        delete[] resF;
        delete[] resG;
        //check for error and call CoMd if necessary 
        //printf("error %lf\n", out->error);
        if(out->error > inp.kr_threshold)
        {
            double stopKr = getUnixTime();
	        out->error = 0.0;
            double startCo = getUnixTime();
		    //Call CoMD
		    //4 strains
		    double strain_xx = in->w.w[0];
		    double strain_xy = in->w.w[1];
		    double strain_yx = in->w.w[2];
		    double strain_yy = in->w.w[3];
		    //2 momentum_fluxes
		    double momentum_x = in->w.w[4];
		    double momentum_y = in->w.w[5];
		    //enery_flux
		    double energy = in->w.w[6];

#ifdef CoMD
		    CoMD_input theInput;

		    ///FIXME: Fix potdir
		    strcpy(theInput.potDir,potDir);
		    strcpy(theInput.potName,"Cu01.eam.alloy");
		    strcpy(theInput.potType,"setfl");
		    theInput.doeam = 1; 
		    theInput.nx = 6;
		    theInput.ny = 6;
		    theInput.nz = 6;
		    theInput.nSteps = 1000;
		    theInput.printRate = 1; 
		    //MUST SPECIFY THE FOLLOWING               
		    theInput.dt = 10.0;
		    theInput.lat = 3.6186;
		    theInput.temperature = 0;
		    theInput.initialDelta = 0.0;
		    theInput.defGrad[0] = strain_xx;
		    theInput.defGrad[1] = strain_xy;
		    theInput.defGrad[2] = strain_yx;
		    theInput.defGrad[3] = strain_yy;
		    theInput.enDens = energy;
		    theInput.momDens[0] = momentum_x;
		    theInput.momDens[1] = momentum_y;
		    theInput.momDens[2] = 0.0;
		    //theInput.rank = rank;
		    //theInput.calls = calls;
		    theInput.rank = 0;
		    theInput.calls = 0;

		    //Call it
		    CoMD_return theRet = CoMD_lib(&theInput);

		    //4 stresses
		    double rho = 1.0;
		    out->f[0] = momentum_x/rho;
		    out->f[1] = momentum_y/rho;
		    out->f[2] = 0.0;
		    out->f[3] = 0.0;
		    out->f[4] = theRet.stressXX;
		    out->f[5] = theRet.stressXY;
		    out->g[0] = 0.0;
		    out->g[1] = 0.0;
		    out->g[2] = momentum_x/rho;
		    out->g[3] = momentum_y/rho;
		    out->g[4] = theRet.stressYX;
		    out->g[5] = theRet.stressYY;
		    out->f[6] = -theRet.energyDensX;
		    out->g[6] = -theRet.energyDensY;
            //delete &theRet;
#elif defined(C_RAND)
            //gaussian noise, gamma = strength
            unsigned int tid = omp_get_thread_num();
            unsigned s_seed = 1103515245 * tid + floor(1000*(*startCo));
            std::srand(s_seed);
            int seed = rand()%10000; 
            //fake comd values
		    //4 stresses
		    double rho = 1.0;
		    //out.f[0] = theRet.momX;
            double *num;
            num = new double[2];
            gaussian_random(&seed, num);
            //printf("num %lf %lf\n", num[0], num[1]);
		    out->f[0] = momentum_x/rho + (inp.noise*num[0]);
		    out->f[1] = momentum_y/rho + (inp.noise*num[1]);
		    out->f[2] = 0.0;
		    //o->t.f[2] = theRet.momY;
		    out->f[3] = 0.0;
            gaussian_random(&seed, num);
		    out->f[4] = strain_xx-1 + 0.75*(strain_yy-1) + (strain_gamma*num[0]);
		    out->f[5] = 1.9*strain_xy + (strain_gamma*num[1]);
		    out->g[0] = 0.0;
		    //o->t.g[1] = theRet.momX;
		    out->g[1] = 0.0;
            gaussian_random(&seed, num);
		    out->g[2] = momentum_x/rho + (inp.noise*num[0]);
		    //o->t.g[3] = theRet.momY;
		    out->g[3] = momentum_y/rho + (inp.noise*num[1]);
            gaussian_random(&seed, num);
		    out->g[4] = 1.9*strain_yx + (strain_gamma*num[0]);
		    out->g[5] = strain_yy-1 + 0.75*(strain_xx-1) + (strain_gamma*num[1]);
		    //f->6] = g[6] = -0.295;
		    out->f[6] = -out->f[0]*sqrt(out->f[4]*out->f[4] + out->f[5]*out->f[5]);
		    out->g[6] = -out->g[3]*sqrt(out->g[4]*out->g[4] + out->g[5]*out->g[5]);
            delete[] num;
#else
            //fake comd values
    		//4 stresses
    		double rho = 1.0;
    		//out.f[0] = theRet.momX;
    		out->f[0] = momentum_x/rho;
    		out->f[1] = momentum_y/rho;
    		out->f[2] = 0.0;
    		//o->t.f[2] = theRet.momY;
    		out->f[3] = 0.0;
    		out->f[4] = strain_xx-1 + 0.75*(strain_yy-1);
    		out->f[5] = 1.9*strain_xy;
    		out->g[0] = 0.0;
    		//o->t.g[1] = theRet.momX;
    		out->g[1] = 0.0;
    		out->g[2] = momentum_x/rho;
    		//o->t.g[3] = theRet.momY;
    		out->g[3] = momentum_y/rho;
    		out->g[4] = 1.9*strain_yx;
            out->g[5] = strain_yy-1 + 0.75*(strain_xx-1);
            //f->6] = g[6] = -0.295;
            out->f[6] = -out->f[0]*sqrt(out->f[4]*out->f[4] + out->f[5]*out->f[5]);
            out->g[6] = -out->g[3]*sqrt(out->g[4]*out->g[4] + out->g[5]*out->g[5]);
#endif//COMD
            //mark as CoMD'd
            in->callCoMD = true;
#ifdef DB
		    //Write result to database
		    putData(in->w.w, out->f, out->g, (char *)"comd", rTask, comdDigits);		
#endif
             double stopCo = getUnixTime();
        }
        else
        {
#ifdef DB
            //Write result to database
            putData(in->w.w, out->f, out->g, (char *)"krig", rTask, krigDigits);		
#endif
            double stopKr = getUnixTime();
        }
	}
    //out->diffKr = stopKr - startKr;

#ifdef DB
	//All pertinent values are set, let's disconneect from Redis and return
	redisFree(rTask);
#endif
#endif
//#endif//CHARM
#ifdef CNC
#if 0
	outVal.index=id;
    for(int i=0;i<7;++i){
      outVal.f[i]=0.0;
	  outVal.g[i]=0.0;
    }
	outVal.error=1.0;
    outVal.diffCo=0.0;
    outVal.diffKr=0.0;
#endif
	outVal.index=id;
    //printf("outval %d\n", outVal.index);
	fluxText.fluxOutp.put(id, outVal);

    return CnC::CNC_Success;
#endif
}
