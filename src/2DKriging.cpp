/** 2D Macrosolver written by
 * Dominic Roehm (dominic.roehm@gmail.com)
 with 
 * Robert Pavel
 * Christoph Junghans
 */
/****************C-STANDARDS**************/
#define _XOPEN_SOURCE 700
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <string.h>
#include <cfloat>
#include <math.h>
#include <string>
#include <algorithm>
#include <iostream>
#include <ctime> 
/****************CPP-STANDARDS************/
#include <map>
#include <vector>
#include <list>
/****************CHARM++******************/
#ifdef CHARM
#include "main_charm.hpp"
#include "main.decl.h"
#include "krigingMod.decl.h"
#include "krigingMod.h"
#define printf CkPrintf 
#endif//CHARM
/****************C-STUFF******************/
#include "types.h"
extern "C"
{
#include <CoMD_lib.h>
#include <hiredis.h>
#ifdef OMP
#include <omp.h>
#endif//OMP
}
/****************CPP-STUFF***************/
#include "2DKriging.hpp"
#include "kriging.hpp"
#include "redisBuckets.hpp"
#include "output.hpp"
#include "flux.hpp"
/****************CNC******************/
#ifdef CNC
#define _CRT_SECURE_NO_DEPRECATE // to keep the VS compiler happy with TBB
#include <cnc/debug.h>
#endif//CNC
#ifdef CIRCLE
#include <sstream>
#include <cstring>
#include <boost/archive/text_oarchive.hpp>
#endif //CIRCLE
/****************FEATURES****************/
/* define specifies if CoMD is used or the linear "analytic" approach */
/* enable redis database */
/*****************OUTPUT****************/
/* use no output for benchmark */
//#define OUTPUT
/* sepcify additional vtk output (otherwise just gnuplot and ps files) */
//#define VTK_FIELDS
//#define VTK_COLORMAP
/* hidden feature */
//#define LOADBAR
/****************************************/
/*****************GLOBALS****************/
/****************************************/
//FIXME get rid of 'em
const int comdDigits = 4;
const int krigDigits = 4;
const double zeroThresh = 0.0001;
int counter = 0;
#ifdef CIRCLE
static std::vector<fluxInput> *_fluxInArgs;
void enqueue_fluxInArgs(CIRCLE_handle *handle) {
  for(std::vector<fluxInput>::iterator iter = _fluxInArgs->begin(); iter != _fluxInArgs->end(); iter++){
    //serialize
    std::ostringstream archive_stream;
    boost::archive::text_oarchive archive(archive_stream);
    archive << *iter;
    char str[archive_stream.str().size()+1];
    strcpy(str,archive_stream.str().c_str());
    handle->enqueue(&str[0]);
  }
}
#endif

/** initialize the struct containing all node quantities
 * @param node       grid_node containing the conserved and the fluxes
 * @param grid_size  number of grid nodes
 * **/
void init_nodes(Node* node, int grid_size){

  for(int i=0; i<grid_size; ++i){
 
    for(int j=0; j<7; ++j){
      node[i].f.f[j] = 0.0;
      node[i].g.f[j] = 0.0;
      node[i].w.w[j] = 0.0;
    }
    node[i].boundary = 0;
    //printf("init node struct test %g\n", node[i].w.momentum[1]);
  }
}

/** min_mod function to compute the discrete slopes, without using jacobian etc.
 * Jiang&Tadmor equ. (3.1'), (3.1`) and def. one on page 1901
 * @param w_plus   conserved or flux of node+1 (inpout) 
 * @param w        conserved or flux of node (inpout) 
 * @param w_minus  conserved or flux of node-1 (inpout) 
 * @param mm       result (spacial derivative of input) (output) 
 * **/
double min_mod(double w_plus, double w, double w_minus){

  double mm = 0.0;
  double w_plus_w = 0.0;
  double w_w = 0.0;
  double w_minus_w = 0.0;
  double theta = 1.0;

  w_plus_w   = theta*(w_plus - w);
  w_w = 0.5*(w_plus - w_minus);
  w_minus_w  = theta*(w - w_minus);

  if(w_plus_w > 0 && w_w > 0 && w_minus_w > 0){
      mm = fmin(w_plus_w, w_w);
      mm = fmin(mm, w_minus_w);
  }else if(w_plus_w < 0 && w_w < 0 && w_minus_w < 0){
      mm = fmax(w_plus_w, w_w);
      mm = fmax(mm, w_minus_w);
  }else{
      mm = 0.0;
  }

  return mm;
}
/** set intial values for strain, momentum- and energy density on grid nodes
 * @param node       grid_node containing the conserved and the fluxes (output)
 * @param grid_size  number of grid nodes (input)
 * @param in          lattice parameters (input)         
 * @param init       initial conserved values (inpout) 
 * **/
void init_conserved_fields(Node* node_a, Input in, int grid_size){

  //init field and threshold values
  Values init;
  init.strain[0]=1.0;
  init.strain[1]=0.0;
  init.strain[2]=0.0;
  init.strain[3]=1.0;
  init.momentum[0]=0.0;
  init.momentum[1]=0.0;
  init.momentum[2]=0.0;
  init.momentum[3]=0.0;
  init.energy=-0.295;
  init.rho=1.0; 

  int x,y;

  for (int i=0; i<grid_size; ++i) {
    index_to_xy(i, in, &x, &y);
    //A0[i] = A[i] = (i < dimX/2) ? 1.01 : 1.0; // small initial step in deformation gradient
    node_a[i].w.w[0] = 1.0;
    node_a[i].w.w[1] = 0.0;
    node_a[i].w.w[2] = 0.0;
    node_a[i].w.w[3] = 1.0;
    node_a[i].w.w[4] = 0.0;
    node_a[i].w.w[5] = 0.0;
    node_a[i].w.w[6] = init.energy;
    node_a[i].w.rho = init.rho;
    }
}
/** test problem 1: flat wave in x-direction
 *
 */
void init_test_problem1(Node* node_a, Input in, int grid_size){
  //xwave
  int x,y;
  for (int i=0; i<grid_size; ++i) {
    index_to_xy(i, in, &x, &y);

    if((x<in.dim_x/2+in.dim_x/10) && (x>=in.dim_x/2-in.dim_x/10)){
       node_a[i].w.w[0] = 1.04;
       //node_a[i].w.w[3] = 1.1;
    }
  }
}

/** test problem 2: circular inital expansion
 *
 */
void init_test_problem2(Node* node_a, Input in, int grid_size){
  int x,y;
  for (int i=0; i<grid_size; ++i) {
    index_to_xy(i, in, &x, &y);
    int r = floor(sqrt((x-in.dim_x/2)*(x-in.dim_x/2) + (y-in.dim_y/2)*(y-in.dim_y/2)));
    if(r <= 3){
      node_a[i].w.w[0] = 1.02;
      node_a[i].w.w[3] = 1.02;
      node_a[i].w.w[1] = 0.02;
      node_a[i].w.w[2] = -0.02;
    }
  }
}

/** test problem 3: thermal expansion
 *
 */
void init_test_problem3(Node* node_a, Input in, int grid_size){
  //FIXME move init struct into an new input struct
  Values init;
  init.energy=-0.295;
  int x,y;
  for (int i=0; i<grid_size; ++i) {
    index_to_xy(i, in, &x, &y);
    int r = floor(sqrt((x-in.dim_x/2)*(x-in.dim_x/2) + (y-in.dim_y/2)*(y-in.dim_y/2)));
    if(r <= 3){
      node_a[i].w.w[6] = init.energy+0.02;
    }
  }
}

/** zero flux boundaries
 *
 */
void boundaries(Input in, int grid_size, Node* node){
  for(int i=0; i<grid_size; ++i){
    if(node[i].boundary == 1){ 
      for(int j=0; j<7; ++j){
        node[i].f.f[j] = 0.0;
        node[i].g.f[j] = 0.0;
        node[i].w.w[0] = 1.0;
        node[i].w.w[1] = 0.0;
        node[i].w.w[2] = 0.0;
        node[i].w.w[3] = 1.0;
        node[i].w.w[4] = 0.0;
        node[i].w.w[5] = 0.0;
      }
    }
  }
}
/** shift back indices, shifts back the inidices with the help of double buffering (raceconditions!) 
 * @param node_b      grid_node containing the conserved and the fluxes (input)
 * @param grid_size   number of grid nodes (input)
 * @param in           lattice parameters (input)         
 * @param node_a      grid_node containing the conserved and the fluxes (output)
 * **/
void shift_back(Node* node_b, int grid_size, Input in, Node* node_a){

  int x = 0;
  int y = 0;
  int xpoypo;
  for(int i=0; i<grid_size; ++i){
    index_to_xy(i, in, &x, &y);
    // temp vars for indexing
    xpoypo = (x+1)%in.dim_x + in.dim_x*((y+1)%in.dim_y);
    node_b[i].f = node_a[xpoypo].f;
    node_b[i].g = node_a[xpoypo].g;
    node_b[i].w = node_a[xpoypo].w;
  }
  for(int i=0; i<grid_size; ++i){
    //index_to_xy(i, l, &x, &y);
    // temp vars for indexing
    node_a[i].f = node_b[i].f;
    node_a[i].g = node_b[i].g;
    node_a[i].w = node_b[i].w;
  }
}

/** checks if database values are within a certain threshold or exact 
 * @param *w0         field that needs to be computed
 * @param *wVec       field which is available from the database
 * @param dbThresh    threshold for the check
 * **/
bool ifConservedFieldsMatch(double * w0, std::vector<double *> * wVec, double dbThresh)
{
	if(wVec != NULL)
	{
		if(wVec->size() != 0)
		{
			double dist = 0.0;
			for(int i = 0; i < 7; i++)
			{
				dist += fabs(w0[i] - (*wVec)[0][i]);
        //printf("w0: %.16f, wIn: %.16f\n", w0[i], (*wVec)[0][i]);
			}
			if(dist <= dbThresh)
			{
        //printf("hit\n");
				return true;
			}
			else
			{
        //printf("false1\n");
				return false;
			}
		}
		else
		{
      //printf("false2\n");
			return false;
		}
	}
	else
	{
    //printf("false3\n");
		return false;
	}
}

/** executes the mapped fields in parallel 
 * @param *in         flux input contains information about field on specific node (input)
 * @return *out       returns flux output information about computed fluxes on specific node (output)
 * @param *dbCache    database cache pointer (input)
 * **/
template <typename T> void doParallelCalls(Node * fields, Node * fluxes, Input in, std::list<gridPoint> * comdTasks, std::list<gridPoint> * krigTasks, std::map<std::string, std::vector<char *> > *dbCache, Calls* ca, Tms *tm, T context)
{
	//fprintf(stdout, "Starting parallel \n");
	//fflush(stdout);

	int dim_x = in.dim_x;
	int dim_y = in.dim_y;

	//For both lists, enqueue unique tasks, map the rest
	std::map<Conserved, std::list<gridPoint>, ConservedComparator_c> taskMap;

	std::vector<fluxInput> fluxInArgs;
    char headNode[1024];
	strcpy(headNode, in.head_node.c_str());

	//Do kriging first so as to prioritize the cheaper solution
	for(std::list<gridPoint>::iterator iter = krigTasks->begin(); iter != krigTasks->end(); iter++)
	{
		int x = iter->x;
		int y = iter->y;
		Conserved *w = &fields[x + dim_x*y].w;
		taskMap[*w].push_back(*iter);
        fields[x + dim_x*y].f.ca = 6;
        ca->kPoints++;
		//Check if it was unique
		if(taskMap[*w].size() == 1)
		{
			//It was, add it to the task queue
#ifdef CNC
			fluxInArgs.push_back(fluxInput(*w, false, headNode, in.kr_threshold));
			int taskID = fluxInArgs.size();
			context->tags.put(taskID);
			context->fluxInp.put(taskID, fluxInput(*w, false, headNode, in.kr_threshold));
#elif CHARM
			fluxInArgs.push_back(fluxInput(*w, false));
#else
			fluxInArgs.push_back(fluxInput(*w, false, headNode, in.kr_threshold));

#endif
            ca->krig++;
            ca->kPoints--;
            fields[x + dim_x*y].f.ca = 5;
		}
	}
    //Collect CoMD tasks
	for(std::list<gridPoint>::iterator iter = comdTasks->begin(); iter != comdTasks->end(); iter++)
	{
		int x = iter->x;
		int y = iter->y;
		Conserved *w = &fields[x+dim_x*y].w;
		taskMap[*w].push_back(*iter);
        fields[x + dim_x*y].f.ca = 2;
        ca->cPoints++;
		//Check if it was unique
		if(taskMap[*w].size() == 1)
		{
			//It was, add it to the task queue
#ifdef CNC
			fluxInArgs.push_back(fluxInput(*w, true, headNode, in.kr_threshold));
			int taskID = fluxInArgs.size();
			context->tags.put(taskID);
			context->fluxInp.put(taskID, fluxInput(*w, true, headNode, in.kr_threshold));
#elif CHARM
			fluxInArgs.push_back(fluxInput(*w, true));
#else
			fluxInArgs.push_back(fluxInput(*w, true, headNode, in.kr_threshold));
#endif
            ca->comd++;
            ca->cPoints--;
            fields[x + dim_x*y].f.ca = 1;
		}
	}
#ifdef CNC
    //***************************************************//
    //CnC::debug::trace(context->steps); 
    //CnC::debug::trace(context->tags);
    //CnC::debug::collect_scheduler_statistics(*context);
    //***************************************************//
	context->wait();
#endif
#ifndef LOADBAR
	printf("About to run parallel kriging & MD: %d args, %d MD, %d kriging\n", (int)fluxInArgs.size(), ca->comd, ca->krig);
#endif
	
	//Empty both task lists
	comdTasks->clear();
	krigTasks->clear();

#ifdef CHARM
    CkReductionMsg *fluxOut;
    context.fluxFn(fluxInArgs, in, CkCallbackResumeThread((void*&)fluxOut));
    
    fluxOutput *fluxOutCharm = (fluxOutput*)fluxOut->getData();

    std::sort(fluxOutCharm, fluxOutCharm + fluxInArgs.size(), fluxOutOrder());
#endif
#ifdef OMP
    fluxOutput * fluxOutOmp = new fluxOutput[fluxInArgs.size()];
#pragma omp parallel for
    for(int i = 0; i < int(fluxInArgs.size()); i++){
        fluxFn(&fluxInArgs[i], &fluxOutOmp[i], dbCache, in);
    }
#endif//OMP
#ifdef CIRCLE
    _fluxInArgs=&fluxInArgs;
    CIRCLE_cb_create(&enqueue_fluxInArgs);
    CIRCLE_begin();
    //redisCommand(context,"sync");
#endif//CIRCLE

	//Process the results of OMP'd tasks
	for(int i = 0; i < int(fluxInArgs.size()); i++)
	{
#ifdef CNC
		fluxOutput fluxOutCnc;
		context->fluxOutp.get(i+1, fluxOutCnc);
        //printf("got val %d from task %d\n\n", fluxOut.index, i+1);
#elif CIRCLE
		std::vector<double *> wVec;
		std::vector<double *> fVec;
		std::vector<double *> gVec;
    		//get results from db
		if (fluxInArgs[i].callCoMD){
		  getSortedSubBucketNearZero(fluxInArgs[i].w.w, (char *)"comd", context, comdDigits, 1, &wVec, &fVec, &gVec, zeroThresh);
		} else {
		  getSortedSubBucketNearZero(fluxInArgs[i].w.w, (char *)"krig", context, comdDigits, 1, &wVec, &fVec, &gVec, zeroThresh);
		  if (! ifConservedFieldsMatch(fluxInArgs[i].w.w, &wVec, 0.0)){
		  getSortedSubBucketNearZero(fluxInArgs[i].w.w, (char *)"comd", context, comdDigits, 1, &wVec, &fVec, &gVec, zeroThresh);
            ca->kFail++;
		  }
		}
		if (! ifConservedFieldsMatch(fluxInArgs[i].w.w, &wVec, 0.0)){
		  printf("Redis error: COMD value not found\n");
		  exit(1);
		}
#endif
		//Write results to fluxes for all duplicates as this is good
		for(std::list<gridPoint>::iterator iter = taskMap[fluxInArgs[i].w].begin(); iter != taskMap[fluxInArgs[i].w].end(); iter++)
		{
			int x = iter->x;
			int y = iter->y;
			Fluxes * f = &fluxes[x + dim_x*y].f;
			Fluxes * g = &fluxes[x + dim_x*y].g;

#ifdef CHARM
			memcpy(f->f, &fluxOutCharm[i].f, sizeof(double)*7);
			memcpy(g->f, &fluxOutCharm[i].g, sizeof(double)*7);
            fluxInArgs[i].callCoMD = fluxOutCharm[i].callCoMD;
#elif CNC
	    	memcpy(f->f, fluxOutCnc.f, sizeof(double)*7);
			memcpy(g->f, fluxOutCnc.g, sizeof(double)*7);
            fluxInArgs[i].callCoMD = fluxOutCnc.callCoMD;
#elif OMP
            memcpy(f->f, &fluxOutOmp[i].f, sizeof(double)*7);
            memcpy(g->f, &fluxOutOmp[i].g, sizeof(double)*7);
            fluxInArgs[i].callCoMD = fluxOutOmp[i].callCoMD;
#elif CIRCLE
            memcpy(f->f, fVec[0], sizeof(double)*7);
            memcpy(g->f, gVec[0], sizeof(double)*7);
            //fluxInArgs[i].callCoMD = fluxOut[i].callCoMD;
#else
#error Something is wrong.
#endif
			if(fluxInArgs[i].callCoMD == true)
            {
                if(fields[x + dim_x*y].f.ca != 1 && fields[x + dim_x*y].f.ca != 2){
                fields[x + dim_x*y].f.ca = 7;
                //fields[x + dim_x*y].f.ca = 1;
                ca->kFail++;
                //printf("kfail++\n\n");
                }
            }
		}
#ifdef CHARM
        //timings
        tm->kr += fluxOutCharm[i].diffKr;
        tm->co += fluxOutCharm[i].diffCo;
#elif CNC
        //timings
        tm->kr += fluxOutCnc.diffKr;
        tm->co += fluxOutCnc.diffCo;
#elif OMP
        //timings
        tm->kr += fluxOutOmp[i].diffKr;
        tm->co += fluxOutOmp[i].diffCo;
#endif
#ifdef CIRCLE
		freeClear(wVec);
		freeClear(fVec);
		freeClear(gVec);
#endif
	}
#ifdef OMP
  delete[] fluxOutOmp;
#endif

#ifdef CHARM
  delete fluxOut;
#endif
  taskMap.clear();
#ifdef CNC
  context->unsafe_reset();
#endif
}

/** checks the gradient of the fields (used to decide whether trying kriging or not) 
 * @param *fields     fields on specific node (input)
 * @param in           lattice information (input)
 * **/
bool checkGradient(int x, int y, Node * fields, Input in)
{
	double dx = in.dx;
	double dy = in.dy;
	int dimX = in.dim_x;
	int dimY = in.dim_y;

	double delta = 0.0;
	int xP = (x + 1 + dimX) % dimX;
	int xM = (x - 1 + dimX) % dimX;
	int yP = (y + 1 + dimY) % dimY;
	int yM = (y - 1 + dimY) % dimY;

	Conserved * w = &fields[x + dimX*y].w;
	Conserved * wXP = &fields[xP + dimX*y].w;
	Conserved * wXM = &fields[xM + dimX*y].w;
	Conserved * wYP = &fields[x + dimX*yP].w;
	Conserved * wYM = &fields[x + dimX*yM].w;
	for(int i = 0; i < 7; i++)
	{
		double gradient[4];
		gradient[0] = fabs(w->w[i] - wXP->w[i] / dx);
		gradient[1] = fabs(w->w[i] - wXM->w[i] / dx);
		gradient[2] = fabs(w->w[i] - wYP->w[i] / dy);
		gradient[3] = fabs(w->w[i] - wYM->w[i] / dy);
		for(int j = 0; j < 4; j++)
		{
			if(gradient[j] > delta)
			{
				delta = gradient[j];
			}
		}
	}
	if(delta <= in.grad_threshold)
	{
		return true;
	}
	else
	{
		return false;
	}
}

/** main call, decides whether kriging is called or comd, maps all tasks for parallel execution 
 * @param *fields     fields  (input)
 * @param *fluxes     fluxes  (input)
 * @param in          lattice information (input)
 * **/
template <typename T> void doFluxes(Node* fields, Node* fluxes, int grid_size, Input* in, T context, redisContext *headRedis)
{
	//fprintf(stdout, "Starting doFluxes \n");
	//fflush(stdout);
	int dim_x = in->dim_x;
	int dim_y = in->dim_y;

    //measurement stuff
    Calls ca;
    ca.comd=0, ca.cPoints=0, ca.db=0, ca.kdb=0, ca.krig=0, ca.kPoints=0, ca.kFail=0;
    Tms tm;
    tm.db=0.0, tm.krDb=0.0, tm.ca=0.0, tm.kr=0.0, tm.co=0.0, tm.kr2=0.0, tm.co2=0.0;
    double startDb, stopDb;
    double startKrDb, stopKrDb;
    double startCa, stopCa;
    double startKr2, stopKr2;
    double startCo2, stopCo2;

	std::list<gridPoint> comdTasks;
	std::list<gridPoint> krigTasks;
	std::map<std::string, std::vector<char *> > dbCache; 

	//Iterate over every point to see who is getting krig'd and who is getting comd'd
	for(int y = 0; y < dim_y; y++)
	{
		for(int x = 0; x < dim_x; x++)
		{
			//Grab comd database
			std::vector<double *> wVec;
			std::vector<double *> fVec;
			std::vector<double *> gVec;

            bool useDB;
            if(in->redis_db == 1){
              startDb = getUnixTime();

		      getCachedSortedSubBucketNearZero(fields[x + dim_x*y].w.w, (char *)"comd", headRedis, comdDigits, 2, &wVec, &fVec, &gVec, zeroThresh, &dbCache);
		      //getSortedSubBucketNearZero(fields[x + dim_x*y].w.w, "comd", headRedis, comdDigits, 2, &wVec, &fVec, &gVec, zeroThresh);

		      //Check for exact value
              //printf("dbthresh %lf\n", dbT);
		      useDB = ifConservedFieldsMatch(fields[x+dim_x*y].w.w, &wVec, in->db_threshold);
		
              stopDb = getUnixTime();
              tm.db += stopDb - startDb;
            }else{
              useDB = false;
            }
			//If we found it
			if(useDB == true)
			{
                ca.db++;
                fields[x + dim_x*y].f.ca = 3;
				//Set the result from the database
				memcpy(fluxes[x+dim_x*y].f.f, fVec[0], sizeof(double)*7);
				memcpy(fluxes[x+dim_x*y].g.f, gVec[0], sizeof(double)*7);
			
            }
			//We did not
			else
			{
				//Check the gradient
                bool smallGradient;
                if(in->kriging == 1){
				  smallGradient = checkGradient(x, y, fields, *in);
                }else{
				  //If gradient is small enough
                  smallGradient = false;
                }
				if(smallGradient == true)
				{
					//Grab krig database
					std::vector<double *> wVecK;
					std::vector<double *> fVecK;
					std::vector<double *> gVecK;

                    if(in->kriging_db == 1){
                      startKrDb = getUnixTime();

			          getCachedSortedSubBucketNearZero(fields[x + dim_x*y].w.w, (char *)"krig", headRedis, krigDigits, 1, &wVecK, &fVecK, &gVecK, zeroThresh, &dbCache);
                      //if(diff > 1.0) printf("get kdb timing %lf\n", diff);
			          //fflush(stdout);
				      //Check for exact value with kriging key
				      bool useDB = ifConservedFieldsMatch(fields[x+dim_x*y].w.w, &wVecK, 0.0);
				      //If we found it
                      //if(int(wVec.size())>0)printf("wvecsize %i\n", int(wVec.size()));
                      stopKrDb = getUnixTime();
                      tm.krDb += stopKrDb - startKrDb;
                    }else{
                      bool useDB = false;
                    }
					if(useDB == true)
					{
                        ca.kdb++;
                        fields[x + dim_x*y].f.ca = 4;
						//Set the result from the database
						memcpy(fluxes[x+dim_x*y].f.f, fVecK[0], sizeof(double)*7);
						memcpy(fluxes[x+dim_x*y].g.f, gVecK[0], sizeof(double)*7);

					}
                    //use at least 2 w's for kriging
					else if(int(wVec.size()) >= 2)
					{
                    //printf("try krig\n");
						//We are gonna krig this
						gridPoint thePoint;
						thePoint.x = x;
						thePoint.y = y;
						krigTasks.push_back(thePoint);
					}
				    //It was not, so mark it for a CoMD
          else
          {
              //Add to CoMD list
              gridPoint thePoint;
              thePoint.x = x;
              thePoint.y = y;
              comdTasks.push_back(thePoint);
          }
          freeClear(wVecK);
          freeClear(fVecK);
          freeClear(gVecK);
				}
				//It was not, so mark it for a CoMD
				else
				{
					//Add to CoMD list
					gridPoint thePoint;
					thePoint.x = x;
					thePoint.y = y;
					comdTasks.push_back(thePoint);
				}
			}
            freeClear(wVec);
            freeClear(fVec);
            freeClear(gVec);
		} //end for loop x
	} // end for loop y

  //fprintf(stdout, "prepare parallel calls \n");
  //fflush(stdout);
  startCa = getUnixTime();
  //At this point, we know who is getting krig'd and who is getting comd'd, so let's do it
  doParallelCalls(fields, fluxes, *in, &comdTasks, &krigTasks, &dbCache, &ca, &tm, context);
  stopCa = getUnixTime();
  tm.ca = stopCa - startCa;

  if(ca.comd)tm.co /= ca.comd;
  if(ca.krig)tm.kr /= ca.krig;
  //print the calls
#ifdef OUTPUT
#ifdef VTK_COLORMAP
  printf_colormap_vtk(counter, fields, *in, grid_size);
#endif
  printf_calls(counter, ca);
  printf_timings(counter, tm);
#endif//OUTPUT
  counter++;

  if(ca.krig != 0 && ca.kFail == 0){
    if(in->grad_threshold < 10) in->grad_threshold *= 2;
  }
  if(ca.krig != 0 && ca.kFail != 0){
    in->grad_threshold *= 0.25;
  }
#ifdef OUTPUT
#ifndef LOADBAR
  printf("kFail %i, grad_threshold %f\n", ca.kFail, in->grad_threshold);
#endif
#endif//OUTPUT

  //clean the cache
  for(std::map<std::string, std::vector<char *> >::iterator iter = dbCache.begin(); iter != dbCache.end(); iter++)
  {  
     for(std::vector<char *>::iterator vecIt = iter->second.begin(); vecIt != iter->second.end(); vecIt++)
     {
         delete[] *vecIt;
     }
  }
  dbCache.clear();
  //fprintf(stdout, "Done doFluxes \n");
  //fflush(stdout);
  //At this point, everything is done
  return;
}

/** executes second order half time step with average conserved value (w^{1/2}_{jk})
 * @param node_a      grid_node containing the conserved and the fluxes (input)
 * @param node_b      grid_node containing the conserved and the fluxes (output)
 * @param grid_size   number of grid nodes (input)
 * @param l           lattice parameters (input)         
 * **/
template <typename T> void half_step_second_order(Node* node_a, Node* node_b, int grid_size, Input* in, T context, redisContext *headRedis)
{

  int dim_x = in->dim_x;
  int dim_y = in->dim_y;
  double df[7] = {0.0};
  double dg[7] = {0.0};
  double mu = in->dt_x/in->dx;
  double nu = in->dt_y/in->dy;
  int xpo, xmo, ypo, ymo, xpoypo;
  int x = 0;
  int y = 0;

  //Jiang&Tadmor '98 equ. (2.15) half step in time
  doFluxes(node_a, node_a, grid_size, in, context, headRedis);

  for(int i=0; i<grid_size; ++i){
    index_to_xy(i, *in, &x, &y);
    // vars for indexing
    xmo = ((x-1+dim_x)%dim_x) + dim_x*y;
    xpo = ((x+1)%dim_x) + dim_x*y; 
    ymo = x + dim_x*((y-1+dim_y)%dim_y);
    ypo = x + dim_x*((y+1)%dim_y);

    // compute w_n^1/2 (half timestep)
    for(int j=0; j<7; ++j){
      // compute derivative with min_mod operation
      df[j] = min_mod(node_a[xpo].f.f[j], node_a[i].f.f[j], node_a[xmo].f.f[j]);
      dg[j] = min_mod(node_a[ypo].g.f[j], node_a[i].g.f[j], node_a[ymo].g.f[j]);
      //compute new conserved values and store it in grid b
      node_b[i].w.w[j] = node_a[i].w.w[j] - 0.5*mu*df[j] - 0.5*nu*dg[j];
    }
  }//loop over grid
  //Calc fluxes for second halfstep (time and space)
  //prepare computation of Jiang&Tadmor '98 equ. (2.16)
  doFluxes(node_b, node_b, grid_size, in, context, headRedis);


  //temp vars
  // w is w' e.g. w_j1k1 = w'_{j+1,k+1} 
  double w_jk   = 0.;
  double w_j1k  = 0.;
  double w_jk1  = 0.;
  double w_j1k1 = 0.;

  // v is w` e.g. v_j1k1 = w`_{j+1,k+1} 
  double v_jk   = 0.;
  double v_j1k  = 0.;
  double v_jk1  = 0.;
  double v_j1k1 = 0.;

  // flux derivative x e.g. f_jk = f(w_{j,k})
  double f_jk   = 0.;
  double f_j1k  = 0.;
  double f_jk1  = 0.;
  double f_j1k1 = 0.;

  // flux derivative y e.g. g_jk = g(w_{j,k})
  double g_jk   = 0.;
  double g_j1k  = 0.;
  double g_jk1  = 0.;
  double g_j1k1 = 0.;

  //Jiang&Tadmor '98 equ. (2.16) full time step, half grid step
  for(int i=0; i<grid_size; ++i){
    index_to_xy(i, *in, &x, &y);
    // temp vars for indexing
    xmo = ((x-1+dim_x)%dim_x) + dim_x*y;
    xpo = ((x+1)%dim_x) + dim_x*y; 
    ymo = x + dim_x*((y-1+dim_y)%dim_y);
    ypo = x + dim_x*((y+1)%dim_y);
    xpoypo = (x+1)%dim_x + dim_x*((y+1)%dim_y);
    //loop over fluxes, new conserved values stored in grid b
    for(int j=0; j<7; ++j){
      //first part of summation
      node_b[xpoypo].w.w[j] = 0.25*(node_a[i].w.w[j] + node_a[xpo].w.w[j] + node_a[ypo].w.w[j] + node_a[xpoypo].w.w[j]);
      //temp vars for summation
      w_jk   = min_mod(node_a[xpo].w.w[j], node_a[i].w.w[j], node_a[xmo].w.w[j]);
      w_j1k  = min_mod(node_a[(x+2)%dim_x + dim_x*y].w.w[j], node_a[xpo].w.w[j], node_a[i].w.w[j]);
      w_jk1  = min_mod(node_a[xpoypo].w.w[j], node_a[ypo].w.w[j], node_a[(dim_x+x-1)%dim_x + (dim_x*((y+1)%dim_y))].w.w[j]);
      w_j1k1 = min_mod(node_a[(x+2)%dim_x + dim_x*((y+1)%dim_y)].w.w[j], node_a[xpoypo].w.w[j], node_a[ypo].w.w[j]);

      v_jk   = min_mod(node_a[ypo].w.w[j], node_a[i].w.w[j], node_a[ymo].w.w[j]);
      v_j1k  = min_mod(node_a[xpoypo].w.w[j], node_a[xpo].w.w[j], node_a[(x+1)%dim_x + dim_x*((dim_y+y-1)%dim_y)].w.w[j]);
      v_jk1  = min_mod(node_a[x + dim_x*((y+2)%dim_y)].w.w[j], node_a[ypo].w.w[j], node_a[i].w.w[j]);
      v_j1k1 = min_mod(node_a[(x+1)%dim_x + dim_x*((y+2)%dim_y)].w.w[j], node_a[xpoypo].w.w[j], node_a[xpo].w.w[j]);

      f_jk   = node_b[i].f.f[j];
      f_j1k  = node_b[xpo].f.f[j];
      f_jk1  = node_b[ypo].f.f[j];
      f_j1k1 = node_b[xpoypo].f.f[j];

      g_jk   = node_b[i].g.f[j];
      g_j1k  = node_b[xpo].g.f[j];
      g_jk1  = node_b[ypo].g.f[j];
      g_j1k1 = node_b[xpoypo].g.f[j];
      //second part of summation in equ (2.16)
      node_b[xpoypo].w.w[j] += 0.0625*(w_jk - w_j1k + w_jk1 - w_j1k1) - 0.5*mu*(f_j1k - f_jk + f_j1k1 - f_jk1);
      node_b[xpoypo].w.w[j] += 0.0625*(v_jk - v_jk1 + v_j1k - v_j1k1) - 0.5*nu*(g_jk1 - g_jk + g_j1k1 - g_j1k);
      //printf("nodeb[%i].w[%i] = %g\n", xpoypo, j,  node_b[xpoypo].w.w[j]);
    }// loop over fluxes 
  }//loop over grid

}

/**********************************************/
/** Main of the 2D macro solver, executes the full timestep loop
@para argc
@para argv
**/
#ifdef CHARM
void main_2DKriging(Input in, App CoMD, CProxy_krigingChare krigingChareProxy)
#else
void main_2DKriging(Input in, App CoMD)
#endif
//template <typename T> void main_2DKriging(Input in, T context)
{

#ifdef CNC
  flux_context context;
#endif//CNC
#ifdef OMP
  //dummy var
  int context;
#endif//OMP

#ifdef OUTPUT
//#######################################//
  //reset calls.dat
  int bufferSize = 256;
  char file_name[bufferSize];
  sprintf(file_name,"calls.dat");
  FILE *fn = fopen(file_name, "w");
  if (fn == NULL) {
    printf("Error writing file< %s >!\n", file_name);
    exit(0);
  }
  fprintf(fn, "#NO      COMD       COMD_P      DB      KR_DB      KR      KR_P      KR_FAIL\n");
  fclose(fn);

  //reset timings.dat
  sprintf(file_name,"timings.dat");
  FILE *fn2 = fopen(file_name, "w");
  if (fn2 == NULL) {
    printf("Error writing file< %s >!\n", file_name);
#ifdef CHARM
    CkExit();
#else
    exit(0);
#endif//Charm
  }
  fprintf(fn2, "#NO   COMD               COMD_P             DB             KR_DB            KR              KR_P          KR_FAIL\n");
  fclose(fn2);
#else
  //init total_time.dat
  int bufferSize = 256;
  char file_name[bufferSize];
  sprintf(file_name,"total_time.dat");
  FILE *fn3 = fopen(file_name, "a");
  if (fn3 == NULL) {
    printf("Error writing file< %s >!\n", file_name);
#ifdef CHARM
    CkExit();
#else
    exit(0);
#endif//Charm
  }

#endif//OUTPUT
  redisContext *headRedis;
  if(in.redis_db == 1){
    // Init Redis
    headRedis = redisConnect(in.head_node.c_str(), 6379);
    if(headRedis == NULL || headRedis->err)
    {
      printf("Redis error: %s\n", headRedis->errstr);
#ifdef CHARM
      CkExit();
#else
      exit(0);
#endif//Charm
    }
    if(in.flush_db == 1){
      redisCommand(headRedis,"flushall");
      redisCommand(headRedis,"flushdb");
      printf("Flushed Redis DB\n");
    }
  }else{
    headRedis = NULL;
  }

  //grid size (1D indexing)
  int grid_size = in.dim_x*in.dim_y;
  //double buffering
  Node *nodes_a = new Node[grid_size]; 
  Node *nodes_b = new Node[grid_size]; 
  //set all values to zero  
  init_nodes(nodes_a, grid_size);
  init_nodes(nodes_b, grid_size);

  //initial fiels on nodes
  init_conserved_fields(nodes_a, in, grid_size);

  if(in.test_problem == 1){
    init_test_problem1(nodes_a, in, grid_size);
    printf("\nExecuting test problem 1 ...\n\n");
  }else if(in.test_problem == 2){
    init_test_problem2(nodes_a, in, grid_size);
    printf("\nExecuting test problem 2 ...\n\n");
  }else if(in.test_problem == 3){
    init_test_problem3(nodes_a, in, grid_size);
    printf("\nExecuting test problem 3 ...\n\n");
  }else{
      printf("Error unknown test problem!\n Set 1,2 or 3 in input.json\n");
#ifdef CHARM
    CkExit();
#else
    exit(0);
#endif//Charm
  }

  //start total time measurement
  double ttime_start = getUnixTime();
  double ttime_stop, stime_start;
  int prev_step = 0;
  if(in.fault_tolerance == 1){
    prev_step = redisRead_fields(nodes_a, &in, headRedis);
  }
#ifdef CIRCLE
  MPI_Bcast(&in.int_steps, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&prev_step, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
  //set_boundaries(l, grid_size, nodes_a, nodes_b);
  for(int i=prev_step; i<in.int_steps; ++i){
#ifdef OUTPUT
    plot_fields(i, nodes_a, in);
#ifdef VTK_FIELDS
    printf_fields_vtk(i, nodes_a, in, grid_size);
#endif//VTK_FIELDS
#endif//OUTPUT
    if(in.fault_tolerance == 1){
      redisWrite_fields(nodes_a, in, headRedis, i);
      redisDel_fields(headRedis, i);
    }
    //main integration loop
    for (int j = 0; j < 1; j++){
#ifdef LOADBAR
        if(i>0 && in.int_steps >=50) loadBar(i,in.int_steps,50,100); 
        else if (in.int_steps <50) printf("Loadbar failed! Too few int_steps\n");
#else
	    fprintf(stdout, "============Integration Step: %d ==============\n", i);
	    fflush(stdout);
        stime_start = getUnixTime();
#endif
        //first half step
#ifdef CHARM
        half_step_second_order(nodes_a, nodes_b, grid_size, &in, krigingChareProxy, headRedis);
        //second half step
        half_step_second_order(nodes_b, nodes_a, grid_size, &in, krigingChareProxy, headRedis);
#elif CIRCLE
        half_step_second_order(nodes_a, nodes_b, grid_size, &in, headRedis, headRedis);
        //second half step
        half_step_second_order(nodes_b, nodes_a, grid_size, &in, headRedis, headRedis);
#else
        half_step_second_order(nodes_a, nodes_b, grid_size, &in, &context, headRedis);
        //second half step
        half_step_second_order(nodes_b, nodes_a, grid_size, &in, &context, headRedis);
#endif
        shift_back(nodes_b, grid_size, in, nodes_a);
        //apply boundary conditions
        //boundaries(l, grid_size, nodes_a);
#ifndef OUTPUT
        //stop total time measurement
        ttime_stop = getUnixTime();
        fprintf(fn3, "%g\n", ttime_stop-ttime_start);
        printf("time: %g time per step: %g\n", ttime_stop-ttime_start, ttime_stop-stime_start);
        fflush(fn3);
#endif//OUTPUT
    }
  }
  //exit output
  fprintf(stdout, "Done all steps \n");
  fprintf(stdout, "D^2KAS finished\n");
  fflush(stdout);
  //cleanup redis
  if(in.redis_db == 1){
    redisFree(headRedis);
  }
#ifndef OUTPUT
  //close total time file
  fclose(fn3);
#endif//OUTPUT
  //free node structs
  delete [] nodes_a;
  delete [] nodes_b;
}
/*
 * instant function templates
 */
#ifdef CHARM
//template void main_2DKriging(Input in, CProxy_krigingChare krigingChareProxy);
template void half_step_second_order(Node* node_a, Node* node_b, int grid_size, Input* in, CProxy_krigingChare krigingChareProxy, redisContext *headRedis);
template void doFluxes(Node* fields, Node* fluxes, int grid_size, Input* in, CProxy_krigingChare krigingChareProxy, redisContext *headRedis);
template void doParallelCalls(Node * fields, Node * fluxes, Input in, std::list<gridPoint> * comdTasks, std::list<gridPoint> * krigTasks, std::map<std::string, std::vector<char *> > *dbCache, Calls* ca, Tms *tm, CProxy_krigingChare krigingChareProxy);
#elif CNC
//template void main_2DKriging(Input in, flux_context* fluxText);
template void half_step_second_order(Node* node_a, Node* node_b, int grid_size, Input* in, flux_context* fluxText, redisContext *headRedis);
template void doFluxes(Node* fields, Node* fluxes, int grid_size, Input* in, flux_context* fluxText, redisContext *headRedis);
template void doParallelCalls(Node * fields, Node * fluxes, Input in, std::list<gridPoint> * comdTasks, std::list<gridPoint> * krigTasks, std::map<std::string, std::vector<char *> > *dbCache, Calls* ca, Tms *tm, flux_context* fluxText);
#else
//template void main_2DKriging(Input in, int * context);
#ifdef CIRCLE
template void doFluxes(Node* fields, Node* fluxes, int grid_size, Input* in, redisContext *context, redisContext *headRedis);
template void half_step_second_order(Node* node_a, Node* node_b, int grid_size, Input* in, redisContext *context, redisContext *headRedis);
template void doParallelCalls(Node * fields, Node * fluxes, Input in, std::list<gridPoint> * comdTasks, std::list<gridPoint> * krigTasks, std::map<std::string, std::vector<char *> > *dbCache, Calls* ca, Tms *tm, redisContext *context);
#else
template void doFluxes(Node* fields, Node* fluxes, int grid_size, Input* in, int* context, redisContext *headRedis);
template void half_step_second_order(Node* node_a, Node* node_b, int grid_size, Input* in, int* context, redisContext *headRedis);
template void doParallelCalls(Node * fields, Node * fluxes, Input in, std::list<gridPoint> * comdTasks, std::list<gridPoint> * krigTasks, std::map<std::string, std::vector<char *> > *dbCache, Calls* ca, Tms *tm, int* context);
#endif
#endif

#ifdef CHARM
#include "krigingMod.def.h"
#endif
