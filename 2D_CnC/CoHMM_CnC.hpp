#ifndef COHMM_CNC_HPP
#define COHMM_CNC_HPP

#include <string>

#ifdef CNC_DIST
#include <cnc/dist_cnc.h>
#else
#include <cnc/cnc.h>
#endif

#include "2DKriging.hpp"

//Forward delcaration
struct CnCaDContext;

//Typedef for tag (step, phase, task)
struct Flux_Tag
{
public:
	//Constructors
	Flux_Tag();
	Flux_Tag(int step, int phase, int task);
	//Attributes
	int step;
	int phase;
	int task;
};
CNC_BITWISE_SERIALIZABLE(Flux_Tag);

typedef std::string Redis_Tag;
CNC_BITWISE_SERIALIZABLE(Redis_Tag);

//The one singleton item
struct Global_Singleton_Item
{
	char redis_host[MAX_HOST_LENGTH];
	bool doKriging;
	bool doCoMD;
};
CNC_BITWISE_SERIALIZABLE(Global_Singleton_Item);

///The rest are defined in 2DKriging.hpp
//The block item
CNC_BITWISE_SERIALIZABLE(Node);

//The result item
CNC_BITWISE_SERIALIZABLE(FluxOut);

//The future item
CNC_BITWISE_SERIALIZABLE(FluxFuture);

//The actual task item
CNC_BITWISE_SERIALIZABLE(FluxIn);


#endif
