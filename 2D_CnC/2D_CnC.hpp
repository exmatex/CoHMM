#ifndef TWOD_CNC_HPP
#define TWOD_CNC_HPP

#include <tuple>

#ifdef CNC_DIST
#include <cnc/dist_cnc.h>
#else
#include <cnc/cnc.h>
#endif

#include "CoHMM_DaD.hpp"

const unsigned int MAX_HOST_LENGTH = 64;

//Forward delcaration
struct DaDContext;

//Typedef for tag (step, phase, task)
typedef std::tuple<int, int, int> Flux_Tag; 

//The one singleton item
struct Flux_Item
{
	char redis_host[MAX_HOST_LENGTH];
	bool doKriging;
	bool doCoMD;
};
CNC_BITWISE_SERIALIZABLE(Flux_Item);

//Functors to act as steps
struct Flux_Task
{
	int execute(const Flux_Tag &tag, DaDContext &c) const;
};

//Tuners for garbage collection
///TODO: Do this

//Context
struct DaDContext : public CnC::context<DaDContext>
{
	//Flux Tags
	CnC::tag_collection<Flux_Tag> fluxTags;

	//Step
	CnC::step_collection<Flux_Task> fluxTask;

	//The loneliest item
	CnC::item_collection<int, Flux_Item> globalItem;

	//No-arg constructor
	DaDContext();
};

#endif

