#ifndef COHMM_CONTEXT_HPP
#define COHMM_CONTEXT_HPP

#ifdef CNC_DIST
#include <cnc/dist_cnc.h>
#else
#include <cnc/cnc.h>
#endif

#include <utility>

#include "2DKriging.hpp"
#include "RedisWrapper.hpp"


//Forward delcaration
struct CoHMMContext;

//Tags
// < <step, phase>, task>
typedef std::pair<std::pair<unsigned int, unsigned int>, unsigned int> Flux_Tag;
CNC_BITWISE_SERIALIZABLE(Flux_Tag);
//typedef unsigned int Flux_Tag[3];

//Items
struct Flux_Item
{
	char redis_host[RedisWrapper::MAX_HOST_LENGTH];
	bool doKriging;
	bool doCoMD;
    FluxIn task;
};
CNC_BITWISE_SERIALIZABLE(Flux_Item);

struct Flux_Result
{
    FluxOut result;
};
CNC_BITWISE_SERIALIZABLE(Flux_Result);

//Steps
struct Flux_Task
{
	int execute(const Flux_Tag &tag, CoHMMContext &c) const;
};

//Context
struct CoHMMContext : public CnC::context<CoHMMContext>
{
	//Flux Tags
	CnC::tag_collection<Flux_Tag> fluxTags;

	//Step
	CnC::step_collection<Flux_Task> fluxTasks;

	//Items
	CnC::item_collection<Flux_Tag, Flux_Item> fluxItems;
    CnC::item_collection<Flux_Tag, Flux_Result> fluxResults;

	//No-arg constructor
	CoHMMContext();
};

#endif
