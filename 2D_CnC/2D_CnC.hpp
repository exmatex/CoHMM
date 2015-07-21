#ifndef TWOD_CNC_HPP
#define TWOD_CNC_HPP

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

struct Retry_Tag
{
public:
	//Constructors
	Retry_Tag();
	Retry_Tag(int step, int phase, int task, int round);
	//Attributes
	int step;
	int phase;
	int task;
	int round;
};
CNC_BITWISE_SERIALIZABLE(Retry_Tag);

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

struct Retry_Task
{
	int execute(const Retry_Tag &tag, DaDContext &c) const;
};

//Tuners for garbage collection
///TODO: Do this

//Context
struct DaDContext : public CnC::context<DaDContext>
{
	//Flux Tags
	CnC::tag_collection<Flux_Tag> fluxTags;
	CnC::tag_collection<Retry_Tag> retryTags;

	//Step
	CnC::step_collection<Flux_Task> fluxTask;
	CnC::step_collection<Retry_Task> retryTask;

	//The loneliest item
	CnC::item_collection<int, Flux_Item> globalItem;

	//No-arg constructor
	DaDContext();
};

#endif
