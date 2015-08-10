#ifndef COHMM_CNC_HPP
#define COHMM_CNC_HPP

#include <string>

#ifdef CNC_DIST
#include <cnc/dist_cnc.h>
#else
#include <cnc/cnc.h>
#endif

#include "2DKriging.hpp"


const unsigned int MAX_HOST_LENGTH = 64;

//Forward delcaration
struct CnCaDContext;

//Typedef for tag (step, phase, task)
///TODO: Probably need tuple
/// < <step, phase>, task>
typedef std::pair<std::pair<unsigned int, unsigned int>, unsigned int> Flux_Tag;
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
typedef NodeBlock Block_Item;
CNC_BITWISE_SERIALIZABLE(Block_Item);

//The result item
typedef FluxOut Result_Item;
CNC_BITWISE_SERIALIZABLE(Result_Item);

//The future item
typedef FluxFuture Future_Item;
CNC_BITWISE_SERIALIZABLE(Future_Item);

//The actual task item
typedef FluxIn Task_Item;
CNC_BITWISE_SERIALIZABLE(Task_Item);

//Functors to act as steps
struct Flux_Task
{
	int execute(const Flux_Tag &tag, CnCaDContext &c) const;
};


//Context
struct CnCaDContext : public CnC::context<CnCaDContext>
{
	//Flux Tags
	CnC::tag_collection<Flux_Tag> fluxTags;

	//Step
	CnC::step_collection<Flux_Task> fluxTask;

	//The loneliest item
	CnC::item_collection<int, Global_Singleton_Item> globalItem;

	//Block DB
	CnC::item_collection<Redis_Tag, Block_Item> blockItems;

	//Task DB
	CnC::item_collection<Redis_Tag, Task_Item> taskItems;

	//Future Items
	CnC::item_collection<Redis_Tag, Future_Item> futureItems;

	//No-arg constructor
	CnCaDContext();
};

#endif
