#ifndef COHMM_HPP
#define COHMM_HPP

#include <map>
#include <vector>
#include <array>

#include <hiredis.h>

#include "2DKriging.hpp"

//Abstract base class for CoHMM system/solver. All implementations implement this
///TODO: Figure out how to integrate with Swift. Singleton? Serializable?
///TODO: Do proper comments...
class CoHMM
{
public:
    ///TODO: Constructor(s). JSON can be one version
    CoHMM(unsigned int dimX, unsigned int dimY, unsigned int nSteps, bool doKriging, bool doCoMD, const char * redisHost = "localhost");
    CoHMM(unsigned int dims[2], double delta[2], double dt[2], double gamma[3], unsigned int nSteps, bool doKriging, bool doCoMD, const char * redisHost = "localhost");
    CoHMM();
    virtual ~CoHMM();

    void lateInit(unsigned int dims[2], double delta[2], double dt[2], double gamma[3], unsigned int nSteps, bool doKriging, bool doCoMD, const char * redisHost = "localhost");
    void lateInit(unsigned int dimX, unsigned int dimY, unsigned int nSteps, bool doKriging, bool doCoMD, const char * redisHost = "localhost");
    //Method to initialize the conserved fields
    virtual void initializeConservedFields(InitialConditions_e initCondition = InitialConditions_e::X);
    //Method to write output to VTK files
    virtual void outputVTK();
    //Methods to perform full step
    virtual void prepFirstFlux();
    virtual void prepSecondFlux();
    virtual void prepThirdFlux();
    virtual void prepLastFlux();
    virtual void finishStep();
    virtual int getNumberOfTasks();
    //Runtime implementation specific method to execute flux tasks and get results to fluxes and/or fields
    virtual void processFluxes() = 0;

protected:
    bool doCoMD;
    bool doKriging;
    unsigned int dims[2];
    double dt[2];
    double delta[2];
    double gamma[3];
    unsigned int nSteps;
    unsigned int curStep;
    unsigned int curPhase;
    //Map used to avoid duplicate flux calls
    std::map<Conserved, int> taskMap;
    //Vector used to collect processed flux data as needed
    std::vector<FluxFuture> futures;
    //Actual data. Each is a 2D array of size [2][dims[0]*dims[1]]
    std::array<std::vector<Node>, 2> fields;

    //Runtime implementation specific method to actually queue a task
    virtual bool queueTask(FluxIn &task, unsigned int taskID) = 0;
    //Method to apply kriging checks, db checks, and duplicate task checks
    virtual int prepTasks(std::vector<Node> &field);
};



#endif
