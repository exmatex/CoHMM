#ifndef KRIGING2D_HPP
#define KRIGING2D_HPP

#include <math.h>
#include <string.h>
#include <string>
#include "types.h"

#ifdef CNC
#ifdef _DIST_
#include <cnc/dist_cnc.h>
#else
#include <cnc/cnc.h>
#endif
//struct containing the flux information
struct flux_context;
#endif//CNC

#ifdef CIRCLE
#include <mpi.h>
#endif

//some extern constants
extern const int comdDigits;
extern const int krigDigits;
extern const double zeroThresh;

/** struct containing grid point location
 * **/
struct gridPoint
{
	int x;
	int y;
};

/** struct with input values for flux computation
 * **/
struct fluxInput
{
#if defined (CNC) || (OMP) || (CIRCLE) 
    char headNode[1024];
    double kr_threshold; 
#endif
    //enable char here or make add collection input
	Conserved w;
	bool callCoMD;

    fluxInput() { }
#if defined (CNC) || (OMP) || (CIRCLE)
    // constructor
    fluxInput(Conserved w_, bool callCoMD_, char* headNode_, double kr_threshold_)
    : w(w_), callCoMD(callCoMD_), kr_threshold(kr_threshold_)
    {
	    strcpy(headNode, headNode_);
    }

#elif CHARM
    fluxInput(Conserved w_, bool callCoMD_)
    : w(w_), callCoMD(callCoMD_)
    {
    }
  void pup(PUP::er &p) {
    p|w;
    p|callCoMD;
  }
#endif//CHARM
#ifdef CIRCLE
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version) {
    ar & headNode;
    ar & w;
    ar & callCoMD;
  }
#endif//CIRCLE
};
#ifdef CNC
CNC_BITWISE_SERIALIZABLE(fluxInput);
#endif

/** struct with output values for flux computation 
 * **/
typedef struct
{
	int index;
    double f[7];
	double g[7];
	double error;
    //timedifference[]
    double diffCo;
    double diffKr;
	bool callCoMD;
} fluxOutput;
#ifdef CNC
CNC_BITWISE_SERIALIZABLE(fluxOutput);

// declaration of compute step class
struct fluxFn
{
    // declaration of execute method goes here
    int execute( const int & t, flux_context & c ) const;
};
#if 1
struct fluxTuner: public CnC::step_tuner<>
{
	template<class dependency_consumer>
		void depends(const int &tag, flux_context &c, dependency_consumer & dC) const
		{
			dC.depends(c.fluxInp, tag);
		}
};
CNC_BITWISE_SERIALIZABLE(fluxTuner);
#endif
// this is our context containing collections and defining their depndencies
struct flux_context : public CnC::context< flux_context > // derive from CnC::context
{
    // the step collection for the instances of the compute-kernel
    //CnC::step_collection< fluxFn > steps;
	CnC::step_collection<fluxFn, fluxTuner> steps;
    // item collection holding the flux number(s)
	CnC::item_collection<int, fluxInput> fluxInp;
	CnC::item_collection<int, fluxOutput> fluxOutp;
    // tag collection to control steps
    CnC::tag_collection<int> tags;

    // constructor
    flux_context(): CnC::context< flux_context >(),
          // pass context to collection constructors
          steps( *this ),
		  fluxInp(*this),
		  fluxOutp(*this),
          tags( *this )

    {
        // prescribe compute steps with this (context) as argument
        tags.prescribes( steps, *this );
        // step consumes fluxInp
		steps.consumes(fluxInp);
        // step produces fluxOutp
		steps.produces(fluxOutp);
    }
};
#endif//CNC
/** bring charm flux output into flux input order
 * **/
struct fluxOutOrder
{
    bool operator()(const fluxOutput &lhs, const fluxOutput &rhs)
    {
        return lhs.index < rhs.index;
    }
};

/** comparing inputs
 * **/
struct ConservedComparator_c
{
	bool operator() (const Conserved& lhs, const Conserved& rhs) const
	{
		//Need to compare 7 different values... somehow
		//Build two "hashes"
		std::string ls((char *)&lhs, sizeof(Conserved));
		std::string rs((char *)&rhs, sizeof(Conserved));
		//Compare them
		return ls < rs;
	}
};

/** free list/map with dynamic content
 * **/
template <class C> void freeClear( C & cntr ) {
    for ( typename C::iterator it = cntr.begin(); it != cntr.end(); ++it ) {
        delete[] *it;
    }
    cntr.clear();
}

/**randomgenerator which generates numbers [0,1]
 * * @param *rn Pointer to randomnumber array of the local node or particle
 * **/
inline float random_01(int *state){

    const float mxi = 1.f/(float)(1ul<<31);
    int curr = *state;

    curr = 1103515245 * curr + 12345;
    float rno = (float)(curr & ((1ul<<31)-1))*mxi;
    *state = curr;
    return (double)rno;

}

/** gaussian random nummber generator for thermalisation
 * * @param *rn Pointer to randomnumber array of the local node node or particle
 * */
inline void gaussian_random(int *state, double* rn){

    double x1, x2;
    double r2, fac;
    /** On every second call two gaussian random numbers are calculated
    * via the Box-Muller transformation.*/
    /** draw two uniform random numbers in the unit circle */
    do {
        x1 = 2.*random_01(state)-1.;
        x2 = 2.*random_01(state)-1.;
        r2 = x1*x1 + x2*x2;
    } while (r2 >= 1. || r2 == 0.);

    /** perform Box-Muller transformation */
    fac = sqrt(-2.*log(r2)/r2);
    rn[0] = x2*fac;
    rn[1] = x1*fac;
}

/** precise timing measurement
 * **/
inline double getUnixTime(void)
{
    struct timespec ts;

    if(clock_gettime(CLOCK_REALTIME, &ts) != 0) return 0;

    return (((double) ts.tv_sec) + (double) (ts.tv_nsec / 1000000000.0));
}

#endif
