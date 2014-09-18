#ifndef KRIGING2D_HPP
#define KRIGING2D_HPP

#include "types.h"

struct gridPoint
{
	int x;
	int y;
};


struct fluxContext;

struct fluxInput
{
	char headNode[1024];
	Conserved *w;
	bool callCoMD;
};

struct fluxOutput
{
	double f[7];
	double g[7];
	double error;
};

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

template <class C> void freeClear( C & cntr ) {
    for ( typename C::iterator it = cntr.begin(); it != cntr.end(); ++it ) {
        delete[] *it;
    }
    cntr.clear();
}

/**randomgenerator which generates numbers [0,1]
 * * @param *rn Pointer to randomnumber array of the local node or particle
 * */
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
#endif

