/** header file for flux.cpp
 * **/
#ifndef FLUX_HPP
#define FLUX_HPP

#ifdef CIRCLE
#include <libcircle.h>
#endif

#ifdef CHARM
void fluxFn(fluxInput *in, fluxOutput *out, Input inp);
#elif defined (SERIAL) || (OMP)
void fluxFn(fluxInput *in, fluxOutput *out, std::map<std::string, std::vector<char *> > *dbCache, Input inp);
#elif CIRCLE
void fluxFn(CIRCLE_handle *handle);
#endif

#endif
