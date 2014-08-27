/** header file for flux.cpp
 * **/
#ifndef FLUX_HPP
#define FLUX_HPP

#ifdef CHARM
void fluxFn(fluxInput *in, fluxOutput *out, Input inp);
#endif
#ifdef OMP
void fluxFn(fluxInput *in, fluxOutput *out, std::map<std::string, std::vector<char *> > *dbCache, double* startKr, double* stopKr, double* startCo, double* stopCo, Input inp);
#endif

#endif
