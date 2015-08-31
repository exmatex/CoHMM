#ifndef FLUXPUPWRAPPER_HPP
#define FLUXPUPWRAPPER_HPP

#include "2DKriging.hpp"

#include "charm++.h"

PUPbytes(FluxIn);
PUPbytes(FluxOut);
PUPbytes(Node);
PUPbytes(FluxFuture);

#endif
