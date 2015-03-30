2D_Kriging
=========

HMM computation of 2D elastodynamics based on G.-S. Jiang, E. Tadmor, Nonoscillatory central schemes for multidimensional hyperbolic conservation laws, SIAM Journal on Scientific Computing 19 (6) (1998) 1892–1917 theory. The [Redis database](http://redis.io) is used via the [hredis](https://github.com/redis/hiredis) interface. The database supports the prediction of CoMD results with the ordinary Kriging method. 

Original authors of the source code include Dominic Roehm and Robert Pavel. CoHMM is now maintained by ExMatEx: [Exascale Co-Design Center for Materials in Extreme Environments](exmatex.org).

Dependencies
------------

1. [redis](http://redis.io)

2. [hredis](https://github.com/redis/hiredis)

3. [mkl](https://software.intel.com/en-us/intel-mkl) or (another implementation of the cblas and lapacke interfaces)

4. [boost](http://www.boost.org/)

###Additional for Charm++:

5. [charm](http://charm.cs.illinois.edu/software)

###Additional for CnC:

6. [cnc](https://software.intel.com/en-us/articles/intel-concurrent-collections-for-cc)

###Additional for libcircle:

7. [libcircle](https://github.com/hpc/libcircle)


Build
-----

For Charm++:
```
make SET=charm
```
For Cnc:
```
make SET=cnc
```
For OpenMP
```
make SET=omp
```
For Serial
```
make SET=serial
```
For libcircle
```
make SET=circle
```

Your might want to add `BOOST_INCLUDES=/path/to/boost/include` and `HIREDIS_INCLUDES=/path/to/hiredis/include` to the make call, e.g.
```
make -j4 SET=omp OPTFLAGS="-O3 -Werror" HIREDIS_INCLUDES=/usr/include BOOST_INCLUDES=/usr/include
```
to system provided boost and hiredis.

Execution
---------

###Starting the redis server

1. start redis server in background (`redis-server &`)

  1. distributed -> start redis using ./start_redis.sh hostfile

  2. hostfile includes the name of the available hosts (e.g. cn40;cn41;...)

###Charm++

1. run 2D_Kriging with './charmrun +p#processes  ++mpiexec ./2D_Kriging input.json +stacksize 512000'

  1. optional +isomalloc_sync

  2. example command line:

     distributed: ./charmrun +p96  ++mpiexec ./2D_Kriging input.json +stacksize 512000
     local: ./charmrun +p8  ++local ./2D_Kriging input.json +stacksize 512000

3. for use of only one local redis master use:

    
     distributed: ./charmrun +p96  ++mpiexec ./2D_Kriging input.json `hostname` +stacksize 512000

###CnC

1. MPI: run 2D_Kriging distributed with 'env DIST_CNC=MPI mpirun -n $(NPROCS x 48) ./2D_Kriging input.json' 

   on Darwin (no infiniband) 'env DIST_CNC=MPI mpirun -n 96 -env I_MPI_FABRICS shm:tcp -hostfile hostfile ~/2014/CoHMM/2D_Kriging input.json'

   hint: use the hostfile generated for Redis as input for Intel MPI, too.

2. SOCKETS: write client nodes in hostfile (eg.g cn30;cn31)

  1. run 'env DIST_CNC=SOCKETS HOST_FILE=hostfile CNC_SOCKET_HOST=./start.sh ./2D_Kriging input.json'

3. SHARED MEM: in Makefile remove -D_DIST_ option in CNC_FLAG and rebuild 

  1. run './2D_Kriging input.json'

###OpenMP

1. set $OMP_NUM_THREADS e.g. to the number of processors

2. run './2D_Kriging input.json'

###Serial

1. run './2D_Kriging input.json'

###libcircle

1. run 'mpirun -np XXX ./2D_Kriging input.json'

##Test Problems

###TP1

1. features: define DB, KRIGING, KR_DB, FLUSHDB

2. in flux.cpp enable define COMD and DB !!! (disable C_RAND)

3. for total execution time measurement only disable define OUTPUT in 2DKriging.cpp

4. inout.json macro_solver:

        {"id": "dim x",                     "value": 100},
        {"id": "dim y",                     "value": 10},
        {"id": "integration steps",         "value": 50},
        {"id": "redis headnode",            "value": "localhost"},
        {"id": "database threshold",        "value": 0.00001},
        {"id": "kriging err. threshold",    "value": 0.001},
        {"id": "Gaussian noise",            "value": 0.0},
        {"id": "dx",                        "value": 1.0},
        {"id": "dy",                        "value": 1.0},
        {"id": "dt_x",                      "value": 0.1},
        {"id": "dt_y",                      "value": 0.1},
        {"id": "test problem",              "value": 1}

###TP2

1. features: define DB, KRIGING, KR_DB, FLUSHDB

2. in flux.cpp enable define COMD and DB !!! (disable C_RAND)

3. for total execution time measurement only disable define OUTPUT in 2DKriging.cpp

4. inout.json macro_solver:

        {"id": "dim x",                     "value": 50},
        {"id": "dim y",                     "value": 50},
        {"id": "integration steps",         "value": 50},
        {"id": "redis headnode",            "value": "localhost"},
        {"id": "database threshold",        "value": 0.00001},
        {"id": "kriging err. threshold",    "value": 0.001},
        {"id": "Gaussian noise",            "value": 0.0},
        {"id": "dx",                        "value": 1.0},
        {"id": "dy",                        "value": 1.0},
        {"id": "dt_x",                      "value": 0.1},
        {"id": "dt_y",                      "value": 0.1},
        {"id": "test problem",              "value": 2}


##Generating a Trace of Database Accesses
---------

1. start redis server in background (`redis-server &`)

2. start `tools/redisTrace.sh` in the background (From the source directory, `./tools/redisTrace.sh &`)

3. run 2DKriging normally

4. when 2DKriging is complete, terminate the redisTrace process and examine the resulting redisTrace.log file

Copyright and license
---------------------

Los Alamos National Security, LLC (LANS) owns the copyright to CoHMM, which it identifies as LA-CC-2012-065 (ExMatEx: Scale-Bridging Materials Evaluation and Test Suite, Version 1). The license is BSD-sh with a "modifications must be indicated" clause.  See LICENSE.md for the full text.

