Database Assisted Distribution
=========

HMM code based on CoHMM-Red with the addition of a distributed shared memory implemented through the use of a database so as to facilitate fault tolerance and an alternative method of utilizing a wide range of
distributed runtime systems.

Based on code previously developed by Dominic Roehm and the students of the 2013 Los Alamos Co-Design Summer School. CoHMM is now maintained by ExMatEx: [Exascale Co-Design Center for Materials in Extreme Environments](exmatex.org).

Dependencies
------------

1. [redis](http://redis.io)

2. [hredis](https://github.com/redis/hiredis)

3. [Intel Math Kernel Library](https://software.intel.com/en-us/intel-mkl) or (another implementation of the cblas and lapacke interfaces)

4. Optionally [SWIG](http://www.swig.org/) if using the [Swift/T](swift-lang.org/Swift-T/) example

5. Optionally [Chunks and Tasks](http://chunks-and-tasks.org)

6. Optionally [Intel Concurrent Collections](https://icnc.github.io/)

7. Optionally [Libcircle](http://hpc.github.io/libcircle/)

8. Optionally [Charm++](http://charm.cs.illinois.edu/research/charm)

9. Optionally [twemproxy](https://github.com/twitter/twemproxy)

Installation
------------

1. `mkdir bld`

2. `cd bld`

3. `cmake ..`
    (As usual, edit CMakeCache.txt or use ccmake to specify non-standard location of libraries)

4. `make`

Execution
---------

1. start redis server in background (`redis-server &`)

2. run benchmark with arguments `<dim_x> <dim_y> <nsteps> <redis_server>`

   example command line with OpenMP driver:

   `./2D_OMP 100 100 1000 localhost`

Swift/T Example
---------
To build and execute the Swift/T example, ensure that the CoHMM library has been built in a subdirectory `bld` and execute `swiftTBuild` to generate .tic files

Chunks and Tasks
---------
If Chunks and Tasks is detected, a driver will be automatically built in the `2D_ChunksAndTasks` subdirectory.

Intel Concurrent Collections
---------
If Intel CnC is detected, a driver will be automatically built in the `2D_CnC` subdirectory.

Libcircle
---------
If Libcircle is detected, a driver will be automatically built in the `2D_Libcircle` subdirectory.

Charm++
---------
If Charm++ is detected, a driver will be automatically built in the `2D_Charm++` subdirectory.

OpenMP
---------
If OpenMP is detected, a driver will be automatically built in the `2D_OMP` subdirectory.

MPI
---------
If MPI is detected, a driver will be automatically built in the `2D_MPI` subdirectory.

Distributed Redis
---------
To run with the redis database distributed across multiple nodes of a cluster, twemproxy is required. From the driver's directory, use python to execute the `buildRedisHostfile.py` file in the `scripts` subdirectory with a command line parameter of the desired number of redis servers. Then source the `startRedisServers.sh` file in the working subdirectory and run the driver with the path to the created `lutFile` in place of the `<redis_server>` parameter. When done, source the `endRedisServers.sh` file to remotely kill the instances of `redis-server` and `nutcracker`.

Finer-Grain Fault Tolerance
---------
By default, CoHMM-DaD provides fault tolerance at the granularity of a timestep. If the simulation fails, CoHMM-DaD will be able to resume from the last complete timestep. If this is not sufficient, the `FINER_GRAIN_FT` flag may be enabled during CMake configuration and each Flux task will be verified and re-executed until it has been successfully completed or the simulation fails at the runtime level.

Copyright and license
---------------------

Los Alamos National Security, LLC (LANS) owns the copyright to CoHMM, which it identifies as LA-CC-2012-065 (ExMatEx: Scale-Bridging Materials Evaluation and Test Suite, Version 1). The license is BSD-sh with a "modifications must be indicated" clause.  See LICENSE.md for the full text.
