Database Assisted Distribution
=========

HMM code based on CoHMM-Red with the addition of a distributed shared memory implemented through the use of a database so as to facilitate fault tolerance and an alternative method of utilizing a wide range of
distributed runtime systems.

Based on code previously developed by Dominic Roehm and the students of the 2013 Los Alamos Co-Design Summer School. CoHMM is now maintained by ExMatEx: [Exascale Co-Design Center for Materials in Extreme Environments](exmatex.org).

Dependencies
------------

1. [redis](http://redis.io)

2. [hredis](https://github.com/redis/hiredis)

3. mkl or (another implementation of the cblas and lapacke interfaces)

4. Optionally [SWIG](http://www.swig.org/) if using the [Swift/T](swift-lang.org/Swift-T/) example

5. Optional [Chunks and Tasks](http://chunks-and-tasks.org)

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

2. run 2D_DaDTest with `./2D_dadTest <dim_x> <dim_y> <nsteps> <redis_server>`

   example command line:

   `./2D_dadTest 100 100 1000 localhost 0.0001 0.0001 0`

Swift/T Example
---------
To build and execute the Swift/T example, ensure that the CoHMM library has been built in a subdirectory `bld` and execute ./swiftTBuild to generate .tic files

Chunks and Tasks
---------
If Chunks and Tasks is detected, it will be automatically built alongside the 2D_DaDTest driver

Copyright and license
---------------------

Los Alamos National Security, LLC (LANS) owns the copyright to CoHMM, which it identifies as LA-CC-2012-065 (ExMatEx: Scale-Bridging Materials Evaluation and Test Suite, Version 1). The license is BSD-sh with a "modifications must be indicated" clause.  See LICENSE.md for the full text.

