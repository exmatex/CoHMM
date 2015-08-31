CnC Assisted Distribution
=========

HMM code based on CoHMM-DaD which replaces the database with Concurrent Collections collections used in place of database.

Based on code previously developed by Dominic Roehm and the students of the 2013 Los Alamos Co-Design Summer School. CoHMM is now maintained by ExMatEx: [Exascale Co-Design Center for Materials in Extreme Environments](exmatex.org).

Dependencies
------------

1. [redis](http://redis.io)

2. [hredis](https://github.com/redis/hiredis)

3. [Intel Math Kernel Library](https://software.intel.com/en-us/intel-mkl) or (another implementation of the cblas and lapacke interfaces)

4. [Intel Concurrent Collections](https://icnc.github.io/)


Installation
------------

1. `mkdir bld`

2. `cd bld`

3. `cmake ..`
    (As usual, edit CMakeCache.txt or use cmake to specify non-standard location of libraries)

4. `make`

Execution
---------

1. start redis server in background (`redis-server &`)

2. run benchmark with arguments `<dim_x> <dim_y> <nsteps> <redis_server>`

   example command line with OpenMP driver:

   `./2D_CnC 100 100 1000 localhost`

Copyright and license
---------------------

Los Alamos National Security, LLC (LANS) owns the copyright to CoHMM, which it identifies as LA-CC-2012-065 (ExMatEx: Scale-Bridging Materials Evaluation and Test Suite, Version 1). The license is BSD-sh with a "modifications must be indicated" clause.  See LICENSE.md for the full text.
