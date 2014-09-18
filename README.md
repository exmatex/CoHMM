2DKriging
=========

HMM computation of 2D elastodynamics based on G.-S. Jiang, E. Tadmor, Nonoscillatory central schemes for multidimensional hyperbolic conservation laws, SIAM Journal on Scientific Computing 19 (6) (1998) 1892–1917 theory. The [Redis database](http://redis.io) is used via the [hredis](https://github.com/redis/hiredis) interface. The database supports the prediction of CoMD results with the ordinary Kriging method. 

Original authors of the source code include Dominic Roehm and the students of the 2013 Los Alamos Co-Design Summer School. CoHMM is now maintained by ExMatEx: [Exascale Co-Design Center for Materials in Extreme Environments](exmatex.org).

Dependencies
------------

1. [redis](http://redis.io)

2. [hredis](https://github.com/redis/hiredis)

3. mkl or (another implementation of the cblas and lapacke interfaces)

Installation
------------

1. `mkdir build`

2. `cd build`

3. `cmake ..`
    (As usual, edit CMakeCache.txt or use ccmake to specify non-standard location of libraries)
    
4. `make`

Execution
---------

1. start redis server in background (`redis-server &`)

2. run 2DKriging with `./2DKriging <dim_x> <dim_y> <nsteps> <redis_server> <database error threshold> <Kriging error threshold> <Gaussian noise strength>`

   example command line:

   `./2DKriging 100 100 1000 localhost 0.0001 0.0001 0`

Generating a Trace of Database Accesses
---------

1. start redis server in background (`redis-server &`)

2. start tools/redisTrace.sh in the background (From the source directory, './tools/redisTrace.sh &')

3. run 2DKriging normally

4. when 2DKriging is complete, terminate the redisTrace process and examine the resulting redisTrace.log file
