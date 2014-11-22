/** file containing the CNC main function
 * calls 2D Kriging main
 * **/
#include "main_generic.hpp"
#include "input.hpp"
#include "types.h"
#include "2DKriging.hpp"

#include <cstring>
#include <stdio.h>
#include <stdlib.h>
#ifdef OMP
#include <omp.h>
#endif//OMP

#ifdef CIRCLE
#include <mpi.h>
#include <libcircle.h>
#include "flux.hpp"
#endif//CIRCLE

int main( int argc, char* argv[] )
{

#ifdef _DIST_
  CnC::dist_cnc_init<flux_context> dinit;
#endif

#ifdef CIRCLE
  int rank = CIRCLE_init(argc, argv,CIRCLE_DEFAULT_FLAGS);
  CIRCLE_cb_process(&fluxFn);
  //default logging is too verbose
  CIRCLE_enable_logging(CIRCLE_LOG_ERR);
  if (rank!=0){
    int steps;
    MPI_Bcast(&steps, 1, MPI_INT, 0, MPI_COMM_WORLD);
    int prev_step;
    MPI_Bcast(&prev_step, 1, MPI_INT, 0, MPI_COMM_WORLD);
    for(int i=prev_step; i<steps; ++i){
      //2 *2 =4 half time steps
      for (int j=0;j<4;j++){ 
        CIRCLE_begin();
      }
    }
    CIRCLE_finalize();
    exit(0);
  }
#endif
// Display some info about this execution
// for the user.
  printf("**************************************************\n");
  printf("**************************************************\n");
  printf("**                                              **\n");
  printf("**                                              **\n");
#ifdef CNC
  printf("**   Running \"CNC 2D Kriging x processors       **\n");
#elif defined OMP
  printf("**   Running \"OpenMP 2D Kriging %d processors    **\n",
         omp_get_max_threads());
#elif CIRCLE
  int size;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  printf("**   Running \"libcircle 2D Kriging %i processors  **\n",
         size);
#else
#error Something is wrong
#endif
  printf("**                                              **\n");
  printf("**                                              **\n");
  printf("**************************************************\n");
  printf("**************************************************\n");

  //get input values from json file
  App CoMD;
  Input in;
  char input_file[1024];

  if(argc < 2) {
    fprintf(stderr,"Missing argument - json file\n");
    exit(1);
  }
  strcpy(input_file, argv[1]);
  parse_input((string)input_file, &in, &CoMD);

  if(argc >2){ 
    char host[1024];
    strcpy(host, argv[2]);
    in.head_node = string(host); 
    printf("Use single redis host on node: %s\n", in.head_node.c_str());
  }

  main_2DKriging(in, CoMD);
#ifdef CIRCLE
  CIRCLE_finalize();
#endif
}
