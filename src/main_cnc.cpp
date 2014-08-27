/** file containing the CNC main function
 * calls 2D Kriging main
 * **/
#include "main_cnc.hpp"
#include "input.hpp"
#include "types.h"
#include "2DKriging.hpp"

#include <cstring>
#include <stdio.h>
#ifdef OMP
#include <omp.h>
#endif//OMP

#ifdef CIRCLE
#include <libcircle.h>
#endif//CIRCLE

int main( int argc, char* argv[] )
{

#ifdef _DIST_
  CnC::dist_cnc_init<flux_context> dinit;
#endif
// Display some info about this execution
// for the user.
  printf("**************************************************\n");
  printf("**************************************************\n");
  printf("**                                              **\n");
  printf("**                                              **\n");
#ifdef CNC
  printf("**   Running \"CNC 2D Kriging x processors    **\n"
#elif defined OMP
  printf("**   Running \"OpenMP 2D Kriging %d processors    **\n",
         omp_get_max_threads());
#elif defined CIRCLE
  printf("**   Running \"libcircle 2D Kriging x processors  **\n"
         );
#else
#error Something is wrong
#endif
  printf("**                                              **\n");
  printf("**                                              **\n");
  printf("**************************************************\n");
  printf("**************************************************\n");

  //get input values from json file
  Input in;
  char input_file[1024];

  strcpy(input_file, argv[1]);
  parse_input((string)input_file, &in);

  main_2DKriging(in);

}
