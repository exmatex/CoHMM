#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <string.h>

#define GNUPLOT "/usr/bin/gnuplot -persist"

int main(int argc, char **argv) {
  
  FILE *gp1;
  gp1 = popen(GNUPLOT,"w");
  if (gp1==NULL) {
    printf("Error opening pipe to GNU plot < %s > t! \n", GNUPLOT);
    exit(0);
  }
  
  fprintf(gp1, "set terminal postscript\n");
  fprintf(gp1, "set output \'totalCoMDcalls_JMAK_res100000.ps\' \n");



  fprintf(gp1, "set xlabel 'Continuum time step'\n");
  fprintf(gp1, "set ylabel 'Number of finer-scale simulations'\n");
  fprintf(gp1, "set key bottom right\n");
  fprintf(gp1, "set logscale y\n");
  fprintf(gp1, "set xrange [-1000:101000]\n");
  fprintf(gp1, "plot 'l_totalCoMDcalls_res100000_timelength100000_threshold9.975e-09.dat' u 1:2 ti 'Double adaptive interpolation' w l lw 4 linecolor rgb \"green\", '' u 1:3 ti 'Coarsening step' w l lw 4 linecolor rgb \"red\", '' u 1:4 ti 'Refinement step' w l lw 4 linecolor rgb \"blue\", '' u 1:(($1)*400000) ti 'Brute force approach' w l lw 4 linecolor rgb \"black\"\n\n");

  fflush(gp1);


  FILE *gp2;
  gp2 = popen(GNUPLOT,"w");
  if (gp2==NULL) {
    printf("Error opening pipe to GNU plot < %s > t! \n", GNUPLOT);
    exit(0);
  }

  fprintf(gp2, "set terminal postscript\n");
  fprintf(gp2, "set output \'percentageOfCoMDcallsPerTimeStep_JMAK_res100000.ps\' \n");



  fprintf(gp2, "set xlabel 'Continuum time step'\n");
  fprintf(gp2, "set ylabel 'Percentage of finer-scale simulations'\n");
  fprintf(gp2, "set key top left\n");
  fprintf(gp2, "set xrange [-1000:101000]\n");
  fprintf(gp2, "plot 'l_percentageOfCoMDcallsPerTimeStep_res100000_timelength100000_threshold9.975e-09.dat' u 1:2 ti 'Double adaptive interpolation' w l lw 2 linecolor rgb \"green\", '' u 1:3 ti 'Coarsening step' w l lw 2 linecolor rgb \"red\", '' u 1:4 ti 'Refinement step' w l lw 2 linecolor rgb \"blue\"\n\n");

  fflush(gp2);

  pclose(gp1);
  pclose(gp2);
  //  pclose(gp3);

}
