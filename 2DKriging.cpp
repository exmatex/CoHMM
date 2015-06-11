/** 2D Macrosolver written by
 * Dominic Roehm
 with the help of
 * Robert Pavel
 * CoDesign summer school 2013
 */
#define _XOPEN_SOURCE 700
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <string.h>
#include <cfloat>
#include <math.h>

#include <string>

extern "C"
{
#include "types.h"
#include <CoMD_lib.h>
#include <omp.h>
#include <hiredis.h>
}


#include "2DKriging.hpp"
#include "kriging.hpp"
#include "redisBuckets.hpp"

#include <map>
#include <vector>
#include <list>
//#include <random>
#include <algorithm>
#include <iostream>
#include <ctime> 

#define GNUPLOT "/usr/bin/gnuplot -persist"

/******************************************/
//define specifies if CoMD is used or the linear "analytic" approach
#define CoMD
//#define C_RAND
#define KRIGING
#define DB
#define KR_DB
//#define GIF
//#define OUTPUT
//#define XWAVE
#define CIRCULAR
/*****************************************/


int comdDigits = 4;
int krigDigits = 8;
double gradThresh = 0.001;
double dbT = -1;
double mom_gamma = 0.0;
double strain_gamma = mom_gamma;
double en_gamma = 0.1*strain_gamma;
//double dbT = 0.0;
double errorThresh = -1;
//double errorThresh = 0.01;
double zeroThresh = 0.0001;
char  potDir[] = "../pots";
//char * potDir = "../pots";
redisContext *headRedis;
char headNode[1024];
int counter = 0;

/** precise timing measurement */
double getUnixTime(void)
{
    struct timespec ts;

    if(clock_gettime(CLOCK_REALTIME, &ts) != 0) return 0;

    return (((double) ts.tv_sec) + (double) (ts.tv_nsec / 1000000000.0));
}
/** initialize the struct containing all node quantities
 * @param node       grid_node containing the conserved and the fluxes
 * @param grid_size  number of grid nodes
 * **/
void init_nodes(Node* node, int grid_size){

  for(int i=0; i<grid_size; ++i){
 
    for(int j=0; j<7; ++j){
      node[i].f.f[j] = 0.0;
      node[i].g.f[j] = 0.0;
      node[i].w.w[j] = 0.0;
    }
    node[i].boundary = 0;

    //printf("init node struct test %g\n", node[i].w.momentum[1]);
 }

}

/** min_mod function to compute the discrete slopes, without using jacobian etc.
 * Jiang&Tadmor equ. (3.1'), (3.1`) and def. one on page 1901
 * @param w_plus   conserved or flux of node+1 (inpout) 
 * @param w        conserved or flux of node (inpout) 
 * @param w_minus  conserved or flux of node-1 (inpout) 
 * @param mm       result (spacial derivative of input) (output) 
 * **/
double min_mod(double w_plus, double w, double w_minus){

  double mm = 0.0;
  double w_plus_w = 0.0;
  double w_w = 0.0;
  double w_minus_w = 0.0;
  double theta = 1.0;

  w_plus_w   = theta*(w_plus - w);
  w_w = 0.5*(w_plus - w_minus);
  w_minus_w  = theta*(w - w_minus);

  if(w_plus_w > 0 && w_w > 0 && w_minus_w > 0){
      mm = fmin(w_plus_w, w_w);
      mm = fmin(mm, w_minus_w);
  }else if(w_plus_w < 0 && w_w < 0 && w_minus_w < 0){
      mm = fmax(w_plus_w, w_w);
      mm = fmax(mm, w_minus_w);
  }else{
      mm = 0.0;
  }

  return mm;
}
/** set intial values for strain, momentum- and energy density on grid nodes
 * @param node       grid_node containing the conserved and the fluxes (output)
 * @param grid_size  number of grid nodes (input)
 * @param l          lattice parameters (input)         
 * @param init       initial conserved values (inpout) 
 * **/
void init_conserved_fields(Node* node_a, Lattice l, int grid_size, Values init){

  int x,y;
  for (int i=0; i<grid_size; ++i) {
    index_to_xy(i, l, &x, &y);
      //A0[i] = A[i] = (i < dimX/2) ? 1.01 : 1.0; // small initial step in deformation gradient
    node_a[i].w.w[0] = 1.0;
    node_a[i].w.w[1] = 0.0;
    node_a[i].w.w[2] = 0.0;
    node_a[i].w.w[3] = 1.0;
    node_a[i].w.w[4] = 0.0;
    node_a[i].w.w[5] = 0.0;
    //node_a[index].w.w[6] = ((i < l.dim_x/2+l.dim_x/10)&&(i >= l.dim_x/2-l.dim_x/10)) ? init.energy+0.05 : init.energy;
    node_a[i].w.w[6] = init.energy;
    node_a[i].w.rho = init.rho;
#if 0
    if((x<l.dim_x/2+l.dim_x/10) && (x>=l.dim_x/2-l.dim_x/10) && (y<l.dim_y/2+l.dim_y/10) && (y>=l.dim_y/2-l.dim_y/10)){
    //if(y==l.dim_y/2){
       node_a[i].w.w[0] = 1.1;
       node_a[i].w.w[3] = 1.1;
    }
#endif
#ifdef XWAVE
//#########################################
    //xwave
    if((x<l.dim_x/2+l.dim_x/10) && (x>=l.dim_x/2-l.dim_x/10)){
    //if(y==l.dim_y/2){
       node_a[i].w.w[0] = 1.04;
       //node_a[i].w.w[3] = 1.1;
    }
//#########################################
#endif//XWAVE
    //if(x==y){
    //if(x==l.dim_x/2 && y==l.dim_y/2){
#if 0
    //set w[2] = -w[1] for symmetric behaviour with xx, yy
    if((x<l.dim_x/2+l.dim_x/10) && (x>=l.dim_x/2-l.dim_x/10) && (y<l.dim_y/2+l.dim_y/10) && (y>=l.dim_y/2-l.dim_y/10)){
    //if(x==y){
       node_a[i].w.w[1] = 0.1;
       node_a[i].w.w[2] = -0.1;
    }
#endif
#ifdef CIRCULAR
//#########################################
    //ciruclar wave
    //if((x<l.dim_x/2+l.dim_x/10) && (x>=l.dim_x/2-l.dim_x/10) && (y<l.dim_y/2+l.dim_y/10) && (y>=l.dim_y/2-l.dim_y/10)){
      int r = floor(sqrt((x-l.dim_x/2)*(x-l.dim_x/2) + (y-l.dim_y/2)*(y-l.dim_y/2)));
      if(r <= 3){
        node_a[i].w.w[0] = 1.02;
        node_a[i].w.w[3] = 1.02;
        node_a[i].w.w[1] = 0.02;
        node_a[i].w.w[2] = -0.02;
      }
    //}
//#########################################
#endif//CIRCULAR    

    //if(y==l.dim_y/2) node_a[i].w.w[3] = 1.02;
  }
  //node_a[4].w.w[0] = 1.02;
    //int half = floor(dimX/2);
}
/** print the calls
 * @param i         number of output file (input)         
 * @param node      grid_node containing the conserved and the fluxes (input)
 * @param l         lattice parameters (input)         
 * @param grid_size  number of grid nodes (input)
 * **/
void printf_calls(int counter, Calls ca){

#if 0
    //correct for double counting
    if(ca.kPoints!=0){
        ca.kPoints -= ca.krig;
    }
    if(ca.cPoints!=0){
        ca.cPoints -= ca.comd;
    }
#endif
    int bufferSize = 256;
    char file_name[bufferSize];
    sprintf(file_name,"calls.dat");
    FILE *fn = fopen(file_name, "a");
    if (fn == NULL) {
        printf("Error writing file< %s >!\n", file_name);
        exit(0);
    }

    fprintf(fn, "%i     %.0f       %.0f      %i      %i      %.0f      %.0f      %i\n", counter, double(ca.comd+(ca.kFail*0.1)), double(ca.cPoints+(ca.kFail*0.9)), ca.db, ca.kdb, double(ca.krig-(ca.kFail*0.1)), double(ca.kPoints-(ca.kFail*0.9)), ca.kFail);
    fclose(fn);
}

/** print the calls
 * @param i         number of output file (input)         
 * @param node      grid_node containing the conserved and the fluxes (input)
 * @param l         lattice parameters (input)         
 * @param grid_size  number of grid nodes (input)
 * **/
void printf_timings(int counter, Tms tm){

  int bufferSize = 256;
  char file_name[bufferSize];
  sprintf(file_name,"timings.dat");
  FILE *fn = fopen(file_name, "a");
  if (fn == NULL) {
    printf("Error writing file< %s >!\n", file_name);
    exit(0);
  }

  fprintf(fn, "%i     %.8f       %.8f      %.8f      %.8f      %.8lf      %.8f     %.8f\n", counter, tm.co, tm.co2, tm.db, tm.krDb, tm.kr, tm.kr2, tm.ca);
  fclose(fn);
}

/** print the colormap output in vtk format
 * @param i         number of output file (input)         
 * @param node      grid_node containing the conserved and the fluxes (input)
 * @param l         lattice parameters (input)         
 * @param grid_size  number of grid nodes (input)
 * **/
void printf_colormap_vtk(int i, Node* fields, Lattice l, int grid_size){

  int bufferSize = 256;
  char file_name[bufferSize];
  sprintf(file_name,"colormap_%i.vtk",i);
  FILE *fn = fopen(file_name, "w");
  if (fn == NULL) {
    printf("Error writing file< %s >!\n", file_name);
    exit(0);
  }

  fprintf(fn, "# vtk DataFile Version 2.0\nstrain\nASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS %i %i 1\nORIGIN 0.0 0.0 0.0\nSPACING %lf %lf 1.0\nPOINT_DATA %i\nSCALARS OutArray floats 1\nLOOKUP_TABLE default\n", l.dim_x, l.dim_y, l.dx, l.dy, grid_size); 
  for(int index=0; index<grid_size; ++index){
    fprintf(fn, "%i\n", fields[index].f.ca);
  }
  fclose(fn);
}
/** print the output in vtk format
 * @param i         number of output file (input)         
 * @param node      grid_node containing the conserved and the fluxes (input)
 * @param l         lattice parameters (input)         
 * @param grid_size  number of grid nodes (input)
 * **/
void printf_fields_vtk(int i, Node* node_a, Lattice l, int grid_size){

  int bufferSize = 256;
  char file_name[bufferSize];
  sprintf(file_name,"strain_%i.vtk",i);
  FILE *fn = fopen(file_name, "w");
  if (fn == NULL) {
    printf("Error writing file< %s >!\n", file_name);
    exit(0);
  }
  char file_name2[bufferSize];
  sprintf(file_name2,"mom_%i.vtk",i);
  FILE *fn2 = fopen(file_name2, "w");
  if (fn2 == NULL) {
    printf("Error writing file< %s >!\n", file_name2);
    exit(0);
  }
  char file_name3[bufferSize];
  sprintf(file_name3,"energy_%i.vtk",i);
  FILE *fn3 = fopen(file_name3, "w");
  if (fn3 == NULL) {
    printf("Error writing file< %s >!\n", file_name3);
    exit(0);
  }
  char file_name4[bufferSize];
  sprintf(file_name4,"mom_flux_%i.vtk",i);
  FILE *fn4 = fopen(file_name4, "w");
  if (fn4 == NULL) {
    printf("Error writing file< %s >!\n", file_name4);
    exit(0);
  }
  char file_name5[bufferSize];
  sprintf(file_name5,"flux_%i.vtk",i);
  FILE *fn5 = fopen(file_name5, "w");
  if (fn5 == NULL) {
    printf("Error writing file< %s >!\n", file_name5);
    exit(0);
  }
  char file_name6[bufferSize];
  sprintf(file_name6,"energy_flux_%i.vtk",i);
  FILE *fn6 = fopen(file_name6, "w");
  if (fn6 == NULL) {
    printf("Error writing file< %s >!\n", file_name6);
    exit(0);
  }

  fprintf(fn, "# vtk DataFile Version 2.0\nstrain\nASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS %i %i 1\nORIGIN 0.0 0.0 0.0\nSPACING %lf %lf 1.0\nPOINT_DATA %i\nSCALARS OutArray floats 4\nLOOKUP_TABLE default\n", l.dim_x, l.dim_y, l.dx, l.dy, grid_size); 
  fprintf(fn2, "# vtk DataFile Version 2.0\nmom\nASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS %i %i 1\nORIGIN 0.0 0.0 0.0\nSPACING %lf %lf 1.0\nPOINT_DATA %i\nSCALARS OutArray floats 2\nLOOKUP_TABLE default\n", l.dim_x, l.dim_y, l.dx, l.dy, grid_size); 
  fprintf(fn3, "# vtk DataFile Version 2.0\nenergy\nASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS %i %i 1\nORIGIN 0.0 0.0 0.0\nSPACING %lf %lf 1.0\nPOINT_DATA %i\nSCALARS OutArray floats 1\nLOOKUP_TABLE default\n", l.dim_x, l.dim_y, l.dx, l.dy, grid_size); 
  fprintf(fn4, "# vtk DataFile Version 2.0\nmom_flux\nASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS %i %i 1\nORIGIN 0.5 0.5 0.5\nSPACING %lf %lf 1.0\nPOINT_DATA %i\nSCALARS OutArray floats 2\nLOOKUP_TABLE default\n", l.dim_x, l.dim_y, l.dx, l.dy, grid_size); 
  fprintf(fn5, "# vtk DataFile Version 2.0\nflux\nASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS %i %i 1\nORIGIN 0.5 0.5 0.5\nSPACING %lf %lf 1.0\nPOINT_DATA %i\nSCALARS OutArray floats 4\nLOOKUP_TABLE default\n", l.dim_x, l.dim_y, l.dx, l.dy, grid_size); 
  fprintf(fn6, "# vtk DataFile Version 2.0\nenergy_flux\nASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS %i %i 1\nORIGIN 0.5 0.5 0.5\nSPACING %lf %lf 1.0\nPOINT_DATA %i\nSCALARS OutArray floats 2\nLOOKUP_TABLE default\n", l.dim_x, l.dim_y, l.dx, l.dy, grid_size); 
  for(int index=0; index<grid_size; ++index){
    fprintf(fn, "%f %f %f %f\n", node_a[index].w.w[0]-1, node_a[index].w.w[1], node_a[index].w.w[2], node_a[index].w.w[3]-1);
    fprintf(fn2, "%f %f \n", node_a[index].w.w[4], node_a[index].w.w[5]);
    fprintf(fn3, "%f \n", node_a[index].w.w[6]);
    fprintf(fn4, "%f %f \n", node_a[index].f.f[0], node_a[index].g.f[3]);
    fprintf(fn5, "%f %f %f %f\n", node_a[index].f.f[4], node_a[index].f.f[5], node_a[index].g.f[4], node_a[index].g.f[5]);
    fprintf(fn6, "%f %f \n", node_a[index].f.f[6], node_a[index].g.f[6]);
  }

  fclose(fn);
  fclose(fn2);
  fclose(fn3);
  fclose(fn4);
  fclose(fn5);
  fclose(fn6);
}
/** plot the strain, stress and energy dens to a gnuplot .gif file
 * @param gp         filepointer for gnuplot output (input)         
 * @param node       grid_node containing the conserved and the fluxes (output)
 * @param l          lattice parameters (input)         
 * **/
void plot_fields(int i, Node* node_a, Lattice l, FILE * gp){

//########################################//
#ifndef GIF
  int bufferSize = 256;
  char file_name[bufferSize];
  sprintf(file_name,"gp_%i.dat",i);
  FILE *gp_out = fopen(file_name, "w");
  if (gp_out == NULL) {
    printf("Error writing file< %s >!\n", file_name);
    exit(0);
  }
  gp = popen(GNUPLOT,"w");
  if (gp==NULL) {
    printf("Error opening pipe to GNU plot < %s > t! \n", GNUPLOT);
    exit(0);
  }
  fprintf(gp, "set terminal postscript portrait size 8.5cm,6.3cm enhanced color dashed lw 1 \"Helvetica,8\"\n");
  fprintf(gp, "set output \'cohmm_%04d.ps\' \n", i);
#endif//GIF
//########################################//
  fprintf(gp, "set multiplot \n");
  fprintf(gp, "set label 1 'y [0..%i]' at -15.0,5.0,0.0 rotate by -35.5\n", l.dim_y);
  fprintf(gp, "set label 2 'x [0..%i]' at 60.0, 22.5,0.0 rotate by 4.5\n", l.dim_x);
  fprintf(gp, "set label 3 'Strain [MPa]' at -110.0,10.0,54\n");
  fprintf(gp, "set label 4 'Calltype' at -110.0,10.0,64\n");
  fprintf(gp, "set xrange [ :%i ] noreverse nowriteback\n", (l.dim_x-1)*2);
  fprintf(gp, "set yrange [ :%i ] noreverse nowriteback\n", (l.dim_y-1)*2);
  fprintf(gp, "set zrange [ 0.0 : 40 ] noreverse nowriteback\n");
  //fprintf(gp, "set cbrange [ 0.0 : 0.02 ] noreverse nowriteback\n");
  //fprintf(gp, "set xtics offset 0,-0.5,0\n");
  //fprintf(gp, "set ytics 0,2,8\n");
  fprintf(gp, "unset xtics\n");
  fprintf(gp, "unset ytics\n");
  fprintf(gp, "set border 4095 front linetype -1 linewidth 1.000\n");
  fprintf(gp, "set ticslevel 0\n");
  fprintf(gp, "set style line 100  linetype 5 linecolor rgb \"#f0e442\"  linewidth 0.500 pointtype 5 pointsize default pointinterval 0\n");
  fprintf(gp, "set view 110, 20, 1, 1\n");
  fprintf(gp, "set samples %i,%i\n", (2*l.dim_x), (2*l.dim_y));
  fprintf(gp, "set isosample %i,%i\n", (2*l.dim_x), (2*l.dim_y));
  fprintf(gp, "set pm3d at s\n");
  fprintf(gp, "set pm3d interpolate 1,1 flush begin noftriangles hidden3d 100 corners2color mean\n");
   //fprintf(gp, "set pm3d \n");
  fprintf(gp, "set colorbox vertical user origin 0.86, .245 size .04,.45\n");
  int index = 0;
  //fprintf(gp, "set xrange [0:%i]\n", l.dim_x);
  //fprintf(gp, "set dgrid3d %i,%i, 1\n", l.dim_y, l.dim_x);
  //fprintf(gp, "splot '-' u 1:2:3 ti 'strain_xx' w l, '-' u 1:2:4 ti 'strain_yy' w l, '-' u ($1+0.5):($2+0.5):5 ti 'stress_xx' w l\n\n");
  fprintf(gp, "splot '-' u 1:2:((sqrt($3**2+$4**2+$5**2+$6**2))*1000) ti '' w l palette\n\n");
  //fprintf(gp, "splot '-' u 1:2:3 ti 'strain' w l, '-' u ($1+0.5):($2+0.5):5 ti 'stress' w l\n\n");
  //fprintf(gp, "splot '-' u 1:2:(sqrt($3**2+$6**2)) ti 'strain' w l\n\n");
  // gnuplot forces us to repeat ourselves
  int x = 0;
  int y = 0;
  for (int k=0; k<1; ++k) {
    for (int i=0; i<(2*l.dim_x); i++) {
        y=0;
      for (int j=0; j<(2*l.dim_y); j+=2) {
        xy_to_index(x, y, &index, l);
        y++;
        //gnuplot gif
        fprintf(gp, "%f %f %f %f %f %f\n", i*l.dx, j*l.dy, node_a[index].w.w[0]-1, node_a[index].w.w[1], node_a[index].w.w[2], node_a[index].w.w[3]-1);
        fprintf(gp, "%f %f %f %f %f %f\n", i*l.dx, (j+1)*l.dy, node_a[index].w.w[0]-1, node_a[index].w.w[1], node_a[index].w.w[2], node_a[index].w.w[3]-1);
        //gnuplot dat file
        fprintf(gp_out, "%f %f %f %f %f %f %i\n", i*l.dx, j*l.dy, node_a[index].w.w[0]-1, node_a[index].w.w[1], node_a[index].w.w[2], node_a[index].w.w[3]-1, node_a[index].f.ca);
        fprintf(gp_out, "%f %f %f %f %f %f %i\n", i*l.dx, (j+1)*l.dy, node_a[index].w.w[0]-1, node_a[index].w.w[1], node_a[index].w.w[2], node_a[index].w.w[3]-1, node_a[index].f.ca);
      }
      if(i%2 == 1) x++;
      fprintf(gp, "\n");
      fprintf(gp_out, "\n");
    }
    fprintf(gp, "e\n");
    fprintf(gp_out, "e\n");
  }
  fflush(gp);
  fprintf(gp, "reset\n");
  fprintf(gp, "unset surface\n");
  fprintf(gp, "set border 4095 front linetype -1 linewidth 1.000\n");
  fprintf(gp, "set style line 100 linetype 5 linecolor rgb \"#f0e442\" linewidth 0.500 pointtype 5 pointsize default pointinterval 0\n");
  fprintf(gp, "set view 110, 20, 1, 1\n");
  fprintf(gp, "set samples %i,%i\n", (2*l.dim_x), (2*l.dim_y));
  fprintf(gp, "set isosample %i,%i\n", (2*l.dim_x), (2*l.dim_y));
  fprintf(gp, "set ticslevel 0\n");
  fprintf(gp, "set palette rgbformulae 33,13,10\n");
  fprintf(gp, "set pal maxcolors 6\n");
  fprintf(gp, "set pm3d implicit at t\n");
  fprintf(gp, "set xrange [ :%i ] noreverse nowriteback\n", (l.dim_x-1)*2);
  fprintf(gp, "set yrange [ :%i ] noreverse nowriteback\n", (l.dim_y-1)*2);
  fprintf(gp, "set cbrange [ 1 : 6 ] noreverse nowriteback\n");
  fprintf(gp, "set palette defined (1 \"blue\", 2 \"light-blue\", 3 \"turquoise\", 4 \"dark-green\", 5 \"red\", 6 \"orange\")\n");
  fprintf(gp, "set cbtics (\"CoMD\" 1, \"C. Dupl.\" 2, \"DB\" 3, \"Kr. DB\" 4, \"Kr.\" 5, \"Kr. Dupl.\" 6)\n");
  fprintf(gp, "unset xtics\n");
  fprintf(gp, "unset ytics\n");
  fprintf(gp, "set colorbox horizontal user origin .155, 0.925 size .7,.04\n");
  fprintf(gp, "splot '-' u 1:2:3 ti '' \n");
  x = 0;
  for (int k=0; k<1; ++k) {
    for (int i=0; i<(2*l.dim_x); ++i) {
        y=0;
      for (int j=0; j<(2*l.dim_y); j+=2) {
        xy_to_index(x, y, &index, l);
        y++;
        fprintf(gp, "%f %f %i \n", i*l.dx, j*l.dy, node_a[index].f.ca);
        fprintf(gp, "%f %f %i \n", i*l.dx, (j+1)*l.dy, node_a[index].f.ca);
      }
      if(i%2 == 1) x++;
      fprintf(gp, "\n");
    }
    fprintf(gp, "e\n");
  }
  fflush(gp);
#ifndef GIF
  fclose(gp);
  fclose(gp_out);
#endif//GIF
}

/** set boundaries
 *
 */
void set_boundaries(Lattice l, int grid_size, Node* node_a, Node* node_b){
  int x,y;
  for(int i=0; i<grid_size; ++i){
    index_to_xy(i, l, &x, &y);
    if(x==0 || x==(l.dim_x-1) || y==0 || y==(l.dim_y-1)){ 
        node_a[i].boundary = 1;
        node_b[i].boundary = 1;
    }
  }

}
/** zero flux boundaries
 *
 */
void boundaries(Lattice l, int grid_size, Node* node){
  for(int i=0; i<grid_size; ++i){
    if(node[i].boundary == 1){ 
      for(int j=0; j<7; ++j){
        node[i].f.f[j] = 0.0;
        node[i].g.f[j] = 0.0;
        node[i].w.w[0] = 1.0;
        node[i].w.w[1] = 0.0;
        node[i].w.w[2] = 0.0;
        node[i].w.w[3] = 1.0;
        node[i].w.w[4] = 0.0;
        node[i].w.w[5] = 0.0;
      }
    }
  }

}
/** shift back indices, shifts back the inidices with the help of double buffering (raceconditions!) 
 * @param node_b      grid_node containing the conserved and the fluxes (input)
 * @param grid_size   number of grid nodes (input)
 * @param l           lattice parameters (input)         
 * @param node_a      grid_node containing the conserved and the fluxes (output)
 * **/
void shift_back(Node* node_b, int grid_size, Lattice l, Node* node_a){

  int x = 0;
  int y = 0;
  int xpoypo;
  for(int i=0; i<grid_size; ++i){
    index_to_xy(i, l, &x, &y);
    // temp vars for indexing
    xpoypo = (x+1)%l.dim_x + l.dim_x*((y+1)%l.dim_y);
    node_b[i].f = node_a[xpoypo].f;
    node_b[i].g = node_a[xpoypo].g;
    node_b[i].w = node_a[xpoypo].w;
  }
  for(int i=0; i<grid_size; ++i){
    //index_to_xy(i, l, &x, &y);
    // temp vars for indexing
    node_a[i].f = node_b[i].f;
    node_a[i].g = node_b[i].g;
    node_a[i].w = node_b[i].w;
  }
}

/** flux computation executed in parallel 
 * @param *in         flux input contains information about field on specific node (input)
 * @return *out       returns flux output information about computed fluxes on specific node (output)
 * @param *dbCache    database cache pointer (input)
 * **/
void fluxFn(fluxInput *in, fluxOutput *out, std::map<std::string, std::vector<char *> > *dbCache, double* startKr, double* stopKr, double* startCo, double* stopCo)
{
	//Prep outputs
	out->error = 0.0;

	//Make redis context
	redisContext * rTask;
	rTask = redisConnect(in->headNode, 6379);
	if(rTask == NULL || rTask->err)
	{
		//printf("Redis error: %s\n", rTask->errstr);
		//printf("Redis error: %s\n", in->headNode);
		printf("Redis connection delay\n");
    //sleep(20);
	  rTask = redisConnect(in->headNode, 6379);
		//return NULL;
	}
    //if(diff1 > 1.0) printf("CONNECTION timing %lf\n", diff1);
	if(rTask == NULL || rTask->err)
	{
		printf("Redis error: %s\n", rTask->errstr);
		printf("Redis error: %s\n", in->headNode);
    exit(0);
		//return NULL;
	}
	//Is this CoMD?
	if(in->callCoMD == true)
	{
    *startCo = getUnixTime();
		//Call CoMD
		//4 strains
		double strain_xx = in->w->w[0];
		double strain_xy = in->w->w[1];
		double strain_yx = in->w->w[2];
		double strain_yy = in->w->w[3];
		//2 momentum_fluxes
		double momentum_x = in->w->w[4];
		double momentum_y = in->w->w[5];
		//enery_flux
		double energy = in->w->w[6];

#ifdef CoMD
		CoMD_input theInput;

		///FIXME: Fix potdir
		strcpy(theInput.potDir,potDir);
		strcpy(theInput.potName,"Cu01.eam.alloy");
		strcpy(theInput.potType,"setfl");
		theInput.doeam = 1; 
		theInput.nx = 6;
		theInput.ny = 6;
		theInput.nz = 6;
		theInput.nSteps = 1000;
		theInput.printRate = 1; 
		//MUST SPECIFY THE FOLLOWING               
		theInput.dt = 10.0;
		theInput.lat = 3.6186;
		theInput.temperature = 0;
		theInput.initialDelta = 0.0;
		theInput.defGrad[0] = strain_xx;
		theInput.defGrad[1] = strain_xy;
		theInput.defGrad[2] = strain_yx;
		theInput.defGrad[3] = strain_yy;
		theInput.enDens = energy;
		theInput.momDens[0] = momentum_x;
		theInput.momDens[1] = momentum_y;
		theInput.momDens[2] = 0.0;
		//theInput.rank = rank;
		//theInput.calls = calls;
		theInput.rank = 0;
		theInput.calls = 0;

		//Call it
		CoMD_return theRet = CoMD_lib(&theInput);

		//4 stresses
		double rho = 1.0;
		out->f[0] = momentum_x/rho;
		out->f[1] = momentum_y/rho;
		out->f[2] = 0.0;
		out->f[3] = 0.0;
		out->f[4] = theRet.stressXX;
		out->f[5] = theRet.stressXY;
		out->g[0] = 0.0;
		out->g[1] = 0.0;
		out->g[2] = momentum_x/rho;
		out->g[3] = momentum_y/rho;
		out->g[4] = theRet.stressYX;
		out->g[5] = theRet.stressYY;
		out->f[6] = -theRet.energyDensX;
		out->g[6] = -theRet.energyDensY;
#elif defined(C_RAND)
    //gaussian noise, gamma = strength
    unsigned int tid = omp_get_thread_num();
    unsigned s_seed = 1103515245 * tid + floor(1000*(*startCo));
    std::srand(s_seed);
    int seed = rand()%10000; 
    //fake comd values
	  //4 stresses
	  double rho = 1.0;
	  //out.f[0] = theRet.momX;
    double *num;
    num = new double[2];
    gaussian_random(&seed, num);
    //printf("num %lf %lf\n", num[0], num[1]);
	  out->f[0] = momentum_x/rho + (mom_gamma*num[0]);
	  out->f[1] = momentum_y/rho + (mom_gamma*num[1]);
	  out->f[2] = 0.0;
	  //o->t.f[2] = theRet.momY;
	  out->f[3] = 0.0;
    gaussian_random(&seed, num);
	  out->f[4] = strain_xx-1 + 0.75*(strain_yy-1) + (strain_gamma*num[0]);
	  out->f[5] = 1.9*strain_xy + (strain_gamma*num[1]);
	  out->g[0] = 0.0;
	  //o->t.g[1] = theRet.momX;
	  out->g[1] = 0.0;
    gaussian_random(&seed, num);
	  out->g[2] = momentum_x/rho + (mom_gamma*num[0]);
	  //o->t.g[3] = theRet.momY;
	  out->g[3] = momentum_y/rho + (mom_gamma*num[1]);
    gaussian_random(&seed, num);
	  out->g[4] = 1.9*strain_yx + (strain_gamma*num[0]);
	  out->g[5] = strain_yy-1 + 0.75*(strain_xx-1) + (strain_gamma*num[1]);
	  //f->6] = g[6] = -0.295;
	  out->f[6] = -out->f[0]*sqrt(out->f[4]*out->f[4] + out->f[5]*out->f[5]);
	  out->g[6] = -out->g[3]*sqrt(out->g[4]*out->g[4] + out->g[5]*out->g[5]);
    delete[] num;
#else
    //fake comd values
		//4 stresses
		double rho = 1.0;
		//out.f[0] = theRet.momX;
		out->f[0] = momentum_x/rho;
		out->f[1] = momentum_y/rho;
		out->f[2] = 0.0;
		//o->t.f[2] = theRet.momY;
		out->f[3] = 0.0;
		out->f[4] = strain_xx-1 + 0.75*(strain_yy-1);
		out->f[5] = 1.9*strain_xy;
		out->g[0] = 0.0;
		//o->t.g[1] = theRet.momX;
		out->g[1] = 0.0;
		out->g[2] = momentum_x/rho;
		//o->t.g[3] = theRet.momY;
		out->g[3] = momentum_y/rho;
		out->g[4] = 1.9*strain_yx;
		out->g[5] = strain_yy-1 + 0.75*(strain_xx-1);
		//f->6] = g[6] = -0.295;
		out->f[6] = -out->f[0]*sqrt(out->f[4]*out->f[4] + out->f[5]*out->f[5]);
		out->g[6] = -out->g[3]*sqrt(out->g[4]*out->g[4] + out->g[5]*out->g[5]);
#endif//COMD
		//Write result to database
		putData(in->w->w, out->f, out->g, (char *)"comd", rTask, comdDigits);		
    *stopCo = getUnixTime();
    //if(diff1 > 1.0) printf("put timing %lf\n", diff1);
	}
	//It is kriging time!
	else
	{
    *startKr = getUnixTime();
    //printf("krigingtime\n");
		//Get data for kriging
		std::vector<double *> oldWs;
		std::vector<double *> oldFs;
		std::vector<double *> oldGs;
		getCachedSortedSubBucketNearZero(in->w->w, (char *)"comd", rTask, comdDigits, 10, &oldWs, &oldFs, &oldGs, zeroThresh, dbCache); 
		//getSortedSubBucketNearZero(in->w->w, "comd", rTask, comdDigits, 10, &oldWs, &oldFs, &oldGs, zeroThresh); 

		//Call Kriging on each point
    double *resF = new double[2];
    double *resG = new double[2];
    int info;
		for(int i = 0; i < 7; i++)
		{
			///TODO: Consider refactoring this in to kriging.cpp
			std::vector<double> oldF;
			std::vector<double> oldG;
			for(int j = 0; j < oldFs.size(); j++)
			{
				oldF.push_back(oldFs[j][i]);
				oldG.push_back(oldGs[j][i]);
			}
			info = kriging(in->w->w, 7, oldWs, oldF, resF);	
			info = kriging(in->w->w, 7, oldWs, oldG, resG);	
			//Write result
			out->f[i] = resF[0];
			out->g[i] = resG[0];
			//Set error
			if(resF[1] > out->error)
			{
				out->error = resF[1];
			}
			if(resG[1] > out->error)
			{
				out->error = resG[1];
			}
		}
    freeClear(oldWs);
    freeClear(oldFs);
    freeClear(oldGs);
    delete[] resF;
    delete[] resG;
    //check for error and call CoMd if necessary 
    //printf("error %lf\n", out->error);
    if(out->error > errorThresh)
    {
        *stopKr = getUnixTime();
	      out->error = 0.0;
        *startCo = getUnixTime();
		    //Call CoMD
		    //4 strains
		    double strain_xx = in->w->w[0];
		    double strain_xy = in->w->w[1];
		    double strain_yx = in->w->w[2];
		    double strain_yy = in->w->w[3];
		    //2 momentum_fluxes
		    double momentum_x = in->w->w[4];
		    double momentum_y = in->w->w[5];
		    //enery_flux
		    double energy = in->w->w[6];

#ifdef CoMD
		    CoMD_input theInput;

		    ///FIXME: Fix potdir
		    strcpy(theInput.potDir,potDir);
		    strcpy(theInput.potName,"Cu01.eam.alloy");
		    strcpy(theInput.potType,"setfl");
		    theInput.doeam = 1; 
		    theInput.nx = 6;
		    theInput.ny = 6;
		    theInput.nz = 6;
		    theInput.nSteps = 1000;
		    theInput.printRate = 1; 
		    //MUST SPECIFY THE FOLLOWING               
		    theInput.dt = 10.0;
		    theInput.lat = 3.6186;
		    theInput.temperature = 0;
		    theInput.initialDelta = 0.0;
		    theInput.defGrad[0] = strain_xx;
		    theInput.defGrad[1] = strain_xy;
		    theInput.defGrad[2] = strain_yx;
		    theInput.defGrad[3] = strain_yy;
		    theInput.enDens = energy;
		    theInput.momDens[0] = momentum_x;
		    theInput.momDens[1] = momentum_y;
		    theInput.momDens[2] = 0.0;
		    //theInput.rank = rank;
		    //theInput.calls = calls;
		    theInput.rank = 0;
		    theInput.calls = 0;

		    //Call it
		    CoMD_return theRet = CoMD_lib(&theInput);

		    //4 stresses
		    double rho = 1.0;
		    out->f[0] = momentum_x/rho;
		    out->f[1] = momentum_y/rho;
		    out->f[2] = 0.0;
		    out->f[3] = 0.0;
		    out->f[4] = theRet.stressXX;
		    out->f[5] = theRet.stressXY;
		    out->g[0] = 0.0;
		    out->g[1] = 0.0;
		    out->g[2] = momentum_x/rho;
		    out->g[3] = momentum_y/rho;
		    out->g[4] = theRet.stressYX;
		    out->g[5] = theRet.stressYY;
		    out->f[6] = -theRet.energyDensX;
		    out->g[6] = -theRet.energyDensY;
        //delete &theRet;
#elif defined(C_RAND)
        //gaussian noise, gamma = strength
        unsigned int tid = omp_get_thread_num();
        unsigned s_seed = 1103515245 * tid + floor(1000*(*startCo));
        std::srand(s_seed);
        int seed = rand()%10000; 
        //fake comd values
		    //4 stresses
		    double rho = 1.0;
		    //out.f[0] = theRet.momX;
        double *num;
        num = new double[2];
        gaussian_random(&seed, num);
        //printf("num %lf %lf\n", num[0], num[1]);
		    out->f[0] = momentum_x/rho + (mom_gamma*num[0]);
		    out->f[1] = momentum_y/rho + (mom_gamma*num[1]);
		    out->f[2] = 0.0;
		    //o->t.f[2] = theRet.momY;
		    out->f[3] = 0.0;
        gaussian_random(&seed, num);
		    out->f[4] = strain_xx-1 + 0.75*(strain_yy-1) + (strain_gamma*num[0]);
		    out->f[5] = 1.9*strain_xy + (strain_gamma*num[1]);
		    out->g[0] = 0.0;
		    //o->t.g[1] = theRet.momX;
		    out->g[1] = 0.0;
        gaussian_random(&seed, num);
		    out->g[2] = momentum_x/rho + (mom_gamma*num[0]);
		    //o->t.g[3] = theRet.momY;
		    out->g[3] = momentum_y/rho + (mom_gamma*num[1]);
        gaussian_random(&seed, num);
		    out->g[4] = 1.9*strain_yx + (strain_gamma*num[0]);
		    out->g[5] = strain_yy-1 + 0.75*(strain_xx-1) + (strain_gamma*num[1]);
		    //f->6] = g[6] = -0.295;
		    out->f[6] = -out->f[0]*sqrt(out->f[4]*out->f[4] + out->f[5]*out->f[5]);
		    out->g[6] = -out->g[3]*sqrt(out->g[4]*out->g[4] + out->g[5]*out->g[5]);
        delete[] num;
#else
        //fake comd values
    		//4 stresses
    		double rho = 1.0;
    		//out.f[0] = theRet.momX;
    		out->f[0] = momentum_x/rho;
    		out->f[1] = momentum_y/rho;
    		out->f[2] = 0.0;
    		//o->t.f[2] = theRet.momY;
    		out->f[3] = 0.0;
    		out->f[4] = strain_xx-1 + 0.75*(strain_yy-1);
    		out->f[5] = 1.9*strain_xy;
    		out->g[0] = 0.0;
    		//o->t.g[1] = theRet.momX;
    		out->g[1] = 0.0;
    		out->g[2] = momentum_x/rho;
    		//o->t.g[3] = theRet.momY;
    		out->g[3] = momentum_y/rho;
    		out->g[4] = 1.9*strain_yx;
        out->g[5] = strain_yy-1 + 0.75*(strain_xx-1);
        //f->6] = g[6] = -0.295;
        out->f[6] = -out->f[0]*sqrt(out->f[4]*out->f[4] + out->f[5]*out->f[5]);
        out->g[6] = -out->g[3]*sqrt(out->g[4]*out->g[4] + out->g[5]*out->g[5]);
#endif//COMD
        //mark as CoMD'd
        in->callCoMD = true;
		    //Write result to database
		    putData(in->w->w, out->f, out->g, (char *)"comd", rTask, comdDigits);		
        *stopCo = getUnixTime();

    }
    else
    {
	    //Write result to database
	    putData(in->w->w, out->f, out->g, (char *)"krig", rTask, krigDigits);		
      *stopKr = getUnixTime();
    }
	}

	//All pertinent values are set, let's disconneect from Redis and return
	redisFree(rTask);
	//return 0;

}

/** checks if database values are within a certain threshold or exact 
 * @param *w0         field that needs to be computed
 * @param *wVec       field which is available from the database
 * @param dbThresh    threshold for the check
 * **/
bool ifConservedFieldsMatch(double * w0, std::vector<double *> * wVec, double dbThresh)
{
	if(wVec != NULL)
	{
		if(wVec->size() != 0)
		{
			double dist = 0.0;
			for(int i = 0; i < 7; i++)
			{
				dist += fabs(w0[i] - (*wVec)[0][i]);
        //printf("w0: %.16f, wIn: %.16f\n", w0[i], (*wVec)[0][i]);
			}
			if(dist <= dbThresh)
			{
        //printf("hit\n");
				return true;
			}
			else
			{
        //printf("false1\n");
				return false;
			}
		}
		else
		{
      //printf("false2\n");
			return false;
		}
	}
	else
	{
    //printf("false3\n");
		return false;
	}
}

/** executes the mapped fields in parallel 
 *FIXME
 * @param *in         flux input contains information about field on specific node (input)
 * @return *out       returns flux output information about computed fluxes on specific node (output)
 * @param *dbCache    database cache pointer (input)
 * **/
void doParallelCalls(Node * fields, Node * fluxes, Lattice l, std::list<gridPoint> * comdTasks, std::list<gridPoint> * krigTasks, std::map<std::string, std::vector<char *> > *dbCache, Calls* ca, Tms *tm)
{
	//fprintf(stdout, "Starting parallel \n");
	//fflush(stdout);

	int dim_x = l.dim_x;
	int dim_y = l.dim_y;

	//For both lists, enqueue unique tasks, map the rest
	std::map<Conserved, std::list<gridPoint>, ConservedComparator_c> taskMap;

	//std::vector<Conserved> taskRefs;

	std::vector<fluxInput *> fluxInArgs;

	//Do kriging first so as to prioritize the cheaper solution
	for(std::list<gridPoint>::iterator iter = krigTasks->begin(); iter != krigTasks->end(); iter++)
	{
		int x = iter->x;
		int y = iter->y;
		Conserved *w = &fields[x + dim_x*y].w;
		taskMap[*w].push_back(*iter);
    fields[x + dim_x*y].f.ca = 6;
    ca->kPoints++;
		//Check if it was unique
		if(taskMap[*w].size() == 1)
		{
			//It was, add it to the task queue
			fluxInput * fluxArg = new fluxInput;
      fluxArg->w = new Conserved;
			memcpy(fluxArg->w, w, sizeof(Conserved));
			fluxArg->callCoMD = false;
			strcpy(fluxArg->headNode, headNode);

			fluxInArgs.push_back(fluxArg);
      ca->krig++;
      ca->kPoints--;
      fields[x + dim_x*y].f.ca = 5;
		}
	}
	for(std::list<gridPoint>::iterator iter = comdTasks->begin(); iter != comdTasks->end(); iter++)
	{
		int x = iter->x;
		int y = iter->y;
		Conserved *w = &fields[x+dim_x*y].w;
		taskMap[*w].push_back(*iter);
    fields[x + dim_x*y].f.ca = 2;
    ca->cPoints++;
		//Check if it was unique
		if(taskMap[*w].size() == 1)
		{
			//It was, add it to the task queue
			fluxInput * fluxArg = new fluxInput;
      fluxArg->w = new Conserved;
			memcpy(fluxArg->w, w, sizeof(Conserved));
			fluxArg->callCoMD = true;
			strcpy(fluxArg->headNode, headNode);

			fluxInArgs.push_back(fluxArg);
      ca->comd++;
      ca->cPoints--;
      fields[x + dim_x*y].f.ca = 1;
		}
	}
	
	//Empty both task lists
	comdTasks->clear();
	krigTasks->clear();

	//Set up outputs
	fluxOutput * fluxOutArgs = new fluxOutput[fluxInArgs.size()];

	//init timers
  double startKr[int(fluxInArgs.size())];
  double stopKr[int(fluxInArgs.size())];
  double startCo[int(fluxInArgs.size())];
  double stopCo[int(fluxInArgs.size())];
	for(int i = 0; i < int(fluxInArgs.size()); i++)
  {
      startKr[i] = 0.0;
      stopKr[i] = 0.0;
      startCo[i] = 0.0;
      stopCo[i] = 0.0;
  }
    //tm->kr = 0.0;
    //tm->co = 0.0;
	//OMP all them tasks
#pragma omp parallel for
	for(int i = 0; i < int(fluxInArgs.size()); i++)
	{
    //printf("openmp task %i\n", i);
		fluxFn(fluxInArgs[i], &fluxOutArgs[i], dbCache, &startKr[i], &stopKr[i], &startCo[i], &stopCo[i]);
	}
	//fprintf(stdout, "finished parallel \n");
	//fflush(stdout);

	for(int i = 0; i < int(fluxInArgs.size()); i++)
  {
      tm->kr += (stopKr[i] - startKr[i]);
      tm->co += (stopCo[i] - startCo[i]);
  }
	//Process the results of OMP'd tasks
	for(int i = 0; i < int(fluxInArgs.size()); i++)
	{
		//Write results to fluxes for all duplicates as this is good
		for(std::list<gridPoint>::iterator iter = taskMap[*fluxInArgs[i]->w].begin(); iter != taskMap[*fluxInArgs[i]->w].end(); iter++)
		{
			int x = iter->x;
			int y = iter->y;
			Fluxes * f = &fluxes[x + dim_x*y].f;
			Fluxes * g = &fluxes[x + dim_x*y].g;

			memcpy(f->f, &fluxOutArgs[i].f, sizeof(double)*7);
			memcpy(g->f, &fluxOutArgs[i].g, sizeof(double)*7);
			if(fluxInArgs[i]->callCoMD == true)
      {
          if(fields[x + dim_x*y].f.ca != 1 && fields[x + dim_x*y].f.ca != 2){
            fields[x + dim_x*y].f.ca = 7;
            //fields[x + dim_x*y].f.ca = 1;
            ca->kFail++;
          }
          //ca->cPoints++;
      }
      //else
      //{
      //    ca->kPoints++;
      //}
		}
	}
  delete[] fluxOutArgs;
  for(std::vector<fluxInput *>::iterator iter = fluxInArgs.begin(); iter != fluxInArgs.end(); iter++)
  {
      delete (*iter)->w;
      delete *iter;
  }
  taskMap.clear();
}

/** checks the gradient of the fields (used to decide whether trying kriging or not) 
 * @param *fields     fields on specific node (input)
 * @param l           lattice information (input)
 * **/
bool checkGradient(int x, int y, Node * fields, Lattice l)
{
	double dx = l.dx;
	double dy = l.dy;
	int dimX = l.dim_x;
	int dimY = l.dim_y;

	double delta = 0.0;
	int xP = (x + 1 + dimX) % dimX;
	int xM = (x - 1 + dimX) % dimX;
	int yP = (y + 1 + dimY) % dimY;
	int yM = (y - 1 + dimY) % dimY;

	Conserved * w = &fields[x + dimX*y].w;
	Conserved * wXP = &fields[xP + dimX*y].w;
	Conserved * wXM = &fields[xM + dimX*y].w;
	Conserved * wYP = &fields[x + dimX*yP].w;
	Conserved * wYM = &fields[x + dimX*yM].w;
	for(int i = 0; i < 7; i++)
	{
		double gradient[4];
		gradient[0] = fabs(w->w[i] - wXP->w[i] / dx);
		gradient[1] = fabs(w->w[i] - wXM->w[i] / dx);
		gradient[2] = fabs(w->w[i] - wYP->w[i] / dy);
		gradient[3] = fabs(w->w[i] - wYM->w[i] / dy);
		for(int j = 0; j < 4; j++)
		{
			if(gradient[j] > delta)
			{
				delta = gradient[j];
			}
		}
	}
	if(delta <= gradThresh)
	{
		return true;
	}
	else
	{
		return false;
	}
}

/** main call, decides whether kriging is called or comd, maps all tasks for parallel execution 
 * @param *fields     fields  (input)
 * @param *fluxes     fluxes  (input)
 * @param l           lattice information (input)
 * **/
void doFluxes(Node* fields, Node* fluxes, int grid_size, Lattice l)
{
	//fprintf(stdout, "Starting doFluxes \n");
	//fflush(stdout);
	int dim_x = l.dim_x;
	int dim_y = l.dim_y;

    //measurement stuff
    Calls ca = { .comd=0, .cPoints=0, .db=0, .kdb=0, .krig=0, .kPoints=0, .kFail=0 };
    Tms tm = { .db=0.0, .krDb=0.0, .ca=0.0, .kr=0.0, .co=0.0, .kr2=0.0, .co2=0.0 };
    double startDb, stopDb;
    double startKrDb, stopKrDb;
    double startCa, stopCa;
    double startKr2, stopKr2;
    double startCo2, stopCo2;

	std::list<gridPoint> comdTasks;
	std::list<gridPoint> krigTasks;
	std::map<std::string, std::vector<char *> > dbCache; 

	//Iterate over every point to see who is getting krig'd and who is getting comd'd
	for(int y = 0; y < dim_y; y++)
	{
		for(int x = 0; x < dim_x; x++)
		{
			//Grab comd database
			std::vector<double *> wVec;
			std::vector<double *> fVec;
			std::vector<double *> gVec;
#ifdef DB
            startDb = getUnixTime();

			getCachedSortedSubBucketNearZero(fields[x + dim_x*y].w.w, (char *)"comd", headRedis, comdDigits, 2, &wVec, &fVec, &gVec, zeroThresh, &dbCache);
			//getSortedSubBucketNearZero(fields[x + dim_x*y].w.w, "comd", headRedis, comdDigits, 2, &wVec, &fVec, &gVec, zeroThresh);
            //if(int(wVec.size())>1)printf("wvecsize %i\n", int(wVec.size()));
            //printf("wvecsize %i\n", int(wVec.size()));

			//Check for exact value
            //printf("dbthresh %lf\n", dbT);
			bool useDB = ifConservedFieldsMatch(fields[x+dim_x*y].w.w, &wVec, dbT);
		
            stopDb = getUnixTime();
            tm.db += stopDb - startDb;
#else
            bool useDB = false;
#endif//DB
			//If we found it
			if(useDB == true)
			{
                ca.db++;
                fields[x + dim_x*y].f.ca = 3;
				//Set the result from the database
				memcpy(fluxes[x+dim_x*y].f.f, fVec[0], sizeof(double)*7);
				memcpy(fluxes[x+dim_x*y].g.f, gVec[0], sizeof(double)*7);
			
            }
			//We did not
			else
			{
				//Check the gradient
#ifdef KRIGING
				bool smallGradient = checkGradient(x, y, fields, l);
#else
				//If gradient is small enough
                bool smallGradient = false;
#endif//KRIGING
				if(smallGradient == true)
				{
					//Grab krig database
					std::vector<double *> wVecK;
					std::vector<double *> fVecK;
					std::vector<double *> gVecK;

#ifdef KR_DB
                    startKrDb = getUnixTime();

			        getCachedSortedSubBucketNearZero(fields[x + dim_x*y].w.w, (char *)"krig", headRedis, krigDigits, 1, &wVecK, &fVecK, &gVecK, zeroThresh, &dbCache);
					//getSortedSubBucketNearZero(fields[x + dim_x*y].w.w, "krig", headRedis, krigDigits, 1, &wVecK, &fVecK, &gVecK, zeroThresh);
                    //if(diff > 1.0) printf("get kdb timing %lf\n", diff);
			        //fflush(stdout);
					//Check for exact value with kriging key
					bool useDB = ifConservedFieldsMatch(fields[x+dim_x*y].w.w, &wVecK, 0.0);
					//If we found it
                    //if(int(wVec.size())>0)printf("wvecsize %i\n", int(wVec.size()));
                    stopKrDb = getUnixTime();
                    tm.krDb += stopKrDb - startKrDb;
#else
                    bool useDB = false;
#endif//KR_DB
					if(useDB == true)
					{
                        ca.kdb++;
                        fields[x + dim_x*y].f.ca = 4;
						//Set the result from the database
						memcpy(fluxes[x+dim_x*y].f.f, fVecK[0], sizeof(double)*7);
						memcpy(fluxes[x+dim_x*y].g.f, gVecK[0], sizeof(double)*7);

					}
                    //use at least 2 w's for kriging
					else if(int(wVec.size()) >= 2)
					{
                    //printf("try krig\n");
						//We are gonna krig this
						gridPoint thePoint;
						thePoint.x = x;
						thePoint.y = y;
						krigTasks.push_back(thePoint);
					}
				    //It was not, so mark it for a CoMD
          else
          {
            //Add to CoMD list
            gridPoint thePoint;
            thePoint.x = x;
            thePoint.y = y;
            comdTasks.push_back(thePoint);
          }
          freeClear(wVecK);
          freeClear(fVecK);
          freeClear(gVecK);
				}
				//It was not, so mark it for a CoMD
				else
				{
					//Add to CoMD list
					gridPoint thePoint;
					thePoint.x = x;
					thePoint.y = y;
					comdTasks.push_back(thePoint);
				}
			}
            freeClear(wVec);
            freeClear(fVec);
            freeClear(gVec);
		} //end for loop x
	} // end for loop y

  //fprintf(stdout, "prepare parallel calls \n");
  //fflush(stdout);
  startCa = getUnixTime();
  //At this point, we know who is getting krig'd and who is getting comd'd, so let's do it
  doParallelCalls(fields, fluxes, l, &comdTasks, &krigTasks, &dbCache, &ca, &tm);
  stopCa = getUnixTime();
  tm.ca = stopCa - startCa;

  if(ca.comd)tm.co /= ca.comd;
  if(ca.krig)tm.kr /= ca.krig;
  //print the calls
#ifdef OUTPUT
  printf_colormap_vtk(counter, fields, l, grid_size);
  printf_calls(counter, ca);
  printf_timings(counter, tm);
#endif//OUTPUT
  counter++;

  if(ca.krig != 0 && ca.kFail == 0){
    if(gradThresh < 10) gradThresh *= 2;
  }
  if(ca.krig != 0 && ca.kFail != 0){
    gradThresh *= 0.25;
  }
#ifdef OUTPUT
  printf("kFail %i, gradThresh %f\n", ca.kFail, gradThresh);
#endif//OUTPUT

  //clean the cache
  for(std::map<std::string, std::vector<char *> >::iterator iter = dbCache.begin(); iter != dbCache.end(); iter++)
  {  
     for(std::vector<char *>::iterator vecIt = iter->second.begin(); vecIt != iter->second.end(); vecIt++)
     {
         delete[] *vecIt;
     }
  }
  //FIXME free maybe every tens timestep
  dbCache.clear();
  //fprintf(stdout, "Done doFluxes \n");
  //fflush(stdout);
  //At this point, everything is done
  return;
}

/** executes second order half time step with average conserved value (w^{1/2}_{jk})
 * @param node_a      grid_node containing the conserved and the fluxes (input)
 * @param node_b      grid_node containing the conserved and the fluxes (output)
 * @param grid_size   number of grid nodes (input)
 * @param l           lattice parameters (input)         
 * **/
void half_step_second_order(Node* node_a, Node* node_b, int grid_size, Lattice l){

  int dim_x = l.dim_x;
  int dim_y = l.dim_y;
  double df[7] = {0.0};
  double dg[7] = {0.0};
  double mu = l.dt_x/l.dx;
  double nu = l.dt_y/l.dy;
  int xpo, xmo, ypo, ymo, xpoypo;
  int x = 0;
  int y = 0;

	int i;
  //Jiang&Tadmor '98 equ. (2.15) half step in time
  doFluxes(node_a, node_a, grid_size, l);

  for(i=0; i<grid_size; ++i){
    index_to_xy(i, l, &x, &y);
    // vars for indexing
    xmo = ((x-1+dim_x)%dim_x) + dim_x*y;
    xpo = ((x+1)%dim_x) + dim_x*y; 
    ymo = x + dim_x*((y-1+dim_y)%dim_y);
    ypo = x + dim_x*((y+1)%dim_y);

    // compute w_n^1/2 (half timestep)
    for(int j=0; j<7; ++j){
      // compute derivative with min_mod operation
      df[j] = min_mod(node_a[xpo].f.f[j], node_a[i].f.f[j], node_a[xmo].f.f[j]);
      dg[j] = min_mod(node_a[ypo].g.f[j], node_a[i].g.f[j], node_a[ymo].g.f[j]);
      //compute new conserved values and store it in grid b
      node_b[i].w.w[j] = node_a[i].w.w[j] - 0.5*mu*df[j] - 0.5*nu*dg[j];
    }
  }//loop over grid
  //Calc fluxes for second halfstep (time and space)
  //prepare computation of Jiang&Tadmor '98 equ. (2.16)

  doFluxes(node_b, node_b, grid_size, l);


  //temp vars
  // w is w' e.g. w_j1k1 = w'_{j+1,k+1} 
  double w_jk   = 0.;
  double w_j1k  = 0.;
  double w_jk1  = 0.;
  double w_j1k1 = 0.;

  // v is w` e.g. v_j1k1 = w`_{j+1,k+1} 
  double v_jk   = 0.;
  double v_j1k  = 0.;
  double v_jk1  = 0.;
  double v_j1k1 = 0.;

  // flux derivative x e.g. f_jk = f(w_{j,k})
  double f_jk   = 0.;
  double f_j1k  = 0.;
  double f_jk1  = 0.;
  double f_j1k1 = 0.;

  // flux derivative y e.g. g_jk = g(w_{j,k})
  double g_jk   = 0.;
  double g_j1k  = 0.;
  double g_jk1  = 0.;
  double g_j1k1 = 0.;

  //Jiang&Tadmor '98 equ. (2.16) full time step, half grid step
  for(int i=0; i<grid_size; ++i){
    index_to_xy(i, l, &x, &y);
    // temp vars for indexing
    xmo = ((x-1+dim_x)%dim_x) + dim_x*y;
    xpo = ((x+1)%dim_x) + dim_x*y; 
    ymo = x + dim_x*((y-1+dim_y)%dim_y);
    ypo = x + dim_x*((y+1)%dim_y);
    xpoypo = (x+1)%dim_x + dim_x*((y+1)%dim_y);
    //loop over fluxes, new conserved values stored in grid b
    for(int j=0; j<7; ++j){
      //first part of summation
      node_b[xpoypo].w.w[j] = 0.25*(node_a[i].w.w[j] + node_a[xpo].w.w[j] + node_a[ypo].w.w[j] + node_a[xpoypo].w.w[j]);
      //temp vars for summation
      w_jk   = min_mod(node_a[xpo].w.w[j], node_a[i].w.w[j], node_a[xmo].w.w[j]);
      w_j1k  = min_mod(node_a[(x+2)%dim_x + dim_x*y].w.w[j], node_a[xpo].w.w[j], node_a[i].w.w[j]);
      w_jk1  = min_mod(node_a[xpoypo].w.w[j], node_a[ypo].w.w[j], node_a[(dim_x+x-1)%dim_x + (dim_x*((y+1)%dim_y))].w.w[j]);
      w_j1k1 = min_mod(node_a[(x+2)%dim_x + dim_x*((y+1)%dim_y)].w.w[j], node_a[xpoypo].w.w[j], node_a[ypo].w.w[j]);

      v_jk   = min_mod(node_a[ypo].w.w[j], node_a[i].w.w[j], node_a[ymo].w.w[j]);
      v_j1k  = min_mod(node_a[xpoypo].w.w[j], node_a[xpo].w.w[j], node_a[(x+1)%dim_x + dim_x*((dim_y+y-1)%dim_y)].w.w[j]);
      v_jk1  = min_mod(node_a[x + dim_x*((y+2)%dim_y)].w.w[j], node_a[ypo].w.w[j], node_a[i].w.w[j]);
      v_j1k1 = min_mod(node_a[(x+1)%dim_x + dim_x*((y+2)%dim_y)].w.w[j], node_a[xpoypo].w.w[j], node_a[xpo].w.w[j]);

      f_jk   = node_b[i].f.f[j];
      f_j1k  = node_b[xpo].f.f[j];
      f_jk1  = node_b[ypo].f.f[j];
      f_j1k1 = node_b[xpoypo].f.f[j];

      g_jk   = node_b[i].g.f[j];
      g_j1k  = node_b[xpo].g.f[j];
      g_jk1  = node_b[ypo].g.f[j];
      g_j1k1 = node_b[xpoypo].g.f[j];
      //second part of summation in equ (2.16)
      node_b[xpoypo].w.w[j] += 0.0625*(w_jk - w_j1k + w_jk1 - w_j1k1) - 0.5*mu*(f_j1k - f_jk + f_j1k1 - f_jk1);
      node_b[xpoypo].w.w[j] += 0.0625*(v_jk - v_jk1 + v_j1k - v_j1k1) - 0.5*nu*(g_jk1 - g_jk + g_j1k1 - g_j1k);
      //printf("nodeb[%i].w[%i] = %g\n", xpoypo, j,  node_b[xpoypo].w.w[j]);
    }// loop over fluxes 
  }//loop over grid

}

/**********************************************/
/** Main of the 2D macro solver, executes the full timestep loop
@para argc
@para argv
**/
int main(int argc, char **argv) {

  //default values
	Lattice l;
	l.dim_x = 10;
	l.dim_y = 10;
	l.dx = 10;
	l.dy = 1.0;
	l.dt_x = 0.1;
	l.dt_y = 0.1;
//  Lattice l = {.dim_x=10, .dim_y=10,
//               .dx=1.0, .dy=1.0,
//               .dt_x=0.1, .dt_y=0.1};
#ifdef OUTPUT
//#######################################//
  // open pipe to gnuplot instance
  FILE * gp;
#ifdef GIF
  gp = popen(GNUPLOT,"w");
  if (gp==NULL) {
    printf("Error opening pipe to GNU plot < %s > t! \n", GNUPLOT);
    exit(0);
  }
   	
  //Use this to see the output/avoid seg faults
  fprintf(gp, "set terminal gif animate delay 10 \n");
  fprintf(gp, "set output \'cohmm.gif\' \n");
#endif//GIF
//#######################################//
  //reset calls.dat
  int bufferSize = 256;
  char file_name[bufferSize];
  sprintf(file_name,"calls.dat");
  FILE *fn = fopen(file_name, "w");
  if (fn == NULL) {
    printf("Error writing file< %s >!\n", file_name);
    exit(0);
  }
  fprintf(fn, "#NO      COMD       COMD_P      DB      KR_DB      KR      KR_P      KR_FAIL\n");
  fclose(fn);

  //reset timings.dat
  sprintf(file_name,"timings.dat");
  FILE *fn2 = fopen(file_name, "w");
  if (fn2 == NULL) {
    printf("Error writing file< %s >!\n", file_name);
    exit(0);
  }
  fprintf(fn2, "#NO      COMD       COMD_P      DB      KR_DB      KR      KR_P      KR_FAIL\n");
  fclose(fn2);
#else
  //init total_time.dat
  int bufferSize = 256;
  char file_name[bufferSize];
  sprintf(file_name,"total_time.dat");
  FILE *fn3 = fopen(file_name, "a");
  if (fn3 == NULL) {
    printf("Error writing file< %s >!\n", file_name);
    exit(0);
  }

  if(argc != 8){
	fprintf(stderr, "%s $(DIMX) $(DIMY) $(NSTEPS) $(REDIS_HOST) $(DBTHRESH) $(KTHRESH) $(GAMMA)\n",argv[0]);
	return 1;
  }
#endif//OUTPUT
	// Init Redis
	headRedis = redisConnect(argv[4], 6379);
	if(headRedis == NULL || headRedis->err)
	{
		printf("Redis error: %s\n", headRedis->errstr);
		return -1;
	}
	strcpy(headNode, argv[4]);

  //init field and threshold values
  l.dim_x = atoi(argv[1]);
  l.dim_y = atoi(argv[2]);
  dbT = atof(argv[5]);
  errorThresh = atof(argv[6]);
  mom_gamma = atof(argv[7]);
  int time_steps = atoi(argv[3]);
  printf("set database Tresh to: %.8f\n", dbT);
  printf("set kriging Tresh to: %.8f\n", errorThresh);
  printf("set gamma to: %.8f\n", mom_gamma);
  Values init;
  init.strain[0]=1.0;
  init.strain[1]=0.0;
  init.strain[2]=0.0;
  init.strain[3]=1.0;
  init.momentum[0]=0.0;
  init.momentum[1]=0.0;
  init.momentum[2]=0.0;
  init.momentum[3]=0.0;
  init.energy=-0.295;
  init.rho=1.0; 
  //lat = atof(argv[4]);

  //grid size (1D indexing)
  int grid_size = l.dim_x*l.dim_y;
  //double buffering
  Node nodes_a[grid_size]; 
  Node nodes_b[grid_size]; 
  //set all values to zero  
  init_nodes(nodes_a, grid_size);
  init_nodes(nodes_b, grid_size);
  //initial fiels on nodes
  init_conserved_fields(nodes_a, l, grid_size, init);
  //define warnings
#ifndef DB
    printf("Database turned off!!!\n Check for defines\n");
#endif//DB
#ifndef KRIGING
    printf("Kriging turned off!!!\n Check for defines\n");
#endif//KRIGING
#ifndef KR_DB
    printf("Kriging database turned off!!!\n Check for defines\n");
#endif//KR_DB
  //start total time measurement
  double ttime_start = getUnixTime();
  double ttime_stop;

  //set_boundaries(l, grid_size, nodes_a, nodes_b);
  for(int i=0; i<time_steps; ++i){
#ifdef OUTPUT
#ifdef GIF
    plot_fields(i, nodes_a, l, gp);
#else
    plot_fields(i, nodes_a, l, NULL);
#endif//GIF
#endif//OUTPUT
    //printf_fields_vtk(i, nodes_a, l, grid_size);
    for (int j = 0; j < 1; j++){
	    fprintf(stdout, "=================== Step %d =================\r", i);
	    fflush(stdout);
        half_step_second_order(nodes_a, nodes_b, grid_size, l);
        half_step_second_order(nodes_b, nodes_a, grid_size, l);
        shift_back(nodes_b, grid_size, l, nodes_a);
        //boundaries(l, grid_size, nodes_a);
        //stop total time measurement
#ifdef TRACE
		writeMessage("STEP_DONE", headRedis);
#endif
#ifndef OUTPUT
        ttime_stop = getUnixTime();
        fprintf(fn3, "%g\n", ttime_stop-ttime_start);
        printf("time: %g\n", ttime_stop-ttime_start);
#endif//OUTPUT
    }
  }
  //exit output
  fprintf(stdout, "Done allSteps \n");
  fflush(stdout);
  //cleanup redis
  redisFree(headRedis);
#ifndef OUTPUT
  //close total time file
  fclose(fn3);
#endif//OUTPUT

 return 0;
}
