/** Output for gnuplot, vtk, terminal
 * **/

#include <stdio.h>
#include <stdlib.h>
#ifdef CHARM
#include "main.decl.h"
#endif
//extern "C"
//{
#include "types.h"
//}
#define GNUPLOT "/usr/bin/gnuplot -persist"
#include <unistd.h>
//#define PS_FILES
/** print the calls
 * @param i         number of output file (input)         
 * @param node      grid_node containing the conserved and the fluxes (input)
 * @param grid_size  number of grid nodes (input)
 * **/
void printf_calls(int counter, Calls ca){

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

/** print the timings
 * @param i         number of output file (input)         
 * @param node      grid_node containing the conserved and the fluxes (input)
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
 * @param in         lattice parameters (input)         
 * @param grid_size  number of grid nodes (input)
 * **/
void printf_colormap_vtk(int i, Node* fields, Input in, int grid_size){

  int bufferSize = 256;
  char file_name[bufferSize];
  sprintf(file_name,"colormap_%i.vtk",i);
  FILE *fn = fopen(file_name, "w");
  if (fn == NULL) {
    printf("Error writing file< %s >!\n", file_name);
    exit(0);
  }

  fprintf(fn, "# vtk DataFile Version 2.0\nstrain\nASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS %i %i 1\nORIGIN 0.0 0.0 0.0\nSPACING %lf %lf 1.0\nPOINT_DATA %i\nSCALARS OutArray floats 1\nLOOKUP_TABLE default\n", in.dim_x, in.dim_y, in.dx, in.dy, grid_size); 
  for(int index=0; index<grid_size; ++index){
    fprintf(fn, "%i\n", fields[index].f.ca);
  }
  fclose(fn);
}
/** print the output in vtk format
 * @param i         number of output file (input)         
 * @param node      grid_node containing the conserved and the fluxes (input)
 * @param in         lattice parameters (input)         
 * @param grid_size  number of grid nodes (input)
 * **/
void printf_fields_vtk(int i, Node* node_a, Input in, int grid_size){

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

  fprintf(fn, "# vtk DataFile Version 2.0\nstrain\nASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS %i %i 1\nORIGIN 0.0 0.0 0.0\nSPACING %lf %lf 1.0\nPOINT_DATA %i\nSCALARS OutArray floats 4\nLOOKUP_TABLE default\n", in.dim_x, in.dim_y, in.dx, in.dy, grid_size); 
  fprintf(fn2, "# vtk DataFile Version 2.0\nmom\nASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS %i %i 1\nORIGIN 0.0 0.0 0.0\nSPACING %lf %lf 1.0\nPOINT_DATA %i\nSCALARS OutArray floats 2\nLOOKUP_TABLE default\n", in.dim_x, in.dim_y, in.dx, in.dy, grid_size); 
  fprintf(fn3, "# vtk DataFile Version 2.0\nenergy\nASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS %i %i 1\nORIGIN 0.0 0.0 0.0\nSPACING %lf %lf 1.0\nPOINT_DATA %i\nSCALARS OutArray floats 1\nLOOKUP_TABLE default\n", in.dim_x, in.dim_y, in.dx, in.dy, grid_size); 
  fprintf(fn4, "# vtk DataFile Version 2.0\nmom_flux\nASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS %i %i 1\nORIGIN 0.5 0.5 0.5\nSPACING %lf %lf 1.0\nPOINT_DATA %i\nSCALARS OutArray floats 2\nLOOKUP_TABLE default\n", in.dim_x, in.dim_y, in.dx, in.dy, grid_size); 
  fprintf(fn5, "# vtk DataFile Version 2.0\nflux\nASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS %i %i 1\nORIGIN 0.5 0.5 0.5\nSPACING %lf %lf 1.0\nPOINT_DATA %i\nSCALARS OutArray floats 4\nLOOKUP_TABLE default\n", in.dim_x, in.dim_y, in.dx, in.dy, grid_size); 
  fprintf(fn6, "# vtk DataFile Version 2.0\nenergy_flux\nASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS %i %i 1\nORIGIN 0.5 0.5 0.5\nSPACING %lf %lf 1.0\nPOINT_DATA %i\nSCALARS OutArray floats 2\nLOOKUP_TABLE default\n", in.dim_x, in.dim_y, in.dx, in.dy, grid_size); 
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
/** plot the strain, stress and energy dens to gnuplot .ps files
 * @param gp         filepointer for gnuplot output (input)         
 * @param node       grid_node containing the conserved and the fluxes (output)
 * @param in          lattice parameters (input)         
 * **/
void plot_fields(int i, Node* node_a, Input in){

//########################################//
  int bufferSize = 256;
  char file_name[bufferSize];
  sprintf(file_name,"gp_%i.dat",i);
  //sleep(2);
  FILE *gp_out = fopen(file_name, "w");
  if (gp_out == NULL) {
    printf("Error writing file< %s >!\n", file_name);
    exit(0);
  }
#ifdef PS_FILES
  FILE* gp = popen(GNUPLOT,"w");
  if (gp==NULL) {
    printf("Error opening pipe to GNU plot < %s > t! \n", GNUPLOT);
    exit(0);
  }
  fprintf(gp, "set terminal postscript portrait size 8.5cm,6.3cm enhanced color dashed lw 1 \"Helvetica,8\"\n");
  fprintf(gp, "set output \'colorplot_%04d.ps\' \n", i);
  fprintf(gp, "set multiplot \n");
  fprintf(gp, "set label 1 'y [0..%i]' at -15.0,5.0,0.0 rotate by -35.5\n", in.dim_y);
  fprintf(gp, "set label 2 'x [0..%i]' at 60.0, 22.5,0.0 rotate by 4.5\n", in.dim_x);
  fprintf(gp, "set label 3 'Strain [MPa]' at -110.0,10.0,47.5\n");
  fprintf(gp, "set label 4 'Calltype' at -110.0,10.0,56.0\n");
  fprintf(gp, "set border 4095 front linetype -1 linewidth 0.500\n");
  fprintf(gp, "set style line 100  linetype 5 linecolor rgb \"#f0e442\"  linewidth 0.500 pointtype 5 pointsize default \n");
  fprintf(gp, "set palette rgbformulae 22,13,-31\n");
  fprintf(gp, "set view 110, 20, 1, 1\n");
  fprintf(gp, "set samples %i,%i\n", (2*in.dim_x), (2*in.dim_y));
  fprintf(gp, "set isosample %i,%i\n", (2*in.dim_x), (2*in.dim_y));
  fprintf(gp, "set ticslevel 0\n");
  fprintf(gp, "set xrange [ :%i ] noreverse nowriteback\n", (in.dim_x-1)*2);
  fprintf(gp, "set yrange [ :%i ] noreverse nowriteback\n", (in.dim_y-1)*2);
  fprintf(gp, "set zrange [ 0.0 : 35 ] noreverse nowriteback\n");
  fprintf(gp, "unset xtics\n");
  fprintf(gp, "unset ytics\n");
  fprintf(gp, "set pm3d implicit at s\n");
  fprintf(gp, "set pm3d interpolate 1,1 flush begin noftriangles hidden3d 100 corners2color mean\n");
   //fprintf(gp, "set pm3d \n");
  fprintf(gp, "set colorbox vertical user origin 0.89, .135 size .04,.7\n");
  fprintf(gp, "splot '-' u 1:2:((sqrt($3**2+$4**2+$5**2+$6**2))*1000) w l lc rgb \"black\" lw 0.1 ti ''\n\n");
#endif
  // gnuplot forces us to repeat ourselves
  int index = 0;
  int x = 0;
  int y = 0;
  for (int k=0; k<1; ++k) {
    for (int i=0; i<(2*in.dim_x); i++) {
        y=0;
      for (int j=0; j<(2*in.dim_y); j+=2) {
        xy_to_index(x, y, &index, in);
        y++;
#ifdef PS_FILES
        //gnuplot ps 
        fprintf(gp, "%f %f %f %f %f %f\n", i*in.dx, j*in.dy, node_a[index].w.w[0]-1, node_a[index].w.w[1], node_a[index].w.w[2], node_a[index].w.w[3]-1);
        fprintf(gp, "%f %f %f %f %f %f\n", i*in.dx, (j+1)*in.dy, node_a[index].w.w[0]-1, node_a[index].w.w[1], node_a[index].w.w[2], node_a[index].w.w[3]-1);
#endif
        //gnuplot dat file
        fprintf(gp_out, "%f %f %f %f %f %f %i\n", i*in.dx, j*in.dy, node_a[index].w.w[0]-1, node_a[index].w.w[1], node_a[index].w.w[2], node_a[index].w.w[3]-1, node_a[index].f.ca);
        fprintf(gp_out, "%f %f %f %f %f %f %i\n", i*in.dx, (j+1)*in.dy, node_a[index].w.w[0]-1, node_a[index].w.w[1], node_a[index].w.w[2], node_a[index].w.w[3]-1, node_a[index].f.ca);
      }
      if(i%2 == 1) x++;
#ifdef PS_FILES
      fprintf(gp, "\n");
#endif
      fprintf(gp_out, "\n");
    }
#ifdef PS_FILES
    fprintf(gp, "e\n");
#endif
    fprintf(gp_out, "e\n");
  }
#ifdef PS_FILES
  fflush(gp);
  fprintf(gp, "reset\n");
  fprintf(gp, "unset surface\n");
  fprintf(gp, "set border 4095 front linetype -1 linewidth 0.500\n");
  fprintf(gp, "set view 110, 20, 1, 1\n");
  fprintf(gp, "set samples %i,%i\n", (2*in.dim_x), (2*in.dim_y));
  fprintf(gp, "set isosample %i,%i\n", (2*in.dim_x), (2*in.dim_y));
  fprintf(gp, "set ticslevel 0\n");
  fprintf(gp, "set autoscale cbfix\n");
  fprintf(gp, "set pm3d implicit at t\n");
  fprintf(gp, "set xrange [ :%i ] noreverse nowriteback\n", (in.dim_x-1)*2);
  fprintf(gp, "set yrange [ :%i ] noreverse nowriteback\n", (in.dim_y-1)*2);
  fprintf(gp, "set pal maxcolors 6\n");
  fprintf(gp, "set cbrange [ 1 : 6 ] noreverse nowriteback\n");
  fprintf(gp, "set palette defined (1 \"blue\", 2 \"light-blue\", 3 \"turquoise\", 4 \"dark-green\", 5 \"red\", 6 \"orange\")\n");
  fprintf(gp, "set cbtics (\"CoMD\" 1, \"C. Dupl.\" 2, \"DB\" 3, \"Kr. DB\" 4, \"Kr.\" 5, \"Kr. Dupl.\" 6)\n");
  fprintf(gp, "unset xtics\n");
  fprintf(gp, "unset ytics\n");
  fprintf(gp, "set colorbox horizontal user origin .155, 0.925 size .7,.04\n");
  fprintf(gp, "splot '-' u 1:2:3 w l lc rgb \"black\" lw 2 ti '' \n");
  x = 0;
  for (int k=0; k<1; ++k) {
    for (int i=0; i<(2*in.dim_x); ++i) {
        y=0;
      for (int j=0; j<(2*in.dim_y); j+=2) {
        xy_to_index(x, y, &index, in);
        y++;
        fprintf(gp, "%f %f %i \n", i*in.dx, j*in.dy, node_a[index].f.ca);
        fprintf(gp, "%f %f %i \n", i*in.dx, (j+1)*in.dy, node_a[index].f.ca);
      }
      if(i%2 == 1) x++;
      fprintf(gp, "\n");
    }
    fprintf(gp, "e\n");
  }
  fflush(gp);
  fclose(gp);
#endif
  fclose(gp_out);
}
