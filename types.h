/** header file for 2D macro solver, inlcudes mainly typedefs
 *
 * **/

#ifndef _TYPES_H
#define _TYPES_H

/** struct containing the fluxes
 * **/
typedef struct{

  //use of the following notation:
  //momentum_dens_x = w[0] //for g=0 
  //momentum_dens_x = w[1] //for f=0
  //momentum_dens_y = w[2] //for g=0
  //momentum_dens_y = w[3] //for f=0
  //stress_fxx/gyx  = w[4]
  //stress_fxy/gyy  = w[5]
  //energy_dens     = w[6]
  double f[7];
  int ca;

} Fluxes;

/** struct containg the conserved quantities
 * **/
typedef struct{

  //use of the following notation:
  //strain_xx  = w[0]
  //strain_yx  = w[1]
  //strain_xy  = w[2]
  //strain_yy  = w[3]
  //momentum_x = w[4]
  //momentum_y = w[5]
  //energy     = w[6]
  double w[7];
  double rho;

} Conserved;

/** struct containing the 2D field 
 * **/
typedef struct{

  int dim_x;
  int dim_y;
  double dx;
  double dy;
  double dt_x;
  double dt_y;

} Lattice;

/** struct containing the node values
 *  **/
typedef struct{

  Fluxes f;
  Fluxes g;
  Conserved w;
  int boundary;

} Node;

/** struct containing the initial values given via command line
 * **/
typedef struct{

  //strain
  double strain[4];
  //momentum density
  double momentum[4];
  //energy density
  double energy;
  //rho0
  double rho;

} Values;


/** struct containing the counts 
 * **/
typedef struct{

  int comd;
  int cPoints;
  int db;
  int kdb;
  int krig;
  int kPoints;
  int kFail;

} Calls;

/** struct containing the timer 
 * **/
typedef struct{

  double db;
  double krDb;
  double ca;
  double kr;
  double co;
  double kr2;
  double co2;

} Tms;

/** calc x,y coordiante out of index
 * @param index  1-dim index of the grid node   (input) 
 * @param l      lattice parameters (input)         
 * @param x      x-coordinate of the grid node (output)
 * @param y      y-coordinate of the grid node (output)
 * **/
inline void index_to_xy(int index, Lattice l, int* x, int* y){

  *x = index%l.dim_x;
  index /= l.dim_x;
  *y = index%l.dim_y;

}

/** calc index out of x,y coordinate
 * @param x      x-coordinate of the grid node (input)
 * @param y      y-coordinate of the grid node (input)
 * @param index  1-dim index of the grid node (output) 
 * @param l      lattice infos
 * **/
inline void xy_to_index(int x, int y, int* index, Lattice l){

  *index = x + y*l.dim_x;

}
#endif
