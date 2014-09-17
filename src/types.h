/** header file for 2D macro solver, inlcudes 
 * typdefs
 * utilities
 * **/

#ifndef _TYPES_H
#define _TYPES_H

#include <string>
#ifdef CHARM
#include <charm++.h>
#endif
#include <stdio.h>

using std::string;
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

#ifdef CHARM
  void pup(PUP::er &p) {
    PUParray(p, w, 7);
    p|rho;
  }
#endif
#ifdef CIRCLE
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version) {
    for(int i=0;i<7;i++){
      ar & w[i];
    }
    ar & rho;
  }
#endif  
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
  int int_steps;
  int redis_db;
  int flush_db;
  double db_threshold;
  int kriging;
  double kr_threshold;
  int kriging_db;
  int gauss_noise;
  double noise; 
  double grad_threshold;
  int test_problem;
  int fault_tolerance;

} Save_Input;
/** struct containing the 2D field 
 * **/
typedef struct : public Save_Input {

  string head_node;

#ifdef CHARM
  void pup(PUP::er &p){
    //pipe arguments
    p|dim_x;
    p|dim_y;
    p|dx;
    p|dy;
    p|dt_x;
    p|dt_y;
    p|int_steps;
    p|redis_db;
    p|flush_db;
    p|db_threshold;
    p|kriging;
    p|kr_threshold;
    p|kriging_db;
    p|gauss_noise;
    p|noise; 
    p|grad_threshold;
    p|head_node;
    p|test_problem;
    p|fault_tolerance;
    p|flush_db;
  }

#endif
} Input;

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
/** struct containing the initial values for the mini_app
 * **/
typedef struct{

  int comd_on;
  string pot_name;
  string pot_type;
  int eam_on;
  int dim_x;
  int dim_y;
  int dim_z;
  int integration_steps;
  double lattice_spacing;
  double dt;

#ifdef CHARM
  void pup(PUP::er &p){
    p|comd_on;
    p|pot_name;
    p|pot_type;
    p|eam_on;
    p|dim_x;
    p|dim_y;
    p|dim_z;
    p|integration_steps;
    p|lattice_spacing;
    p|dt;
  }
#endif

} App;

/** calc x,y coordiante out of index
 * @param index  1-dim index of the grid node   (input) 
 * @param l      lattice parameters (input)         
 * @param x      x-coordinate of the grid node (output)
 * @param y      y-coordinate of the grid node (output)
 * **/
inline void index_to_xy(int index, Input in, int* x, int* y){

  *x = index%in.dim_x;
  index /= in.dim_x;
  *y = index%in.dim_y;

}

/** calc index out of x,y coordinate
 * @param x      x-coordinate of the grid node (input)
 * @param y      y-coordinate of the grid node (input)
 * @param index  1-dim index of the grid node (output) 
 * @param l      lattice infos
 * **/
inline void xy_to_index(int x, int y, int* index, Input in){

  *index = x + y*in.dim_x;

}
static inline void loadBar(int x, int n, int r, int w)
{
  // Only update r times.
  if ( x % (n/r) != 0 ) return;
  // // Calculuate the ratio of complete-to-incomplete.
  float ratio = x/(float)n;
  int c = ratio * w;
  // //Custom 2013 summer school loading bar, essential to the code
  if((int)(ratio*100)==2)   printf("\n\n:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n");
  if((int)(ratio*100)==4)   printf(":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n");
  if((int)(ratio*100)==8)   printf(":::::::::::::::::::::::::::::::::::::::::::::-'    `-::::::::::::::::::\n");
  if((int)(ratio*100)==12)  printf("::::::::::::::::::::::::::::::::::::::::::-'          `::::::::::::::::\n");
  if((int)(ratio*100)==16)  printf(":::::::::::::::::::::::::::::::::::::::-  '   /(_M_)\\  `:::::::::::::::\n");
  if((int)(ratio*100)==20)  printf(":::::::::::::::::::::::::::::::::::-'        |       |  :::::::::::::::\n");
  if((int)(ratio*100)==24)  printf("::::::::::::::::::::::::::::::::-         .   \\/~V~\\/  ,:::::::::::::::\n");
  if((int)(ratio*100)==28)  printf("::::::::::::::::::::::::::::-'             .          ,::::::::::::::::\n");
  if((int)(ratio*100)==32)  printf(":::::::::::::::::::::::::-'                 `-.    .-::::::::::::::::::\n");
  if((int)(ratio*100)==36)  printf(":::::::::::::::::::::-'                  _,,-::::::::::::::::::::::::::\n");
  if((int)(ratio*100)==40)  printf("::::::::::::::::::-'                _,--:::::::::::::::::::::::::::::::\n");
  if((int)(ratio*100)==44)  printf("::::::::::::::-'               _.--:::::::::::::::::::::::#####::::::::\n");
  if((int)(ratio*100)==48)  printf("::::::::-'    ##     ###.-:::::###::::::::::::::::::::::::#####::::####\n");
  if((int)(ratio*100)==52)  printf("::::-'       ###_.::######:::::###::::::::::::::#####:##########:::####\n");
  if((int)(ratio*100)==56)  printf(":'         .-###::########:::::###::::::::::::::#####:##########:::####\n");
  if((int)(ratio*100)==60)  printf("      ..--:::###::########:::::###:::::######:::#####:##########:::####\n");
  if((int)(ratio*100)==64)  printf(" _.-':::##:::###:#########:::::###:::::######:::#####:#################\n");
  if((int)(ratio*100)==68)  printf("'#########:::###:#########::#########::######:::#####:#################\n");
  if((int)(ratio*100)==72)  printf(":#########:::#############::#########::######:::#######################\n");
  if((int)(ratio*100)==76)  printf("##########:::########################::################################\n");
  if((int)(ratio*100)==80)  printf("##########:::##########################################################\n");
  if((int)(ratio*100)==84)  printf("##########:::##########################################################\n");
  if((int)(ratio*100)==88)  printf("#######################################################################\n");
  if((int)(ratio*100)==92)  printf("#######################################################################\n");
  if((int)(ratio*100)==96)  printf("##################################################### DGR && RSP ######\n");
  if((int)(ratio*100)==98)  printf("#######################################################################\n\n\n");
}
#endif
