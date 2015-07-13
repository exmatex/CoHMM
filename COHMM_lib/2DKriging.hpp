#ifndef KRIGING2D_HPP
#define KRIGING2D_HPP

#include <string>
#include <cmath>
#include <vector>

const int comdDigits = 4;
const int krigDigits = 8;
const double dbT = 0.0001;
const double zeroThresh = 0.0001;
const double gradThresh = 0.001;
const double errorThresh = 0.01;

/** struct containing the fluxes
 * **/
struct Fluxes
{
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
};

/** struct containg the conserved quantities
 * **/
struct Conserved
{
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

	bool operator<(const Conserved &other) const
	{
		//Need to compare 7 different values... somehow
		//Build two "hashes"
		std::string ls((char *)this, sizeof(Conserved));
		std::string rs((char *)&other, sizeof(Conserved));
		//Compare them
		return ls < rs;
	}
};

/** struct containing the 2D field
 * **/
struct Save_Input
{
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
};

/** struct containing the 2D field
 * **/
struct Lattice : public Save_Input
{
	int dim_x;
	int dim_y;
	double dx;
	double dy;
	double dt_x;
	double dt_y;
};

typedef Lattice Input;

/** struct containing the node values
 *  **/
struct Node
{
	Fluxes f;
	Fluxes g;
	Conserved w;
	int boundary;
};

/** struct containing the initial values given via command line
 * **/
struct Values
{
	//strain
	double strain[4];
	//momentum density
	double momentum[4];
	//energy density
	double energy;
	//rho0
	double rho;
};

/** struct containing the counts
 * **/
struct Calls
{
	int comd;
	int cPoints;
	int db;
	int kdb;
	int krig;
	int kPoints;
	int kFail;
};

/** struct containing the timer
 * **/
struct Tms
{
	double db;
	double krDb;
	double ca;
	double kr;
	double co;
	double kr2;
	double co2;
};

struct gridPoint
{
	int x;
	int y;
};

struct ConservedComparator_c
{
	bool operator() (const Conserved& lhs, const Conserved& rhs) const
	{
		//Need to compare 7 different values... somehow
		//Build two "hashes"
		std::string ls((char *)&lhs, sizeof(Conserved));
		std::string rs((char *)&rhs, sizeof(Conserved));
		//Compare them
		return ls < rs;
	}
};

struct FluxIn
{
	///TODO: Do we need anything else? Rest of data should be in the function call
	Conserved fields;
	bool tryKriging;
};

struct FluxOut
{
	///TODO: Decide if we keep w[7] here. More storage, but fewer calls
	Conserved fields;
	double f[7];
	double g[7];
	double error;
};

struct FluxFuture
{
	Fluxes f;
	Fluxes g;
	unsigned int taskID;
	bool alreadyComputed;
};

struct RetryRedirect
{
	unsigned int realTaskID;
};

void init_conserved_fields(Node* node_a, int * dims, int grid_size);
void shift_back(Node* node_b, int grid_size, int * dims, Node* node_a);
double min_mod(double w_plus, double w, double w_minus);
void shift_back(Node* node_b, int grid_size, int * dims, Node* node_a);
bool checkGradient(int x, int y, Node * fields, int dims[2], double deltas[2]);
void wNSqrt(Node * field, int * dims, double *dt, double * delta);
void wSummation(Node * node_a, Node * node_b, int * dims, double * dt, double * delta);
bool ifConservedFieldsMatch(double * w0, std::vector<double *> * wVec, double dbThresh);
void printf_fields_vtk(int i, Node* node_a, Lattice l, int grid_size);

template <class C> void freeClear( C & cntr )
{
    for ( typename C::iterator it = cntr.begin(); it != cntr.end(); ++it )
	{
		delete[] *it;
	}
	cntr.clear();
}

#endif
