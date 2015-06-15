/** 2D Macrosolver written by
 * Dominic Roehm
 with the help of
 * Robert Pavel
 * CoDesign summer school 2013
 */

#include <cstring>
#include <cstdio>
#include <cmath>
#include <vector>

#include "2DKriging.hpp"



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
			}
			if(dist <= dbThresh)
			{
				return true;
			}
			else
			{
				return false;
			}
		}
		else
		{
			return false;
		}
	}
	else
	{
		return false;
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

/** checks the gradient of the fields (used to decide whether trying kriging or not) 
 * @param *fields     fields on specific node (input)
 * @param l           lattice information (input)
 * **/
bool checkGradient(int x, int y, Node * fields, int dims[2], double deltas[2])
{
	double dx = deltas[0];
	double dy = deltas[1];
	int dimX = dims[0];
	int dimY = dims[1];

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

void wSummation(Node * node_a, Node * node_b, int * dims, double * dt, double * delta)
{
	int dim_x = dims[0];
	double mu = dt[0] / delta[0];
	double nu = dt[1] / delta[1];
	int dim_y = dims[1];
	for(int y = 0; y < dims[1]; y++)
	{
		for(int x = 0; x < dims[0]; x++)
		{
			int i = x + dim_x * y;
			int xmo = ((x-1+dim_x)%dim_x) + dim_x*y;
			int xpo = ((x+1)%dim_x) + dim_x*y; 
			int ymo = x + dim_x*((y-1+dim_y)%dim_y);
			int ypo = x + dim_x*((y+1)%dim_y);
			int xpoypo = (x+1)%dim_x + dim_x*((y+1)%dim_y);
			for(int j = 0; j < 7; j++)
			{
				node_b[xpoypo].w.w[j] = 0.25*(node_a[i].w.w[j] + node_a[xpo].w.w[j] + node_a[ypo].w.w[j] + node_a[xpoypo].w.w[j]);
				//temp vars for summation
				double w_jk   = min_mod(node_a[xpo].w.w[j], node_a[i].w.w[j], node_a[xmo].w.w[j]);
				double w_j1k  = min_mod(node_a[(x+2)%dim_x + dim_x*y].w.w[j], node_a[xpo].w.w[j], node_a[i].w.w[j]);
				double w_jk1  = min_mod(node_a[xpoypo].w.w[j], node_a[ypo].w.w[j], node_a[(dim_x+x-1)%dim_x + (dim_x*((y+1)%dim_y))].w.w[j]);
				double w_j1k1 = min_mod(node_a[(x+2)%dim_x + dim_x*((y+1)%dim_y)].w.w[j], node_a[xpoypo].w.w[j], node_a[ypo].w.w[j]);

				double v_jk   = min_mod(node_a[ypo].w.w[j], node_a[i].w.w[j], node_a[ymo].w.w[j]);
				double v_j1k  = min_mod(node_a[xpoypo].w.w[j], node_a[xpo].w.w[j], node_a[(x+1)%dim_x + dim_x*((dim_y+y-1)%dim_y)].w.w[j]);
				double v_jk1  = min_mod(node_a[x + dim_x*((y+2)%dim_y)].w.w[j], node_a[ypo].w.w[j], node_a[i].w.w[j]);
				double v_j1k1 = min_mod(node_a[(x+1)%dim_x + dim_x*((y+2)%dim_y)].w.w[j], node_a[xpoypo].w.w[j], node_a[xpo].w.w[j]);

				double f_jk   = node_b[i].f.f[j];
				double f_j1k  = node_b[xpo].f.f[j];
				double f_jk1  = node_b[ypo].f.f[j];
				double f_j1k1 = node_b[xpoypo].f.f[j];

				double g_jk   = node_b[i].g.f[j];
				double g_j1k  = node_b[xpo].g.f[j];
				double g_jk1  = node_b[ypo].g.f[j];
				double g_j1k1 = node_b[xpoypo].g.f[j];

				//second part of summation in equ (2.16)
				node_b[xpoypo].w.w[j] += 0.0625*(w_jk - w_j1k + w_jk1 - w_j1k1) - 0.5*mu*(f_j1k - f_jk + f_j1k1 - f_jk1);
				node_b[xpoypo].w.w[j] += 0.0625*(v_jk - v_jk1 + v_j1k - v_j1k1) - 0.5*nu*(g_jk1 - g_jk + g_j1k1 - g_j1k);
			}
		}
	}
}

void wNSqrt(Node * field, int * dims, double *dt, double * delta)
{
	//Just do an in-place update as we don't save the fields from this phase, we just build tasks
	//Compute w_n^(1/2), which is a half timestep
	double mu = dt[0] / delta[0];
	double nu = dt[1] / delta[1];
	for(unsigned int y = 0; y < dims[1]; y++)
	{
		for(unsigned int x = 0; x < dims[0]; x++)
		{
			//Easy index vars
			int xmo = ( (x+dims[0]-1)%dims[0]) + dims[0]*y;
			int xpo = ( (x+1)%dims[0]) + dims[0]*y;
			int ymo =  x + dims[0]* ( (y+dims[1]-1)%dims[1] );
			int ypo =  x + dims[0]* ( (y+1)%dims[1] );
			int i = x + dims[0]*y;
			for(unsigned int j = 0; j < 7; j++)
			{
				double df = min_mod(field[xpo].f.f[j], field[i].f.f[j], field[xmo].f.f[j]);
				double dg = min_mod(field[ypo].g.f[j], field[i].g.f[j], field[ymo].g.f[j]);
				field[i].w.w[j] = field[i].w.w[j] - 0.5 * mu * df - 0.5 * nu * dg;
			}
		}
	}
}

/** set intial values for strain, momentum- and energy density on grid nodes
 * @param node       grid_node containing the conserved and the fluxes (output)
 * @param grid_size  number of grid nodes (input)
 * **/
void init_conserved_fields(Node* node_a, int * dims, int grid_size)
{
	int dimX = dims[0];
	int dimY = dims[1];
	for(int y = 0; y < dimY; y++)
	{
		for(int x = 0; x < dimX; x++)
		{
			int i = x + y*dimX;
			node_a[i].w.w[0] = 1.0;
			node_a[i].w.w[1] = 0.0;
			node_a[i].w.w[2] = 0.0;
			node_a[i].w.w[3] = 1.0;
			node_a[i].w.w[4] = 0.0;
			node_a[i].w.w[5] = 0.0;
			node_a[i].w.w[6] = -0.295;
			if( ( x < (dimX/2 + dimX/10) )  && ( x >= (dimX/2 - dimX/10) ) )
			{
				node_a[i].w.w[0] = 1.04;
				node_a[i].w.w[6] = -0.295;
			}
		}
	}
}

/** print the output in vtk format
 * @param i         number of output file (input)         
 * @param node      grid_node containing the conserved and the fluxes (input)
 * @param l         lattice parameters (input)         
 * @param grid_size  number of grid nodes (input)
 * **/
void printf_fields_vtk(int i, Node* node_a, Lattice l, int grid_size)
{
	int bufferSize = 256;
	char file_name[bufferSize];
	sprintf(file_name,"strain_%i.vtk",i);
	FILE *fn = fopen(file_name, "w");
	if (fn == NULL)
	{
		printf("Error writing file< %s >!\n", file_name);
	}
	char file_name2[bufferSize];
	sprintf(file_name2,"mom_%i.vtk",i);
	FILE *fn2 = fopen(file_name2, "w");
	if (fn2 == NULL)
	{
		printf("Error writing file< %s >!\n", file_name2);
	}
	char file_name3[bufferSize];
	sprintf(file_name3,"energy_%i.vtk",i);
	FILE *fn3 = fopen(file_name3, "w");
	if (fn3 == NULL)
	{
		printf("Error writing file< %s >!\n", file_name3);
	}
	char file_name4[bufferSize];
	sprintf(file_name4,"mom_flux_%i.vtk",i);
	FILE *fn4 = fopen(file_name4, "w");
	if (fn4 == NULL)
	{
		printf("Error writing file< %s >!\n", file_name4);
	}
	char file_name5[bufferSize];
	sprintf(file_name5,"flux_%i.vtk",i);
	FILE *fn5 = fopen(file_name5, "w");
	if (fn5 == NULL)
	{
		printf("Error writing file< %s >!\n", file_name5);
	}
	char file_name6[bufferSize];
	sprintf(file_name6,"energy_flux_%i.vtk",i);
	FILE *fn6 = fopen(file_name6, "w");
	if (fn6 == NULL)
	{
		printf("Error writing file< %s >!\n", file_name6);
	}

	fprintf(fn, "# vtk DataFile Version 2.0\nstrain\nASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS %i %i 1\nORIGIN 0.0 0.0 0.0\nSPACING %lf %lf 1.0\nPOINT_DATA %i\nSCALARS OutArray floats 4\nLOOKUP_TABLE default\n", l.dim_x, l.dim_y, l.dx, l.dy, grid_size); 
	fprintf(fn2, "# vtk DataFile Version 2.0\nmom\nASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS %i %i 1\nORIGIN 0.0 0.0 0.0\nSPACING %lf %lf 1.0\nPOINT_DATA %i\nSCALARS OutArray floats 2\nLOOKUP_TABLE default\n", l.dim_x, l.dim_y, l.dx, l.dy, grid_size); 
	fprintf(fn3, "# vtk DataFile Version 2.0\nenergy\nASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS %i %i 1\nORIGIN 0.0 0.0 0.0\nSPACING %lf %lf 1.0\nPOINT_DATA %i\nSCALARS OutArray floats 1\nLOOKUP_TABLE default\n", l.dim_x, l.dim_y, l.dx, l.dy, grid_size); 
	fprintf(fn4, "# vtk DataFile Version 2.0\nmom_flux\nASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS %i %i 1\nORIGIN 0.5 0.5 0.5\nSPACING %lf %lf 1.0\nPOINT_DATA %i\nSCALARS OutArray floats 2\nLOOKUP_TABLE default\n", l.dim_x, l.dim_y, l.dx, l.dy, grid_size); 
	fprintf(fn5, "# vtk DataFile Version 2.0\nflux\nASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS %i %i 1\nORIGIN 0.5 0.5 0.5\nSPACING %lf %lf 1.0\nPOINT_DATA %i\nSCALARS OutArray floats 4\nLOOKUP_TABLE default\n", l.dim_x, l.dim_y, l.dx, l.dy, grid_size); 
	fprintf(fn6, "# vtk DataFile Version 2.0\nenergy_flux\nASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS %i %i 1\nORIGIN 0.5 0.5 0.5\nSPACING %lf %lf 1.0\nPOINT_DATA %i\nSCALARS OutArray floats 2\nLOOKUP_TABLE default\n", l.dim_x, l.dim_y, l.dx, l.dy, grid_size); 
	for(int index=0; index<grid_size; ++index)
	{
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

void shift_back(Node* node_b, int grid_size, int * dims, Node* node_a)
{
	int x = 0;
	int y = 0;
	int xpoypo;
	for(int i=0; i<grid_size; ++i)
	{
		int x = i % dims[0];
		int y = i / dims[0];
		// temp vars for indexing
		xpoypo = (x+1)%dims[0] + dims[0]*((y+1)%dims[1]);
		memcpy(&node_b[i].f,& node_a[xpoypo].f, sizeof(Fluxes));
		memcpy(&node_b[i].g,& node_a[xpoypo].g, sizeof(Fluxes));
		memcpy(&node_b[i].w,& node_a[xpoypo].w, sizeof(Conserved));
	}
}

