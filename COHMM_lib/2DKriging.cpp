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

#include <iostream>

#include "2DKriging.hpp"
#include "kriging.hpp"
#include "redisBuckets.hpp"

extern "C"
{
#include "CoMD_lib.h"
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
bool checkGradient(int x, int y, std::vector<Node> & fields, unsigned int dims[2], double deltas[2])
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

void wSummation(std::vector<Node> & node_a, std::vector<Node> & node_b, unsigned int * dims, double * dt, double * delta)
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

void wNSqrt(std::vector<Node> & inField, std::vector<Node> & outField, unsigned int * dims, double *dt, double * delta)
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
				double df = min_mod(inField[xpo].f.f[j], inField[i].f.f[j], inField[xmo].f.f[j]);
				double dg = min_mod(inField[ypo].g.f[j], inField[i].g.f[j], inField[ymo].g.f[j]);
				outField[i].w.w[j] = inField[i].w.w[j] - 0.5 * mu * df - 0.5 * nu * dg;
			}
		}
	}
}

void init_conserved_fields(std::vector<Node> & node_a, unsigned int * dims, int grid_size)
{
	return init_conserved_fields(node_a, dims, grid_size, InitialConditions_e::X);
}


/** set intial values for strain, momentum- and energy density on grid nodes
 * @param node       grid_node containing the conserved and the fluxes (output)
 * @param grid_size  number of grid nodes (input)
 * **/
void init_conserved_fields(std::vector<Node> & node_a, unsigned int * dims, int grid_size, InitialConditions_e prob_type)
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
			//Initial stimulus
			switch(prob_type)
			{
				case InitialConditions_e::LINE:
					if( ( x < (dimX/2 + dimX/10) )  and ( x >= (dimX/2 - dimX/10) ) )
					{
						node_a[i].w.w[0] = 1.04;
						node_a[i].w.w[6] = -0.295;
					}
				break;
				case InitialConditions_e::CENTRALIZED:
					if(
						( x < (dimX/2 + dimX/10) )  and ( x >= (dimX/2 - dimX/10) )
							and
						( y < (dimY/2 + dimY/10) )  and ( y >= (dimY/2 - dimY/10) )
					)
					{
						node_a[i].w.w[0] = 1.04;
						node_a[i].w.w[6] = -0.295;
					}
				break;
				case InitialConditions_e::X:
					//Just doing a plus as they are close enough
					if(
						( x < (dimX/2 + dimX/10) )  and ( x >= (dimX/2 - dimX/10) )
							or
						( y < (dimY/2 + dimY/10) )  and ( y >= (dimY/2 - dimY/10) )
					)
					{
						node_a[i].w.w[0] = 1.04;
						node_a[i].w.w[6] = -0.295;
					}
				break;
				default:
					//No initial conditions, so no stimulus
				break;
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
void printf_fields_vtk(int i, std::vector<Node> & node_a, Lattice l, int grid_size)
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

void inPlace_shift_back(Node * arr, unsigned int * dims)
{
	//Save the top row
	Node topRow[dims[0]];
	memcpy(topRow, &arr[0 + dims[0] * (dims[1] - 1)], sizeof(Node) * dims[0]);
	//Do all but the top row
	for(int j = 0; j < dims[1] - 1; j++)
	{
		int yM = (j + dims[0] - 1) % dims[1];
		for(int i = 0; i < dims[0]; i++)
		{
			int xM = (i + dims[0] - 1) % dims[0];
			memcpy(&arr[xM + yM * dims[0]], &arr[i + j * dims[0]], sizeof(Node));
		}
	}
	//Do the top row
	{
		int j = dims[1] - 1;
		int yM = (j + dims[0] - 1) % dims[1];
		for(int i = 0; i < dims[0]; i++)
		{
			int xM = (i + dims[0] - 1) % dims[0];
			memcpy(&arr[xM + yM * dims[0]], &topRow[i], sizeof(Node));
		}
	}
	//Done?
}

void shift_back(Node* node_b, int grid_size, unsigned int * dims, Node* node_a)
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

FluxOut fluxFn(bool doKriging, bool doCoMD, FluxIn * input)
{
	FluxOut output;

	//See if we are kriging or solving
	bool doSolver = (not doKriging) or (not input->tryKriging);
	bool tryKriging = not doSolver;
	//Do we bother kriging?
	if(tryKriging == true)
	{
		//We do
		//Get data for kriging
		std::vector<double *> oldWs;
		std::vector<double *> oldFs;
		std::vector<double *> oldGs;
		getSortedSubBucketNearZero(input->fields.w, (char *)"comd", comdDigits, 10, &oldWs, &oldFs, &oldGs, zeroThresh);
		//Call Kriging on each point
		double resF[2];
		double resG[2];
		double error;
		for(int i = 0; i < 7; i++)
		{
			int info;
			///TODO: Consider refactoring this to kriging.cpp
			std::vector<double> oldF;
			std::vector<double> oldG;
			for(int j = 0; j < oldFs.size(); j++)
			{
				oldF.push_back(oldFs[j][i]);
				oldG.push_back(oldGs[j][i]);
			}
			info = kriging(input->fields.w, 7, oldWs, oldF, resF);
			info = kriging(input->fields.w, 7, oldWs, oldG, resG);
			//Write result
			output.f[i] = resF[0];
			output.g[i] = resG[0];
			//Set error
			if(resF[1] > error)
			{
				error = resF[1];
			}
			if(resG[1] > error)
			{
				error = resG[1];
			}
		}
		freeClear(oldWs);
		freeClear(oldFs);
		freeClear(oldGs);
		//Was error too high?
		if(error > errorThresh)
		{
			//It was, so fall back to comd
			doSolver = true;
		}
		else
		{
			//It was not, so put it to the kriging db
			putData(input->fields.w, output.f, output.g, (char *)"krig", krigDigits);
		}
	}
	//Either kriging failed or we never tried
	if(doSolver == true)
	{
		//Call CoMD
		//4 strains
		double strain_xx = input->fields.w[0];
		double strain_xy = input->fields.w[1];
		double strain_yx = input->fields.w[2];
		double strain_yy = input->fields.w[3];
		//2 momentum_fluxes
		double momentum_x = input->fields.w[4];
		double momentum_y = input->fields.w[5];
		//enery_flux
		double energy = input->fields.w[6];
		//CoMD or approximation?
		if(doCoMD == true)
		{
			CoMD_input theInput;
			///FIXME: Fix potdir
			strcpy(theInput.potDir,"../pots");
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
			output.f[0] = momentum_x/rho;
			output.f[1] = momentum_y/rho;
			output.f[2] = 0.0;
			output.f[3] = 0.0;
			output.f[4] = theRet.stressXX;
			output.f[5] = theRet.stressXY;
			output.g[0] = 0.0;
			output.g[1] = 0.0;
			output.g[2] = momentum_x/rho;
			output.g[3] = momentum_y/rho;
			output.g[4] = theRet.stressYX;
			output.g[5] = theRet.stressYY;
			output.f[6] = -theRet.energyDensX;
			output.g[6] = -theRet.energyDensY;
		}
		else
		{
			//4 stresses
			double rho = 1.0;
			//out.f[0] = theRet.momX;
			output.f[0] = momentum_x/rho;
			output.f[1] = momentum_y/rho;
			output.f[2] = 0.0;
			//o->t.f[2] = theRet.momY;
			output.f[3] = 0.0;
			output.f[4] = strain_xx-1 + 0.75*(strain_yy-1);
			output.f[5] = 1.9*strain_xy;
			output.g[0] = 0.0;
			//o->t.g[1] = theRet.momX;
			output.g[1] = 0.0;
			output.g[2] = momentum_x/rho;
			//o->t.g[3] = theRet.momY;
			output.g[3] = momentum_y/rho;
			output.g[4] = 1.9*strain_yx;
			output.g[5] = strain_yy-1 + 0.75*(strain_xx-1);
			//f->6] = g[6] = -0.295;
			output.f[6] = -output.f[0]*sqrt(output.f[4]*output.f[4] + output.f[5]*output.f[5]);
			output.g[6] = -output.g[3]*sqrt(output.g[4]*output.g[4] + output.g[5]*output.g[5]);
		}
		//Put result to DB for future use: Warning, flush if we switch to comd as this is horrible
		putData(input->fields.w, output.f, output.g, (char *)"comd", comdDigits);
	}
	return output;
}
