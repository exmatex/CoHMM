#ifndef COMD_LIB
#define COMD_LIB

typedef struct CoMD_return_t
{
	int simIter;
	double stressXX;
	double stressXY;
	double stressXZ;
	double stressYX;
	double stressYY;
	double stressYZ;
	double stressZX;
	double stressZY;
	double stressZZ;
	double simTime;
	double eTotal;
	double eK;
	double eU;
	double temp;
	double energyDensX;
	double energyDensY;
	double energyDensZ;
	double momentumDensX;
	double momentumDensY;
	int numAtoms;
} CoMD_return;

typedef struct CoMD_input_t
{
	char potDir[1024];  //!< the directory where EAM potentials reside
	char potName[1024]; //!< the name of the potential
	char potType[1024]; //!< the type of the potential (funcfl or setfl)
	int doeam;          //!< a flag to determine whether we're running EAM potentials
	int nx;             //!< number of unit cells in x
	int ny;             //!< number of unit cells in y
	int nz;             //!< number of unit cells in z
	int nSteps;         //!< number of time steps to run
	int printRate;      //!< number of steps between output
	double dt;          //!< time step (in femtoseconds)
	double lat;         //!< lattice constant (in Angstroms)
	double temperature; //!< simulation initiag temperature (in Kelvin)
	double initialDelta;//!< magnitude of initial displacement from lattice (in Angstroms)
	double defGrad[4];   //!< deformation gradient == 1,1 component of the strain tensor
	double enDens;
	double momDens[3];
    int rank;
	int calls;
} CoMD_input;







CoMD_return CoMD_lib(CoMD_input *inputStruct);

#endif

