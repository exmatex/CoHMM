#include "parseInputs.h"
#include <math.h>


#include <string.h>
#include <stdlib.h>

#include <stdio.h>

Command parseInputStruct(CoMD_input * inputStruct)
{

	Command cmd;
	
	memset(cmd.potDir, 0, 1024);
	memset(cmd.potName, 0, 1024);
	memset(cmd.potType, 0, 1024);
	if(inputStruct->potName[0] != 0)
	{
		strcpy(cmd.potName, inputStruct->potName);
	}
	else
	{
		strcpy(cmd.potName, "\0"); // default depends on potType
	}
	if(inputStruct->potDir[0] != 0)
	{
		strcpy(cmd.potDir,  inputStruct->potDir);
	}
	else
	{
		strcpy(cmd.potDir,  "pots");
	}
	if(inputStruct->potType[0] != 0)
	{
		strcpy(cmd.potType, inputStruct->potType);
	}
	else
	{
		strcpy(cmd.potType, "funcfl");
	}

	if(inputStruct->doeam != -1)
	{
		cmd.doeam = inputStruct->doeam;
	}
	else
	{
		cmd.doeam = 0;
	}
	
	
	if(inputStruct->nx != -1)
	{
		cmd.nx = inputStruct->nx;
	}
	else
	{
		cmd.nx = 20;
	}
	if(inputStruct->ny != -1)
	{
		cmd.ny = inputStruct->ny;
	}
	else
	{
		cmd.ny = 20;
	}
	if(inputStruct->nz != -1)
	{
		cmd.nz = inputStruct->nz;
	}
	else
	{
		cmd.nz = 20;
	}

	//We don't care about these, leave as default
	cmd.xproc = 1;
	cmd.yproc = 1;
	cmd.zproc = 1;

	if(inputStruct->nSteps != -1)
	{
		cmd.nSteps = inputStruct->nSteps;
	}
	else
	{
		cmd.nSteps = 100;
	}
	
	if(inputStruct->printRate != -1)
	{
		cmd.printRate = inputStruct->printRate;
	}
	else
	{
		cmd.printRate = 10;
	}

	//Need to specify dt

	cmd.dt = inputStruct->dt;
	cmd.lat = inputStruct->lat;
	cmd.energy = inputStruct->energy;
	cmd.temperature = inputStruct->temperature;
	cmd.initialDelta = inputStruct->initialDelta;
	cmd.defGrad = inputStruct->defGrad;
	

	
	memset(cmd.stressSuffix, 0, 1024);

	//If users didn't specify potName
	if (strlen(cmd.potName) == 0) 
	{		
		if (strcmp(cmd.potType, "setfl" ) == 0)
		{
			strcpy(cmd.potName, "Cu01.eam.alloy");
		}
		if (strcmp(cmd.potType, "funcfl") == 0)
		{
			strcpy(cmd.potName, "Cu_u6.eam");
		}
	}

	return cmd;


}



