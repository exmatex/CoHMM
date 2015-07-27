import blob;
import io;
import cohmm_swiftt;

main
{
	//Initialize program constraints and variables
	boolean doKriging = true;
	boolean doCoMD = false;
	int dims[];
	dims[0] = 50;
	dims[1] = 50;
	float dt[];
	dt[0] = 0.1;
	dt[1] = 0.1;
	float delta[];
	delta[0] = 1.0;
	delta[1] = 1.0;
	float gamma[];
	gamma[0] = 0;
	gamma[1] = gamma[0];
	gamma[2] = 0.1 * gamma[1];
	int numSteps = 10;

	//Build blobs of all arrays
	blob bDims = blob_from_ints(dims);
	blob bDt = blob_from_floats(dt);
	blob bDelta = blob_from_floats(delta);
	blob bGamma = blob_from_floats(gamma);

	//Initialize the data
	boolean rInit = initWrapper(doKriging, doCoMD, bDims, bDt, bDelta, bGamma);
	//Run through timesteps after waiting on initialization
	boolean doneRun;
	wait(rInit)
	{
		//WARNING: Need this loop to execute sequentially
		int j[];
		for(int t = 0; t < numSteps; t = j[t])
		{
			wait(t)
			{
				boolean doneFirst[];
				boolean doneSecond[];
				boolean doneThird[];
				boolean doneLast[];
				prepVTK(doKriging, doCoMD, bDims, bDt, bDelta, bGamma, t) =>
				printf("%d: Starting", t) =>
				int firstTasks = prepFirstWrapper(doKriging, doCoMD, bDims, bDt, bDelta, bGamma, t) =>
				printf("%d: First Tasks: %d", t, firstTasks) =>
				foreach i in [0:(firstTasks-1):1]
				{
					doneFirst[i] = cloudFluxWrapper(doKriging, doCoMD, t, 0, i);
				}
				int secondTasks = prepSecondWrapper(doKriging, doCoMD, bDims, bDt, bDelta, bGamma, t, doneFirst) =>
				printf("%d: Second Tasks: %d", t, secondTasks) =>
				foreach i in [0:(secondTasks-1):1]
				{
					doneSecond[i] = cloudFluxWrapper(doKriging, doCoMD, t, 1, i);
				}
				int thirdTasks = prepThirdWrapper(doKriging, doCoMD, bDims, bDt, bDelta, bGamma, t, doneSecond) =>
				printf("%d: Third Tasks: %d", t, thirdTasks) =>
				foreach i in [0:(thirdTasks-1):1]
				{
					doneThird[i] = cloudFluxWrapper(doKriging, doCoMD, t, 2, i);
				}
				int lastTasks = prepLastWrapper(doKriging, doCoMD, bDims, bDt, bDelta, bGamma, t, doneThird) =>
				printf("%d: Fourth Tasks: %d", t, lastTasks) =>
				foreach i in [0:(lastTasks-1):1]
				{
					doneLast[i] = cloudFluxWrapper(doKriging, doCoMD, t, 3, i);
				}
				int finalTasks = finishStepWrapper(doKriging, doCoMD, bDims, bDt, bDelta, bGamma, t, doneLast) =>
				printf("%d: Final Tasks: %d", t, finalTasks) =>
				if(finalTasks == 0)
				{
					j[t] = t+1;
					if(j[t] == numSteps)
					{
						doneRun = true;
					}
				}
			}
		}
	}
	//Do final vis of the data when all else is done
	wait(doneRun)
	{
		prepVTK(doKriging, doCoMD, bDims, bDt, bDelta, bGamma, numSteps);
	}
	//No clean-up because Swift/T
}

