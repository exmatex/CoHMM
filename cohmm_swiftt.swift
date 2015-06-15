(boolean done) initWrapper (boolean doKriging, boolean doCoMD, blob bDims, blob bDt, blob bDelta, blob bGamma) "cohmm_swiftt" "1.0"
[
	"set <<done>> [ cohmm_swiftt::initWrapper <<doKriging>> <<doCoMD>> <<bDims>> <<bDt>> <<bDelta>> <<bGamma>> ]"
];

(int nTasks) prepFirstWrapper (boolean doKriging, boolean doCoMD, blob bDims, blob bDt, blob bDelta, blob bGamma, int step) "cohmm_swiftt" "1.0"
[
	"set <<nTasks>> [ cohmm_swiftt::prepFirstWrapper <<doKriging>> <<doCoMD>> <<bDims>> <<bDt>> <<bDelta>> <<bGamma>> <<step>> ]"
];

(int nTasks) prepSecondWrapper (boolean doKriging, boolean doCoMD, blob bDims, blob bDt, blob bDelta, blob bGamma, int step, boolean dummy[]) "cohmm_swiftt" "1.0"
[
	"set <<nTasks>> [ cohmm_swiftt::prepSecondWrapper <<doKriging>> <<doCoMD>> <<bDims>> <<bDt>> <<bDelta>> <<bGamma>> <<step>> ]"
];

(int nTasks) prepThirdWrapper (boolean doKriging, boolean doCoMD, blob bDims, blob bDt, blob bDelta, blob bGamma, int step, boolean dummy[]) "cohmm_swiftt" "1.0"
[
	"set <<nTasks>> [ cohmm_swiftt::prepThirdWrapper <<doKriging>> <<doCoMD>> <<bDims>> <<bDt>> <<bDelta>> <<bGamma>> <<step>> ]"
];

(int nTasks) prepLastWrapper (boolean doKriging, boolean doCoMD, blob bDims, blob bDt, blob bDelta, blob bGamma, int step, boolean dummy[]) "cohmm_swiftt" "1.0"
[
	"set <<nTasks>> [ cohmm_swiftt::prepLastWrapper <<doKriging>> <<doCoMD>> <<bDims>> <<bDt>> <<bDelta>> <<bGamma>> <<step>> ]"
];

(int nTasks) finishStepWrapper (boolean doKriging, boolean doCoMD, blob bDims, blob bDt, blob bDelta, blob bGamma, int step, boolean dummy[]) "cohmm_swiftt" "1.0"
[
	"set <<nTasks>> [ cohmm_swiftt::finishStepWrapper <<doKriging>> <<doCoMD>> <<bDims>> <<bDt>> <<bDelta>> <<bGamma>> <<step>> ]"
];

(boolean done) prepVTK (boolean doKriging, boolean doCoMD, blob bDims, blob bDt, blob bDelta, blob bGamma, int step) "cohmm_swiftt" "1.0"
[
	"set <<done>> [ cohmm_swiftt::prepVTK <<doKriging>> <<doCoMD>> <<bDims>> <<bDt>> <<bDelta>> <<bGamma>> <<step>> ]"
];

(boolean done) cloudFluxWrapper (boolean doKriging, boolean doCoMD, int step, int phase, int task) "cohmm_swiftt" "1.0"
[
	"set <<done>> [ cohmm_swiftt::cloudFluxWrapper <<doKriging>> <<doCoMD>> <<step>> <<phase>> <<task>> ]"
];


