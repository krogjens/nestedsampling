#include "mex.h"

void mexFunction( int nlhs, mxArray *plhs[],
		  int nrhs, const mxArray *prhs[] )
{

	double *obs;
	double *theta;
	double *M;
	
	double *outMatrix;

	mwSize dim;
	mwSize num_obs;

	/*checks for input type and number */
	if(nrhs != 3) {
	mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","Inputs: obs [dim x #steps] theta [1x5] M [1x4]");
	}
	obs = mxGetPr(prhs[0]);
	theta = mxGetPr(prhs[1]);
	M = mxGetPr(prhs[2]);

	dim = mxGetM(prhs[0]);
	num_obs = mxGetN(prhs[0]);

	if( !mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1]) || !mxIsDouble(prhs[2]) ) {
	mexErrMsgIdAndTxt("MyToolBox:arrayProduct:notScalar","All inputs should be of type double.");
	}

	/* create output matrix */
	plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);

	/* assign pointer to output*/
	outMatrix = mxGetPr(plhs[0]);

	/* run function */

	fbm_logl(obs,theta,M,
		dim, num_obs, outMatrix);

}
