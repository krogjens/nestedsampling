#include "mex.h"

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{

   double *obs;
   double *theta;

   double *outMatrix;

   mwSize num_obs;

   /*Declare inputs*/
   obs = mxGetPr(prhs[0]);
   theta = mxGetPr(prhs[1]);
   
   /*Fetch data set length*/
   num_obs = mxGetM(prhs[0]);
   if(num_obs == 1){
      num_obs = mxGetN(prhs[0]);
   }

   /*Create output matrix*/
   plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);

   /*Assign pointer to output*/
   outMatrix = mxGetPr(plhs[0]);

   /*run routine*/

   SDE_logl(obs,theta,
           num_obs, outMatrix);
}
