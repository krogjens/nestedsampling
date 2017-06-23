#define _USE_MATH_DEFINES
#include <math.h>
#include "mex.h"

void SDE_logl(double *obs, double *theta,
      mwSize num_obs,
      double *logprob)
{

   mwSize obs_idx;

   *logprob = 0;

   double mob;
   double sgn;
   double var;
   double var_d;
   double mean_x;
   double F; 
   
   mob = theta[0] * pow(fabs(obs[0]),theta[1]);
   sgn = (obs[0] > 0) - (obs[0] < 0);
   mean_x = theta[2] * mob + sgn * theta[1] * theta[0] * pow(fabs(obs[0]),theta[1] - 1.0);
   F = obs[1] - obs[0] - mean_x;
   var_d = 2.0 * mob;
   var = var_d + 2.0 * theta[3];
   *logprob += - pow(F,2.0) / (2 * var) - log(sqrt(2 * M_PI * var));

   for(obs_idx = 2; obs_idx < num_obs; obs_idx++){
      mob = theta[0] * pow(fabs(obs[obs_idx - 1]),theta[1]);
      var_d = 2.0 * mob;
      sgn = (obs[obs_idx - 1] > 0) - (obs[obs_idx - 1] < 0);
      mean_x = theta[2] * mob + sgn * theta[1] * theta[0] * pow(fabs(obs[obs_idx - 1]),theta[1] - 1.0);
      F = obs[obs_idx] - obs[obs_idx - 1] - mean_x + theta[3] / var * F;
      var = var_d + 2.0 * theta[3] - pow(theta[3],2) / var;
      *logprob += - pow(F,2) / (2.0 * var) - log(sqrt(2.0 * M_PI * var));
   }
}
