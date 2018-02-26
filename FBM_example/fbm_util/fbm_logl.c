#define _USE_MATH_DEFINES
#include <math.h>
#include "mex.h"

/* Likelihood calculation routine */

void fbm_logl(double *obs, double *theta, double *M,
	mwSize dim, mwSize num_obs,
	double *logprob)
{

   mwSize obs_idx;
   mwSize dim_idx;


   *logprob = 0; 
   
   if (M[3] == 0) {       	
	/*Calculation with conditional variance from noise for pure Brownian motion*/

	double obs_sum = 0;
	
	double var = 0;
	double var_prev = 0;
	double var_0; 
	double* F = (double*)mxCalloc(dim,sizeof(double));
	double* F_prev = (double*)mxCalloc(dim,sizeof(double));

	for(dim_idx = 0; dim_idx<dim; dim_idx++){
		F_prev[dim_idx] = obs[dim_idx] - theta[dim_idx + 1];
		obs_sum += pow(F_prev[dim_idx],2);
	}
	var_0 = pow(theta[0],2) + 2 * pow(theta[3],2);
	var_prev = var_0;
	*logprob += - obs_sum / (2 * var_0) - log(pow(2*M_PI,(double)dim/2) * pow(var_0,(double)dim/2));

	for (obs_idx = 1; obs_idx<num_obs; obs_idx++) {

		double obs_sum = 0;

		var = var_0 - pow(theta[3],4) / var_prev;
		for (dim_idx = 0; dim_idx < dim; dim_idx++) {
			F[dim_idx] = obs[obs_idx*dim + dim_idx] - theta[1 + dim_idx] + pow(theta[3],2) / var_prev * F_prev[dim_idx];
			obs_sum += pow(F[dim_idx],2);
			var_prev = var;
			F_prev[dim_idx]=F[dim_idx];
		}
		
		*logprob += -obs_sum / (2 * var) - log(pow(2*M_PI,(double)dim/2) * pow(var,(double)dim/2));
	}
	mxFree(F_prev);
	mxFree(F);
   }
   else { 

	mwSize obs_idx_1;
	
	double* gamma = (double*)mxCalloc(num_obs+1,sizeof(double));
	double* phi_last = (double*)mxCalloc(num_obs,sizeof(double));
	double* phi_next = (double*)mxCalloc(num_obs,sizeof(double));
	double* mu_cond = (double*)mxCalloc(dim*num_obs,sizeof(double)); 

	double var = 0;
	double var_prev = 0;
	double gamma_sum = 0;
	double* phi_sum = (double*)mxCalloc(dim,sizeof(double));
	double obs_sum = 0;
	
	gamma[1] = 0.5 * pow(theta[0],2) * (pow(2,2*theta[4]) - 2 * pow(1,2*theta[4]) + pow(0,2*theta[4])) - pow(theta[3],2);
	for (obs_idx = 2; obs_idx<num_obs; obs_idx++) {
		gamma[obs_idx] = 0.5 * pow(theta[0],2) * (pow(obs_idx + 1,2*theta[4]) - 2 * pow(obs_idx,2*theta[4]) + pow(obs_idx - 1,2*theta[4]));
	}
	var_prev = pow(theta[0],2) + 2 * pow(theta[3],2);
	phi_last[0]=gamma[1]/var_prev;
	for(dim_idx = 0; dim_idx<dim; dim_idx++){
		obs_sum += pow(obs[dim_idx] - theta[dim_idx + 1],2);
	}

	*logprob += -obs_sum / (2 * var_prev) - log(pow(2*M_PI,(double)dim/2) * pow(var_prev,(double)dim/2));

	for(obs_idx = 1; obs_idx<num_obs; obs_idx++){
		
		obs_sum = 0; 
		gamma_sum = 0;
	
		var = var_prev * (1 - pow(phi_last[obs_idx - 1],2));
		var_prev = var;		
		
		for(dim_idx = 0; dim_idx<dim; dim_idx++){
			phi_sum[dim_idx] = 0;
			for(obs_idx_1 = 0; obs_idx_1<obs_idx; obs_idx_1++){
				phi_sum[dim_idx] += phi_last[obs_idx_1] * (obs[dim * (obs_idx - obs_idx_1 - 1) + dim_idx] - theta[dim_idx + 1]);
			}
			mu_cond[dim * obs_idx + dim_idx] = theta[1 + dim_idx] + phi_sum[dim_idx]; 
			obs_sum += pow(obs[dim * obs_idx + dim_idx]-mu_cond[dim * obs_idx + dim_idx],2);
		}
		
		*logprob += - obs_sum / (2 * var) - log(pow(2*M_PI,(double)dim/2) * pow(var,(double)dim/2));

		for(obs_idx_1 = 0; obs_idx_1<obs_idx; obs_idx_1++){
			gamma_sum += phi_last[obs_idx_1] * gamma[obs_idx - obs_idx_1];
		}
		phi_next[obs_idx] = (gamma[obs_idx+1] - gamma_sum) / var;

		for(obs_idx_1 = 0; obs_idx_1<obs_idx; obs_idx_1++){
			phi_next[obs_idx_1] = phi_last[obs_idx_1] - phi_last[obs_idx - obs_idx_1 - 1] * phi_next[obs_idx];
		}
	
		for(obs_idx_1 = 0; obs_idx_1 < obs_idx + 1 ; obs_idx_1++){
			phi_last[obs_idx_1] = phi_next[obs_idx_1];
		}
	}
   }
}
