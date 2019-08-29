/*
* function to add two logarithms*/


#define _USE_MATH_DEFINES
#include <math.h>
#include "mex.h"

	
	double bik(double B, double *x, double *y, double ux,double uy, int t)
	{
		double b;
			b=((x[t]-ux)*(ux-x[t])+(y[t]-uy)*(uy-y[t]))/(B*B*2.0);
			return b;
	}
	
	
		double dbik(double B, double x, double y, double *u,int J,int j)
	{	
			double b=((x-u[j*2])*(u[j*2]-x)+(y-u[j*2+1])*(u[j*2+1]-y))/(B*B*2.0)-2.0*log(B)-log(2.0*M_PI);
			return b;
	}
  
 double logplus(double x, double y)
 {
	 double z;
    if ( x>y ) {
        /* Execute these statements if TRUE */
        z = x+log(1.0+exp(y-x));
    }
    else {
        /* Execute these statements if FALSE */
        z = y+log(1.0+exp(x-y));
    }
	return z;
 }
 
 void forward(double *p,double *logA,double *B, double *x, double *y, double *u, int T,int J, double *logP)
	{
			
	if (J==1){	
		logP[0]=0;
		for (int t =0; t<T;t++){
			logP[0]=logP[0]+bik(B[0],x,y,u[0],u[1],t);
		}
		logP[0]=logP[0]-T*(2.0*log(B[0])+log(2.0*M_PI)); 
		}
	else{
		 
		 double loga[J];
		 for (int j=0;j<J;j++){
			  loga[j]=log(p[j])+dbik(B[j],x[0],y[0],u,J,j); /*initial emission probabilities for each state*/
		 }
		
		 double gamma;
		 for(int t = 1; t<T;t++){
			 
			for(int j = 0; j<J;j++){
				gamma = logplus(loga[0]+logA[j],loga[1]+logA[j+J]); /*first state only I MAY HAVE TO CHANGE THE ORDER AROUND BECAUSE MATRIX IS NOT TRANSPOSED ANYMORE*/
				
				for (int i=2;i<J;i++){  /*rest of states*/
					gamma = logplus(gamma,loga[i]+logA[i*J+j]); 
				}
				loga[j]=gamma+dbik(B[j],x[t],y[t],u,J,j);
			}
		 }
		 logP[0]=logplus(loga[0],loga[1]);
		 for (int i=2;i<J;i++){
			 logP[0]=logplus(logP[0],loga[i]);

		 }

		 /*implementation to handle missing data*/
		/*-(constant) /*could be calculated and inputted*/
		/* should not loop over 2 but length of missing data*/
		/*will contain missing data*/
			/*mangler missing B parametre */
			/*
		 if (M>0){
			 if (uniprior==1){
				
				 logP[0]=logP[0];
			 }
			 else{		

				for (int m =0; m<M;m++){ 
					
					logP[0]=logP[0]-bik(missingB,xm,ym,0.0,0.0,m); 
				}		
			 }
		 }
		 */
		}
		
	}
 
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
                 
 /*ForwardAlgorithm(obs, p, logA,B, u nstates,ndrifts,T)*/

 {
	
	double *obs =mxGetPr(prhs[0]); /*observations */
	double *p =mxGetPr(prhs[1]);
	double *logA =mxGetPr(prhs[2]);
	double *B =mxGetPr(prhs[3]);
	double *u =mxGetPr(prhs[4]);
	
/*	double *theta =mxGetPr(prhs[1]);*/
	int nstates = mxGetScalar(prhs[5]);; /*nstates*/
	int ndrifts =mxGetScalar(prhs[6]);; /* ndrifts */
	int T =mxGetScalar(prhs[7]); /* number of observations */
/*	double *B; /* diffusion parameters */
/*	double *u; /* drift parameters */
/*	double *logA; /* log(Transition probability matrix) */
/*	double *logP; /* log(probability) */
/*	double *p; /*log(starting probabilities)  */
	
	double *x;
	double *y;
    double *outMatrix; 

	/* Now to assign the inputs to their proper variables */
	x=&obs[0]; /*obs should be stores column wise, hence we set the pointer to x at 0 */
	y=&obs[T]; /* and the pointer to y at T */
	
   /* B = &theta[0]; /*the first parameters are diffusion */
	
	/* u has  ndrift values + nstates-ndrift value that are 0*/

/*	double arr[10]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}; /*variable sized array was not allowed, but i make one large enough */
/*	for(int i = 0; i<ndrifts;i++){
		arr[i]=theta[nstates+i];
	}
	u = &arr[0];	
	
	/*logA = &theta[nstates+ndrifts]; /*calculate llog of this* loga=transpose(expm(SPar.dt*newAP))/
	
	p=&theta[nstates*2+ndrifts];
	
	
	/**	/*mxArray *M, *uniprior,*missB;
	M= mxGetFieldByNumber(prhs[8], 0,5); det er  antallet af missing steps
	int M;
	M= 0; 
	 noget til udregning hvis der er missing data
	int unip;
	uniprior=mxGetFieldByNumber(prhs[8], 0,0); boolean for uniprior
	unip=*mxGetPr(uniprior);

    mwSize N;
     mere til missing data
	double missingB;
	missB=mxGetFieldByNumber(prhs[8], 0,4);
	missingB=*mxGetPr(missB);
	*/
	
	
	
    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
    

    
    /* get a pointer to the real data in the output matrix */
    outMatrix = mxGetPr(plhs[0]);
    
    forward(p,logA,B,x,y,u,T,nstates,outMatrix);
 }



