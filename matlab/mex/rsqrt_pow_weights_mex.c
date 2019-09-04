#include "mex.h"
#include "matrix.h"

#include "linequad.h"

/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    // Read input
    double* tj = mxGetPr(prhs[0]);  
    int n = mxGetNumberOfElements(prhs[0]);
    double* troot_re_p = mxGetPr(prhs[1]);    
    double* troot_im_p = mxGetPi(prhs[1]);
    int N = mxGetNumberOfElements(prhs[1]);   
    // Prepare output
    plhs[0] = mxCreateDoubleMatrix(n, N, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(n, N, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(n, N, mxREAL);
    double* w1_all = mxGetPr(plhs[0]);
    double* w3_all = mxGetPr(plhs[1]);
    double* w5_all = mxGetPr(plhs[2]);
    // Call function
    int pmax;
    if (nlhs==1)
      pmax = 1;
    if (nlhs==2)
      pmax = 3;
    if (nlhs==3)
      pmax = 5;

#pragma omp parallel for schedule(guided)
    for (int i=0; i<N; i++)
      {
	double complex troot = troot_re_p[i] + I*troot_im_p[i];
	double* w1 = w1_all + i*n;
	double* w3 = w3_all + i*n;
	double* w5 = w5_all + i*n;	
	rsqrt_pow_weights(tj, troot, n, pmax, w1, w3, w5);
      }
}
