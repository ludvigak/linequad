#include "mex.h"
#include "matrix.h"

#include "linequad.h"

/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  // Read input
  int n = mxGetNumberOfElements(prhs[0]);
  double* xhat = mxGetPr(prhs[0]);
  double* yhat = mxGetPr(prhs[1]);
  double* zhat = mxGetPr(prhs[2]);
  int N = mxGetNumberOfElements(prhs[3]);    
  double* x0 = mxGetPr(prhs[3]);
  double* y0 = mxGetPr(prhs[4]);
  double* z0 = mxGetPr(prhs[5]);
  double* tinit_re_p = mxGetPr(prhs[6]);    
  double* tinit_im_p = mxGetPi(prhs[6]);    
  // Prepare output
  plhs[0] = mxCreateDoubleMatrix(N, 1, mxCOMPLEX);
  double* troot_re_p = mxGetPr(plhs[0]);    
  double* troot_im_p = mxGetPi(plhs[0]);        
  plhs[1] = mxCreateLogicalMatrix(N, 1);
  mxLogical* converged_p = mxGetLogicals(plhs[1]);
  // Call function for all inputs
#pragma omp parallel for schedule(guided)
  for (int i=0; i<N; i++)
    {
      double complex troot = 0;
      int converged = 0;
      double complex tinit = tinit_re_p[i] + I*tinit_im_p[i];
      converged = rootfinder(xhat, yhat, zhat, n, x0[i], y0[i], z0[i], tinit, &troot);
      troot_re_p[i] = creal(troot);
      troot_im_p[i] = cimag(troot);
      converged_p[i] = converged;
    }
}
