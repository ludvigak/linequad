#include "math.h"
#include "mex.h"
#include "matrix.h"

#include "linequad.h"

/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  // Read input
  int n = mxGetNumberOfElements(prhs[0]);
  double* tj = mxGetPr(prhs[0]);  
  double* xj = mxGetPr(prhs[1]);
  double* yj = mxGetPr(prhs[2]);
  double* zj = mxGetPr(prhs[3]);
  int N = mxGetNumberOfElements(prhs[4]);    
  double* x0 = mxGetPr(prhs[4]);
  double* y0 = mxGetPr(prhs[5]);
  double* z0 = mxGetPr(prhs[6]);
  // Prepare output
  plhs[0] = mxCreateDoubleMatrix(N, 1, mxCOMPLEX);
  double* tinit_re_p = mxGetPr(plhs[0]);    
  double* tinit_im_p = mxGetPi(plhs[0]);
  // Call function for all inputs
#pragma omp parallel for schedule(guided)
  for (int i=0; i<N; i++)
    {
      double complex tinit = rootfinder_initial_guess(tj, xj, yj, zj, n, x0[i], y0[i], z0[i]);
      tinit_re_p[i] = creal(tinit);
      tinit_im_p[i] = cimag(tinit);
    }
}
