#include "mex.h"
#include "matrix.h"

#include "linequad.h"


void bclag_interp_weights(const double* restrict x, double* restrict w, int n)
{
  for (int j=0; j<n; j++)
    {
      double wj = 1;
      for (int i=0; i<n; i++)
	{
	  if (i != j)
	    wj = wj*(x[j]-x[i]);
	}
      w[j] = 1/wj;
    }
}

void bclag_interp_matrix(const double* restrict x, const double* restrict xi,
			 double* restrict B, double* restrict w, double* restrict denom,
			 int* restrict exact,
			 int n, int N)
{
  // Setup
  for (int i=0; i<N; i++) exact[i] = -1;
  // Compute weights
  bclag_interp_weights(x, w, n);
  // Setup matrix elements
  for (int i=0; i<N; i++)  
    {
      for (int j=0; j<n; j++)
	{
	  double xdiff = xi[i]-x[j];
	  double temp = w[j]/xdiff;
	  B[i + N*j] = temp;
	  denom[i] = denom[i] + temp;
	  if (xdiff==0)
	    exact[i] = j;
	}
    }
  for (int j=0; j<n; j++)
    {
      for (int i=0; i<N; i++)
	{
	  B[i + N*j] = B[i + N*j] / denom[i];
	}
    }
  // Fix exact points
  for (int i=0; i<N; i++)
    {
      if (exact[i] >= 0)
	{
	  for (int j=0; j<n; j++)
            B[i + j*N] = 0;
	  B[i + N*exact[i]] = 1;
        }
    }  
}

/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  // Read input
  double* x = mxGetPr(prhs[0]);  
  int n = mxGetNumberOfElements(prhs[0]);
  double* xi = mxGetPr(prhs[1]);  
  int N = mxGetNumberOfElements(prhs[1]);   
  // Prepare output
  plhs[0] = mxCreateDoubleMatrix(N, n, mxREAL);
  double* B = mxGetPr(plhs[0]);
  // Algorithm
  double* w = mxCalloc(n,sizeof(double));
  double* denom = mxCalloc(N,sizeof(double));
  int* exact = mxCalloc(N,sizeof(int));
  bclag_interp_matrix(x, xi, B, w, denom, exact, n, N);
    
  mxFree(exact);
  mxFree(denom);
  mxFree(w);
}
