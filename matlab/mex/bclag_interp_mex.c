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

void bclag_interp_eval(const double* restrict x, const double* restrict xi,
		       const double* restrict f, double* restrict fi, const double* restrict w,
		       int n, int N)
{
  for (int i=0; i<N; i++)
    {
      double numer = 0;
      double denom = 0;
      int exact = -1;
      for (int j=0; j<n; j++)
	{
	  double xdiff = xi[i]-x[j];
	  if (xdiff==0)
	    {
	      exact = j;
	      break;
	    }
	  double tmp = w[j]/xdiff;
	  numer += f[j]*tmp;
	  denom += tmp;
	}
      fi[i] = (exact >= 0) ? f[exact] : numer/denom;
    }
}

/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  // Read input
  double* x = mxGetPr(prhs[0]);  
  int n = mxGetNumberOfElements(prhs[0]);
  double* f = mxGetPr(prhs[1]);    
  double* xi = mxGetPr(prhs[2]);  
  int N = mxGetNumberOfElements(prhs[2]);   
  // Prepare output
  plhs[0] = mxCreateDoubleMatrix(N, 1, mxREAL);
  double* fi = mxGetPr(plhs[0]);
  // Algorithm
  double* w = mxCalloc(n,sizeof(double));
  bclag_interp_weights(x, w, n);
  bclag_interp_eval(x, xi, f, fi, w, n, N);
      
  mxFree(w);
}
