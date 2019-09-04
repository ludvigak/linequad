static inline void pvand(int n, const double* restrict alpha, double* restrict x, const double* restrict b)
{
  // x = pvand(n, alpha, x, b)
  //
  // Solves system A*x = b
  // A is Vandermonde matrix, with nonstandard definition
  // A(i,j) = alpha(j)^i
  //
  // Algorithm by Bjorck & Pereyra
  // Mathematics of Computation, Vol. 24, No. 112 (1970), pp. 893-903
  // https://doi.org/10.2307/2004623    
  //
  for (int i=0; i<n; i++)
    x[i] = b[i];
  for (int k=0; k<n; k++)
    for (int j=n-1; j>k; j--)
      x[j] = x[j]-alpha[k]*x[j-1];      
  for (int k=n-1; k>0; k--)
    {
      for (int j=k; j<n; j++)
	x[j] = x[j]/(alpha[j]-alpha[j-k]);	
      for (int j=k-1; j<n-1; j++)
	x[j] = x[j]-x[j+1];
    }
}
