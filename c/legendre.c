static inline void legendre_deriv_scalar(int n, double complex x,
					 double complex * restrict P,
					 double complex * restrict D)
{
  P[0] = 1.0; // l=0
  D[0] = 0.0;
  P[1] = x; // l=1
  D[1] = 1.0;
  for (int l=1; l<n; l++)
    {
      // Compute l+1
      P[l+1] = ( (2*l+1)*x * P[l] - l*P[l-1] ) / (l+1);
      D[l+1] = ( (2*l+1)*(P[l] + x * D[l]) - l*D[l-1] ) / (l+1);
    }
}

static inline void legendre_scalar(int n, double complex x,
				   double complex * restrict P)
{
  P[0] = 1.0; // l=0
  P[1] = x; // l=1
  for (int l=1; l<n; l++)
    {
      // Compute l+1
      P[l+1] = ( (2*l+1)*x * P[l] - l*P[l-1] ) / (l+1);
    }
}

static inline double complex eval_legendre_expa(const double complex* restrict P,
						const double* restrict fhat,
						int n)
{
  double complex f = 0.0;
  for (int i=0; i<n; i++)
    f += P[i]*fhat[i];
  return f;
}
