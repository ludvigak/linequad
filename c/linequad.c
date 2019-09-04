#include "stdio.h"
#include "math.h"
#include "stdbool.h"

#include "linequad.h"

#include "vandermonde.c"
#include "legendre.c"
#include "coeffs.c"

#ifdef VERBOSE
#define VERBINFO(...) printf(__VA_ARGS__)
#else
#define VERBINFO(...) /* nothing */
#endif

void rsqrt_pow_weights(const double* restrict tj,
		       double complex troot,
		       int n, int pmax,
		       double* restrict w1,
		       double* restrict w3,
		       double* restrict w5)
{
  double I1[MAX_EXPANSION_ORDER], I3[MAX_EXPANSION_ORDER], I5[MAX_EXPANSION_ORDER];    
  rsqrt_pow_integrals(troot, n, pmax, I1, I3, I5);
  pvand(n, tj, w1, I1);
  if (pmax >= 3)
    pvand(n, tj, w3, I3);
  if (pmax >= 5)
    pvand(n, tj, w5, I5);    
  for (int i=0; i<n; i++)
    {
      double a = cabs(tj[i]-troot);
      double a2 = a*a;
      w1[i] *= a;
      w3[i] *= a*a2;
      w5[i] *= a*a2*a2;
    }
}

void rsqrt_pow_integrals(double complex z, int N, int pmax,
			 double* restrict I1,
			 double* restrict I3,
			 double* restrict I5)
{
  double zr = creal(z);
  double zi = cimag(z);
  // (t-zr)^2+zi^2 = t^2-2*zr*t+zr^2+zi^2 = t^2 + b*t + c
  double b = -2*zr;
  double c = zr*zr+zi*zi;
  double d = zi*zi; // d = c - b^2/4;

  //u1 = sqrt(1 - b + c);
  //u2 = sqrt(1 + b + c);
  // The form below gets *much* better accuracy than ui = sqrt(ti^2 + b*ti + c)
  double u1 = sqrt((1+zr)*(1+zr) + zi*zi);
  double u2 = sqrt((1-zr)*(1-zr) + zi*zi);

  // Compute I1
  // k=1
  // Series expansion inside rhombus
  if (4*fabs(zi) < 1-fabs(zr))
  {
    double arg1 = (1-fabs(zr))*eval_series(coeffs_I1, 1-fabs(zr), zi, 11);
    double arg2 = 1+fabs(zr) + sqrt((1+fabs(zr))*(1+fabs(zr)) + zi*zi);
    I1[0] = log(arg2)-log(arg1);   
  }
  else
  {
    // Evaluate after substitution zr -> -|zr|
    double arg1 = -1+fabs(zr) + sqrt((-1+fabs(zr))*(-1+fabs(zr)) + zi*zi);        
    double arg2 =  1+fabs(zr) + sqrt((1+fabs(zr))*(1+fabs(zr)) + zi*zi);
    I1[0] = log(arg2)-log(arg1);   
  } 
  // k>1
  if (N>1)
    I1[1] = u2-u1 - b/2*I1[0];    
  double s = 1.0;
  for (int n=2; n<N; n++) {
    s = -s; // (-1)^(n-1)
    I1[n] = (u2-s*u1 + (1-2*n)*b/2*I1[n-1] - (n-1)*c*I1[n-2])/n;
  }
  if (pmax < 3) return;

  // Compute I3
  // Series is needed in cones extending around real axis from
  // interval endpoints
  double w = fmin(fabs(1+zr), fabs(1-zr)); // distance interval-singularity
  zi = fabs(zi);
  double zi_over_w = zi/w;
  bool in_cone = (zi_over_w < 0.6);
  bool outside_interval = (fabs(zr)>1);
  bool use_series = (outside_interval & in_cone);
  if (!use_series)
    {
      I3[0] = (b+2)/(2*d*u2) - (b-2)/(2*d*u1);
    }
  else
    {
      // Power series for shifted integral
      // pick reasonable number of terms
      int Ns;
      if (zi_over_w < 0.01)
	Ns = 4;
      else if (zi_over_w < 0.1)
	Ns = 10;
      else if (zi_over_w < 0.2)
	Ns = 15;
      else // zi_over_w < 0.6
	Ns = 30;
      double x1 = -1-zr;
      double x2 = 1-zr;
      double Fs1 = fabs(x1)/(x1*x1*x1) * (-0.5 + eval_series(coeffs_I3, x1, zi, Ns));
      double Fs2 = fabs(x2)/(x2*x2*x2) * (-0.5 + eval_series(coeffs_I3, x2, zi, Ns));	
      I3[0] = Fs2 - Fs1;
    }
  if (N>1)
    I3[1] = 1/u1-1/u2 - b/2*I3[0];
  for (int n=2; n<N; n++)
    I3[n] = I1[n-2] - b*I3[n-1] - c*I3[n-2];
  if (pmax < 5) return;

  // Compute I5
  // Here too use power series for first integral, in cone around real axis
  in_cone = (zi_over_w < 0.7);
  use_series = (outside_interval && in_cone);
  if (!use_series)
    {
      I5[0] = (2+b)/(6*d*u2*u2*u2) - (-2+b)/(6*d*u1*u1*u1) + 2/(3*d)*I3[0];
    }
  else
    {
      // Power series for shifted integral
      int Ns;
      if (zi_over_w < 0.01)
	Ns = 4;
      else if (zi_over_w < 0.2)
	Ns = 10;
      else if (zi_over_w < 0.5)
	Ns = 24;        
      else if (zi_over_w < 0.6)
	Ns = 35;                
      else // zi/w < 0.7
	Ns = 50;
      double x1 = -1-zr;
      double x2 = 1-zr;      
      double Fs1 = 1/(x1*x1*x1*fabs(x1)) *(-0.25 + eval_series(coeffs_I5, x1, zi, Ns));    
      double Fs2 = 1/(x2*x2*x2*fabs(x2)) *(-0.25 + eval_series(coeffs_I5, x2, zi, Ns));    
      I5[0] = Fs2 - Fs1;
    }      

  if (N>1)
    {
      // Second integral computed using shifted version, and then corrected by first
      // integral (which is correctly evaluated using power series)
      // This is analogous to the formula for I3(1), but was overlooked by Tornberg & Gustavsson
      I5[1] = 1/(3*u1*u1*u1) - 1/(3*u2*u2*u2)  - b/2*I5[0];
    }

  for (int n=2; n<N; n++)
    I5[n] = I3[n-2] - b*I5[n-1] - c*I5[n-2];
}

int rootfinder(const double* restrict xhat,
	       const double* restrict yhat,
	       const double* restrict zhat,
	       int n,
	       double x0,
	       double y0,
	       double z0,
	       double complex tinit,
	       double complex* troot)
{
  // Find roots using Newton and Muller
  double complex t = tinit;
  double tol = 1e-14;
  int maxiter_newton = 20;
  int maxiter_muller = 20;    
  // === Newton
  // Setup history variables (needed in Muller)
  double complex Fp, tp, Fpp, tpp;
  double complex dt, F, Fprime;
  int converged = 0;
  int iter;
  double absres;
  double complex P[MAX_EXPANSION_ORDER];
  double complex D[MAX_EXPANSION_ORDER];
    
  for (iter=0; iter<maxiter_newton; iter++)
    {
      double complex x, y, z, xp, yp, zp;
      double complex dx, dy, dz;
      // Legendre eval
      legendre_deriv_scalar(n-1, t, P, D);	
      x = eval_legendre_expa(P, xhat, n);
      y = eval_legendre_expa(P, yhat, n);
      z = eval_legendre_expa(P, zhat, n);
      xp = eval_legendre_expa(D, xhat, n);
      yp = eval_legendre_expa(D, yhat, n);
      zp = eval_legendre_expa(D, zhat, n);
      // Compute F and F'
      dx = x-x0;
      dy = y-y0;
      dz = z-z0;
      F = dx*dx + dy*dy + dz*dz;
      Fprime = 2*(dx*xp + dy*yp + dz*zp);	
      dt = -F/Fprime;
      // Update history
      tpp = tp;
      Fpp = Fp;
      Fp = F;
      tp = t;
      // Update root
      t = t+dt;
      absres = cabs(dt);
      if (absres < tol)
	{
	  converged = 1;
	  break;
	}    
    }
  if (converged==1)
    {
      VERBINFO("Newton converged in %d iterations.\n", iter);
      *troot = t;
      return converged;
    }
  // === Muller
  VERBINFO("Newton did not converge after %d iterations (abs(dt)=%g), switching to Muller\n", iter, absres);
  converged = 0;
  for (iter=0; iter<maxiter_muller; iter++)
    {
      double complex x, y, z;
      double complex dx, dy, dz;	
      double complex q, A, B, C, d1, d2;
      // Compute
      legendre_scalar(n-1, t, P);
      x = eval_legendre_expa(P, xhat, n);
      y = eval_legendre_expa(P, yhat, n);
      z = eval_legendre_expa(P, zhat, n);
      dx = x-x0;
      dy = y-y0;
      dz = z-z0;
      F = dx*dx + dy*dy + dz*dz;
      // Mullers method
      q = (t-tp)/(tp - tpp);
      A = q*F - q*(q+1)*Fp + q*q*Fpp;
      B = (2*q+1)*F - (1+q)*(1+q)*Fp + q*q*Fpp;
      C =(1+q)*F;
      d1 = B + csqrt(B*B-4*A*C);
      d2 = B - csqrt(B*B-4*A*C);
      if (cabs(d1) > cabs(d2))
	dt = -(t-tp)*2*C/d1;
      else
	dt = -(t-tp)*2*C/d2;
      // Update history
      tpp = tp;
      Fpp = Fp;
      Fp = F;
      tp = t;
      // Update root
      t = t+dt;
      absres = cabs(dt);
      if (absres < tol)
	{
	  converged = 1;
	  break;
        }
    }    
  if (converged)
    VERBINFO("Muller converged in %d iterations.\n", iter);
  else
    VERBINFO("Muller did not converge after %d iterations. abs(dt)=%g",
	     iter, absres);
  *troot = t;    
  return converged;
}

double complex rootfinder_initial_guess(const double* restrict tj,
					const double* restrict xj,
					const double* restrict yj,
					const double* restrict zj,
					int n, double x0, double y0, double z0)
{
  double complex tinit;
  // Find two closest point  
  double Rsqmin1 = INFINITY;
  double Rsqmin2 = INFINITY;
  int imin1, imin2;
  for (int i=0; i<n; i++)
    {
      double dx = xj[i]-x0;
      double dy = yj[i]-y0;
      double dz = zj[i]-z0;
      double Rsq = dx*dx + dy*dy + dz*dz;
      if (Rsq < Rsqmin1)
	{
	  Rsqmin2 = Rsqmin1;
	  imin2 = imin1;	  	  
	  Rsqmin1 = Rsq;
	  imin1 = i;
	}
      else if (Rsq < Rsqmin2)
	{
	  Rsqmin2 = Rsq;
	  imin2 = i;	  
	}
    }
  // Now compute initial guess
  int i1 = imin1;
  int i2 = imin2;
  double t1 = tj[i1];
  double t2 = tj[i2];
  double p[3];
  p[1] = xj[i1]-xj[i2];
  p[2] = yj[i1]-yj[i2];
  p[3] = zj[i1]-zj[i2];
  double pnorm = sqrt(p[1]*p[1] + p[2]*p[2] + p[3]*p[3]);
  double r[3];
  r[1] = x0-xj[i1];
  r[2] = y0-yj[i1];
  r[3] = z0-zj[i1];  
  double rnormsq = r[1]*r[1] + r[2]*r[2] + r[3]*r[3];  
  double rdotp = r[1]*p[1] + r[2]*p[2] + r[3]*p[3];
  double a = (t1-t2)*rdotp/(pnorm*pnorm);
  double b = sqrt(rnormsq-rdotp*rdotp/(pnorm*pnorm)) * (t1-t2)/pnorm;
  tinit = t1 + a + I*b;
  return tinit;
}
