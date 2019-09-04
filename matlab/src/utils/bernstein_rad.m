function rho = bernstein_rad(t)

rho1 = abs(t + sqrt(t.^2-1));
rho2 = abs(t - sqrt(t.^2-1));
rho = max(rho1, rho2);
