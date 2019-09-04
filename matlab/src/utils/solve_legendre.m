function [t, relres, iter] = solve_legendre(c, zr, t0, maxiter, tol)

% These should both be arguments
maxiter1 = maxiter;
maxiter2 = maxiter;
% Break if diverging
Rbreak = 100;
% Setup and run Newton
maxiters = maxiter1+maxiter2;
t = t0;
lmax = numel(c)-1;
relres = tol;
thist = zeros(maxiters,1);
fhist = zeros(maxiters,1);
for iter=1:maxiter1
    [P, D] = legendre_deriv_scalar(lmax, t);
    f = P*c - zr;
    thist(iter) = t;
    fhist(iter) = f;
    fp = D*c;
    dt = -f / fp;
    t = t + dt;
    relres = abs(dt/t);
    if relres < tol || abs(t) > Rbreak
        break
    end
end

if relres < tol || abs(t) > Rbreak
    return
end

for iter=1:maxiter2
    [P, ~] = legendre_deriv_scalar(lmax, t);
    F = P*c - zr;
    idx = 2+iter; % Restart after first Newton step
    thist(idx) = t;
    fhist(idx) = F;

    % Mullers method
    tp = thist(idx-1);
    tpp = thist(idx-2);

    Fp = fhist(idx-1);
    Fpp = fhist(idx-2);
    
    q = (t-tp)/(tp - tpp);
    A = q*F - q*(q+1)*Fp + q*q*Fpp;
    B = (2*q+1)*F - (1+q)*(1+q)*Fp + q*q*Fpp;
    C =(1+q)*F;
    d1 = B + sqrt(B*B-4*A*C);
    d2 = B - sqrt(B*B-4*A*C);
    if (abs(d1) > abs(d2))
        dt = -(t-tp)*2*C/d1;
    else
        dt = -(t-tp)*2*C/d2;
    end        

    t = t + dt;
    relres = abs(dt/t);
    if relres < tol || abs(t) > Rbreak
        break
    end
end

return % Disable warnings
if relres > tol 
    warning('Muller failed')
end

