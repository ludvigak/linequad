clear

p = 16
k = 0.9

g = @(t) t + 1i*k*t.^2;
freq = 0.1; f = @(t) sin(1+freq*t);
tj = lgwt(p,-1,1);
fj = f(tj);
tauj = g(tj);

Ac = ones(p,p); for j=2:p, Ac(:,j) = Ac(:,j-1).*tauj; end
Ar = ones(p,p); for j=2:p, Ar(:,j) = Ar(:,j-1).*tj; end

cc = Ac\fj;
cr = Ar\fj;

% Construct matrix M that maps expansion coeffs
% from complex monomial basis to real
q = [0 1 1i*k];
c = 1;
M = zeros(2*p-1, p);
M(1) = 1;
for i=2:p
    % c is poly coeffs for (t+1i*k*t^2)^(p-1)
    c = conv(c, q);
    M(1:2*i-1, i) = c;
end
Mtrunc = M(1:p, 1:p); % Truncate to square

condM = cond(Mtrunc)
err = norm(Ar*(Mtrunc*cc) - fj) / norm(fj)

norm(Mtrunc, inf)