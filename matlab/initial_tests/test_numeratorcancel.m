% test numerator catastrophic cancellation fixes for 3D lines.
% Aborted.
% Barnett 9/22/19

clear    % let's do the t-plane, integral [-1,1]. test m=3 case, q=2

d = 1e-5; %1e-4;  % dist
s0 = 0.0;    % nr param pt (try 0 for a hope of fixing)
t0 = s0 + 1i*d;   % aka troot.  s0 must = Re(t0)

%f = @(t) sin(t+1);   % generic "density"
q=2;     % known numerator power in...
f = @(t) sin(t+1).*(t-s0+0.6*d).^q; % "density", RR^T style vanishing near s0

R = @(t) abs(t-t0);   % dist func in C
integrand3 = @(t) f(t)./R(t).^3;
I3e = integral(integrand3, -1, 1, 'reltol', eps, 'abstol', eps);
fprintf('exact       \tI=%.15g\n',I3e)

% quad
n = 16;
[tj,wj] = legendre.gauss(n);

% direct... (just for debug tj,wj)
fj = f(tj);  % samples
I3d = sum(fj./R(tj).^3.*wj);
%fprintf('n=%d direct\tI=%.15g \trel err %.3g\n',n,I3d,abs((I3d-I3e)/I3e))

shift = 0;       % whether to use shifted monomial basis for c
if ~shift, sh = 0;
  [p1, p3, p5] = rsqrt_pow_integrals(t0, n);    % recurrences
else, sh = s0;  % note I3(3) only is not acc here...
  [p1, p3, p5] = rsqrt_pow_integrals_shift(t0, n);    % recurrences for shift s0
end

% plain non-adj method... (note wei w3 include |t0-tj|^m factors)
V = ones(n,n); for k=2:n, V(:,k) = V(:,k-1).*(tj-sh); end   % Vandermonde not transp, and poss shifted
c = V\fj(:);           % c_k coeffs of f
fprintf('rel err on c(1): %.3g\n',abs((f(sh)-c(1))/f(sh)))
c(1) = f(sh);   % try repair c(1) ! -> works (up to p3 err), shows the culprit
I3 = sum(c.*p3);
fprintf('recur, non-adj\tI=%.15g \trel err %.3g\n',I3,abs((I3-I3e)/I3e))

% adjoint meth: get lambda wei first
lam = V'\p3;          % as in rsqrt_pow_weights, but without the |tj-t0|^3 fac
I3a = sum(lam.*fj);
fprintf('recur, adj   \tI=%.15g \trel err %.3g\n',I3a,abs((I3a-I3e)/I3e))

% following fail... due to not having right q-order-root location in f.
% (F is in fact not smooth)

% repaired non-adj:
F = @(t) f(t)./(t-t0).^2;   % density w/ sim-to-known-factor killed
Fj = F(tj);   % assume still relative acc for each eval
d = V\Fj(:);           % d_k coeffs of F
F(sh)
fprintf('rel err on d(1): %.3g\n',abs((F(sh)-d(1))/F(sh)))
M = toeplitz([t0^2,-2*t0,1,zeros(1,n-3)],[t0^2,zeros(1,n-1)]); % q=2 only
c = M*d;
fprintf('rel err on repaired c(1): %.3g\n',abs((f(sh)-c(1))/f(sh)))
I3r = sum(c.*p3);
fprintf('repair non-adj\tI=%.15g \trel err %.3g\n',I3r,abs((I3r-I3e)/I3e))

% repaired adj meth:
lam = V'\(M'*p3);
I3ra = sum(lam.*Fj);
fprintf('repair adj\tI=%.15g \trel err %.3g\n',I3ra,abs((I3ra-I3e)/I3e))
