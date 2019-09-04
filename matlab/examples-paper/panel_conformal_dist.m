function u = panel_conformal_dist(g,z,n)
% PANEL_CONFORMAL_DIST.  conformal distance of points from analytic panel arc
%
% u = panel_conformal_dist(g,z)
%  returns the conformal distances u of each point z in the neighborhood of
%  a panel parameterized by g : [-1,1] -> C, in the complex plane C.
%
% u = panel_conformal_dist(g,z,n) controls the accuracy, taking O(n^3) time.
%
% Uses potential theory: solves u=0 on curve, with u(r) ~ log(r) as r -> infty.
% This is solved via crude SLP on curve, plus constant, with nodes of Cheby
% density, and plain charges (no quadr wei) at each, so that charge sum = 1.
%
% The singular quadrature is O(1/n) for now, using int(ln x) = x(ln x - 1)
% applied from 0 to h/2, (h = local grid spacing), then with overall scaling
% of 1/h which cancels out the factor h out front. A hack, but n=100 gives 1e-3.
% Can't quite get O(1/n^2), but prefactor in 1/n is small.
%
% Without arguments, does a self-test.
%
% Needs: diagind.

% Barnett 3/20/19

if nargin==0, test_panel_conformal_dist; return; end
if nargin<3, n = 100; end        % default # nodes

tj = cos(pi*(0:n-1)/(n-1));   % row vec of Clenshaw-Curtis nodes on [-1,1]
zj = g(tj);                   % nodes in C plane
A = log(abs(zj - zj.'));      % potentials between all pairs of pts (SLP eval)
%A(diagind(A)) = 0;            % no sing quadr rule at all! :)
% eff quadr node box sizes:
hj = cos(pi*(-0.5:n-0.5)/(n-1)); hj(1) = 1; hj(n) = -1; hj = abs(diff(hj));
%fprintf('err in sum of grid spacings hj = %.3g\n',sum(hj)-2)  % O(1/n^2)
A(diagind(A)) = log(hj/2)-1;  % 1st-order sing quadr term (integr ln x)
%figure; imagesc(A)                 % check diag matches near-diag
A = [A, ones(n,1)];           % plus unknown const affecting u everywhere
A = [A; [ones(1,n) 0]];       % augment lin sys w/ constraint total charge = 1
b = [zeros(n,1); 1];
q = A\b;                      % solve for charges
%fprintf('resid nrm = %.3g, sum(q)-1 = %.3g\n',norm(A*q-b),sum(q(1:n))-1)
u = q(n+1) + log(abs(zj - z(:))) * q(1:n);  % eval potential (const + SLP)
u = reshape(u,size(z));

%%%%%
function test_panel_conformal_dist
pan = @(t) t;         % flat panel [-1,1]
g = -2:0.02:2;
[xx yy] = meshgrid(g,g); zz = xx + 1i*yy;
u = panel_conformal_dist(pan,zz);
figure; subplot(1,2,1);
contourf(g,g,u); colorbar; title('potential u (should be 0 on line)');
axis equal tight;
uex = log(abs(zz+sqrt(zz-1).*sqrt(zz+1)));   % Trefethen Spec Meth book Ch.5
subplot(1,2,2); imagesc(g,g,u-uex); colorbar; title('error in potential');
axis equal tight xy;
i = round(numel(g)*0.75); fprintf('u err @ pt: %.3g\n',u(i,i)-uex(i,i))

k = -1.0; pan = @(t) t + 1i*k*t.^2;   % parabola panel
u = panel_conformal_dist(pan,zz);
figure; contourf(g,g,u); colorbar; title('potential u (should be 0 on panel)');
axis equal tight;
u2 = panel_conformal_dist(pan,zz,200);
fprintf('parabola: u @ pt = %.3g, diff: %.3g\n',u(i,i), u(i,i)-u2(i,i))
hold on; plot(pan(linspace(-1,1,50)),'w-');
