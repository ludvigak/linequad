% test close eval to fiber for 3D Laplace SLP line source.
% fiber is [-1,1] on x-axis, target coords (a,b). P is 1/r Laplace kernel.
% Q is integral w/ kernel r, P is w/ kernel 1/r (the 3D Laplace).
% Ideas of Ludvig af Klinteberg and Alex Barnett 5/1/18
% Needs: gauss, clencurt.

clear; verb = 1;
a=-0.7; b = 1e-3; % targ loc (b = min dist if a in [-1,1])
dens = @(t) cos(2*t + 0.8);
%dens = @(t) 2.1+0.7*t;  % tests q_0 and q_1 in isolation
rint = @(t) dens(t) .* sqrt((t-a).^2+b^2);     % integrand dens.r  (for Q)
irint = @(t) dens(t) ./ sqrt((t-a).^2+b^2);     % integrand dens/r  (for P)

p=16;    % # density rep nodes, sets maxdegree = p-1
tj = gauss(p);  % dens nodes
fj = dens(tj);  % dens data

t = linspace(-1,1,1e3)';  % fine grid for tests, col vec
if verb, figure; plot(t,dens(t),'b-'); hold on;
  plot(t,rint(t),'r-'); plot(t,irint(t),'g-');
  plot(tj,fj,'bo'); plot(tj,rint(tj),'ro'); plot(tj,irint(tj),'go');
  axis([-1 1 -2 2]); legend('dens','rint','irint'); title('integrands'); end

ne = ceil(20/b);   % "exact" integrals of dens*r, dens/r
[te we]=clencurt(ne);
Qe = sum(we.'.*rint(te)); fprintf('Q exact (ne=%d):\t %.16g\n',ne,Qe)
Pe = sum(we.'.*irint(te)); fprintf('P exact (ne=%d):\t %.16g\n',ne,Pe)

% check (shifted?) poly rep is good for dens...
T = a;  % Taylor center for poly
V = ones(p); for j=2:p, V(:,j) = V(:,j-1).*(tj-T); end % Vandermonde on tj
%c = V\fj;   % poly coeffs in O(p^3), stable
[L,U]=lu(V);  % computed once and for all
c = U\(L\fj);  % O(p^2), stable
Vf = ones(numel(t),p); for j=2:p, Vf(:,j) = Vf(:,j-1).*(t-T); end % Vander eval
err = Vf*c - dens(t);
fprintf('dens poly interp sup err = %.3g\n',max(abs(err)))
if verb, figure; plot(t,err,'-'); hold on; plot(tj,0*tj,'o');
  title('dens poly approx error func'); end

q = nan(p,1);   % r-kernel: q_k array for this targ
fq0 = @(s) (s*sqrt(s^2+b^2) + b^2*log(s+sqrt(s^2+b^2)))/2; % antideriv
q(1) = fq0(1-a)-fq0(-1-a);  % q_0  (note 1-offset index from k)
fq1 = @(s) (s^2+b^2)^(3/2)/3;
q(2) = fq1(1-a)-fq1(-1-a);  % q_1
for k=0:p-3
  q(k+3) = (-(k+1)*b^2*q(k+1) + (1-a)^(k+1)*((1-a)^2+b^2)^(3/2) - (-1-a)^(k+1)*((-1-a)^2+b^2)^(3/2))/(4+k);
end
%q  % note q_k grows for a near ends

Q = sum(q.*c);  % interpolatory quadr
fprintf('Q by recurrence:\t %.16g\t\trel err = %.3g\n',Q,(Q-Qe)/Qe)

% now 1/r kernel: p_k array for this targ
pk = nan(p,1);
fp0 = @(s) log(s+sqrt(s^2+b^2));
pk(1) = fp0(1-a)-fp0(-1-a);
pk(2) = sqrt((1-a)^2+b^2) - sqrt((-1-a)^2+b^2);
for k=0:p-3, pk(k+3) = q(k+1) - b^2*pk(k+1); end   % p recurrence

P = sum(pk.*c);  % interpolatory quadr
fprintf('P by recurrence:\t %.16g\t\trel err = %.3g\n',P,(P-Pe)/Pe)

