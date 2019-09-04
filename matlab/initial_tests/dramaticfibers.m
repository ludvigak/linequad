% ideas for dramatic 3d fiber demos. needs gauss.m
% Barnett 5/23/19

clear

% helix from t in [0,1], param so could use uniform panelization of [0,1]
w = 100; % freq
a = 2;   % growth rate
e = @(t) exp(a*t);
x = @(t) e(t);
y = @(t) e(t).*cos(w*t);
z = @(t) e(t).*sin(w*t);
npan = 50; figure(1); showfiber(x,y,z,npan); title('geometric helix')

% random squiggle loop in t in [0,1]
K = 20;   % max Fourier mode in squiggle
K0 = 5;   % mode beyond which decay kicks in
k = -K:K;
nk = numel(k);
rng(0);
ampl = (K0+abs(k)).^(-1);   % Sobolev-type decay
c = (randn(3,nk)+1i*randn(3,nk)) .* ampl;    % cmplx F coeffs in each coord
x = @(t) real(c(1,:)*exp(2i*pi*k'*t(:)'));   % outer prod to eval the F series
y = @(t) real(c(2,:)*exp(2i*pi*k'*t(:)'));
z = @(t) real(c(3,:)*exp(2i*pi*k'*t(:)'));
npan = 100; figure(2); showfiber(x,y,z,npan); title('squiggle')

function showfiber(x,y,z,npan)  % plots unif and G-L panel versions
clf
subplot(1,2,1);
n=3e3; t=(0.5:n-0.5)/n;
plot3(x(t),y(t),z(t),'-'); axis equal tight vis3d
subplot(1,2,2);
p = 16;
[X W] = gauss(p); X = (X+1)/2;   % on [0,1]
t = bsxfun(@plus, X/npan, (0:npan-1)/npan);
t = t(:);             % G-L panelization of param
plot3(x(t),y(t),z(t),'.-');
te = (0:npan)/npan;     % panel param ends
hold on; plot3(x(te),y(te),z(te),'k.','markersize',10);
axis equal tight vis3d
end
