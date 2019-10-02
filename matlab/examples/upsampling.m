% Illustration of why upsampling is needed sometimes (2D)
% =======================================================
%
% integrand here is 
% i*g'(t)/(g(t0)-g(t))
% = i(1+2ikt) / (t0-t+ik(t0^2-t^2))
%
% regularized integrand is
% i*g'(t)*(t0-t)/(g(t0)-g(t))
% = i(1+2ikt) / (1+ik(t0+t))
%
% so smoothness of regularized integrand is limited by 
% location of secondary root, t1=i/k-t0

clear

% Parabolic panel
k = 0.4;
g = @(t) t + 1i*k*t.^2 -1i*k;
gp = @(t) 1 + 1i*2*k*t;    

% Nearby target point
t0 = 0.123+0.1i;
z = g(t0);

% Simplest possible density
density = @(t) t.^0; % sigma(t) = 1

% Laplace double layer kernel
kernel = @(x, y, ny) ( ny ./ (x-y) );    
integrand = @(t,z) kernel(z, g(t), 1i*gp(t)) .* density(t);

% Gauss-Legendre quadrature nodes and weights
[t16,w16] = lgwt2(16,-1,1);
[t32,w32] = lgwt2(32,-1,1);

% Sample
g16 = g(t16);
n16 = 1i*gp(t16);
sigma16 = density(t16);

% Upsample
M = bclag_interp_matrix(t16, bclag_interp_weights(t16), t32);
g32 = M*g16;
n32 = M*n16;
sigma32 = M*sigma16;

% Compute SSQ
wssq16 = real_quad_2d(t16, t0);        
Qssq16 = real(sum(kernel(z, g16, n16) .* sigma16 .* wssq16));
wssq32 = real_quad_2d(t32, t0);        
Qssq32 = real(sum(kernel(z, g32, n32) .* sigma32 .* wssq32));

% Compute errors
I = real(integral(@(t) integrand(t, z),  -1,1,'abstol',eps,'reltol',eps))
Essq16 = abs(Qssq16-I) 
Essq32 = abs(Qssq32-I)

% Plot regularized integrand
t = linspace(-1,1,300);
sfigure(1);clf
plot(t, real( integrand(t,z)), '-', 'DisplayName', 'original')
hold on
plot(t, real( integrand(t,z) .* (t-t0)), '-', 'DisplayName', 'regularized')
xlabel('t')
title('integrand')
legend('location','best')
grid on

% This is equal to real( integrand(t,z) .* (t-t0))
regularized = @(t) -real( 1i*(1+2i*k*t) .* 1 ./ (1+1i*k*(t0+t)));
t1 = 1i/k-t0;

% Plot poly coeffs for original and regularized integrands
sfigure(2);clf
L32 = legendre.matrix(32);
c32 = L32 * (real( integrand(t32,z)));
creg32 = L32 * (real( integrand(t32,z) .* (t32-t0)));
semilogy(abs(c32),'.-', 'DisplayName', 'original')
hold on
semilogy(abs(creg32),'.-', 'DisplayName', 'regularized')
grid on
title('Legendre coeffs')
legend('Location','best')
ylabel('|c_k|')
xlabel('k')

% Plot abs. of complex continuation
tr = linspace(-1,1);
ti = linspace(-3,3);
[Tr, Ti] = meshgrid(tr,ti);
T = Tr + 1i*Ti;

sfigure(3);clf
pcolor(Tr, Ti, abs((integrand(T,z)))), shading interp
title('abs. of integrand')
hold on
plot(t0, 'or', 'MarkerSize',10)
plot(t1, 'or', 'MarkerSize',10)
text(real(t0)+0.1, imag(t0), 't0', 'color','w')
text(real(t1)+0.1, imag(t1), 't1', 'color','w')
xlabel('Re t')
ylabel('Im t')

sfigure(4);clf
pcolor(Tr, Ti, abs((integrand(T,z).*(T-t0)))), shading interp
title('abs. of regularized integrand')
xlabel('Re t')
ylabel('Im t')
hold on
plot(t1, 'or', 'MarkerSize',10)
text(real(t1)+0.1, imag(t1), 't1', 'color','w')
