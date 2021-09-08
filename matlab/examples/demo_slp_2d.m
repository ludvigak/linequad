function demo_slp_2d()
% Demo for how to evaluate the 2D Laplace single layer potential on one
% panel

% Parametrization of panel: parabola
% k=Curving of panel, higher is harder to compute       
% try range 0.2--0.4
k = 0.3; 
g = @(t) t + 1i*k*t.^2 -1i*k;
gp = @(t) 1 + 1i*2*k*t;    

% Example density to use
density_xy = @(t) real(g(t)).*imag(g(t)); % sigma(t) = x(t)*y(t)
density_t = @(t) t; % sigma(t) = t
density_1 = @(t) t.^0; % sigma(t) = 1, (trivial case for complex interpolatory)

density = density_xy; % try xy/t/1

slp_integrand = @(t, z) density(t) .* log(abs(z-g(t))) .* abs(gp(t))
reference = @(z) integral(@(t) slp_integrand(t, z), -1, 1, 'abstol', eps, 'reltol', eps);

% Gauss-Legendre quadrature nodes and weights
[t16,w16] = lgwt2(16,-1,1);

% Create a grid
Ngrid = 40;
t1 = linspace(-2.9, 2.9, Ngrid);
t2 = linspace(-1.8, 1.8, Ngrid);

% Create Legendre expansion of panel (in rescaled coordinates)
za = g(-1);
zb = g(1);
L = legendre.matrix(numel(t16));
coeffs = L*rotate_and_scale(za, zb, g(t16));

% Create a grid around curve segment by complexifying parameter
[T1, T2] = ndgrid(t1, t2);
T = T1 + 1i*T2;
Z = g(T);

% Matrix used for computing weights
A = fliplr(vander(t16));

% Compute values on grid
[I, Q16, SS16] = deal(zeros(size(Z)));
tic
parfor i=1:numel(Z)
        z = Z(i);
        % === Reference
        I(i) = reference(z);
        % === Direct
        Q16(i) = sum(slp_integrand(t16, z) .* w16);    
        % === Singularity swap
        % * Find root
        zr = rotate_and_scale(za, zb, z);
        t0 = zr;
        t = legendre.newton(coeffs, zr, t0, 20, 1e-13);
        % * Compute complex integrals for t
        [pL, p1] = complex_integrals(t, 16);
        % * Compute weights
        wL = A.' \  real(pL); % Can use LU here since A is constant
        wLcorr = wL - log(abs(t-t16)) .* w16; % Log correction weights
        % * Evaluate singularity swap
        SS16(i) = Q16(i) + sum(density(t16) .* abs(gp(t16)) .* wLcorr);
end
toc

% Compute rel. errors
Inorm = norm(I(:), inf);
E16 = abs(Q16-I) / Inorm + eps;
ESS16 = abs(SS16-I) / Inorm + eps;

sfigure(1); clf
pcolor(real(Z), imag(Z), log10(E16));    
title(['Direct, 16 pts, k=' num2str(k)])
hold on

sfigure(2); clf
pcolor(real(Z), imag(Z), log10(ESS16));    
title(['Singularity swap, 16 pts, k=' num2str(k)])
hold on

for i=1:2
    sfigure(i)
    plot(complex(g(t16)), '.-r')            
    ax = [-2.4, 2.4, -1.9, 1.4];
    axis equal
    axis(ax)
    %axis image
    shading interp
    caxis([-16 -6])
    axis off        
end