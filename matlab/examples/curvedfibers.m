function curvedfibers()
% Test close eval of kernels 1/R, 1/R^3, and 1/R^5 for single point close to 3D fiber
% Fiber parametrized in t from -1 to 1
% Ludvig af Klinteberg and Alex Barnett, May 2018   
    
    % Parametrization of fibre (part of helix)
    r = 1/2;
    k = pi/2;
    c = 1;
    scale = 100; % To check invariance
    x = @(t) r*cos(k*t)*scale;
    y = @(t) r*sin(k*t)*scale;
    z = @(t) c*t*scale;       
     
    % Target point near bdry:
    s0 = 0.2; d = .01;
    
    % Target point on extension:
    %s0 = 1.01; d = 0;
    
    distance = d*scale
    x0 = x(s0);
    y0 = y(s0)+distance;
    z0 = z(s0);
    
    % Density to integrate
    density = @(t) exp(t/2);    
    
    % === END PARAMS
    
    % - Integrands with 1/R^p kernels, p=1,3,5
    % - These are the kernels needed for the stokeslet and the doublet
    % of slender body theory.
    % - This is actually missing line element! (not important here)
    rx = @(t) x(t)-x0;
    ry = @(t) y(t)-y0;
    rz = @(t) z(t)-z0;
    R = @(t) sqrt(rx(t).^2 + ry(t).^2 + rz(t).^2);
    integrand1 = @(t) density(t)./ R(t);
    integrand3 = @(t) density(t)./ R(t).^3;
    integrand5 = @(t) density(t)./ R(t).^5;
    
    % Compute reference
    I1 = integral(integrand1, -1, 1, 'reltol', eps, 'abstol', eps);
    I3 = integral(integrand3, -1, 1, 'reltol', eps, 'abstol', eps);
    I5 = integral(integrand5, -1, 1, 'reltol', eps, 'abstol', eps);
    
    % Discretization
    n = 16;
    [tj,wj] = legendre.gauss(n);
    xj = x(tj);
    yj = y(tj);
    zj = z(tj);
    f1j = integrand1(tj);
    f3j = integrand3(tj);
    f5j = integrand5(tj);
        
    % Direct quadrature
    Q1n = wj'*f1j;
    Q3n = wj'*f3j;
    Q5n = wj'*f5j;
    
    % Legendre expansions of discretization
    L = legendre.matrix(n); % O(n^2) per panel
    xhat = L*xj;
    yhat = L*yj;
    zhat = L*zj;
       
    % Find root for this target point
    tinit = rootfinder_initial_guess(tj, xj, yj, zj, x0, y0, z0); % O(n) per point
    troot = rootfinder(xhat, yhat, zhat, x0, y0, z0, tinit); % O(n) per point
    
    % * This is the place to decide if special quad is needed,
    % based on radius of Bernstein ellipse of troot
    % e.g.
    bernstein_radius = @(z) abs(z + sign(real(z)).*sqrt(z.^2 - 1));
    rho = bernstein_radius(troot);
    % rho > 2-3 usually gives full accuarcy for n=16
    % In general error is prop to rho^(-2*n-1)
    if rho < 2.5^(16/n)
        fprintf('Bernstein radius: %g\n', rho)
        disp('-- Special quadrature probably needed')
    end
    
    % Compute special quadrature weights
    [w1, w3, w5] = rsqrt_pow_weights(tj, troot); % O(n^2) per point
    
    % Eval special quadrature
    Q1spec = sum(w1.*f1j);
    Q3spec = sum(w3.*f3j);
    Q5spec = sum(w5.*f5j);        
    
    % Compute errors
    R1n = abs(I1-Q1n);
    R3n = abs(I3-Q3n);
    R5n = abs(I5-Q5n);
    R1spec = abs(I1-Q1spec);
    R3spec = abs(I3-Q3spec);
    R5spec = abs(I5-Q5spec);
    
    d1.n = n;
    d1.troot = troot;
    d1.I1= I1;
    d1.R1_direct = R1n / abs(I1);
    d1.R1_spec = R1spec / abs(I1);

    d3.n = n;
    d3.troot = troot;
    d3.I3= I3;
    d3.R3_direct = R3n / abs(I3);
    d3.R3_spec = R3spec / abs(I3);

    d5.n = n;
    d5.troot = troot;
    d5.I5= I5;    
    d5.R5_direct = R5n  / abs(I5);
    d5.R5_spec = R5spec / abs(I5);

    disp(' ')
    disp(struct2table(d1))
    fprintf('\n\n')
    disp(struct2table(d3))
    fprintf('\n\n')    
    disp(struct2table(d5))
    
    % Plot problem
    clf
    t = linspace(-1,1);
    plot3(x(t), y(t), z(t))
    hold on 
    plot3(x(tj), y(tj), z(tj),'.b')   
    plot3(x0, y0, z0, '*r')
    axis equal
end

