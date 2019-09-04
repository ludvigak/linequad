function curvedfibers_field()
% Test close eval of kernels 1/R, 1/R^3, and 1/R^5 close to 3D fiber
% Fiber parametrized in t from -1 to 1
% Ludvig af Klinteberg and Alex Barnett, May 2018   
    
    % Number of points on panel
    nquad = 16;
    % Critical Bernstein radius (sufficent for R5)
    rho = 4^(16/nquad);
    
    %curve = 'helix'
    curve = 'parabola'
       
    % Parametrization of fibre
    switch curve
      case 'helix'
        % Part of helix
        r = 1;
        k = pi/4;
        c = 2;
        scale = 1; % To check invariance
        x = @(t) r*cos(k*t)*scale;
        y = @(t) r*sin(k*t)*scale;
        z = @(t) c*t*scale;       
        xp = @(t) -k*r*sin(k*t)*scale;
        yp = @(t) k*r*cos(k*t)*scale;
        zp = @(t) c*scale*t.^0;  
      case 'parabola'            
        k = 0.25;
        x = @(t) t;
        y = @(t) k*t.^2;
        z = @(t) t;
        xp = @(t) t.^0;
        yp = @(t) 2*k*t;
        zp = @(t) t.^0;
    end    
    dL = @(t) sqrt(xp(t).^2 + yp(t).^2 + zp(t).^2);    
    
    % Grid to compute potentials on
    Ngrid = 50;
    dgrid = linspace(-2, 2, Ngrid);
    tgrid = linspace(-2,2, Ngrid);
    [T,D] = meshgrid(tgrid,dgrid);
    X = x(T);
    Y = y(T) + D.*zp(T);
    Z = z(T) - D.*yp(T);

    % Density to integrate
    density = @(t) x(t) + 10;
    
    % === END PARAMS
    
    % - Integrands with 1/R^p kernels, p=1,3,5
    % - These are the kernels needed for the stokeslet and the doublet
    % of slender body theory.
    %
    % Very convoluted way of expressing things, but flexible when testing.
    rx = @(t,x0) x(t)-x0;
    ry = @(t,y0) y(t)-y0;
    rz = @(t,z0) z(t)-z0;
    R = @(t,x0,y0,z0) sqrt(rx(t,x0).^2 + ry(t,y0).^2 + rz(t,z0).^2);
    f = @(t) density(t).*dL(t);
    integrand1 = @(t,x0,y0,z0) f(t)./ R(t,x0,y0,z0);
    integrand3 = @(t,x0,y0,z0) f(t)./ R(t,x0,y0,z0).^3;
    integrand5 = @(t,x0,y0,z0) f(t)./ R(t,x0,y0,z0).^5;           

    % Discretization of panel
    [tj,wj] = legendre.gauss(nquad);
    xj = x(tj);
    yj = y(tj);
    zj = z(tj);
    
    % Storage for results
    [I1, I3, I5, Q1, Q3, Q5] = deal(zeros(size(X)));
    
    % Compute and apply weights
    fprintf('* Computing quadrature weights: ')
    tic    
    [all_w1, all_w3, all_w5] = line3_near_weights(tj, wj, xj, yj, zj, X, Y, Z, rho);
    toc
    for i=1:numel(X)    
        x0 = X(i);
        y0 = Y(i);
        z0 = Z(i);
        f1j = integrand1(tj,x0,y0,z0);
        f3j = integrand3(tj,x0,y0,z0);
        f5j = integrand5(tj,x0,y0,z0);
        w1 = all_w1(:,i);
        w3 = all_w3(:,i);
        w5 = all_w5(:,i);            
        Q1(i) = sum(w1.*f1j);
        Q3(i) = sum(w3.*f3j);
        Q5(i) = sum(w5.*f5j);
    end 
    
    % Compute reference with adaptive quadrature    
    fprintf('* Computing reference: ')
    tic
    parfor i=1:numel(X)    
        x0 = X(i);
        y0 = Y(i);
        z0 = Z(i);
        warning off
        I1(i) = integral(@(t) integrand1(t,x0,y0,z0), -1, 1, 'reltol', eps, 'abstol', eps);
        I3(i) = integral(@(t) integrand3(t,x0,y0,z0), -1, 1, 'reltol', eps, 'abstol', eps);
        I5(i) = integral(@(t) integrand5(t,x0,y0,z0), -1, 1, 'reltol', eps, 'abstol', eps);
        warning on
    end
    toc        
    
    % Compute errors
    R1 = abs(I1-Q1) ./ abs(I1) + eps();
    R3 = abs(I3-Q3) ./ abs(I3) + eps();
    R5 = abs(I5-Q5) ./ abs(I5) + eps();
    
    % Print
    rms = @(v) sqrt( mean(v(:).^2) );    
    fprintf('R1 relerr max=%g, rms=%g\n', norm(R1, inf), rms(R1));
    fprintf('R3 relerr max=%g, rms=%g\n', norm(R3, inf), rms(R3));
    fprintf('R5 relerr max=%g, rms=%g\n', norm(R5, inf), rms(R5));
    
    % Plot
    clf();
    function ploterr(R)    
        surf(X, Y, Z, log10(R));
        hold on;
        plot3(xj, yj, zj, '.-k')
        axis off
        axis image
        hcb = colorbar()    ;
        xlabel(hcb, 'log10(err)')
        shading interp
        %caxis([-15 -12])
        switch curve
          case 'helix'        
            view(-145,-1.8)
          case 'parabola'
            view(-55, 30)
        end
    end    
    subplot(1, 3, 1);ploterr(R1); title('p=1'); 
    subplot(1, 3, 2);ploterr(R3); title('p=3'); 
    subplot(1, 3, 3);ploterr(R5); title('p=5'); 
end

