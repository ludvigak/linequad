function test_mex()
% Test that MEX and MATLAB codes do the same thing
%
% Note that quadrature weights cannot be compared directly, since they are unstable,
% and the MEX code is compiled with full optimizations. Instead the quadrature results
% are compared.
%
    start = tic();
    rng(1)
    straight_line()   
    curve()
    disp('')
    disp('* ALL TESTS PASSED')    
    toc(start)
end

function straight_line()    
% Compare quadratures on straight line, essentially comparing recursions
    n = 16;
    % In zones of power series evaluations
    err        = compare_quadratures(0.1, 1e-5, n);
    err(end+1) = compare_quadratures(1.2, 1e-5, n);
    err(end+1) = compare_quadratures(1.2, 1e-4, n);
    err(end+1) = compare_quadratures(1.2, 1e-3, n);
    % Special cases
    err(end+1) = compare_quadratures(0, 1, n);        
    err(end+1) = compare_quadratures(0, 0.5, n);        
    err(end+1) = compare_quadratures(0, 1e-10, n);    
    err(end+1) = compare_quadratures(-1, 1e-16, n);
    err(end+1) = compare_quadratures(1, 1e-16, n);
    err(end+1) = compare_quadratures(1.01, 0, n);    
    err(end+1) = compare_quadratures(1.1, 0, n);    
    err(end+1) = compare_quadratures(1.2, 0, n);        
    % Random points
    for i=1:10
        zr = 1.5*(1-2*rand());
        zi = 1.0*(1-2*rand());        
        err(end+1) = compare_quadratures(zr, zi, n);
    end    
    summary.reldiff_1 = max([err.reldiff_1]);
    summary.reldiff_3 = max([err.reldiff_3]);
    summary.reldiff_5 = max([err.reldiff_5]);
    disp('* Comparison MATLAB/MEX weights')    
    format short g
    disp(struct2table(err))
    disp(' maximums:')
    disp(struct2table(summary))
    assert(summary.reldiff_1 < 1e-12)
    assert(summary.reldiff_3 < 1e-12)
    assert(summary.reldiff_5 < 1e-12)
    fprintf('=passed\n\n')

    function err = compare_quadratures(zr, zi, n)
        z = complex(zr+1i*zi);
        tj = legendre.gauss(n);
        % Random polynomials with prescribed pole
        coeff = rand(1,n);
        f1 = polyval(coeff, tj) ./ abs(tj-z);
        f3 = polyval(coeff, tj) ./ abs(tj-z).^3;
        f5 = polyval(coeff, tj) ./ abs(tj-z).^5;
        % Compute quad weights using different codes
        [W1, W3, W5] = rsqrt_pow_weights(tj, z);
        [w1, w3, w5] = rsqrt_pow_weights_mex(tj, z);
        % Compare quadrature results
        quaderr = @(f,w,W) abs( sum(f.*(w-W)) / sum(f.*W) );
        err.z = z;
        err.reldiff_1 = quaderr(f1,w1,W1);
        err.reldiff_3 = quaderr(f3,w3,W3);
        err.reldiff_5 = quaderr(f5,w5,W5);
    end
end

function curve()
% Compare codes on helix shaped curve
    r = 1;
    k = pi/3;
    c = 2;
    x = @(t) r*cos(k*t);
    y = @(t) r*sin(k*t);
    z = @(t) c*t;       
    % Discretization
    n = 16;
    tj = legendre.gauss(n);
    xj = x(tj);
    yj = y(tj);
    zj = z(tj);        
    % Legendre expansions of discretization
    Legmat = legendre.matrix(n);
    xhat = Legmat*xj;
    yhat = Legmat*yj;
    zhat = Legmat*zj;
    
    % Cloud of points around curve
    Npt = 10000;
    trand = 1-2*rand(Npt,1);
    dist = 3;
    X = dist*(1-2*rand(Npt,1)) + x(trand);
    Y = dist*(1-2*rand(Npt,1)) + y(trand);
    Z = dist*(1-2*rand(Npt,1)) + z(trand);
    %clf; plot3(xj,yj,zj, 'k',X,Y,Z, '.'); axis equal
    
    %% Root finding test
    disp(' * Test root finding')    
    % Run root finding
    roots_converged = false(size(X));
    all_roots = zeros(size(X));    
    all_tinits = complex(zeros(size(X))); 
    tic
    for i=1:numel(X)    
        x0 = X(i);
        y0 = Y(i);
        z0 = Z(i);
        % Find root for this target point
        tinit = rootfinder_initial_guess(tj, xj, yj, zj, x0, y0, z0); % O(n) per point
        [troot, converged] = rootfinder(xhat, yhat, zhat, x0, y0, z0, tinit); % O(n) per point
        all_tinits(i) = tinit;
        all_roots(i) = troot;
        roots_converged(i) = converged;        
    end
    time_roots = toc();
    % Mex interface for roots
    tic
    all_tinits_mex = rootfinder_initial_guess_mex(tj, xj, yj, zj, X, Y, Z);
    [all_roots_mex, roots_converged_mex] = rootfinder_mex(xhat, yhat, zhat, X, Y, Z, all_tinits_mex);    
    time_roots_mex = toc();
    % Print some runtimes    
    fprintf('MATLAB roots per second: %g\n', Npt/time_roots);
    fprintf('MEX    roots per second: %g\n', Npt/time_roots_mex);                
    % Compare roots
    all_passed = true;
    max_root_diff = 0;
    for i=1:numel(X)  
        root_diff = abs(all_roots(i) - all_roots_mex(i))/abs(all_roots(i));
        if roots_converged(i)
            max_root_diff = max(max_root_diff, root_diff);
        end
        if (roots_converged(i) ~= roots_converged_mex(i)) || (roots_converged(i) && root_diff > 1e-14)
            all_passed = false;
            fprintf('Rootfinder difference at point (%g, %g, %g)\n', X(i), Y(i), Z(i));
            fprintf(' * MATLAB root=%g,\t converged=%d\n', all_roots(i), roots_converged(i));
            fprintf(' * MEX    root=%g,\t converged=%d\n', all_roots_mex(i), roots_converged_mex(i));
            fprintf(' root_diff=%e\n', abs(all_roots_mex(i)-all_roots(i)))            
        end        
    end
    fprintf('Max root diff=%g\n', max_root_diff);
    assert(all_passed, 'Different results from root finders')
    fprintf('=passed\n\n')
    
    %% Quadrature test
    disp(' * Test weights')
    % Compute weights and compare for all converged roots
    idx = find(roots_converged);
    tic
    [w1_mex, w3_mex, w5_mex] = rsqrt_pow_weights_mex(tj, all_roots(roots_converged));
    time_weights_mex = toc();
    [w1, w3, w5] = deal(zeros(size(w1_mex)));
    tic
    for i=1:length(idx)
        [w1(:, i), w3(:, i), w5(:, i)] = rsqrt_pow_weights(tj, all_roots(idx(i)));
    end
    time_weights = toc();
    % Print metrics
    weights_computed = length(idx);
    fprintf('MATLAB weight sets per second: %g\n', weights_computed/time_weights);
    fprintf('MEX    weight sets per second: %g\n\n', weights_computed/time_weights_mex);    
    fprintf('MATLAB targets per second: %g\n', weights_computed/(time_roots+time_weights));    
    fprintf('MEX    targets per second: %g\n\n', weights_computed/(time_roots_mex+time_weights_mex));    
    
    % Compare relative differences in quadrature results
    bernstein_radius = @(z) abs(z + sign(real(z)).*sqrt(z.^2 - 1));
    all_passed = true;
    results.max_reldiff1 = 0;
    results.max_reldiff3 = 0;
    results.max_reldiff5 = 0;    
    for j=1:length(idx)
        i = idx(j);
        troot = all_roots(i);
        rho = bernstein_radius(troot);
        % Only compare quadratures withing a certain Bernstein radius
        if rho < 4
            Rj = sqrt( (xj-X(i)).^2 + (yj-Y(i)).^2 + (zj-Z(i)).^2 );        
            D1 = sum((w1(:, j)-w1_mex(:, j))./Rj.^1) / sum(w1(:, j)./Rj.^1);
            D3 = sum((w3(:, j)-w3_mex(:, j))./Rj.^3) / sum(w3(:, j)./Rj.^3);
            D5 = sum((w5(:, j)-w5_mex(:, j))./Rj.^5) / sum(w5(:, j)./Rj.^5);
            results.max_reldiff1 = max(results.max_reldiff1, D1);
            results.max_reldiff3 = max(results.max_reldiff3, D3);
            results.max_reldiff5 = max(results.max_reldiff5, D5);
            if ~(D1 < 1e-14 && D3 < 1e-14 && D5 < 1e-14)
                disp('FAIL')
                i
                rho
                D1
                D3
                D5            
                all_passed = false;
            end
        end
    end
    assert(all_passed, 'Quadratures differed')
    disp('Relative differences between quadrature results')
    disp(struct2table(results))
    fprintf('=passed\n\n')
end