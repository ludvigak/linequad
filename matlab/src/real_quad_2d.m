function w = real_quad_2d(tleg, z)
% Interpolatory quadrature using real basis
    
    N = numel(tleg);
    A = fliplr(vander(tleg));
    [~, p] = complex_integrals(z, N);
    warning('off', 'MATLAB:nearlySingularMatrix')
    w = (A.'\p) .* (tleg-z);
    warning('on', 'MATLAB:nearlySingularMatrix')                    
end
