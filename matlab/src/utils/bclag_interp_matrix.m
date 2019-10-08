function B = bclag_interp_matrix(x, xi)
% B = bclag_interp_matrix(x, xi)
%
% Barycentric Lagrange Interpolation.
% Berrut, J.-P., & Trefethen, L. N. (2004).
% SIAM Review, 46(3), 501–517. doi:10.1137/S0036144502417715
    assert(size(x,2)==1)
    assert(size(xi,2)==1)
    n = numel(x);
    N = numel(xi);
    w = bclag_interp_weights(x);
    B = zeros(N,n);
    [denom, exact] = deal(zeros(size(xi)));
    for j=1:n
        for i=1:N
            xdiff = xi(i)-x(j);
            temp = w(j)/xdiff;
            B(i,j) = temp;
            denom(i) = denom(i) + temp;
            if xdiff==0
                exact(i) = j;
            end
        end
    end
    B = bsxfun(@rdivide,B,denom);
    for i=1:N
        if exact(i)
            B(i,:) = 0;
            B(i + N*(exact(i)-1)) = 1;
        end
    end
end