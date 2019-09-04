function B = bclag_interp_matrix(x, w, xi)

% Barycentric Lagrange Interpolation.
% Berrut, J.-P., & Trefethen, L. N. (2004).
% SIAM Review, 46(3), 501–517. doi:10.1137/S0036144502417715

assert(size(x,2)==1)
assert(size(xi,2)==1)
assert(all(size(x)==size(w)));

n = numel(x);
N = numel(xi);

B = zeros(N,n);
[denom exact] = deal(zeros(size(xi)));

for j=1:n
    xdiff = xi-x(j);
    temp = w(j)./xdiff;
    B(:,j) = temp;
    denom = denom + temp;
    exact(xdiff==0) = j;
end

B = bsxfun(@rdivide,B,denom);
jj = find(exact);
B(jj,:) = 0;
B(jj + N*(exact(jj)-1)) = 1;
