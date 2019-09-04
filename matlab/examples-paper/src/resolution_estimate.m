function est = resolution_estimate(fj)

n = numel(fj);
L = legendre_matrix(n);
c = abs(L*fj);

est = max(c(end-1), c(end)) / max(abs(fj));
