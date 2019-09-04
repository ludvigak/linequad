function [w1, w3, w5, specquad_needed] = line3_near_weights(tj, wj, xj, yj, zj, X, Y, Z, varargin)
%
% [w1,w3,w5] = line3_near_weights(tj, zj, yj, zj, X, Y, Z, [rho])
%
% Near eval quadrature weights for kernel 1/R^p, p=1,3,5
%
% INPUT:
% tj: Nodes in parametrization space [-1, 1] (probably Gauss-Legendre points)
% wj: Quadrature weights accompanying tj
% xj, yj, zj: Points on curve, corresponding to tj
% X, Y, Z: List of target points
% rho: Critical Bernstein radius, weights computed for roots inside it.
%      Defaults to 4^(16/n).
%
% n = length of tj,wj,xj,yj,zj
%
% OUTPUT:
% w1,w3,w5: Vectors of size (n, numel(X)) with quadratures weights,
%           replacing wj if needed
% specquad_needed: List of bool values indicating if special quadrature weights
%                  were computed

USE_MEX_CODE = true; % Switch off to debug codes
CHECK_SCHWARZ = false; % Additional checks if nearby Schwarz point
                       % (experimental)

% Reshape all inputs
tj = tj(:);
wj = wj(:);
xj = xj(:);
yj = yj(:);
zj = zj(:);

% Legendre expansions of discretization
n = length(xj);
Legmat = legendre.matrix(n); % O(n^2) per panel
xhat = Legmat*xj;
yhat = Legmat*yj;
zhat = Legmat*zj;

% Cap expansion at 16 coefficients.
% This seems to be more stable near boundary
if n > 16
    xhat = xhat(1:16);
    yhat = yhat(1:16);
    zhat = zhat(1:16);
end
n_expa = numel(xhat);

% Set rho if not passed
if isempty(varargin)
    rho = 4^(16/n); % Sufficient for R5
else
    rho = varargin{1};
end

if CHECK_SCHWARZ
    % EXPERIMENTAL: Check where nearest Schwarz singularity is
    % Unclear whether this is actually needed
    [~, Dmat] = legendre.deriv_vec(n-1, tj);
    diffmat = Dmat*Legmat;
    % Expansion of |g'|^2
    gp2hat = Legmat*( (diffmat*xj).^2 + (diffmat*yj).^2 + (diffmat*zj).^2 );
    gp2hat = gp2hat(1:n_expa);
    % Go all in: Comrade matrix + Newton polish (this is overkill)
    ss_all = legendre.roots(gp2hat);
    ss = ss_all(1);
    ss = rootfinder(legendre.diff(xhat), legendre.diff(yhat), legendre.diff(zhat),0,0,0, ss);
    % % Sanity check
    % [~, Dss] = legendre.deriv_scalar(n_expa-1, ss);
    % check = (Dss*xhat)^2 + (Dss*yhat)^2 + (Dss*zhat)^2
    nearby_schwarz = bernstein_radius(ss) < 1.1*rho;
end

%% Rootfinding: initial guesses
if USE_MEX_CODE
    all_tinits = rootfinder_initial_guess_mex(tj, xj, yj, zj, X, Y, Z); % O(n) per point
else
    all_tinits = complex(zeros(size(X)));
    % Standard guess
    for i=1:numel(X)    
        tinit = rootfinder_initial_guess(tj, xj, yj, zj, X(i), Y(i), Z(i)); % O(n) per point
        all_tinits(i) = tinit;
    end
end

% First filter: Don't check points whose initial guess is far away
cp = (bernstein_radius(all_tinits) < 1.5*rho); % cp: Compute Points

if CHECK_SCHWARZ && nearby_schwarz
    % Let initial guess be results of matrix-based rootfinding
    % Newton will polish these roots by a digit or two
    % TODO: If keeping this feature, implement in C (with call to LAPACK/DHSEQR for eigenvalues).
    % TODO: add cp filtering here
    g2hat = Legmat*(xj.^2 + yj.^2 + zj.^2); % Expansion of |g|^2
    g2hat = g2hat(1:n_expa);
    for i=1:numel(X)    
        % Form expansion of |g-x|^2 and compute roots
        coeff = g2hat - 2*(X(i)*xhat+Y(i)*yhat+Z(i)*zhat);
        coeff(1) = coeff(1) + (X(i)^2+Y(i)^2+Z(i)^2);
        tall = legendre.roots(coeff(1:n_expa)); % (O(n^3) per point)
        all_tinits(i) = tall(1);
    end    
end

%% Rootfinding: Run
all_roots = deal(complex(zeros(size(X))));
rootfinder_converged = false(size(X));
if USE_MEX_CODE
    [all_roots(cp), rootfinder_converged(cp)] = rootfinder_mex(xhat, yhat, zhat, X(cp), Y(cp), Z(cp), all_tinits(cp));
else
    for i=find(cp(:)')
        tinit = all_tinits(i);
        [troot, converged] = rootfinder(xhat, yhat, zhat,  X(i), Y(i), Z(i), tinit); % O(n) per point        
        all_roots(i) = troot;
        rootfinder_converged(i) = converged;    
    end
end

% Check which points need special attention
all_bernstein_radii = bernstein_radius(all_roots);
specquad_needed = rootfinder_converged & (all_bernstein_radii < rho);

% Default is to return regular weights
[w1, w3, w5] = deal(repmat(wj(:), 1, numel(X)));

%% Compute special weights where needed
if USE_MEX_CODE
    [tmp1, tmp3, tmp5] = rsqrt_pow_weights_mex(tj, all_roots(specquad_needed));
    w1(:, specquad_needed) = tmp1;
    w3(:, specquad_needed) = tmp3;
    w5(:, specquad_needed) = tmp5;
else
    for i=1:numel(X)
        if specquad_needed(i)
            [w1(:,i), w3(:,i), w5(:,i)] = ...
                rsqrt_pow_weights(tj, all_roots(i)); 
        end
    end
end
