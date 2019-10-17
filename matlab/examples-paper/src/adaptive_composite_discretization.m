function g = adaptive_composite_discretization(curve, opt)

Nmax = 1e4;
order = opt.order; 
tol = opt.tol;

if isfield(opt, 'R')
    R = opt.R;
else
    R = 0;
end

[tgl, wgl] = lgwt2(order,-1,1);
t_edges = [0]; % Grows dynamically
dt = 2*pi/2; % Start with ridiculously large element
L = legendre_matrix(order);
t_edges = [0; 2*pi]; % Grows dynamically
all_pass = false;
while ~all_pass
    all_pass = true;
    for i=1:numel(t_edges)-1
        % Check interval
        ta = t_edges(i);
        tb = t_edges(i+1);
        dt = tb-ta;
        % Discretize and check if panel resolved
        tj = ta + dt/2*(tgl+1);
        zpj = curve.dtau(tj);
        % est = resolution_estimate(abs(zpj));
        est = resolution_estimate(zpj);
        resolution_ok = (est < tol);
        % Check length compared to neighbors
        if i==1
            ldt = t_edges(end)-t_edges(end-1);
        else
            ldt = ta - t_edges(i-1);
        end
        if i==numel(t_edges)-1
            rdt = t_edges(2)-t_edges(1);
        else
            rdt = t_edges(i+2)-tb;
        end
        length_ok = (dt <= ldt*2.01 && dt <= rdt*2.01);
        if resolution_ok && length_ok
            % Pass
        else
            all_pass = false;
            % Bisect element
            new_edge = ta+dt/2;
            t_edges = [t_edges(1:i); new_edge; t_edges(i+1:end)]; 
            break
        end    
    end
    if numel(t_edges)>Nmax
        error('Discretization did not converge')
    end
end

numpanels = length(t_edges)-1;
N = numpanels*order;
[t,w] = deal(zeros(N, 1));
nearby_schwarz = false(numpanels, 1);
for i=1:numpanels
    ta = t_edges(i);
    tb = t_edges(i+1);
    j = (1:order) + order*(i-1);
    tj = ta + (tb-ta)/2*(tgl+1);
    t(j) = tj;
    w(j) = (tb-ta)/2*wgl;    
    % Check if there is a nearby Schwarz singularity
    % on the business side of the panel
    zpj = curve.dtau(tj);    
    coeff = L*zpj;
    s = solve_legendre(coeff, 0, 0, 20, 1e-13);
    rhos = bernstein_rad(s);
    nearby_schwarz(i) = rhos < R && imag(s) > 0;
end

g.t_edges = t_edges;
g.z_edges = curve.tau(t_edges);
g.n_edges = curve.normal(t_edges);

g.t = t;
g.z = curve.tau(t);
g.zp = curve.dtau(t);
g.zpp = curve.d2tau(t);
g.n = curve.normal(t);
g.w = w;
g.dz = g.zp.*w;
g.ds = abs(g.dz);

g.N = N;
g.order = order;
g.numpanels = numpanels;

g.nearby_schwarz = nearby_schwarz;
