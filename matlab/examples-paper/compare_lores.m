% Show difference between 2D quadratures for layer potential in starfish,
% when just resolving enough.
%
% Slow code, very demo-specific
function compare_lores(varargin)

if nargin==0
    save_paper_plots = false;
else
    save_paper_plots = varargin{1};
end

% Load test problem
[curve, fexact] = test_problem();

% Boundary discretization
opt.order = 16;
opt.tol = 1e-6; % Resolution
opt.R = 3; % Bernstein radius for eval

g = adaptive_composite_discretization(curve, opt);
disp(['Num panels: ', num2str(numel(g.t_edges)-1)])

if save_paper_plots
    Ngrid = 300;
else
    Ngrid = 100;
end

% Volume grid
R = max(abs(g.z));
x = linspace(-R,R,Ngrid);
[X,Y] = meshgrid(x,x);
Z = X+1i*Y;

% Assemble matrix and solve directly
tic
A = full_matrix(g);
rhs = fexact(g.z);
den = A\rhs;
toc

tic
REF = fexact(Z);
[Urealn, Ureal2n, Ucpxn, Ucpx2n, exterior] = dlp_eval(g, den, Z, opt);
toc

% Hide points in exterior
REF(exterior) = NaN;
int = ~exterior;

% Comparison
results = {Urealn, Ucpxn, Ureal2n, Ucpx2n};
names = {'Real, 16 pts', 'Complex, 16 pts', 'Real, 32 pts', 'Complex, 32 pts'};
s = struct();
for i=1:4
    sfigure(i); clf(); publication_fig();
    U = results{i};
    pcolor(X, Y, log10(abs(REF-U)/norm(REF(int), inf)));
    hold on
    plot(g.z, '-k')
    plot_panel_edges(g, 0.05, '-k')
    axis off
    axis image
    shading flat
    xlabel(colorbar(), 'log_{10} E_{rel}')
    caxis([-16, -8])    
    title(names{i})    
    maxrelerr = norm(REF(int)-U(int), inf)/norm(REF(int),inf);    
    text(min(real(g.z)), min(imag(g.z))-0.2, ['$\max E_{rel}=$', sprintf('%.1e',maxrelerr)],'interpreter','latex')
    
    %xlabel(['$\max E_{rel}$'])
    
    s(i).method = names{i};
    s(i).err = maxrelerr;
    drawnow
    
end
disp(struct2table(s))


function saveplot(num, filename)
    sfigure(num);
    publication_fig();
    path = ['../fig/', filename];
    print('-dpng', '-r200', path)
    system(['mogrify -trim ' path])    
    disp(['Wrote ' path])
end

if save_paper_plots
    saveplot(1, 'laplace_lores_real16.png')
    saveplot(2, 'laplace_lores_complex16.png')
    saveplot(3, 'laplace_lores_real32.png')
    saveplot(4, 'laplace_lores_complex32.png')    
end
end