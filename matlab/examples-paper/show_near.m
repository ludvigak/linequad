% Show full accuracy results close up for layer potential in starfish
%
% Slow code, very demo-specific
function show_near(varargin)

if nargin==0
    save_paper_plots = false;
else
    save_paper_plots = varargin{1};
end

% Load test problem
[curve, fexact] = test_problem();

% Boundary discretization
opt.order = 16;
opt.tol = 1e-14;
opt.R = 3;
g = adaptive_composite_discretization(curve, opt);
disp(['Num panels: ', num2str(numel(g.t_edges)-1)])


% Assemble matrix and solve directly
tic
A = full_matrix(g);
rhs = fexact(g.z);
den = A\rhs;
toc

if save_paper_plots
    Ngrid = 250;
else
    Ngrid = 50;
end

% Volume grid
R = max(abs(g.z));
x = linspace(-R,R,Ngrid);
[X,Y] = meshgrid(x,x);
Z = X+1i*Y;

% Slice grid
dt = 0.1*pi;
t0 = pi*1.66;
tslice = linspace(t0, t0+dt, Ngrid);
dslice = linspace(0.001, 0.15, Ngrid);
[T,D] = meshgrid(tslice, dslice);
Zslice = curve.tau(T+1i*D);
slicebdry = [flipud(Zslice(:, 1)); Zslice(1,:).' ; Zslice(:, end); flipud(Zslice(end,:).')];

% Log-scaled slice grid
tslice = linspace(t0, t0+dt, Ngrid);
dslice = logspace(-8, log10(dslice(end)), Ngrid);
[Tlog,Dlog] = meshgrid(tslice, dslice);
Zlog = curve.tau(Tlog+1i*Dlog);
    

%% Plot full grid

REF = fexact(Z);
tic
[Urealn, Ureal2n, Ucpxn, Ucpx2n, exterior] = dlp_eval(g, den, Z, opt);
toc
% Hide points in exterior
REF(exterior) = NaN;
int = ~exterior;
refnorm = norm(REF(int), inf);

% Comparison
results = {Urealn, Ucpxn, Ureal2n, Ucpx2n};
names = {'SSQ, 16 pts', 'HO, 16 pts', 'SSQ, 32 pts', 'HO, 32 pts'};
s = struct();
for i=1:4
    sfigure(i); clf; publication_fig();
    U = results{i};
    maxrelerr = norm(REF(int)-U(int), inf)/refnorm;        
    pcolor(X, Y, log10(abs(REF-U)/refnorm));
    hold on
    plot(g.z, '-k')    
    plot(slicebdry, '-r', 'LineWidth', 1)    
    
    axis off
    axis image
    shading flat
    xlabel(colorbar(), 'log_{10} E_{rel}')
    caxis([-16, -14])
    title(names{i})    

    text(min(real(g.z)), min(imag(g.z))-0.2, ['$\max E_{rel}=$', sprintf('%.1e',maxrelerr)],'interpreter','latex')    
    s(i).method = names{i};
    s(i).err = maxrelerr;
    drawnow
end
disp('FULL')
disp(struct2table(s))

%% Plot Slice grid
Z = Zslice;
X = real(Zslice);
Y = imag(Zslice);
REF = fexact(Z);
tic
[Urealn, Ureal2n, Ucpxn, Ucpx2n, exterior] = dlp_eval(g, den, Z, opt);
toc
results = {Urealn, Ucpxn, Ureal2n, Ucpx2n};
s = struct();
for i=1:4
    sfigure(10+i); clf; publication_fig();
    U = results{i};
    maxrelerr = norm(REF-U, inf)/refnorm;    
    s(i).method = names{i};
    s(i).err = maxrelerr;
    
    pcolor(X, Y, log10(abs(REF-U)/refnorm));
    hold on
    plot(slicebdry, '-r', 'linewidth',1)
    plot(g.z, '-k')
    plot_panel_edges(g, 0.01, '-k')    
    text(min(X(:)), min(Y(:))-0.06, ['$\max E_{rel}=$', sprintf('%.1e',maxrelerr)],'interpreter','latex')        
    axis off
    axis equal
    
    ax = [min(X(:)) max(X(:)) min(Y(:)) max(Y(:))];
    ax = ax + 0.02*[-1 1 -1 1].*abs(ax);
    axis(ax)
    
    shading flat
    xlabel(colorbar(), 'log_{10} E_{rel}')
    caxis([-16, -12])

    title(names{i})    
    drawnow
end
disp('SLICE')
disp(struct2table(s))

%% Plot log-scaled slice grid
Z = Zlog;
REF = fexact(Z);
tic
[Urealn, Ureal2n, Ucpxn, Ucpx2n, exterior] = dlp_eval(g, den, Z, opt);
toc
results = {Urealn, Ucpxn, Ureal2n, Ucpx2n};
s = struct();
for i=1:4
    sfigure(20+i); clf; publication_fig();
    U = results{i};
    maxrelerr = norm(REF-U, inf)/refnorm;        
    pcolor(Tlog, log10(Dlog), log10(abs(REF-U)/refnorm));
    hold on
    xmin = t0;
    xmax = t0+dt;
    ymin = min(log10(Dlog(:)));
    ymax = max(log10(Dlog(:)));    
    e = ones(size(g.t_edges));
    plot( [g.t_edges g.t_edges]', [e*ymin e*(ymin+1/2)]', '-k')
    plot([xmin xmax xmax xmin xmin], ...
         [ymin ymin ymax ymax ymin], '-r')
    
    xlim([xmin, xmax])    
    xlabel('Re t')
    
    ylabel('log_{10}(Im t)')
    shading flat
    xlabel(colorbar(), 'log_{10} E_{rel}')
    caxis([-16, -12]);
    title(names{i})    
    s(i).method = names{i};
    s(i).err = maxrelerr;

    publication_fig
    text(t0+dt*0.75, ymin-0.9, ['$\max E_{rel}=$', sprintf('%.1e',maxrelerr)],'interpreter','latex')    
    
    
    drawnow
end
disp('LOG-SLICE')
disp(struct2table(s))

function saveplot(num, filename)
    sfigure(num);
    publication_fig
    path = ['../fig/', filename];
    print('-dpng', '-r200', path)
    system(['mogrify -trim ' path])
    disp(['Wrote ' path])
end

if save_paper_plots
    saveplot(1, 'laplace_domain_real16.png')
    saveplot(11, 'laplace_slice_real16.png')
    saveplot(21, 'laplace_zoom_real16.png')
        
    saveplot(2, 'laplace_domain_complex16.png')
    saveplot(12, 'laplace_slice_complex16.png')
    saveplot(22, 'laplace_zoom_complex16.png')

    saveplot(3, 'laplace_domain_real32.png')
    saveplot(13, 'laplace_slice_real32.png')
    saveplot(23, 'laplace_zoom_real32.png')
        
    saveplot(4, 'laplace_domain_complex32.png')
    saveplot(14, 'laplace_slice_complex32.png')
    saveplot(24, 'laplace_zoom_complex32.png')
    
    
end


end