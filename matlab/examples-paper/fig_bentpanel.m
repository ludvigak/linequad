% test err for poly interp on complex curved fiber in 2D,
% to disentangle that error from the close-eval error.
% Notation: t = param in [-1,1]; tau = gamma(t) = complex location.
% Barnett 6/4/19. based on err_cmplx_interp.m, March 2019.

clear; init
freq = 3; dens = @(t) sin(1+freq*t);  % smooth func of param to interp; density

sfigure(1); clf
k = 0.6;       % makes schw hit u contour 0.2 (0.19915..)
g = @(t) t + 1i*k*t.^2;        % panel parameterization gamma
gp = @(t) 1 + 1i*2*k*t;        % gamma'
schw = 1i/(4*k);    % z-plane Schwarz singularity of unscaled parabola panel
tschw = 1i/(2*k);    % t-plane Schwarz singularity of unscaled parabola panel
rschw = exp(asinh(tschw/1i))   % = exp(panel_conformal_dist(@(t) t,tschw))   :)

p = 16;                 % number of nodes on panel, ie n.
[tj,wj] = lgwt(p,-1,1); tj=tj(:);   % G-L: col of interp pts in param
[X Y] = meshgrid(-2:0.2:2,-2:.2:1.0); Z = X+1i*Y;
plot(Z,'-','color',.8*[1 1 1]); hold on;
plot(Z.','-','color',.8*[1 1 1]);
plot([-1 1],[0 0],'b-');
plot(tj,0*tj,'k.','markersize',10);
ttarg = -0.5 + 0.3i; ztarg = g(ttarg);   % example targ to show
rho = abs(ttarg+sqrt(ttarg-1).*sqrt(ttarg+1));   % Trefethen Spec Meth book Ch.5
N=1e2; tbern = cos(1i*log(rho) + 2*pi*(0:N)/N);     % Bernstein ellipse pts
plot(tbern,':','linewidth',1,'color',[.7 0 0]);
plot(ttarg,'.','markersize',15,'color',[0 .6 0]);    % pre-image of target
plot(tschw,'r*','markersize',10);
axis equal tight; axis([-1.2 1.2 -.6 0.9]);
text(-1.5,1.0,'(a) complex parameter plane and preimages','interpreter','latex')
set(gcf,'paperposition',[0 0 3 2]);
printpng(1,'../fig/tplane')


sfigure(2); clf
% reference eval (from compare_2d_complex_real.m)
kernel_cpx = @(x, y, ny) ( ny ./ (x-y) );    
kernel = @(x, y, ny) real( kernel_cpx(x, y, ny)  );
integrand = @(t,z) kernel(z, g(t), 1i*gp(t)./abs(gp(t))) .* dens(t) .* ...
    abs(gp(t));
eps=1e-14;
reference = @(z) integral(@(t) integrand(t, z), -1, 1, 'abstol', eps, 'reltol', eps);
dx=0.03; gx =-1.7:dx:1.7; gy = -1:dx:1.2;  % plot grid
[xx yy] = meshgrid(gx,gy); zz = xx + 1i*yy;
u = nan*zz; uex = nan*zz;
for i=1:numel(zz)
  uex(i) = reference(zz(i));
  u(i) = sum(wj.*integrand(tj,zz(i)));   % G-L quad
end
imagesc(gx,gy,log10(abs(u-uex))); caxis([-14 0]); colorbar;
axis xy equal tight; hold on;
plot(g(Z),'-','color',1*[1 1 1]);  % image of t complex grid lines
plot(g(Z).','-','color',1*[1 1 1]);
plot(schw,'r*','markersize',10);
plot(g(tbern),':','linewidth',1,'color',[.6 0 0]);     % image of Bernstein ellipse
plot(g(-1:1e-2:1),'b-','linewidth',1);
plot(g(tj),'k.','markersize',15);
plot(ztarg,'.','markersize',18,'color',[0 .7 0]);    % target
axis([min(gx) max(gx) min(gy) max(gy)]);
text(-2,1.3,'log$_{10}$ native Gauss--Legendre error near $\Gamma$','interpreter','latex')
set(gcf,'paperposition',[0 0 4 3]);
printpng(2,'../fig/GLerr')

sfigure(3); clf
gx = -1.6:0.02:1.6;
[xx yy] = meshgrid(gx,gx); zz = xx + 1i*yy;
u = panel_conformal_dist(g,zz);
rcschw = exp(panel_conformal_dist(g,schw))      % rho for poly approx on Gamma
c = colormap(jet(256)); colormap((1+c)/2);  % grey it out
contourf(gx,gx,u,0:0.2:2); colorbar('location','southoutside');
axis equal tight; hold on;
plot(g(-1:1e-2:1),'b-','linewidth',2);
plot(g(tj),'k.','markersize',10);
plot(schw,'r*','markersize',10);
text(-1.8,1.8,'(b) conformal equipotential curves for $\Gamma$','interpreter','latex')
set(gcf,'paperposition',[0 0 3 3.5]);
printpng(3,'../fig/confdist')

% convergence plot
ks = 0:0.05:0.8;  % range of curvatures parabola panel: k = curv/2 = 1/(2*(rad of curv))
ps = 4:2:44;      % range of orders of G-L nodes

nt=400; tt = linspace(-1,1,nt); tt=tt(:);  % col of nt test pts unif in param

err = nan(numel(ps),numel(ks)); con = err; rhos = nan(1,numel(ks));

for i=1:numel(ps); p = ps(i);  % -----------
  [tj,~] = lgwt(p,-1,1); tj=tj(:);   % G-L: col of interp pts in param
  %tj = cos(pi*(0:p-1)'/(p-1));       % Cheby
  fj = dens(tj);   % sample vals
  
  for l=1:numel(ks); k = ks(l);   % ........

    g = @(t) t + 1i*k*t.^2;        % panel parameterization gamma
    gp = @(t) 1 + 1i*2*k*t;        % gamma'
    schw = 1i/(4*k);          % Schwarz singularity of unscaled parabola panel
    % estimate convergence rate (C approx th)...
    rhos(l) = exp(panel_conformal_dist(g,schw)); % rate = exp(phi(sing)), Walsh
  
    zsc = 2/(g(1)-g(-1));             % Helsing geom rescale fac
    tauj = -1 + zsc*(g(tj)-g(-1));    % rescaled panel nodes, for ends at +-1
    taut = -1 + zsc*(g(tt)-g(-1));    % rescaled test pts
    tausing = -1 + zsc*(schw-g(-1));  % rescaled Schw singularity
    
    A = ones(p,p); for j=2:p, A(:,j) = A(:,j-1).*tauj; end  % cmplx Vandermonde
    c = A\fj;    % get monomial coeffs
    
    ft = c(1)*ones(nt,1); zt = ones(nt,1);     % eval the poly at test pts tt
    for j=2:p, zt = zt.*taut; ft = ft + c(j)*zt; end
    
    ftt = dens(tt);     % exact test func vals
    ferr = ft-ftt;   % err on test pts
    
    err(i,l) = max(abs(ferr));
    con(i,l) = cond(A);
    fprintf('p=%d\t kappa=%.3g:\t sup err=%.3g\t cond(A)=%.3g\n',p,k,err(i,l),con(i,l))
  end                               % .......
end                           % -----------

sfigure(4); clf
poff = 6;           % empirical offset in order p
errpred = max(1e-15, rhos.^-(ps-poff)');
ll = [1 3 6 9 13]; % 1:4:numel(ks);
col = 'bgrcm';
for l=1:numel(ll),
  h(l) = semilogy(ps,err(:,ll(l)),[col(l) '+-']); hold on;
  if ks(l)>0, plot(ps,errpred(:,ll(l)),[col(l) '.--']); end
end
axis tight; xlabel('$n$','interpreter','latex');
ylabel('sup approximation error of $f$ on $\Gamma$','interpreter','latex');
v = axis; v(3:4)=[1e-15 1]; axis(v);
legend(h,num2cellstr(ks(ll),2,'k = '),'location','northeast');
text(4,1e1,'(c) best poly. approx. err. (curvature $k$)','interpreter','latex');
set(gcf,'paperposition',[0 0 3.3 3.3]);
printpng(4,'../fig/walsh')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function printpng(num, name)
  sfigure(num);
  fname = [name, '.png'];
  print('-dpng', '-r200', fname)
  disp(['Wrote ' fname])
  system(['mogrify -trim ' fname]);
end    
    