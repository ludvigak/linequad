% test err for poly interp on complex curved fiber in 2D,
% to disentangle that error from the close-eval error.
% Notation: t = param, tau = gamma(t) = complex location.
% Barnett 3/17/19. Complex approx theory pred, 3/21/19.

clear; expt = 'b';   % 's' single p and curvature k; 'p' err vs p (fixed k);
                     % 'k' err vs k (fixed p); 'b' err vs p and k
% choose fixed p, k here:
ks = .7;  % one curvature of parabola panel: k = curv/2 = 1/(2*(rad of curv))
ps = 16;   % one order of G-L nodes
if expt=='p'
  ps = 4:2:44;      % range of orders of G-L nodes  
elseif expt=='k'
  ks = 0:0.05:0.8;  % range of curvatures
elseif expt=='b'
  ks = 0:0.05:0.8;  % range of curvatures
  ps = 4:2:44;      % range of orders of G-L nodes
end

freq = 3; f = @(t) sin(1+freq*t);  % smooth func of param to interp, eg density
%f = @(t) 0.4*exp(t);

nt=400; tt = linspace(-1,1,nt); tt=tt(:);  % col of nt test pts unif in param

err = nan(numel(ps),numel(ks)); con = err; rhos = nan(1,numel(ks));

for i=1:numel(ps); p = ps(i);  % -----------
  [tj,~] = lgwt(p,-1,1); tj=tj(:);   % G-L: col of interp pts in param
  %tj = cos(pi*(0:p-1)'/(p-1));       % Cheby
  fj = f(tj);   % sample vals

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
    
    if expt=='s', figure; subplot(2,1,1);
      plot(real(tauj),imag(tauj),'+'); hold on; plot(real(taut),imag(taut),'.');
      plot([-1 1],[0 0],'ko'); plot(real(tausing),imag(tausing),'*');
      axis equal tight; title(sprintf('rescaled panel: p=%d, k=%.3g',p,k));
    end
    
    A = ones(p,p); for j=2:p, A(:,j) = A(:,j-1).*tauj; end  % cmplx Vandermonde
    c = A\fj;    % get monomial coeffs
  
    ft = c(1)*ones(nt,1); zt = ones(nt,1);     % eval the poly at test pts tt
    for j=2:p, zt = zt.*taut; ft = ft + c(j)*zt; end
    
    ftt = f(tt);     % exact test func vals
    ferr = ft-ftt;   % err on test pts
    if expt=='s', subplot(2,1,2);
      plot(real(tauj),0*tauj,'+'); hold on; plot(real(taut),real(ferr),'.-');
      axis tight; xlabel('rescaled Re \tau'); ylabel('Re interp err');
    end
  
    err(i,l) = max(abs(ferr));
    con(i,l) = cond(A);
    if expt~='b', fprintf('p=%d\t kappa=%.3g:\t sup err=%.3g\t cond(A)=%.3g\n',p,k,err(i,l),con(i,l)), end
  end                              % .......
end                           % -----------

if expt=='p', figure;  % convergence plot
  fprintf('conv rate rho = %.3g\n',rhos(1))
  semilogy(ps,err,'+-'); hold on; plot(ps,rhos(1).^-ps,'r--');  % uses rate rho
  plot(ps,con,'o-');
  xlabel('p'); legend('sup interp err','Schwarz-Walsh','cond(A)'); axis tight;
elseif expt=='k', figure;  % convergence plot
  semilogy(ks,err,'+-'); hold on; plot(ks,con,'o-');
  xlabel('k (curvature)'); legend('sup interp err','cond(A)'); axis tight;
elseif expt=='b', figure;   % pair of images
  subplot(1,2,1); imagesc(ps,ks,log10(err')); colorbar; axis xy;
  xlabel('p'); ylabel('k (curvature)'); title('log_{10}(sup interp err)');
  subplot(1,2,2); imagesc(ps,ks,log10(con')); colorbar; axis xy;
  xlabel('p'); ylabel('k (curvature)'); title('log_{10} cond(A)');
  print -dpng err_cmplx_interp_vs_p_and_k.png
  sc = 1; errpred = max(1e-15, sc*rhos.^-ps');
  figure; ll = 3:4:numel(ks); col = 'bgrcm';
  for l=1:numel(ll),
    h(l) = semilogy(ps,err(:,ll(l)),[col(l) '+-']); hold on; plot(ps,errpred(:,ll(l)),[col(l) '.--']);
  end
  axis tight; xlabel('p'); ylabel('err');
  legend(h,num2cellstr(ks(ll),[],'k = '),'location','southwest');
  
  if 0, figure; subplot(1,2,1);   % not so useful
    title('compare sup interp err to Schwarz-Walsh pred');
    mesh(ps,ks,log10(err')); colorbar; axis vis3d tight
    hold on; mesh(ps,ks,log10(errpred'));
    v = axis; v(5:6)=[-15 0]; axis(v);
    subplot(1,2,2);
    imagesc(ps,ks,log10(err' ./ errpred')); caxis([-2 2]); colorbar; axis xy
  end
end


