% Generate plots illustrating relation preimages / target points

clear

save_paper_plots = true;

addpath(genpath('../../matlab-lib'))
addpath('../../chebfun')
addpath(genpath('src'))

% Notes:
%
% If we want to use Legendre points instead:
% leg2cheb(chebfun.idlt(y(legpts(n)))) == chebcoeffs(Py) 


% Set up curve
syms t
%xsym = cos(0.8*(t-1.5));
%ysym = sin(0.8*(t-1.5)).*(1+0.5*sin(t*pi/8));
xsym = t-0.2
ysym = (t-0.2)^2
clear t

x = matlabFunction(xsym);
y = matlabFunction(ysym);
xderiv = matlabFunction(diff(xsym));
yderiv = matlabFunction(diff(ysym));

z = @(t) x(t) + 1i*y(t);
W = @(t) sqrt(xderiv(t).^2 + yderiv(t).^2);

drawgrid = true;

% Point on single valued region
t01 = 0.3 - 0.2i;
z01 = z(t01);
x01 = real(z01);
y01 = imag(z01);

% Pick a target point in concave region (mapping not single valued)
t0 = -0.6 + 0.34i;
z0 = z(t0);
x0 = real(z0);
y0 = imag(z0);

% Discretize bdry and form cheb representation
n = 16;
tn = chebpts(n);
xn = x(tn);
yn = y(tn);
Px = chebfun(xn);
Py = chebfun(yn);
Pdx = Px-x0;
Pdy = Py-y0;
Qfun = @(t) Pdx(t) + 1i*Pdy(t);
% Form denominator polynomial and get roots
Qn = Pdx + 1i*Pdy;
rall = roots(Qn, 'all');
% Get their Bernstein radii, which poles to include should depend o
rho = abs(rall + sign(real(rall)).*sqrt(rall.^2-1));
r = rall(rho<3); % lets pick 3
% Check accuracy of poles
Qnr=Qn(r);
Qfunr=Qfun(r);
% Find singularity of parametrization
zp = diff(Px+1i*Py);
tsing = roots(zp, 'all');
tsing = tsing(1);
Rb = bernstein_rad(tsing)*0.95
th = linspace(0, 2*pi);
ellipse = (Rb*exp(1i*th) + 1./(Rb*exp(1i*th)))/2;

h = 2/10;
%t1 = -1:h:1;
%t2 = linspace(0, imag(tsing), 5);
t1 = linspace(-(Rb+1/Rb)/2, (Rb+1/Rb)/2, 10);
%t2 = linspace(-(Rb-1/Rb)/2, (Rb-1/Rb)/2, 9);
t2 = linspace(-(Rb-1/Rb)/2, (Rb-1/Rb)/2, 5);
[T1,T2] = ndgrid(t1, t2);
T = T1 + 1i*T2;

% When plotting, use legpts
tn = legpts(n);
xn = x(tn);
yn = y(tn);


% setup figs
if save_paper_plots
    sfigure(1);clf;publication_fig;
    sfigure(2);clf;publication_fig;
end

sfigure(1);
clf
C = 0.8; % greyness
if drawgrid 
    plot(T.','Color', C*[1 1 1])
    hold on
    plot(T, 'Color', C*[1 1 1])
end
plot(tn, tn*0, '.-r')
hold on
plot(tsing, 'pk','markerfacecolor','k')
plot(ellipse, '-k','color', 0.5*[0 1 0])
axis equal

sfigure(2);
clf
Z = x(T) + 1i*y(T);

% Draw smooth mappings of grid
t1f = linspace(t1(1), t1(end));
t2f = linspace(t2(1), t2(end));
if drawgrid
    for i=1:numel(t1)
        plot(complex(z(t1(i)+1i*t2f)),'-', 'Color', C*[1 1 1]);
        hold on
    end
    for j=1:numel(t2)
        plot(complex(z(t1f+1i*t2(j))), 'Color', C*[1 1 1]);
        hold on
    end
end
plot(xn, yn, '.r')
hold on
plot(complex(z(linspace(-1,1))), '-r')

axis equal

% Plot roots and their path from real axis

r = [r; t01];
for i=1:numel(r)
    tr = real(r(i));
    ti = imag(r(i));
    sfigure(2);
    plot(z(tr+1i*ti), 'bo')
    plot(z(tr + 1i*linspace(0, ti)), '-b')
    sfigure(1);
    plot((tr+1i*ti), 'bo')
    plot((tr + 1i*linspace(0, ti)), '-b')
    
end

sfigure(2);
axis tight
axis off
%plot(complex(z0), '.k')
%plot(complex(z01), '.k')
plot(z(ellipse), '-k', 'color', 0.5*[0 1 0])
text(real(z(1))-0.3,imag(z(1))+1,'$\gamma(t) \in \mathbf C$','interpreter','latex')
text(real(z0), imag(z0)+0.2, '$\zeta = \gamma(z_1) = \gamma(z_2)$','interpreter','latex')

zlist = r(imag(r)>0);
zconj = zlist';
plot(z(tsing), 'pk','markerfacecolor','k')
text(real(z(tsing))+0.1, imag(z(tsing))+0.05, '$\gamma(t_{\ast})$', 'interpret', 'latex')


text(real(z01)+0.1, imag(z01), '$\xi=\gamma(w_1)$','interpreter','latex')


sfigure(1);
axis tight
axis off
text(-1,-0.1, '-1')
text(1,-0.1, '1')
text(real(zlist(1))+0.05, imag(zlist(1)), '$ z_1$','interpreter','latex')
text(real(zlist(2))+0.05, imag(zlist(2)), '$ z_2$','interpreter','latex')

text(real(t01)+0.05, imag(t01), '$w_1$','interpreter','latex')
text(-1,1,'$t \in \mathbf{C}$', 'interpreter', 'latex')

text(real(tsing)+0.05, imag(tsing)+0.01, '$t_{\ast}$', 'interpret', 'latex')

if ~save_paper_plots
    return
end

sfigure(1)
publication_fig()
filename = '../fig/all_poles.png';
print('-dpng','-r600',filename)
system(['mogrify -trim ' filename]);
disp(['Wrote ', filename])

% matlab2tikz('../fig/all_poles.tikz', ...
%             'width','\fwidth', ...
%             'extraAxisOptions','enlargelimits=false')

sfigure(2)
publication_fig()
filename = '../fig/all_poles_mapped.png';
print('-dpng','-r600',filename)
system(['mogrify -trim ' filename]);
disp(['Wrote ', filename])

% matlab2tikz('../fig/all_poles_mapped.tikz', ...
%             'width','\fwidth', ...
%             'extraAxisOptions','enlargelimits=false')




% Should get tikz output working for this...
%write_tikz(1, '../fig/all_poles')
%write_tikz(2, '../fig/all_poles_mapped')