function show_3d_mapping()
% Generate plot show circles corresponding to complex roots

addpath(genpath('src'))

t = 0;
% Set up curve
syms t
xsym = cos(t-1);
ysym = sin(t-1).*(1+0.5*sin(t*pi/8));
zsym = 0.5*t;
clear t

x = matlabFunction(xsym);
y = matlabFunction(ysym);
z = matlabFunction(zsym);

% Discretize bdry and form cheb representation
n = 16;
Pg = @(t) [x(t), y(t), z(t)];

function plot_circle(t0, style)
% Root value
    gr = real(Pg(t0));
    gi = imag(Pg(t0));
    % Find vectors spanning plane
    v = cross(gi, rand(1,3));
    v = v / norm(v);
    w = cross(v, gi);
    w = w / norm(w);
    R = norm(gi);
    th = linspace(0, 2*pi);
    c = gr + R*(cos(th')*v + sin(th')*w);
    cx = c(:, 1);
    cy = c(:, 2);
    cz = c(:, 3);
    plot3(cx, cy, cz, style)    
end

% Plot with GL points
tn = lgwt2(n, -1, 1);
xn = x(tn);
yn = y(tn);
zn = z(tn);


sfigure(1);
clf
plot(tn, tn*0, '.-r', 'markersize',10)
hold on
axis equal
%axis off
grid
xlabel('Re t')
ylabel('Im t')
xlim([-1.2, 1.2])
%plot(t0, 'bo')
%plot(t0', 'bo')

sfigure(2);
clf
plot3(xn, yn, zn, '.r', 'markersize',10)
hold on
tfine = linspace(-1,1); % smooth curve
plot3(x(tfine), y(tfine), z(tfine), 'r-')
axis equal
axis tight
grid
xlabel('x')
ylabel('y')
zlabel('z')
%axis off

% Circle along curve
t0list = linspace(-1,1,11) + 0.2i;
for i=1:numel(t0list)
    sfigure(1);
    plot(t0list(i), '+b')
    plot(t0list(i)', '+b')
    sfigure(2);
    plot_circle(t0list(i), '-b');
end
for b=linspace(0.1, 1.0, 5)
    %plot_circle(0.2 + b*1i, '-k');
end

sfigure(3);
clf
plot3(xn, yn, zn, '.r', 'markersize',10)
hold on
tfine = linspace(-1,1); % smooth curve
plot3(x(tfine), y(tfine), z(tfine), 'r-')
axis equal
axis tight
grid
xlabel('x')
ylabel('y')
zlabel('z')

% Concentric circles
t0list = 0.5 + 1i*linspace(0.05, 0.8, 8);
for i=1:numel(t0list)
    sfigure(1);
    plot(t0list(i), 'xk')
    plot(t0list(i)', 'xk')
    sfigure(3);
    plot_circle(t0list(i), '-k');
end

% Output
function makepng(fignum, filename)
    sfigure(fignum);
    publication_fig()
    box on
    print('-dpng','-r600',filename);
    system(['mogrify -trim ', filename])    ;
    disp(['Wrote ', filename]);
end

makepng(1, '../fig/poles_3d.png')
makepng(2, '../fig/poles_3d_mapped1.png')
makepng(3, '../fig/poles_3d_mapped2.png')

end