function show_recurs_in_plane()
%
% Show error plots for all the 3D recursions in the complex plane
%
    
% What to show
zoom = false; % Logarithmic y-scaling
klist = 0:1; % Show recursions for t^k
shifted = false; % Shifted basis (not used)

% How to show it
vpa_reference = true; % Use VPA to compute reference above [-1,1] (only in zoom mode)
if zoom
    % Semilog grid
    x = linspace(-1.5, 1.5, 40+1); % Add one to capture any funny stuff on imag axis
    y = logspace(-10, 0, 40);    
    yscale = 'log';
else
    % Regular grid
    x = linspace(-1.5, 1.5, 70+1);
    y = linspace(-1, 1, 70);
    yscale = 'linear';
end

[X, Y] = meshgrid(x, y);
Z = X + 1i*Y;
Nk = numel(klist);
s = size(Z);
[I1ref,I3ref,I5ref,I1,I3,I5] = deal(zeros([Nk s]));
tic
parfor i=1:numel(Z)
    z = Z(i);
    kmax = max(klist)+1;
    if shifted
        [p1,p3,p5] = rsqrt_pow_integrals(z, kmax);           
    else
        [p1,p3,p5] = rsqrt_pow_integrals_noshift(z, kmax);           
    end
    for j=1:Nk
        k = klist(j);
        I1(j,i) = p1(k+1);
        I3(j,i) = p3(k+1);
        I5(j,i) = p5(k+1);
        if vpa_reference && abs(real(z)) <= 1 && zoom
            % VPA ref is better above [-1,1]
            [p1,p3,p5] = rsqrt_pow_integrals_noshift(vpa(z), kmax);
            I1ref(j,i) = p1(k+1);
            I3ref(j,i) = p3(k+1);
            I5ref(j,i) = p5(k+1);
        else
            % Try to numerically compute references
            ref = cell(1,5);
            warning off
            for p=[1 3 5]
                if abs(real(z)) <= 1
                    % real(z) is inside interval, chop up into two integrals
                    if shifted
                        integrand = @(t) t.^k./abs(t-imag(z)*1i).^p;
                        ref{p} = integral(integrand, -1-real(z), 0, 'abstol',0,'reltol',0) ...
                                 + integral(integrand, 0, 1-real(z), 'abstol',0,'reltol',0);                          
                    else
                        integrand = @(t) t.^k./abs(t-z).^p;
                        ref{p} = integral(integrand, -1, real(z), 'abstol',0,'reltol',0) ...
                                 + integral(integrand, real(z), 1, 'abstol',0,'reltol',0);             
                    end
                else
                    if shifted
                        shift = real(z);
                    else
                        shift = 0;
                    end                
                    integrand = @(t) (t-shift).^k./abs(t-z).^p;                
                    ref{p} = integral(integrand, -1, 1, 'abstol',0,'reltol',0);   
                end
            end
            warning on
            I1ref(j,i) = ref{1};
            I3ref(j,i) = ref{3};
            I5ref(j,i) = ref{5};            
        end                
    end
end
toc

function errplot(I, Iref, k, p, Inorm)
    I = reshape(I, size(Z));
    Iref = reshape(Iref, size(Z));
    Inorm = reshape(Inorm, size(Z));
    % Plot relative error in computed value
    % Probably need a better normalization here
    E = abs(I-Iref) ./ abs(Inorm);
    logE = log10(E );% + min(E(E>0)));
    pcolor(real(Z), imag(Z), logE);
    xlabel('Re t')
    ylabel('Im t')
    %axis image
    set(gca,'yscale',yscale)
    %shading interp
    shading flat
    h = colorbar();
    xlabel(h, 'log10(err)')
    hold on
    plot([-1 1],[0 0],'.-k', 'LineWidth',2,'MarkerSize',10)
    if shifted
        title(sprintf('$\\int (t-a)^{%d}/|t-z|^%d$', k, p), 'interpreter', 'latex')
    else
        title(sprintf('$\\int t^{%d}/|t-z|^%d$', k, p), 'interpreter', 'latex')    
    end
    %caxis([-16 -14])
end

clf
n = 1;
for i=1:Nk    
    k = klist(i);
    subplot(Nk,3,n); errplot(I1(i,:), I1ref(i,:), k, 1, I1ref(1,:)); n = n+1;
    subplot(Nk,3,n); errplot(I3(i,:), I3ref(i,:), k, 3, I3ref(1,:)); n = n+1;
    subplot(Nk,3,n); errplot(I5(i,:), I5ref(i,:), k, 5, I5ref(1,:)); n = n+1;
end

end