function compare_2d_complex_real(varargin)
% 
% To generate figures in paper, run
%
% for k=[0.1, 0.25, 0.4]; compare_2d_complex_real(k); end
%
    
    save_paper_plots = false;
    hires = false; % Slow, but good for paper output
    
    % How to find roots
    rootfinding = 'matrix'; % matrix/newton
    
    shape = 'parabola'
    switch shape
        case 'parabola'
          % Parametrization of panel: parabola
          % k=Curving of panel, higher is harder to compute       
          % try range 0.2--0.4
          k = 0.3; 
          if ~isempty(varargin); k=varargin{1}; end
          g = @(t) t + 1i*k*t.^2 -1i*k;
          gp = @(t) 1 + 1i*2*k*t;    
          gs = 1i/(2*k); % Schwarz singularity
          assert(abs(gp(gs))<1e-14)
      case 'circle'
        % Parametrization of panel: circle segment
        % Try range k=0.1--0.6
        k = 0.5;    
        if ~isempty(varargin); k=varargin{1}; end
        gs = @(t) exp(1i*(pi/2 + k*t));
        g = @(t) -1 + (gs(t) - gs(-1)) * 2 / (gs(1) - gs(-1));
        gp = @(t) 1i*k*gs(t) * 2 / (gs(1) - gs(-1));
    end    
    
    % Density to use
    density_xy = @(t) real(g(t)).*imag(g(t)); % sigma(t) = x(t)*y(t)
    density_t = @(t) t; % sigma(t) = t
    density_1 = @(t) t.^0; % sigma(t) = 1, (trivial case for complex interpolatory)
    
    density = density_xy; % try xy/t/1
        
    % verynear = use log spacing in imag(t)
    verynear = false;
    
    % % Disable subplotting (complete hack)
    % for i=1:9; figure(i); clf(); hold on; end
    % function subplot(i,j,k); figure(k); hold on; publication_fig(); end
    %% END PARAMS
    
    % Laplace double layer kernel
    kernel_cpx = @(x, y, ny) ( ny ./ (x-y) );    
    kernel = @(x, y, ny) real( kernel_cpx(x, y, ny)  );
    integrand = @(t,z) kernel(z, g(t), 1i*gp(t)./abs(gp(t))) .* density(t) .* ...
        abs(gp(t));    
    integrand_cpx = @(t,z) kernel_cpx(z, g(t), 1i*gp(t)) .* density(t);
    
    % Gauss-Legendre quadrature nodes and weights
    [t16,w16] = lgwt2(16,-1,1);
    [t32,w32] = lgwt2(32,-1,1);
    M = bclag_interp_matrix(t16, bclag_interp_weights(t16), t32); % Upsampling matrix t16 -> t32

    % Integrand after interpolating from 16 to 32 points
    integrand_interp32 = @(z) kernel(z, M*g(t16), 1i*M*gp(t16)) .* (M*density(t16));
    integrand_cpx_interp32 = @(z) kernel_cpx(z, M*g(t16), 1i*M*gp(t16)) .* (M*density(t16));
    
    
    % Reference computation and grid setup
    reference = @(z) integral(@(t) integrand(t, z), -1, 1, 'abstol', eps, 'reltol', eps);    
    if verynear
        if strcmp( func2str(density), func2str(density_1) )
            % for density_1
            reference = @(z) imag(log((g(1)-z)/(g(-1)-z)));
            t1 = linspace(-1.5, 1.5, 50);
            t1 = linspace(-1e-10, 1e-10, 50);
            t2 = -logspace(-10, 0, 50);        
        else
            % for density_xy/t (more expensive)            
            t1 = linspace(-1.5, 1.5, 30);
            t2 = -logspace(-6, 0, 30);        
        end    
    else
        if hires
            Ngrid = 400;
        else
            Ngrid = 40;
        end
        t1 = linspace(-2.9, 2.9, Ngrid);
        t2 = linspace(-1.8, imag(gs), Ngrid);
    end
    
    % Create Legendre expansion of panel (in rescaled coordinates)
    za = g(-1);
    zb = g(1);
    L = legendre.matrix(numel(t16));
    coeffs = L*rotate_and_scale(za, zb, g(t16));
    
    % Create a grid around curve segment by complexifying parameter
    [T1, T2] = ndgrid(t1, t2);
    T = T1 + 1i*T2;
    Z = g(T);
    
    % Compute values on grid
    [I, Q16, Q32, Qcpx16, Qcpx32, Qcpxmod16, Qcpxmod32] = deal(zeros(size(Z)));
    tic
    parfor i=1:numel(Z)
        z = Z(i);
        % === Reference
        I(i) = reference(z);
        % === Direct
        Q16(i) = sum(integrand(t16, z) .* w16);
        Q32(i) = sum(integrand_interp32(z) .* w32);
        % === Complex interpolatory with real basis
        % * Find root
        zr = rotate_and_scale(za, zb, z);
        switch rootfinding
          case 'newton'
            t0 = zr;
            t = legendre.newton(coeffs, zr, t0, 20, 1e-13);
          case 'matrix'            
            coeffs_z = coeffs;
            coeffs_z(1) = coeffs_z(1)-zr;
            t_all = legendre.roots(coeffs_z);
            t = t_all(1);
          otherwise
            error
        end
        % t should now be close to T(i)
        % * Compute weights
        wcpx16full = real_quad_2d(t16, t);        
        Qcpxmod16(i) = real(sum(integrand_cpx(t16, z) .* wcpx16full));
        % ... upsampled to 32 points        
        wcpx32full = real_quad_2d(t32, t);        
        Qcpxmod32(i) = real(sum(integrand_cpx_interp32(z) .* wcpx32full));     
        % === Complex interpolatory
        sgn = -imag(T(i)); % Cheating so that we can look at both sides of panel
        wcpx16 = complex_quad_2d(za, zb, z, g(t16), sgn);
        Qcpx16(i) = sum(wcpx16 .* density(t16));
        % ... upsampled to 32 points        
        wcpx32 = complex_quad_2d(za, zb, z, g(M*t16), sgn);
        Qcpx32(i) = sum(wcpx32 .* (M*density(t16)) );                
    end
    toc
    % Compute rel. errors
    Inorm = norm(I(:), inf);
    E16 = abs(Q16-I) / Inorm + eps;
    E32 = abs(Q32-I) / Inorm + eps;    
    Ecpx16 = abs(Qcpx16-I) / Inorm + eps;
    Ecpxmod16 = abs(Qcpxmod16-I) / Inorm + eps;    
    Ecpx32 = abs(Qcpx32-I) / Inorm + eps;
    Ecpxmod32 = abs(Qcpxmod32-I) / Inorm + eps;
    
    % Plot
    if verynear
        clf()
        subplot(3,3,1)
        plot(Z, '-k')
        hold on
        plot(Z.', '-k')    
        plot(complex(g(t16)), '.-r')    
        axis equal

        subplot(3,3,2)
        pcolor(T1, T2, log10(E16));    
        title('Direct quadrature (16 pts)')

        subplot(3,3,3)
        pcolor(T1, T2, log10(E32));    
        title('Direct quadrature (32 pts)')        
        
        subplot(3,3,4)
        pcolor(T1, T2, log10(Ecpx16));    
        title('Complex basis (16 pts)')    

        subplot(3,3,7)
        pcolor(T1, T2, log10(Ecpx32));    
        title('Complex basis (interp to 32)')        
                
        subplot(3,3,6)
        pcolor(T1, T2, log10(Ecpxmod16));    
        title('Real basis (16 pts)')    

        subplot(3,3,9)
        pcolor(T1, T2, log10(Ecpxmod32));
        title('Real basis (interp to 32)')        
        
        for i=2:9
            subplot(3,3,i);
            colorbar
            shading interp
            caxis([-15 -5])
            hold on
            plot(complex(t16), '.-r')    
            set(gca,'yscale','log')
        end    
    else
        
        for i=1:6
            sfigure(i);
            clf
            publication_fig
        end
        
        sfigure(1);
        pcolor(real(Z), imag(Z), log10(E16));    
        title(['Direct, 16 pts, k=' num2str(k)])
        hold on
        plot(complex(g(t16)), '.-r')            

        sfigure(2);
        pcolor(real(Z), imag(Z), log10(E32));    
        title(['Direct, 32 pts, k=' num2str(k)])    
        hold on
        plot(complex(g(t32)), '.-r')                    
        
        sfigure(3);
        pcolor(real(Z), imag(Z), log10(Ecpx16));    
        title(['Complex, 16 pts, k=' num2str(k)])    
        hold on
        plot(complex(g(t16)), '.-r')            

        sfigure(4);
        pcolor(real(Z), imag(Z), log10(Ecpx32));    
        title(['Complex, 32 pts, k=' num2str(k)])        
        hold on
        plot(complex(g(t32)), '.-r')                    
        
        sfigure(5);
        pcolor(real(Z), imag(Z), log10(Ecpxmod16));    
        title(['Real, 16 pts, k=' num2str(k)])    
        hold on
        plot(complex(g(t16)), '.-r')            

        sfigure(6);
        pcolor(real(Z), imag(Z), log10(Ecpxmod32));    
        title(['Real, 32 pts, k=' num2str(k)])        
        hold on
        plot(complex(g(t32)), '.-r')                    
    end    
    
    % Figure out how to frame
    % dom = Z(E16 > 1e-15);
    % xmin = min(real(dom));
    % xmax = max(real(dom));
    % ymin = min(imag(dom));
    % ymax = max(imag(dom));   
    % ax = 1.1*[xmin, xmax, ymin, ymax];
    
    % Fixed frame
    %ax = [-2.0, 2.0, -1.6, 1.0];
    ax = [-2.4, 2.4, -1.9, 1.4];
    
    for i=1:6
        sfigure(i)
        plot(complex(g(gs)), 'pr', 'MarkerFaceColor','r')
        axis equal
        axis(ax)
        %axis image
        shading interp
        caxis([-16 -6])
        axis off        
        publication_fig        
    end
    
    function printpng(num, name)
        sfigure(num);
        fname = [name, '.png'];
        print('-dpng', '-r200', fname)
        disp(['Wrote ' fname])
        system(['mogrify -trim ' fname]);
    end    
    
    if ~save_paper_plots
        return
    end
    printpng(1, sprintf('../fig/compare_2d_direct_k%g_n16',k))
    printpng(2, sprintf('../fig/compare_2d_direct_k%g_n32',k))
    printpng(3, sprintf('../fig/compare_2d_complex_k%g_n16',k))
    printpng(4, sprintf('../fig/compare_2d_complex_k%g_n32',k))
    printpng(5, sprintf('../fig/compare_2d_real_k%g_n16',k))
    printpng(6, sprintf('../fig/compare_2d_real_k%g_n32',k))
    
    % Save colorbar as separate fig
    sfigure(99);
    publication_fig
    clf
    caxis([-16 -6])
    axis off 
    %axis image    
    cb = colorbar('FontSize',15);
    pos = get(cb,'position');
    pos(1) = 1/2;
    xlabel(cb, 'log_{10} E_{rel}')               
    set(cb,'position',pos)
    printpng(99, '../fig/colorbar_16_6')
end
