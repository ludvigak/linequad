function [Urealn, Ureal2n, Ucpxn, Ucpx2n, exterior] = dlp_eval(g, den, Z, opt)

% FMM eval in volume
Ufmm = laplace_dlp_fmm_targ(g, den, Z);
% Exterior marker
marker = laplace_dlp_fmm_targ(g, ones(size(den)), Z);
exterior = abs(marker) < 1e-8;

% Precompute mappings
L = legendre.matrix(g.order);
% zhat = L*reshape(g.z, g.order, []); % No pre-rotation
zhat = complex(zeros(g.order, g.numpanels));
for panel_idx=1:g.numpanels
    j = (1:g.order) + g.order*(panel_idx-1);
    za = g.z_edges(panel_idx);
    zb = g.z_edges(panel_idx+1);                                
    zrs = rotate_and_scale(za, zb, g.z(j));
    zhat(:, panel_idx) = L*zrs;
end


[tn,wn] = lgwt2(g.order,-1,1);
[t2n,w2n] = lgwt2(2*g.order,-1,1);
bcw = bclag_interp_weights(tn);
Mbc = bclag_interp_matrix(tn, bcw, t2n);

% Corrections
[Urealn, Ureal2n, Ucpxn, Ucpx2n] = deal(Ufmm);

% Test disable add/subtract through FMM
nofmm = false;
if nofmm
    [Urealn, Ureal2n, Ucpxn, Ucpx2n] = deal(0*Ufmm);
end

parfor i=1:numel(Z)
    if exterior(i)
        continue
    end    
    zt = Z(i);
    % Find closest point
    [~, imin] = min(abs(zt-g.z));
    % Translate to panel idx
    nearest_panel_idx = ceil(imin / g.order);
    % Check vs neighbor panels
    %for panel_idx = mod((nearest_panel_idx+(-1:1))-1, g.numpanels)+1
    % Check vs all panels
    for panel_idx=1:g.numpanels
        % j=indices into panel
        j = (1:g.order) + g.order*(panel_idx-1);
        h = sum(g.ds(j)); % Panel length
        d = min(abs(zt-g.z(j)));                
        did_specquad = false;
        if d < h % Coarse filter
                 % We seem to be close: Find root
            za = g.z_edges(panel_idx);
            zb = g.z_edges(panel_idx+1);                            
            zr = rotate_and_scale(za, zb, zt); % Use pre-rotation
            % The following part (rootfinding) is the bottleneck for this Matlab code
            if g.nearby_schwarz(panel_idx)
                % Matrix:
                coeffs = zhat(:,panel_idx);
                coeffs(1) = coeffs(1)-zr;
                t_all = legendre.roots(coeffs);
                t = t_all(1);
            else
                % Newton:
                t0 = zr;
                [t, relres, iter] = solve_legendre(zhat(:,panel_idx), zr, t0, 20, 1e-13);
            end
            if imag(t)<0 && panel_idx==nearest_panel_idx
                % Deemed outside by nearest panel
                exterior(i) = true;
                continue
            end                                                      
            if bernstein_rad(t) < opt.R
                % We're within critical Bernstein                
                % Direct interaction for subtraction
                if nofmm
                    udirect = 0;
                else
                    udirect = sum(imag( den(j).*g.dz(j)./(g.z(j)-zt) ))/(2*pi);
                end
                % Upsample panel
                z = Mbc*g.z(j);
                d = Mbc*den(j);
                zp = (Mbc*g.zp(j));
                dz = zp.*w2n * g.w(j(1))/wn(1);
                % Upsampled interaction
                uupsamp = imag(sum( d.*dz./(z-zt) ))./(2*pi);
                % Real-var specquad
                wrealn = real_quad_2d(tn, t) * g.w(j(1))/wn(1);
                urealn = imag(sum( den(j).*g.zp(j).*wrealn./(g.z(j)-zt) ))/(2*pi);
                wreal2n = real_quad_2d(t2n, t) * g.w(j(1))/wn(1);
                ureal2n = imag(sum( d.*zp.*wreal2n./(z-zt) ))/(2*pi);
                % Complex-var specquad
                wcpxn = complex_quad_2d(za, zb, zt, g.z(j), -1);
                ucpxn = sum(wcpxn.*den(j))/(2*pi);
                wcpx2n = complex_quad_2d(za, zb, zt, z, -1);
                ucpx2n = sum(wcpx2n.*d)/(2*pi);                
                % Complex 2n specquad only works inside of upsampling boundary              
                if bernstein_rad(t) > sqrt(opt.R)
                    ucpx2n = uupsamp;
                    ucpxn = uupsamp;
                end                
                % Apply
                Urealn(i) = Urealn(i) - udirect + urealn;
                Ureal2n(i) = Ureal2n(i) - udirect + ureal2n;                
                Ucpxn(i) = Ucpxn(i) - udirect + ucpxn;
                Ucpx2n(i) = Ucpx2n(i) - udirect + ucpx2n;
                did_specquad = true;
            end
        end
        if nofmm && did_specquad==false
            udirect = sum(imag( den(j).*g.dz(j)./(g.z(j)-zt) ))/(2*pi);
            Urealn(i) = Urealn(i) + udirect;
            Ureal2n(i) = Ureal2n(i) + udirect;
            Ucpxn(i) = Ucpxn(i) + udirect;
            Ucpx2n(i) = Ucpx2n(i) + udirect;
        end                        
    end
end
