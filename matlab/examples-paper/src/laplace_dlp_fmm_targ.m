function Ufmm = laplace_dlp_fmm_targ(grid, density, z, varargin)

% 1 => tolerance =.5d-3  =>  3 digits
% 2 => tolerance =.5d-6  =>  6 digits
% 3 => tolerance =.5d-9  =>  9 digits
% 4 => tolerance =.5d-12 => 12 digits
% 5 => tolerance =.5d-15 => 15 digits

if isempty(varargin)
    opt = struct();
else
    opt = varargin{1};
end

if isfield(opt, 'TOL')
    TOL = opt.TOL;
    if TOL > .5d-3
        iprec = 1;
    elseif TOL > .5d-6
        iprec = 2;
    elseif TOL > .5d-9
        iprec = 3;
    elseif TOL > .5d-12
        iprec = 4;
    else 
        iprec = 5;
    end
else
    iprec=5;
end

nsource = grid.N;
source = [real(grid.z(:)) imag(grid.z(:))].';
dipstr = density(:) .* grid.zp(:) .* grid.w(:);
ntarget = numel(z);
target = [real(z(:)) imag(z(:))].';
ifpot = 0;
ifgrad = 0;
ifhess = 0;
ifpottarg = 1;
ifgradtarg = 0;
ifhesstarg = 0;
U = zfmm2dpart(iprec,nsource,source,dipstr,ifpot,ifgrad,ifhess,ntarget,target, ...
               ifpottarg,ifgradtarg,ifhesstarg);
Ufmm = -imag( reshape(U.pottarg, size(z)) )/(2*pi);
