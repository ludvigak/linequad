% Setup compiler and optimization flags
cc = 'gcc';
%cc = 'gcc-8';

openmp = true;
verbflag = false;

my_dir = pwd();
c_dir = [pwd '/../c'];
build_dir = [pwd '/../build'];

cflags = ['-std=c99 -fPIC -march=native -fopt-info-vec -I' c_dir];
ldflags = ' -lm ';

% -Ofast makes code 2x faster than -O3, but causes the results of the unstable parts of
% the computation to differ between MEX and MATLAB (i.e. weights will differ, though
% quadrature results will not)
coptimflags = '-Ofast'; % -Wall 
ldoptimflags = '-Ofast';

if openmp
    coptimflags = [coptimflags ' -fopenmp '];
    ldoptimflags = [ldoptimflags ' -fopenmp '];
end
LD = '';

CC = ['CC=''' cc ''''];
CFLAGS = [' CFLAGS=''' cflags ''''];
LDFLAGS = [' LDFLAGS="\$LDFLAGS ' ldflags '" '];
OPTIMFLAGS = [' COPTIMFLAGS=''' coptimflags '''' ' LDOPTIMFLAGS=''' ldoptimflags ''''];
VERBOSE = '';

if verbflag
    VERBOSE = ' -DVERBOSE ';
end

mex_string = ['mex ' CC LD CFLAGS LDFLAGS OPTIMFLAGS VERBOSE ' -outdir bin/'];

% Compile C code
c_files = [' ' c_dir '/linequad.c '];
obj_file = [' ' build_dir '/linequad.o'];
c_call = [cc ' ' cflags ' ' coptimflags ' -c ' c_files '-o' obj_file ' ' ldflags];
disp(c_call)
assert(system(c_call)==0)

% Compile MEX files
eval([mex_string obj_file ' mex/rootfinder_initial_guess_mex.c -output rootfinder_initial_guess_mex'])
eval([mex_string obj_file ' mex/rootfinder_mex.c -output rootfinder_mex'])
eval([mex_string obj_file ' mex/rsqrt_pow_weights_mex.c -output rsqrt_pow_weights_mex'])
eval([mex_string obj_file ' mex/bclag_interp_matrix_mex.c -output bclag_interp_matrix_mex'])

