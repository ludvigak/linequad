function [curve, fexact] = test_problem()

% Problem: starfish with distant singularity
curve = diff_curve(starfish(0.3, 5));
zsing = 3+3i;
fexact = @(z) log(abs(z-zsing));
