function x = pvand(alpha, b)
% x = pvand(alpha, b)
%
% Solves system A*x = b
% A is Vandermonde matrix, with nonstandard definition
% A(i,j) = alpha(j)^i
%
% Algorithm by Bjorck & Pereyra
% Mathematics of Computation, Vol. 24, No. 112 (1970), pp. 893-903
% https://doi.org/10.2307/2004623    
%
    n = numel(alpha);
    x = b;
    for k=1:n
        for j=n:-1:k+1
            x(j) = x(j)-alpha(k)*x(j-1);
        end
    end
    for k=n-1:-1:1
        for j=k+1:n
            x(j) = x(j)/(alpha(j)-alpha(j-k));
        end
        for j=k:n-1
            x(j) = x(j)-x(j+1);
        end
    end
end
