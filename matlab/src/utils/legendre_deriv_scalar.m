function [P, D] = legendre_deriv_scalar(n, x)
    P = complex(zeros(1, n+1));
    D = complex(zeros(1, n+1));    
    P(1) = 1; % l=0
    D(1) = 0;
    P(2) = x; % l=1
    D(2) = 1;
    for l=1:n-1
        % Compute l+1
        P(l+1+1) = ( (2*l+1)*x * P(l+1) - l*P(l-1+1) ) / (l+1);
        D(l+1+1) = ( (2*l+1)*(P(l+1) + x * D(l+1)) - l*D(l-1+1) ) / (l+1);    
    end
end
