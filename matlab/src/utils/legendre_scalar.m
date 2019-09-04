function P = legendre_scalar(n, x)
% P = legendre_scalar(n, x)
% => P(l+1) = P_l(x)
% 0 <= l <= n
% 
% Recurrence relation: [http://dlmf.nist.gov/18.9]
% (l+1)P_(l+1)(x)-(2l+1)xP_l(x)+lP_(l-1)(x)=0 	
%
    P = complex(zeros(1, n+1));
    P(1) = 1; % l=0
    P(2) = x; % l=1
    for l=1:n-1
        % Compute l+1
        P(l+1+1) = ( (2*l+1)*x * P(l+1) - l*P(l-1+1) ) / (l+1);
    end
end
