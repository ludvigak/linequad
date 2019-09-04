function pk = int_recur_rsqrt(z, p)
    a = real(z);
    b = imag(z);
    q = nan(p,1);   % r-kernel: q_k array for this targ
    fq0 = @(s) (s*sqrt(s^2+b^2) + b^2*log(s+sqrt(s^2+b^2)))/2; % antideriv
    q(1) = fq0(1-a)-fq0(-1-a);  % q_0  (note 1-offset index from k)
    fq1 = @(s) (s^2+b^2)^(3/2)/3;
    q(2) = fq1(1-a)-fq1(-1-a);  % q_1
    for k=0:p-3
        q(k+3) = (-(k+1)*b^2*q(k+1) + (1-a)^(k+1)*((1-a)^2+b^2)^(3/2) - (-1-a)^(k+1)*((-1-a)^2+b^2)^(3/2))/(4+k);
    end
    %q  % note q_k grows for a near ends
    % now 1/r kernel: p_k array for this targ
    pk = nan(p,1);
    fp0 = @(s) log(s+sqrt(s^2+b^2));
    pk(1) = fp0(1-a)-fp0(-1-a);
    pk(2) = sqrt((1-a)^2+b^2) - sqrt((-1-a)^2+b^2);
    for k=0:p-3, pk(k+3) = q(k+1) - b^2*pk(k+1); end   % p recurrence
end    
