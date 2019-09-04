function w = complex_quad_2d(a,b,ztg,zsc,sgn)
% Helsing & Ojala quadrature
    sgn = sgn/abs(sgn);
    N = numel(zsc);
    cc=(b-a)/2;
    ztgtr=(ztg-(b+a)/2)/cc;
    zsctr=(zsc-(b+a)/2)/cc;
    A=complex(ones(N));
    for k=2:N
        A(:,k)=zsctr.*A(:,k-1);
    end
    
    p=complex(zeros(N,1));
    c=(1-(-1).^(1:N))./(1:N);
    p(1)=log((b-ztg)/(a-ztg));
    if sgn*imag(ztgtr) > 0 && abs(real(ztgtr)) < 1
        p(1) = p(1) - 2i*pi*sgn;
    end

    for k=1:N-1
        p(k+1)=ztgtr*p(k)+c(k);
    end

    warning('off', 'MATLAB:nearlySingularMatrix')
    x = A.'\p;
    warning('on', 'MATLAB:nearlySingularMatrix')                    
    w = imag(x);
end
