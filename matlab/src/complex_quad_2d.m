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
    
    if sgn*imag(ztgtr) > 0 
        % Target on wrong side of real axis
        % Check if target in box containing curve (Helsing2008)
        minx = min(-1, min(real(zsctr)));
        maxx = max( 1, max(real(zsctr)));
        miny = min(0, min(imag(zsctr)));
        maxy = max(0, max(imag(zsctr)));       
        x = real(ztgtr);
        y = imag(ztgtr);
        inbox = minx <= x && x <= maxx && ...
                miny <= y && y <= maxy;
        % Add optional condition |Re z| <= 1
        if inbox || abs(real(ztgtr)) <= 1
            % Point is between panel and real axis
            p(1) = p(1) - 2i*pi*sgn;
        end
    end
    

    for k=1:N-1
        p(k+1)=ztgtr*p(k)+c(k);
    end

    warning('off', 'MATLAB:nearlySingularMatrix')
    x = A.'\p;
    warning('on', 'MATLAB:nearlySingularMatrix')                    
    w = imag(x);
end
