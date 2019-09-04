function [troot, converged] = rootfinder(xhat, yhat, zhat, x0, y0, z0, tinit)
% Find roots using Newton and Muller
    
    VERBOSE = false;
    
    t = tinit;
    n = numel(xhat);
    tol = 1e-14;
    maxiter_newton = 20;
    maxiter_muller = 20;    
    % === Newton
    % Setup history variables (needed in Muller)
    Fp = 0; tp = 0; Fpp = 0; tpp = 0;
    converged = false;
    for iter=1:maxiter_newton      
        % Compute
        [P, D] = legendre_deriv_scalar(n-1, t);
        F = (P*xhat-x0)^2 + (P*yhat-y0)^2 + (P*zhat-z0)^2;
        Fprime = 2*(P*xhat-x0)*(D*xhat) + ...
                 2*(P*yhat-y0)*(D*yhat) + ...
                 2*(P*zhat-z0)*(D*zhat);
        dt = -F/Fprime;                        
        % Update history
        tpp = tp;
        Fpp = Fp;        
        Fp = F;
        tp = t;       
        % Update root
        t = t+dt;
        absres = abs(dt);
        if absres < tol
            converged = true;
            break
        end
    end
    if converged
        if VERBOSE
            fprintf('Newton converged in %d iterations.\n', iter);
        end
        troot = t;
        return;
    end    
    % === Mulleri
    if VERBOSE
        fprintf('Newton did not converge after %d iterations (abs(dt)=%g), switching to Muller\n',...
                iter, absres);
    end
    converged = false;
    for iter=1:maxiter_muller
        % Compute
        P = legendre_scalar(n-1, t);
        F = (P*xhat-x0)^2 + (P*yhat-y0)^2 + (P*zhat-z0)^2;            
        % Mullers method
        q = (t-tp)/(tp - tpp);
        A = q*F - q*(q+1)*Fp + q^2*Fpp;
        B = (2*q+1)*F - (1+q)^2*Fp + q^2*Fpp;
        C =(1+q)*F;            
        denoms = [B+sqrt(B^2-4*A*C), B-sqrt(B^2-4*A*C)];
        [~,idx] = max(abs(denoms));            
        dt = -(t-tp)*2*C/denoms(idx);
        % Update history
        tpp = tp;
        Fpp = Fp;        
        Fp = F;
        tp = t;       
        % Update root
        t = t+dt;
        absres = abs(dt);            
        if absres < tol
            converged = true;
            break
        end
    end
    if VERBOSE    
        if converged

            fprintf('Muller converged in %d iterations.\n', iter);
        else
            warning('Muller did not converge after %d iterations. abs(dt)=%g', iter, absres);
            target_point = [x0 y0 z0]
        end            
    end
    troot=t;
end
