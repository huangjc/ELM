function [ p, iter, exitflag ] = trsub_truncCG( B, g, delta, tol, maxit )
% Algorithm based on 
% THE CONJUGATE GRADIENT METHOD AND TRUST REGIONS IN LARGE SCALE OPTIMIZATION 
% by TROND STEIHAUG
% and
% Algorithm 7.2 on Numerical Optimization, Nocedal and Wright


if nargin < 4
    tol = 1e-5;
end

% rho_hist = zeros(k_max,1);

iter = 0;
z = zeros(size(g));
r = g;
d = -g; % thus the 1st step generates the Cauchy point, subsequent steps are no worse than it

rho = r'*r;
if sqrt(rho) < tol
    p = z; % = zeros(size(g)); good enough, no need for further search
    exitflag = '0';
else
    while 1
        iter = iter + 1;
        
        w   = B(d);
        dBd = d'*w;
        if dBd <= 0
            dd = d'*d;
            dz = d'*z;
            zz = z'*z;
            tau = ( -dz + sqrt( dz^2-dd*(zz-delta^2) ))/dd; % positive root of ||z + tau*d||_2 = delta
            p = z + tau*d;
            exitflag = 'n';
            break;
        end
        
        alpha = rho/dBd;
        oldz  = z;
        z     = z + alpha*d;
        if norm(z) >= delta % note that norm(z) is monotonically increasing, by the property of CG
            dd = d'*d;
            dz = d'*oldz;
            zz = oldz'*oldz;
            tau = ( -dz + sqrt( dz^2-dd*(zz-delta^2) ))/dd; % positive root of ||z + tau*d||_2 = delta
            p = oldz + tau*d;
            break;
        end
        
        r = r + alpha*w;
        oldrho = rho;
        rho = r'*r;
        if sqrt(rho) < tol || iter >= maxit
            p = z;
            exitflag = 't';
            break;
        end
        
        beta = rho/oldrho;
        d    = -r + beta*d;
    end
end



