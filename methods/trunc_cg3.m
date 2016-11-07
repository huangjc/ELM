function [ p, info ] = trunc_cg3( B, g, tol, maxit )
%% with function handle B, note that it accepts g, not -g
% subroutine for Newton-CG, solving 
%           B*x = -g (Newton's equation \nabla^2 f * x = - \nabla f),
% where B is a function handle B(x)=B*x, B can be negative.
% In this subroutine, x0 is specifically set to zero.
% adapted from the CG-steihaug algorithm in section 7.2 on book 
% Nocedal and Wright)

if ~isa(B,'function_handle')
    B = @(x) B*x;
end

if nargin < 4
    maxit = 100;
    if nargin < 3
        tol = 1e-5;
    end
end

rho_hist = zeros(maxit+1,1);

iter = 0;
p = zeros( size(g) ); % x0 = 0;
z = zeros( size(g) );
r = g; % r=A(x0)-b=-b;
d = -g; % thus the 1st step generates the Cauchy point, subsequent steps are no worse than it

% exitflag = 0; % status: normal
tau_fix = 0; % no fix needed

rho = r'*r;
rho_hist(1) = rho; % r0
if sqrt(rho) < tol
    p = z; % = zeros(size(g)); good enough, no need for further search
    exitflag = '0';
else
    while 1
        iter = iter + 1;
        
        w   = B(d);
        dBd = d'*w;
        if dBd <= 0
            
            if iter == 1
                p = d; % the same as p = delta/sqrt(d'*d) *d; when delta = norm(g)
            else 
                p = z;
            end
            tau_fix = -dBd/(d'*d);
            
            exitflag = 'n';
            break;
        end
        alpha = rho/dBd;
        oldz  = z;
        z     = z + alpha*d;
        
        r = r + alpha*w;
        oldrho = rho;
        rho = r'*r;
        if sqrt(rho) < tol || iter >= maxit
            p = z;
            if iter >= maxit
                exitflag = 'i';
            else
                exitflag = 't';
            end
            break;
        end
        
        beta = rho/oldrho;
        d    = -r + beta*d;
    end
end

info.res      = sqrt(rho);
info.iter     = iter;
info.rho_hist = rho_hist(1:iter+1); % r0~rk
info.exitflag = exitflag;
info.tau_fix   = tau_fix;

