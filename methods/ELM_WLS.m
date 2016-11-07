function [ x, info ] = ELM_wls( res, rho, x0, pars )
%% Extended Levenberg-Marquardt method using weighted-least-square-like Hessian approximations
%% Thus nabla^2 rho is of diagonal form
% Input:
%     res        [ r, r_g ] = res( x ), inner function
%     rho        [ f, g, Bx ] = rho( r ), outer function with r as variable, 
%                where Bx is a function handle of approximate Hessian
%     x0         starting point
%     pars       parameters of the algorithm
% 
% Output:
%     x          solution x^*

fid     = 1;
eta     = 0.0001;
p1      = 0.25;
p2      = 0.75;
mu_0    = 1e-5;
mu_min  = 1e-8;
mu_max  = 1e5;
maxit   = 1000;
tol     = 1e-9;
cg_maxit = 50;
%%% Initialization %%%
if nargin < 4
    pars = [];
end

% if some fields are not set
if isfield(pars, 'fid')
    fid = pars.fid;
end
if isfield(pars, 'eta')
    eta = pars.eta;
end
if isfield(pars, 'p1')
    p1 = pars.p1;
end
if isfield(pars, 'p2')
    p2 = pars.p2;
end
if isfield(pars, 'mu_0')
    mu_0 = pars.mu_0;
end
if isfield(pars, 'mu_min')
    mu_min = pars.mu_min;
end
if isfield(pars, 'maxit')
    maxit = pars.maxit;
end
if isfield(pars, 'tol')
    tol = pars.tol;
end
if isfield(pars, 'cg_maxit')
    cg_maxit = pars.cg_maxit;
end


%%% End initialization %%%

t0=tic;
       
x = x0;
mu = mu_0;  

% mu_hist    = zeros(maxiter+1,1);
% ratio_hist = zeros(maxiter,1);
f_hist      = zeros(maxit+1,1);
gnorm_hist  = zeros(maxit+1,1);
t_hist      = zeros(maxit+1,1);
% mu_hist(1) = mu;

[ r, J, Jt ] = res(x); 
[ phi, rho_g ] = rho( r );

thres = eps;
% parameter = mu*norm( rho_g ); % mu*F;
% thres = min( 1e-4, parameter );
w = update_W( rho_g, r, thres );

phi_g = Jt(rho_g);
gnorm = norm(phi_g,'inf');
f_hist(1)     = gather(phi);
gnorm_hist(1) = gnorm;
t_hist(1)     = toc(t0);

nfeval = 1;
nJeval = 1;

iter = 0;

fprintf(fid,'%4s \t %4s \t %8s \t %8s \t %8s \t %8s \t %8s \n', 'iter','nfeval','phi','gnorm','dnorm','ratio','cg_iter');
fprintf(fid,'%4d \t %4d \t %8.4e \t %8.4e \t %8.4e \t %8.4e \t %8.4e \n',...
            iter, nfeval, gather(phi),gnorm,0,1,0);

thres = 1;

cg_iter = 0;

[ breakflag, exitflag ] = termination_test( phi, gnorm, iter, 1, 1, toc(t0), cg_iter, pars );
while ~breakflag

    iter = iter + 1;
    
    parameter = mu*norm( rho_g ); % mu*F;

    fh = @(x) JBJp( x, parameter, J, Jt, w );

    % use cg 
    % forcing term
    switch( pars.cg_tol_opt )
        case 1
            cg_tol = pars.cg_tol;
        case 2
            cg_tol = min( 0.5, sqrt(gnorm) ); % superlinear, Nocedal & Wright
        case 3
            cg_tol = min( 0.5, gnorm ); % quadratic, Nocedal & Wright
    end
     
    [ d, info ] = trunc_cg3( fh, phi_g, cg_tol, cg_maxit );
    k = info.iter;
    cg_res = info.res;
    if k == cg_maxit
        cg_tol = cg_res;
    end
    
    oldx = x;
    oldr = r;
    oldJ = J;
    oldJt = Jt;
    oldphi = phi; 
    oldw = w;
    % oldPrho = P_rho;

    x = x+d;

    [ r, J, Jt ] = res(x); 

    [ phi, rho_g ] = rho( r );
    
    nfeval = nfeval + 1;
    nJeval = nJeval + 1;

    
    ared = oldphi - phi;
    pred = -phi_g'*d - 1/2*d'*JBJ(d,oldJ,oldJt,oldw);
    % pred = -phi_g'*d - 1/2*d'*JBJ(d,oldJ,oldJt,oldw) + cg_tol*norm(d); % no need for the last term

    ratio = ared/pred;
    % ratio_hist(iter) = ratio;

    if ratio < p1
        mu = min( 4*mu, mu_max );
    elseif ratio > p2 % mu = mu; in [p1,p2]
        mu = max( mu/4, mu_min );
    end
    % mu_hist(iter+1) = mu;
    
    if ratio <= eta
        x = oldx;
        r = oldr;
        J = oldJ;
        Jt = oldJt;
        phi = oldphi;
%         w = oldw;
        % P_rho = oldPrho;
    else
        % thres = eps;
        % parameter = mu*norm( rho_g ); % mu*F;
        % thres = min( 1, parameter );
        if norm(x - oldx)/norm(oldx) < sqrt(thres)/100;
            thres = thres/10;
        end
        w = update_W( rho_g, r, thres );
    end
    
    % when r>p0, x has changed, so does fvec;
    % when r<p2, fjac has changed (see notes above)
    if any( x ~= oldx )
        phi_g = Jt(rho_g);
        gnorm = norm(phi_g,'inf');
    end
    
    f_hist(iter+1)     = gather(phi); 
    gnorm_hist(iter+1) = gnorm;
    t_hist(iter+1)     = toc(t0);
    
    fprintf(fid,'%4d \t %4d \t %8.4e \t %8.4e \t %8.4e \t %8.4e \t %d \n',...
                iter,nfeval,gather(phi),gnorm,norm(d),gather(ratio),k);
    
    [ breakflag, exitflag ] = termination_test( phi, gnorm, iter, ared, d, toc(t0), cg_iter, pars );
    if breakflag
        break;
    end
end 

time = toc(t0)

info.nF         = nfeval;
info.nJ         = nJeval;
info.res        = phi;
info.gnorm      = gnorm;
info.exit       = exitflag;
info.time       = time;
info.iter       = iter;
info.pars       = pars;
% info.r_hist    = ratio_hist(1:iter);
% info.mu_hist   = mu_hist(1:iter+1);
info.f_hist     = f_hist(1:iter+1);
info.gnorm_hist = gnorm_hist(1:iter+1);
info.t_hist     = t_hist(1:iter+1);
    
end



%% sub-functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = JBJp(x,p,J,Jt,Bx)
    y = JBJ(x,J,Jt,Bx) + p*x;
end

function y = JBJ(x,J,Jt,Bx)
%    B - function handle provides B(x)=B*x
    w = J(x);
    y = Jt(Bx(w));
end

function w = update_W( rho_g, r, thres )
% modification to avoid division over 0
w = @(x) ( rho_g./(sgn(r).*max(abs(r),thres)) ).*x; 
end

function s = sgn(x)
% slightly different from default sign() function.
% zero has sign 1
    s = ones(size(x));
    s(x<0) = -1;
end

