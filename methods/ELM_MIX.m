function [ x, info ] = ELM_tcg_wls6( res, rho, x0, pars )
%% use Hessian_rho-vector product to compute exact Hessian in cg
% Note: You need \nabla \rho(r) to decide whether to use weighted least square or not.
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
% mu_max  = 1e5;
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

if ~isfield(pars, 'delta')
    pars.delta = 0.1;
end
if ~isfield(pars, 'l2_penal')
    pars.l2_penal = 0;
end

delta = pars.delta;
penal = pars.l2_penal;

%%% End initialization %%%

t0=tic;
       
x = x0;
mu = mu_0;  

mu_hist    = zeros(maxit+1,1);
% ratio_hist = zeros(maxit,1);
f_hist      = zeros(maxit+1,1);
gnorm_hist  = zeros(maxit+1,1);
t_hist      = zeros(maxit+1,1);
mu_hist(1) = mu;

[ r, J, Jt ] = res(x); 
if min(abs(r)) < delta
    flag = 'ELM';
    [ rho_f, rho_g, rho_Bx ] = rho( r ); % rho_f is also phi itself
else
    flag = 'WLS';
    [ rho_f, rho_g ] = rho( r );
    rho_Bx = @(x) ( rho_g./r ).*x;
end
f = fp(rho_f,penal,x);
g = gp(Jt,rho_g,penal,x);
B = @(x) Bpx(x,J,Jt,rho_Bx,penal);

gnorm = norm(g,'inf');
f_hist(1)     = f;
gnorm_hist(1) = gnorm;
t_hist(1)     = toc(t0);

nfeval = 1;
nJeval = 1;

iter = 0;

fprintf(fid,'%4s \t %8s \t %8s \t %8s \t %8s \t %8s \t %8s \t %8s \t %8s \n',...
            'iter','phi','gnorm','dnorm','ratio','flag','mu','exitflag','cg_iter');
% fprintf(fid,'%4d \t %8.4e \t %8.4e \t %8.4e \t %8.4e \n',iter,norm(r),gnorm,0,1);
fprintf(fid,'%4d \t %8.4e \t %8.4e \t %8.4e \t %8.4e \t %8s \t %8.4e \t %8s \t %8d \n',...
            iter, f,gnorm,0,1, flag, mu_0, ' ', 0);

t_cg = 0;
t_fJ = 0;
t_res = 0;
info.tau_fix = 0;
info.exitflag = 0;

cg_iter = 0;
[ breakflag, exitflag ] = termination_test( f, gnorm, iter, 1, 1, toc(t0), cg_iter, pars );
while ~breakflag

    iter = iter + 1;
    
    tau = mu*norm( rho_g ); % mu*F;

    Btau = @(x) Btaux( x, B, tau );

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
    
    % trunc_cg
    [ d, info ] = trunc_cg3( Btau, g, cg_tol, cg_maxit );
    
    if strcmp(info.exitflag, '0') && gnorm > pars.tolG
        f_hist(iter+1)     = f; 
        gnorm_hist(iter+1) = gnorm;
        t_hist(iter+1)     = toc(t0);
        
        gnorm = 0.1*gnorm;
        continue;
    end
    
    if info.tau_fix ~= 0;
        mu = mu + info.tau_fix/norm( rho_g ) + mu_min; % mu_fix is the LM parameter of the corresponding TR radius
    end
    k = info.iter;  
   
    cg_iter = cg_iter + k;

    oldx = x;
    oldf = f; 
    oldg = g;
    oldB = B;

    x = x+d;

    tic;
    r = res(x); 
    
    if min(abs(r)) < delta
        flag = 'ELM';
        [ rho_f, rho_g, rho_Bx ] = rho( r );
    else
        flag = 'WLS';
        [ rho_f, rho_g ] = rho( r );
        rho_Bx = @(x) ( rho_g./r ).*x;
    end
    nfeval = nfeval + 1;
    
    f = fp(rho_f,penal,x);

    ared = oldf - f; % 
    pred = -oldg'*d - 1/2*(d'*oldB(d)); 

    ratio = ared/pred;

    if ratio < p1
        mu = 4*mu; % min( 4*mu, mu_max );
        % delta = delta/4;
    elseif ratio > p2 % mu = mu; in [p1,p2]
        mu = max( mu/4, mu_min );
        % delta = delta*4;
    end
    mu_hist(iter+1) = mu;

    if ratio <= eta
        x = oldx;
        f = oldf;
        g = oldg;
        B = oldB;
        % P_rho = oldPrho;
    end

    if ratio > eta % x has changed
        [ ~, J, Jt ] = res(x); % evaluation of J suppose to be much more expensive than r
        nJeval = nJeval + 1;
    end

    if ratio > eta
        g = gp(Jt,rho_g,penal,x);
        B = @(x) Bpx(x,J,Jt,rho_Bx,penal);
        gnorm = norm(g,'inf');
    end

    f_hist(iter+1)     = f; 
    gnorm_hist(iter+1) = gnorm;
    t_hist(iter+1)     = toc(t0);

%     fprintf(fid,'%4d \t %8.4e \t %8.4e \t %8.4e \t %8.4e \n',iter,norm(r),gnorm,norm(d),ratio);
    fprintf(fid,'%4d \t %8.4e \t %8.4e \t %8.4e \t %8.4e \t %4s \t %8.4e \t %s \t %d \n',...
                iter,f,gnorm,norm(d),ratio,flag,mu,info.exitflag,k);

    [ breakflag, exitflag ] = termination_test( f, gnorm, iter, ared, d, toc(t0), cg_iter, pars );
    if breakflag
        break;
    end
end 

time = toc(t0)
info.nF         = nfeval;
info.nJ         = nJeval;
info.res        = f;
info.gnorm      = gnorm;
info.exit       = exitflag;
info.time       = time;
info.iter       = iter;
info.cg_iter    = cg_iter;
info.pars       = pars;
% info.r_hist    = ratio_hist(1:iter);
% info.mu_hist   = mu_hist(1:iter+1);
info.f_hist     = f_hist(1:iter+1);
info.gnorm_hist = gnorm_hist(1:iter+1);
info.t_hist     = t_hist(1:iter+1);
    
end



%% sub-functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = fp(rho_f,penal,x)
    f = rho_f + penal*(x'*x);
end

function g = gp(Jt,rho_g,penal,x)
    g = Jt(rho_g) + 2*penal*x;
end

function y = Bpx(x,J,Jt,Bx,penal)
%    B - function handle provides B(x)=B*x
    w = J(x);
    y = Jt(Bx(w)) + 2*penal*x;
end

function y = Btaux(x,B,tau)
    y = B(x) + tau*x;
end

function s = sgn(x)
% zero has sign 1
    s = ones(size(x));
    s(x<0) = -1;
end