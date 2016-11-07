function [ x, info ] = trust_region( funObj, x0, pars )
% trust region algorithm based on Algorithm 2.1 in 
% A review of trust region algorithms for optimization by Ya-Xiang Yuan
% 
% Input:
%     funObj:     objective function [f,J] = funObj( x )
%     x0:         starting point
%     pars:       parameters of the algorithm
% 
% Output:
%     x:          solution x^*
%     fvec:       residual vector at x^*
%
%% Initialization

if nargin < 3
    pars = [];
end

delta0  = pars.delta0;
% tol     = pars.tol;
maxit   = pars.maxit;
cg_maxit= pars.cg_maxit;
fid     = 1;
eta0    = 1e-4;
eta1    = 0.25;
eta2    = 0.75;
sigma1  = 0.25;
sigma2  = 0.5;
sigma3  = 4;
% delta_max = ;

%%% End initialization %%%%%%%%%%%%%%%%%%%%

t0=tic;

x = x0;

gnorm_hist = zeros(maxit, 1);
f_hist = zeros(maxit, 1);
t_hist = zeros(maxit+1,1);

[f, g, B] = funObj(x);  % we update B here
nfeval = 1;
ngeval = 1;

f_hist(1) = f;
gnorm = norm(g,'inf');
gnorm_hist(1) = gnorm;
t_hist(1)     = toc(t0);
% delta = gnorm;  
delta = delta0;

iter = 0;
cg_iter = 0;

fprintf(fid,'%4s \t %8s \t %8s \t %8s \t %8s \t %8s \n', 'iter', 'Fnorm', 'gnorm', 'norm(d)', 'ratio','cg_iter');
fprintf(fid,'%4d \t %8.4e \t %8.4e \t %8.4e \t %8.4e \t %4d\n', iter, f, gnorm, 0, 0, 0);

[ breakflag, exitflag ] = termination_test( f, gnorm, iter, 1, 1, toc(t0), cg_iter, pars );
while ~breakflag

    iter = iter + 1;

    sub_tol = min( 1/2, gnorm ); % tol; 
    [ d, sub_iter ] = trsub_truncCG( B, g, delta, sub_tol, cg_maxit );
    cg_iter = cg_iter + sub_iter;
    
    if sub_iter == 0 && gnorm > pars.tolG
        f_hist(iter+1)     = f; 
        gnorm_hist(iter+1) = gnorm;
        t_hist(iter+1)     = toc(t0);
        
        gnorm = 0.1*gnorm;
        continue;
    end
    
    oldx = x;
    oldf = f;
    oldg = g;
    oldB = B;
    
    x=x+d;
    
    [f, g, B] = funObj(x);    % we update B here
    nfeval = nfeval+1;
    ngeval = ngeval + 1; 
   

    ared = oldf - f;
    pred = -oldg'*d - d'*oldB(d)/2; 
      
    r = ared/pred;
    
    % update rule based on Trust Region Newton Method for Large-Scale 
    % Logistic Regression by
    % Lin,Weng and Keerthi
    if r < eta1
        delta = sigma1*min(norm(d),delta);
    elseif r > eta2 % mu = mu; in [p1,p2]
        delta = sigma3*delta;
    end
%     % based on Algorithm 4.1 in Numerical Optimization, Nocedal and Wright
%     if r < 1/4
%         delta = delta/4;
%     else
%         if r > 3/4 && (norm(d)-delta) < 1e-6 % norm(d)==delta
%             delta = min(2*delta,delta_max);
%         % else delta = delta
%         end
%     end

    if r <= eta0
        x = oldx;
        f = oldf;
    end
    
    f_hist(iter+1) = f;
    gnorm = norm(g,'inf');
    gnorm_hist(iter+1) = gnorm;
    t_hist(iter+1)     = toc(t0);
    
    fprintf(fid,'%4d \t %8.4e \t %8.4e \t %8.4e \t %8.4e \t %4d \n', iter, f, gnorm, norm(d),r, sub_iter);

    [ breakflag, exitflag ] = termination_test( f, gnorm, iter, ared, d, toc(t0), cg_iter, pars );
    if breakflag
        break;
    end
end 

time = toc(t0)


info.exit = exitflag;
info.res = f;
info.gnorm = gnorm;
info.nF = nfeval;
info.nJ = ngeval;
info.time = time;
info.iter = iter;
info.f_hist     = f_hist(1:iter+1);
info.gnorm_hist = gnorm_hist(1:iter+1);
info.t_hist     = t_hist(1:iter+1);
info.pars = pars;
    
    
