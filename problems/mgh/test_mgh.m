function test_mgh( problems, rd, methods, options )

nprob = length(problems);
nmthd = length(methods);

rec = cell(1, nmthd);

pars = options.pars;
savedata = options.savedata;

for kprob = 1 : nprob
    pid     = problems{kprob}.pid;
    n       = problems{kprob}.n; % intended n
    factor  = problems{kprob}.factor;
    name    = problems{kprob}.name;
    rho_model = problems{kprob}.rho_model;
    noise_model = problems{kprob}.noise_model;
    
    % pars.maxiter = 100*(n+1);
    
    %====== Initialize problem and parameters ========================    
    fprintf('**************************************************\n');
    fprintf('Loading problem %s: \n', name);

    [x0, solution, fmin, m ] = MGH81problems(pid, n);
%     [x0, solution, fmin, m ] = MGH81problems(pid, n, [], 1);
    
    x0 = factor*x0;
    if ( factor ~= 1 )
        fprintf('Factor=%d. Different starting point from default. May not achieve the same minimum.\n', factor);
    end

    % rank deficient
    switch (rd)
        case 1
            soljac = vecjac( solution, pid,m );
            A = ones( size( soljac, 2 ), 1);
            trans = soljac*A*( (A'*A)\A' );
        case 2
            soljac = vecjac( solution, pid,m );
            A = ones( size( soljac, 2 ), 2);
            A(2*(1:floor( end/2)),2) = -1;
            trans = soljac*A*( (A'*A)\A' );
    end
    
    if rd == 0 % full rank
        fh = @(x) vecfcn( x, pid, m );
        Jh = @(x) vecjac( x, pid, m );
    else % rank deficient
        fh = @(x) ( vecfcn( x, pid, m )-trans*(x-solution) );
        Jh = @(x) ( vecjac( x, pid, m )-trans );
    end
    
    % assume the nonlinear equation is perturbed by a certain noise
    noise_pars.sigma = 1; % Gaussian
    noise_pars.b = 1; % Laplace
    noise_pars.nu = 1; % t
    noise_pars.lambda = 1; % Poisson

    res = @(x) funObj(x, fh, Jh);
    
    r0 = res(x0);
    noise = 0.01*max( [abs(r0); 1] )*genNoise( size(r0,1), 1, noise_model, noise_pars );
        
    %==================================================

    % run tests    
    for kmthd = 1:nmthd
%         % [ x{kmthd}, info{kmthd} ] = algo_batching( method{kmthd}, fh, Jh, x0, options.pars );
%         fJ = @(x) funObj(x, fh, Jh);
%         [ x{kmthd}, info{kmthd} ] = algo_batching( method{kmthd}, fJ, x0, pars );
%             rho = @(x) rhoSum( x, 'gauss' );
%             [ x{kmthd}, info{kmthd} ] = algo_batching( methods{kmthd}, res, rho, x0, pars );

            rho = @(x) rhoSum( x, rho_model );
            res_n = @(x) res_noisy( x, noise, res );
            [ x{kmthd}, info{kmthd} ] = algo_batching( methods{kmthd}, res_n, rho, x0, pars );
    end
        
    % save data
    for kmthd = 1:nmthd
        rec{kmthd}.name     = name;
        rec{kmthd}.pid      = pid;
        rec{kmthd}.n        = length(x0);
        rec{kmthd}.m        = m;
        rec{kmthd}.factor   = factor;
        rec{kmthd}.method   = methods{kmthd};
        rec{kmthd}.x0       = x0;
        rec{kmthd}.r0       = r0;
        rec{kmthd}.solution = solution;
        rec{kmthd}.fmin     = fmin;
        rec{kmthd}.x        = x{kmthd};
        % rec{kmthd}.fvec     = fvec{Nmthd(kmthd)};
        rec{kmthd}.info     = info{kmthd};
    end
    
    if savedata == 1
        save(num2str(kprob),'rec');
    end
 
end

% function [f,J] = funObj(x, fh, Jh)
%     f = fh(x);
%     if nargout > 1
%         J = Jh(x);
%     end

function [ r, J, Jt ] = funObj(x, fh, Jh)
if nargout == 1
    r = fh(x);
elseif nargout == 2
    r = fh(x);
    J = Jh(x);
elseif nargout == 3
    r = fh(x);
    Jm = Jh(x);
    J = @(x) Jm*x;
    Jt= @(x) Jm'*x;
end

function [ r, J, Jt ] = res_noisy( x, noise, res )
% add noise
if nargout == 1
    r = res( x );
elseif nargout == 2
    [ r, J ] = res( x );
elseif nargout == 3
    [ r, J, Jt ] = res( x );
end
r = r + noise;


