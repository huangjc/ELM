function [ x, info ] = algo_batching( method, res, rho, x0, pars )

if nargin < 5
    pars = [];    
end

model = @(x) phi( x, res, rho );

switch ( method )
    case 'ELM_WLS'
        [ x, info ] = ELM_WLS( res, rho, x0, pars );
    case 'ELM_MIX'
        [ x, info ] = ELM_MIX( res, rho, x0, pars );
    case 'ELM_MIX_exact'
        pars.cg_tol_opt = 1; % 1: pars.cg_tol; 2: min( 0.5, gnorm ); 
        pars.cg_tol = 1e-16;
        [ x, info ] = ELM_MIX( res, rho, x0, pars );
    case 'ELM_MIX_inexact'
        pars.cg_tol_opt = 3;
        [ x, info ] = ELM_MIX( res, rho, x0, pars );
    case 'ELM_MIX_approxJ_exact'
        pars.cg_tol_opt = 1; % 1: pars.cg_tol; 2: min( 0.5, gnorm ); 
        pars.cg_tol = 1e-16;
        [ x, info ] = ELM_MIX( res, rho, x0, pars );
    case 'ELM_MIX_approxJ_inexact'
        pars.cg_tol_opt = 3;
        [ x, info ] = ELM_MIX( res, rho, x0, pars );
    case 'ELM'
        [ x, info ] = ELM( res, rho, x0, pars );
    case 'ELM_FIM'
        [ x, info ] = ELM_FIM( res, rho, x0, pars );
    case 'LBFGS'
        addpath(genpath('/home/hjc/Public/Codes/OthersCode/minFunc_2012'));
        tic;
        options = [];
        options.optTol = pars.tolG;
        options.progTol = pars.tolFun;
        options.maxIter = pars.maxit*10;
        options.maxFunEvals = pars.maxit*10;
        options.Method = 'lbfgs';
        options.corrections = pars.M; % lbfgs steps
        options.walltime = pars.walltime;
        [x,f,exitflag,output] = minFunc(model,x0,options);
        t = toc;
        switch( exitflag)
            case 0
                exit = 'i';
            case 1
                exit = 'g';
            case 2
                exit = 'f';
            otherwise
                exit = '-'; % line search failed or stopped by output function
        end
        
        info = struct('res', f,...
                      'iter', output.iterations,...
                      'nF', output.funcCount,...
                      'f_hist', output.trace.fval,...
                      'gnorm_hist', output.trace.optCond,...
                      't_hist', output.trace.t_hist,...
                      'exit',exit,...
                      'time',t,...
                      'pars',options );
end

