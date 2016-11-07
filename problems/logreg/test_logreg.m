function test_logreg(problems, method, options)

nprob = length(problems);
nmthd = length(method);

pars = options.pars;
savedata = options.savedata;

for kprob = 1 : nprob
    X      = problems{kprob}.X;
    b      = problems{kprob}.b;
    lambda = problems{kprob}.lambda;
    w0     = problems{kprob}.w0;
    % n_class= problems{kprob}.n_class;
    
    rec = cell(1, nmthd);

    for kmthd = 1:nmthd
        res = @(x) res_x(x);
        rho = @(w) penalizedL2( w, @phi_logistic, lambda, X, b );
        
        pars.X = X;
        pars.b = b;
        pars.lambda = lambda;
        [ x{kmthd}, info{kmthd} ] = algo_batching( method{kmthd}, res, rho, w0, pars );
    end

    % save data
    for kmthd = 1:nmthd
        rec{kmthd}.n        = size(X,2);
        rec{kmthd}.m        = size(X,1);
        rec{kmthd}.method   = method{kmthd};
        rec{kmthd}.x0       = w0;
        rec{kmthd}.x        = x{kmthd};
        rec{kmthd}.info     = info{kmthd};
    end

    if savedata == 1
        save(num2str(kprob),'rec');
    end
end
