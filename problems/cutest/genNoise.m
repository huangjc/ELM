function ns = genNoise( m, p, type, pars )
% generate zero-mean mxp noise vector
% m: number of samples
% p: number of variables/dimensions

switch( type )
    case {'0', 'zero'} % no noise
        ns = zeros(m,p);
    case 'Uniform'
        ns = rand( m,p );
    case 'Gaussian'
        ns = randn( m,p );
    case 'Laplace'
        b = pars.b;
        U = rand(m,p)-1/2;
        ns = - b*(U./abs(U)).*log(1-2*abs(U));
    case 't'
        nu = pars.nu; 
        sigma = pars.sigma;
        if ( length(sigma) == 1 )
            sigma = sigma*eye(p);
        else
            if (length(sigma) ~= p)
                error('Dimension mismatch.');
            end
        end
        ns = mvtrnd( sigma, nu, m );
    case 'Poisson'
        lambda = pars.lambda;
        ns = poissrnd(lambda,m,p);
end


