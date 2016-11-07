function [ f, g, Bx ] = rho_logistic( r )
% NOT TESTED YET

sig = 1./(1+exp(-r));

f = - sum( log( sig ) );
    
if nargout > 1
        g = -(1-sig);
end

if nargout > 2
    Bx = @(x) (sig.*(1-sig)).*x;
end