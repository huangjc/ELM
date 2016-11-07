function [ f, g, Bx, d ] = phi_logistic( w, X, y )
% X(samples, feature)
% y(samples, 1)
% w(feature, 1)
% option = 'f','g','H','All'

% [n,m] = size(X);

% not recommended to formulate from rho_logistic() and res_logistic(), too inefficient.

Xw = X*w;
% Xw = Xw + mean(Xw); % normalized 
yXw = y.*Xw;

f = sum( log( 1+exp(-yXw)));
    
if nargout > 1
    if nargout > 2
        sig = 1./(1+exp(-yXw));
        g = -X.'*(y.*(1-sig));
    else
        g = -(X.'*(y./(1+exp(yXw))));
    end
end

if nargout > 2
    % B = X.'*diag(sparse(sig.*(1-sig)))*X; 
    % Bx = @(x) B*x;
    Bx = @(x) X.'*( (sig.*(1-sig)).* (X*x) );
    if nargout > 3 
        % compute the diagonals of B
        % d = diag(B);
        d = sum( ((X.').^2)*sparse(sig.*(1-sig)), 2 );
    end
end