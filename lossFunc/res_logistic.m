function [ r, J, Jt ] = res_logistic( w, X, y )
% NOT TESTED YET

Xw = X*w;
% Xw = Xw + mean(Xw); % normalized 
r = y.*Xw;

if nargout > 1
    J = @(w) y.*(X*w);
    Jt = @(w) X'*(y.*w);
end