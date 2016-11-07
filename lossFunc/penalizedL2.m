function [nll,g,H,D] = penalizedL2(w,lossFunc,lambda,varargin)
% Adds L2-penalization to a loss function
% (you can use this instead of always adding it to the loss function code)

if nargout <= 1
    [nll] = lossFunc(w,varargin{:});
elseif nargout == 2
    [nll,g] = lossFunc(w,varargin{:});
elseif nargout == 3
    [nll,g,H] = lossFunc(w,varargin{:});
else
    [nll,g,H,D] = lossFunc(w,varargin{:});
end

nll = nll+sum(lambda.*(w.^2));

if nargout > 1
    g = g + 2*lambda.*w;
end

if nargout > 2
    if isscalar(lambda)
        H = @(x) H(x) + 2*lambda*x;
    else
        H = @(x) H(x) + diag(2*lambda)*x;
    end
end

if nargout > 3
    if isscalar(lambda)
        D = D + 2*lambda*ones(length(w),1);
    else
        D = D + 2*lambda;
    end
end