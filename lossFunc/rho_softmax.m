function [ f, g, Hx, d ] = rho_softmax( w, X, y, k )
% w(feature*class,1)
% X(instance,feature)
% y(instance,1)
[n,p] = size(X);
w = reshape(w,[p k]);

Xw = X*w;
Z = sum(exp(Xw),2);
%nll = -sum((sum(X.*w(:,y).',2) - log(Z)));
ind = sub2ind([n k],[1:n]',y);
f = -sum(Xw(ind)-log(Z));

if nargout > 1
   g = zeros(p,k);
   for c = 1:k
      g(:,c) = X'*(exp(Xw(:,c))./Z-(y == c));
   end
   g = reshape(g,[p*k 1]);
end

if nargout > 2
    H = zeros(p*k);
    SM = exp(X*w)./repmat(Z,[1 k]);
    parfor c1 = 1:k
        for c2 = 1:k
            D = SM(:,c1).*((c1==c2)-SM(:,c2));
            H((p*(c1-1)+1):p*c1,(p*(c2-1)+1):p*c2) = X'*diag(sparse(D))*X;
        end
    end
    Hx = @(x) H*x;
    if nargout > 3
        d = diag(H);
    end
end
