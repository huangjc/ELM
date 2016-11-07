function [ r, J, Jt ] = res_Ax( x, A, b )

r = b - A*x;

if nargout == 2 % return matrix
    J = -A;
elseif nargout == 3 % return function handle
    J = @(x) -A*x;
    Jt = @(x) -A'*x;
end