function [ r, J, JT ] = res_x( x )

r = x;

if nargout > 1
    J  = @(x) x;
    JT = @(x) x;
end