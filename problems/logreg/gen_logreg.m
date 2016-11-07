function [ X,b ] = gen_logreg( n, m )
% X(samples, feature)
% b(samples, 1)

if mod( m, 2 ) ~= 0
    error('m must be even.');
end

b = ones(m,1);
b( randperm( m, m/2 ) ) = -1; % choose half index to be -1

X = sparse( ( b*ones(1,n) + randn(m,n) ).*( rand(m,n) < 10/n ) );