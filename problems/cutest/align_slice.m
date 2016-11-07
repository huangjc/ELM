function M = align_slice( Cube, dim )
% get_slice rearrange the 3D array Cube into matrix M along given dimension
% e.g. M = get_slice( Cube, 2 ) results to 
%   M = [M1,M2,...]
% where Mi = Cube(:,i,:);

if nargin < 2 
    dim = 3;
end

[n,m,p] = size(Cube);
switch( dim )
    case 1
        M = zeros(m,p*n);
    case 2
        M = zeros(n,m*p);
    case 3
        M = zeros(p,n*m);
end
if ~all(isnumeric( Cube(:) ))
    M = num2cell(M);
end

for i = 1:size(Cube,dim)
    switch( dim )
        case 1
            M(:,(i-1)*p+(1:p)) = Cube(i,:,:); % m x p
        case 2
            M(:,(i-1)*p+(1:p)) = Cube(:,i,:); % n x p
        case 3
            M(:,(i-1)*m+(1:m)) = Cube(:,:,i); % n x m
    end
end