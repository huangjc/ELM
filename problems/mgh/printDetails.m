function printDetails2( fid, T, rformat, ttl )
% turn cell matrix T into simple latex table
% Input: 
% fid:          file id in fprintf
% T:            matrix or cell array, each cell contain one element only.
% rformat:      the format string that defines the data format in fprintf,
%               e.g. fprintf('%d \t %8.4e \n', k, rho); then '%d \t %8.4e \n' is the rformat
% ttl:          table title in the format {'label1','label2',...}
% Input rformat are optional, and can be set to ''

if nargin < 3; 
    rformat = '';
    if nargin < 4
        ttl = [];
    end
end

[ height, width ] = size(T);

if ~iscell(T) % ismatrix(num2cell(A)) will be 1, though cell is expected...
    T = num2cell(T);
end

if isempty( rformat )
    rformat = '';
    for i = 1:width-1
        if ischar( T{1,i} )
            rformat = strcat( rformat, '%s \t ');
        elseif isinteger( T{1,i} )
            rformat = strcat( rformat, '%d \t ');
        else % assume to be digits
            rformat = strcat( rformat, '%f \t ');
        end
    end
    if ischar( T{1,end} )
        rformat = strcat( rformat, '%s');
    elseif isinteger( T{1,i} )
        rformat = strcat( rformat, '%d');
    else % assume to be digits
        rformat = strcat( rformat, '%f');
    end
end

if ~strcmp( rformat(end-1:end), '\n' )
    rformat = strcat( rformat, '\n' );
end

if ~isempty(ttl)
    fprintf( fid, [repmat('%s \t ',1,width-1),'%s \n'], ttl{:} );
end

for i = 1:height
    fprintf( fid, rformat, T{ i, : } );
end

