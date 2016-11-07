function make_table( fid, T, rsep, csep, rformat )
% turn cell matrix T into simple latex table
% Input: 
% fid:          file id in fprintf
% T:            matrix or cell array, each cell contain one element only.
% rsep, csep:   vectors that store the format of the table, e.g. resp = '-r-rr-'; csep = '|c|cc|';
%               ( rsep: consists of '-', '=' and ' '(white space)
%                       '-'     single border '\hline', 
%                       '='     double border '\hline\hline',
%                       ' '     nothing
%                 csep: the same as the format string in latex )
% rformat:      the format string that defines the data format in fprintf,
%               e.g. fprintf('%d \t %8.4e \n', k, rho); then '%d \t %8.4e \n' is the rformat
%
% Input rsep, csep, rformat are optional, and can be set to ''


[ height, width ] = size(T);

if ~iscell(T) % ismatrix(num2cell(A)) will be 1, though cell is expected...
    T = num2cell(T);
end

if nargin == 2; 
    rsep = '';
    csep = '';
    rformat = '';
end

if isempty( rsep )
    rsep = strcat( repmat('-r',1,height), '-' );
end
if isempty( csep )
    csep = strcat( repmat('|c',1,width), '|' );
end
if isempty( rformat )
    rformat = '';
    for i = 1:width-1
        if ischar( T{1,i} )
            rformat = strcat( rformat, '%s & ');
        elseif isinteger( T{1,i} )
            rformat = strcat( rformat, '%d & ');
        else % assume to be digits
            rformat = strcat( rformat, '%f & ');
        end
    end
    if ischar( T{1,end} )
        rformat = strcat( rformat, '%s');
    elseif isinteger( T{1,i} )
        rformat = strcat( rformat, '%d');
    else % assume to be digits
        rformat = strcat( rformat, '%f');
    end
    rformat = strcat( rformat, '\\\\ \n' ); % avoid invalid call rformat(end-6:end) if rformat is too short
end

if ~strcmp( rformat(end-6:end), '\\\\ \n' )
    rformat = strcat( rformat, '\\\\ \n' );
end

nrow = length( find( rsep == 'r' ));
ncol = length( find( csep == 'c' ));

if ( nrow ~= height ) || ( ncol ~= width )
    error('Table dimension mismatched.');
end



fprintf( fid, strcat(...
    '\\documentclass[english]{article}\n',...
    '\\begin{document}\n',...
    '\\begin{table}\n',...
    '\\begin{tabular}{%s}\n'), csep );
ind = 1;
for i = 1:length( rsep )
    switch( rsep(i) )
        case '-'
            fprintf( fid, '\\hline\n' );
        case '='
            fprintf( fid, '\\hline\\hline\n' );
        case 'r'
            fprintf( fid, rformat, T{ ind, : } );
            ind = ind + 1;
    end
end
fprintf( fid, strcat(...
    '\\end{tabular}\n',...
    '\\end{table}\n',...
    '\\end{document}\n'));


