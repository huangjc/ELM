function output = read_fileContent( file )
% read content from structured file
% the file has content like this:
%         10FOLDTR   NOR2-AN-        V-        0
% output stores them separately.
% output    mxn cells containing single word

fid = fopen( file );

line = 1;
tline = fgetl(fid);
while ischar(tline)
    C = strsplit(tline);
    while( isempty(C{1}) && length(C)>1 ); C = C(2:end); end % remove empty preceding whitespace
    output(line, 1:length(C)) = C(1:end);
    line = line + 1;
    tline = fgetl(fid);
end

fclose(fid);