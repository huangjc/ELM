function Cube = getCube( record, list )
% Cube      np x nm x length(list)
% list      cell of field names, e.g. {'name','n','m'}; etc

[np,nm] = size( record );
Cube = cell(np,nm,length(list));
out = getFieldCell( record, 0, list{:} );
for i = 1:length(list)
    Cube(:,:,i) = out{i};
end