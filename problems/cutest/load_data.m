function record = load_data( DIR, pid, mid )
% pid: np x 1 numeric vector, selected problem id
% mid: nm x 1 numeric vector, selected methods id
% data is assumed to be stored in rec{1:nm} cells

np = length( pid );
nm = length( mid );

record = cell( np,nm );
for i = 1:np
    id = pid(i);
    load( strcat( DIR, num2str(id) )  ); % load rec: nm x 1 cell
    for j = 1:nm
        record{ i, j } = rec{mid(j)};
    end
end