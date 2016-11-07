function show_results( PATH, pid, methods )

cd(PATH); % direct to PATH to save figures

np = length(pid);
nm = length(methods);
mid = 1:nm;
% record = load_data( PATH, 1:np, 1:nm ); % np x nm cell
record = load_data( PATH, pid, mid ); % np x nm cell

% Res = getFieldCell( record, [], 'info.res' );
% fprintf( [ repmat('%e \t', 1, nm ), '\n'], Res' );

fid = 1; % fopen('result.txt','w+');

list_info = {'name','n','m','info.res','info.time','info.exit','info.nF','info.nJ','info.iter'};
% list_info = {'info.nF', 'info.nFac', 'info.time'};
out2 = getFieldCell( record, 0, list_info{:} );
Cube = cell(np+1,nm,length(list_info));
for i = 1:length(list_info)
        Cube(:,:,i) = [ out2{i}; repmat({''},1,size(out2{i},2)) ]; % create an empty line
end
details = reshape( Cube, size(Cube,1)*size(Cube,2), size(Cube,3) );
printDetails( fid, details, '%s \t %d \t %d \t %e \t %e \t %s \t %d \t %d \t %d ', list_info );


