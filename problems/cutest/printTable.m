function printTable(record)

[np,nm] = size(record);

fid = 1; % fopen('result.txt','w+');
list1 = {'name','n','m'};
Cube1 = getCube( record, list1 );
% Cube1(:,:,1) = getFieldCell( record, 0, 'name' );
% Cube1(:,:,2) = num2cell(num2str(getFieldCell( record, 1, 'n' )));
% Cube1(:,:,2) = num2cell(num2str(getFieldCell( record, 1, 'm' )));
comp1 = squeeze( Cube1(:,1,:) );
    
% % add empty rows (How to print digits and char/string simultaneously?)
% Emptyline = num2cell( repmat(' ',np,size(Cube1,3)) );
% % Emptyline = num2cell( double.empty(np,size(Cube1,3)) );
% for i = 1:nm
%     Cube1n(:,:,i) = Emptyline;
% end
% Cube1n(:,:,floor(nm/2)) = comp1;
% comp1 = align_slice( Cube1n, 1 );
% comp1 = comp1';


% list2 = {'info.exit'};
% out2 = getFieldCell( record, 0, list2{:} );
% Cube = cell(np,nm,length(list2));
% for i = 1:length(list2)
% %     if ~iscell(out2{i})
%         Cube(:,:,i) = num2cell(out2{i});
% %     else
% %         out_w_break = [ out2{i}; repmat({''},1,size(out2{i},2)) ]; % create an empty line
% %         Cube(:,:,i) = out_w_break;
% %     end
% end
% details = reshape( Cube, size(Cube,1)*size(Cube,2), size(Cube,3) );
% printDetails2( fid, details, '%s \t %d \t %d \t %e \t %e \t %s \t %d \t %d \t %d \t %d', list_info );
% % printDetails2( fid, details, '%d/%d/%6.2f', list_info );
% 
% Table = [ Res(:), T(:) ];
% % Table = [ NF(:), NJ(:), T(:), PASS(:) ];
% C = mat2cell( Table, repmat(np, 1, nm) );
% Table = cell2mat( C' );
% Table = num2cell( Table );
% fid = fopen('table.tex','w+');
% % make_table( fid, Table, '', '', repmat( '%d/%d/%4.2f&%d&', 1, nm ) );
% make_table( fid, Table, '', '', repmat( '%6.4e & %6.4f & ', 1, nm ) );
% 

% % CUTEST
% fields = {'name','n','m','info.res','info.exit'};
% 
% Table = cell(np*nm,length(fields));
% for j = 1:length(fields)
%     F = getFieldCell( record, fields{j}, 0 );
%     if any( [13:17] == j )
%         F = cell2mat( F ); 
%         F = F./T*100; % time percentage
%         F = num2cell( F );
%     end
%     F = F(:);
%         
%     for i = 1:np*nm
%         Table{i,j} = F{i};
%     end
% end

Iter = getFieldCell( record, 0, 'info.iter' );
Nfe = getFieldCell( record, 0, 'info.nF' );
Res = getFieldCell( record, 0, 'info.res' );
% Gnorm = getFieldCell( record, 0, 'info.gnorm' );
Time = getFieldCell( record, 0, 'info.time' );
Exit = getFieldCell( record, 0, 'info.exit' );
Emptyline = num2cell( repmat(' ',np,nm) );

Cube2(:,1,:) = Emptyline;
Cube2(:,2,:) = Iter;
Cube2(:,3,:) = Nfe;
Cube2(:,4,:) = Res;
Cube2(:,5,:) = Time;
Cube2(:,6,:) = Exit;


comp2 = align_slice( Cube2, 1 );
comp2 = comp2';

Table = [comp1];
fid = fopen('leftTable.tex','w+');
rsep = repmat('r',1,np);
csep = repmat('c',1,size(Table,2));
rformat = '%8s & %4d & %4d ';
make_table( fid, Table, rsep, csep, rformat );

Table = [comp2];
% Table = [Name(:,1),Emptyline,Res(:,1),Exit(:,1),Res(:,2),Exit(:,2)];
rsep = repmat('r',1,np*nm);
csep = repmat('c',1,size(Table,2));
rformat = '%s & %4d & %4d & %6.4e & %6.4e & %2s';
% rformat = ['%8s & %4d & %4d ',repmat('& %1s & %6.4e ', 1, nm ) ];

fid = fopen('rightTable.tex','w+');
make_table( fid, Table, rsep, csep, rformat );
% make_table( fid, Table );
% 
% % plotFigs( record, 0 );
