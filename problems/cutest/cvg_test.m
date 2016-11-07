function PASS = cvg_test( record )

Res = getFieldCell( record, [], 'info.res' );
T = getFieldCell( record, [], 'info.time' );
[np,nm] = size( Res );
tau = 1e-2;

% % print results
% pid = getFieldCell( record, [], 'pid' );
% pid = pid(:,1);
% fprintf( [ '%d \t', repmat('%6.4e \t', 1, nm ), repmat('%6.2f \t', 1, nm ),...
%         '\n'], [pid,Res,T]' );

%--------------------------------------------------------------------------
% convergence test

% residual
if ~isfield( record, 'fmin' )
    minVal = min( Res, [], 2 );
else
    minVal = getFieldCell( record, [], 'fmin' );
    minVal = minVal(:,1);
end

% if isscalar(minVal)
%     minVal = minVal*ones(np,1); % usually all 0s
% end

for i = 1:np
    for j = 1:nm        
        if ( Res(i, j) - minVal(i) )/max( 1,  minVal(i) ) <= tau
%         if ( Res(i, j) ) < tau
            PASS(i,j) = 1;
        else
            PASS(i,j) = 0;
        end
    end
end

% % Time less than 30 seconds
% PASS(:,:,2) = T<30;
% 
% PASS = all(PASS,3);

