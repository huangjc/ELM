function select_cutest(  )
% classf is a nx4 cell array in the following structure:
%         {name}, {XXXr-XX}, {n-}, {m}
% See CUTE paper for details

CUTEST = '/home/hjc/.linuxbrew/opt/cutest/libexec';
MASTSIF = '/home/hjc/.linuxbrew/opt/mastsif/share/mastsif';
myDataDir = '/home/hjc/Projects/Extended_LM/codes/data';

addpath( strcat(CUTEST,'/src/matlab') );
classf = read_fileContent( strcat(MASTSIF,'/CLASSF.DB') );

% PATH = '/home/hjc/Documents/CUTEst_nonlinear_eqns/CUTEst_default/';

% we select the feasibility instances, i.e, problems without objective
% function, with nonlinear equality constraints and without fixed variables 
% for which the Jacobian matrix is available.
feasible_pid = [];
largeProb_id = [];
smallProb_id = [];
undecided_id = [];
pid = [];

for i = 1:size( classf, 1 )
    % fprintf('%s   %s   %s   %s\n', classf{i,:} );
    % skip problems that are not ( no objective &&  not constant objective
    % function )
    if ~strcmp( classf{i,2}(1), 'N' ) && ~strcmp( classf{i,2}(1), 'C' )
        continue;
    else
        % skip unconstrained, fixed constraints only and bounded
        % constraints only
        if strcmp( classf{i,2}(2), 'U' ) || strcmp( classf{i,2}(2), 'X' ) || strcmp( classf{i,2}(2), 'B' ) 
            continue;
        else
            % skip irregular problem (first and second derivative does not exist or discontinuous )
            if ~strcmp( classf{i,2}(3), 'R' )
                % fprintf('%s   %s', classf{i,1}, classf{i,2} );
                continue;
            else
                % degree of the highest analytic derivatives (0,1,or 2)
                if strcmp( classf{i,2}(4), '0' )
                    % fprintf('%s   %s', classf{i,1}, classf{i,2} );
                    continue;
                else
                    feasible_pid = [ feasible_pid, i ];
                    
%                     % check some properties
%                     cd(strcat( PATH, classf{i,1} ));
%                     prob = cutest_setup();
% 
%                     if any(prob.equatn ~= 1)
%                         all_equatn = '.'; %0
%                     else
%                         all_equatn = '*'; %1
%                     end
% 
%                     if any(prob.linear == 1)
%                         any_linear = '*'; %1
%                     else
%                         any_linear = '.'; %0
%                     end
% 
%                     if all(prob.linear == 1)
%                         all_linear = '*'; %1
%                     else
%                         all_linear = '.'; %0
%                     end
% 
%                     fprintf( '%-8s \t %s \t %s \t %s \t %s \t %s \t %s \n', ...
%                             classf{i,1}, classf{i,2}, classf{i,3}, classf{i,4}, all_equatn, any_linear, all_linear );
% 
%                     cutest_terminate();
                        
                    % check problem size
                    str3 = strrep( classf{i,3}, '-', '' );
                    % V: user-specified, or given number of constraints (fixed not included)
                    if strcmp( str3, 'V' ) || ( str2double( str3 ) >= 1000 && str2num( str3 ) <= Inf )
                        largeProb_id = [ largeProb_id, i ];
%                         copyfile(strcat(MASTSIF,'/',classf{i,1},'.SIF'),...
%                             '/home/hjc/Documents/CUTEst_nonlinear_eqns/CUTEst_nonlinear_eqns_default_1000/00_SIF_data/');
% 
%                         copyfile(strcat('/home/hjc/Documents/CUTEst_nonlinear_eqns/CUTEst_default/',classf{i,1}),...
%                             strcat('/home/hjc/Documents/CUTEst_nonlinear_eqns/CUTEst_nonlinear_eqns_default_1000/',classf{i,1}));
                    end
                    if strcmp( str3, 'V' ) || ( str2double( str3 ) <= 1000 )
                        smallProb_id = [ smallProb_id, i ];
                    end
                end
            end
        end
    end
end

fprintf('Feasible instances: %d\n', length(feasible_pid) );
for l = 1:length(feasible_pid)
    fprintf( '%-8s\n', classf{feasible_pid(l),1} );
end


fprintf('Large (>=1000) problems: %d\n', length(largeProb_id) );
for l = 1:length(largeProb_id)
    fprintf( '%-8s\n', classf{largeProb_id(l),1} );
end

fprintf('Small (<=1000) problems: %d\n', length(smallProb_id) );
for l = 1:length(smallProb_id)
    fprintf( '%-8s\n', classf{smallProb_id(l),1} );
end

% check whether all filename have been selected
filename = read_folderContent( strcat(myDataDir,'/CUTEst_nonlinear_eqns_small'), '.SIF' );
check = zeros(length(filename),1); % 53 + 00_SIF_data

for l = 1:length(feasible_pid)
    i = feasible_pid(l);
    for j = 1:length(filename)
        if strcmp(filename{j}, classf{i,1})
            check(j) = 1;
            pid = [ pid, i ];
            break;
        end
    end
end
                    
fprintf('Problems not included in feasible instances:\n')
q = find(check==0);
for i = 1:length(q)
    disp(filename{q(i)});
end

fprintf('Extra problems:\n') % not in that 53 problems
extras = setdiff( feasible_pid, pid );
for i = 1:length(extras)
    fprintf( '%-8s \t %s \t %s \t %s\n', ...
            classf{extras(i),1}, classf{extras(i),2}, classf{extras(i),3}, classf{extras(i),4} );
end




