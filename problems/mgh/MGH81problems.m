function [x0, solu, fmin, m, name ] = MGH81problems(nprob, n, m, findSolution )
% fmin = sum(f(:).^2)
%% size of the problem:
%       x	f
%       n	m
% 1     2	2
% 2     2	2
% 3     2	2
% 4     2	3
% 5     2	3
% 6     2	>=
% 7     3	3
% 8     3	15
% 9     3	15
% 10	3	16
% 11	3	3~100
% 12	3	>=
% 13	4	4
% 14	4	6
% 15	4	11
% 16	4	>=
% 17	5	33
% 18	6	>=
% 19	11	65
% 20	2~31	31
% 21	even	=
% 22	mod4	=
% 23	any	n+1
% 24	any	2n
% 25	any	n+2
% 26	any	=
% 27	any	=
% 28	any	=
% 29	any	=
% 30	any	=
% 31	any	=
% 32	any	>=
% 33	any	>=
% 34	any	>=
% 35	any	>=



%% %%%%%% Setting parameters %%%%%%

% predefined/suggested n, m
switch nprob
    case {1,2,3,4,5,6}
        n0 = 2;
    case {7,8,9,10,11,12}
        n0 = 3;
    case {13,14,15,16}
        n0 = 4;
    case 17 
        n0 = 5;
    case 18
        n0 = 6;
    case 19
        n0 = 11;
    case 20 % n = 2~31
        n0 = 12; % has fmin
    case 21 % even
        n0 = 2;
    case 22 % mod(n,4)=0
        n0 = 4;
    case {23,24} % any
        n0 = 4; % has fmin
    case 35 % any        
        n0 = 3; % has fmin
    otherwise
        n0 = 2; % any
end

if nargin < 2 || isempty(n)
    n = n0;
else
    if nprob < 20 % fixed
        if ( n ~= n0 )
            warning('n is fixed. Set to fixed value.');
        end
        n = n0;
    else
        % use given value with check
        if nprob == 20
            if ~(n>=2&&n<=31);      error('n is not allowed.');     end
        elseif nprob == 21
            if mod( n,2 ) ~= 0;     error('n is not allowed.');     end
        elseif nprob == 22
            if mod( n,4 ) ~= 0      error('n is not allowed.');     end
        end
    end
end

switch nprob
    case {1,2,3,7,13,21,22,26,27,28,29,30,31} % m=n
        m0 = n;
    case {4,5}
        m0 = 3;
    case 6 % m >= n
        m0 = 10; % has fmin
    case {8,9}
        m0 = 15;
    case 10
        m0 = 16;
    case 11 % 3~100
        m0 = n;% arbitrary
    case {12,18,32,33,34} % m >= n
        m0 = n; % arbitrary
    case 14
        m0 = 6;
    case 15
        m0 = 11;
    case 16 % m >= n
        m0 = 20; % has fmin
    case 17
        m0 = 33;
    case 19 
        m0 = 65;
    case 20
        m0 = 31;
    case 23
        m0 = n+1;
    case 24
        m0 = 2*n;
    case 25
        m0 = n+2;
    case 35 % m >= n
        m0 = n; % has fmin 
end

if nargin < 3 || isempty(m)
    m = m0;
else
    % use fixed/given value with check
    if any( nprob == [1:5,7:10,13:15,17,19:20] ) % fixed
        if ( m ~= m0 )
            warning('m is fixed. Set to fixed value.');
        end
        m = m0;
    elseif any( nprob == [6,12,16,18,32:35] ) % m >= n
        if ( m < n )
            error('m is not allowed.');
        end
    elseif nprob == 11
        if m <3 || m > 100
            error('m is not allowed.');
        end
    else % m=n(=m0) or other fixed value
        if ( m ~= m0 )
            error('m is not allowed.');
        end
    end
end

if nargin < 4
    findSolution = 0;
end


%% %%%%%% Load problems %%%%%%

% fprintf('**************************************************\n');
% fprintf('Loading problem %d: \n', nprob);

switch nprob
    % case 1-35 are examples from more81 
        
    case 1 % Rosenbrock function (more81-1)
        name = 'Rosenbrock';
        % n = 2; 
        % m = 2;
        x0 = [-1.2, 1];
        solu = [1, 1];
        fmin = 0;
    
    case 2 % Freudenstein and Roth function 
        name = 'Freudenstein and Roth';
        % n = 2; 
        % m = 2;
        x0 = [0.5, -2];
        solu = [5, 4];
        fmin = 0;
    
    case 3 % Powell badly scaled function (more81-3)
        name = 'Powell badly scaled';
        % n = 2; 
        % m = 2;
        x0 = [0, 1];
        solu = [1.098e-5, 9.106];
        fmin = 0;
        
    case 4 % Brown badly scaled function 
        name = 'Brown badly scaled';
        % n = 2; 
        % m = 3;
        x0 = [1, 1];
        solu = [1e6, 2e-6];
        fmin = 0;
        
    case 5 % Beale function
        name = 'Beale';
        % n = 2; 
        % m = 3;
        x0 = [1, 1];
        solu = [3, 0.5];
        fmin = 0;
        
    case 6 % Jenrich and Sampson function
        name = 'Jenrich and Sampson';
        % n = 2; 
        x0 = [0.3, 0.4];
        if m == 10;
            solu = [0.2578, 0.2578];
            fmin = 124.362;
        else
            solu = [];
            fmin = [];
        end
        
    case 7 %  Helical valley function (more81-7)
        name = 'Helical valley';
        % n = 3; 
        % m = 3;
        x0 = [-1, 0, 0];
        solu = [1, 0, 0];
        fmin = 0;
        
    case 8 % Bard function
        name = 'Bard';
        % n = 3; 
        % m = 15;
        x0 = ones(1,3);
        solu = [];
        fmin = 8.21487e-3;
        
    case 9 % Gaussian function
        name = 'Gaussian';
        % n = 3; 
        % m = 15;
        x0 = [0.4, 1, 0];
        solu = [];
        fmin = 1.12793e-8;
        
    case 10 % Meyer function
        name = 'Meyer';
        % n = 3; 
        % m = 16;
        x0 = [0.02, 4000, 250];
        solu = [];
        fmin = 87.9458;
        
    case 11 % Gulf research and development function
        name = 'Gulf research and development';
        % n = 3; 
        x0 = [5, 2.5, 0.15];
        solu = [50, 25, 1.5];
        fmin = 0;
        
    case 12 % Box three-dimensional function
        name = 'Box three-dimensional';
        % n = 3; 
        warning('Infinite solution.\nx=[1, 10, 1],[10, 1, -1], and x1=x2 and x3=0.\n');
        x0 = [0, 10, 20];
        solu = [1, 10, 1];
        fmin = 0;
        
    case 13 %  Powell singular function (more81-13)
        name = 'Powell singular';
        % n = 4; 
        % m = 4;
        x0 = [3, -1, 0, 1];
        solu = zeros(1,4);
        fmin = 0;
        
    case 14 %  Wood function (more81-14)
        name = 'Wood';
        % n = 4; 
        % m = 6;
        x0 = [-3, -1, -3, -1];
        solu = ones(1,4);
        fmin = 0;

    case 15 % Kowalik and Osborne function
        name = 'Kowalik and Osborne';
        % n = 4; 
        % m = 11;
        x0 = [0.25, 0.39, 0.415, 0.39];
        solu = [];
        fmin = 3.07505e-4;
      
    case 16 % Brown and Dennis function
        name = 'Brown and Dennis';
        % n = 4; 
        x0 = [25, 5, -5, -1];
        solu = [];
        if m == 20;
            fmin = 85822.2;
        else
            fmin = [];
        end
        
    case 17 % Osborne 1 function
        name = 'Osborne 1';
        % n = 5; 
        % m = 33;
        x0 = [0.5, 1.5, -1, 0.01, 0.02];
        solu = [];
        fmin = 5.46489e-5;
        
    case 18 % Biggs EXP6 function
        name = 'Biggs EXP6';
        % n = 6; 
        % f = 5.65565e-3 if m == 13 ??
        x0 = [1, 2, 1, 1, 1, 1];
        solu = [1, 10, 1, 5, 4, 3];
        fmin = 0;
        
    case 19 % Osborne 2 function 
        name = 'Osborne 2';
        % n = 11; 
        % m = 65;
        x0 = [1.3, 0.65, 0.65, 0.7, 0.6, 3, 5, 7, 2, 4.5, 5.5];
        solu = [];
        fmin = 4.01377e-2;

    case 20 %  Watson function (more81-20)
        name = 'Watson';
        % m = 31;
        % if 2 > n || n > 31
        %     error('n is from 2 to 31.');
        % end
        x0(1:n) = 0;
        solu = [];
        if n == 6
            fmin = 2.28767e-3;
        elseif n == 9 
            fmin = 1.39976e-6;
        elseif n == 12
            fmin = 4.72238e-10;
        else 
            fmin = [];
        end
        
    case 21 % Extended Rosenbrock function
        name = 'Extended Rosenbrock';
        % Not sure what x_(i-1) is when i = 1, assume it is equivalent to
        % Rosenbrock function when n=2
        % so f(2*i-1) = 10*( x(2*i) - x(2*i-1)^2 )?
        % if mod(n,2) ~= 0
        %     error('n should be even.');
        % end
        % m = n;
        
        x0 = ones(1,n);
        for i = 1:n
            if mod(i,2) == 1
                x0(i) = -1.2;
            end
        end
        solu = ones(1,n);
        fmin = 0;
        
    case 22 % Extended Powell singular function
        name = 'Extended Powell singular';
        % m = n;
        % if mod(n,4) ~= 0
        %     error('n is not divided by 4');
        % end
        x0 = repmat( [3, -1, 0, 1], 1, n/4 );
        solu = zeros(1,n);
        fmin = 0;
        
    case 23 % Penalty function I
        name = 'Penalty I';
        % m = n + 1;
        x0 = [1:n];
        solu = [];
        if n == 4
            fmin = 2.24997e-5;
        elseif n == 10 
            fmin = 7.08765e-5;
        else 
            fmin = [];
        end
        
        
    case 24 % Penalty function II
        name = 'Penalty II';
        % m = 2*n;
        x0 = ones(1,n)/2;
        solu = [];
        if n == 4
            fmin = 9.37629e-6;
        elseif n == 10 
            fmin = 2.93660e-4;
        else 
            fmin = [];
        end
         
    case 25 %  Variably dimensioned function (more81-25)
        name = 'Variably dimensioned';
        % m = n + 2;
        x0 = 1 - [1:n]/n;
        solu = ones(1,n);
        fmin = 0;
        
    case 26 %  Trigonometric function (more81-26)
        name = 'Trigonometric';
        % m = n;
        x0 = ones(1,n)/n;
        solu = [];
        fmin = 0;
        
    case 27 %  Brown almost-linear function (more81-27)
        name = 'Brown almost-linear';
        % m = n;
        x0 = ones(1,n)/2;
        warning('Multiple solutions: f=0 at (a,...,a,a^(1-n)) with n*a^n - (n+1)*a^(n-1) + 1=0; in particular, a=1.\n f=1 at (0,...,0,n+1). \n');
        solu = zeros(1,n);
        solu(n) = n+1;
        fmin = 1;

    case 28 %  Discrete boundary value function (more81-28)
        name = 'Discrete boundary value';
        % m = n;
        t = [1:n]/(n+1);
        x0 = t.*(t-1);
        solu = [];
        fmin = 0;

    case 29 %  Discrete integral equation function (more81-29)
        name = 'Discrete integral equation';
        % m = n;
        t = [1:n]/(n+1);
        x0 = t.*(t-1);
        solu = [];
        fmin = 0;

    case 30 %  Broyden tridiagonal function (more81-30)
        name = 'Broyden tridiagonal';
        % m = n;
        x0 = -ones(1,n);
        solu = [];
        fmin = 0;

    case 31 %  Broyden banded function (more81-31)
        name = 'Broyden banded';
        % m = n;
        x0 = -ones(1,n);
        solu = [];
        fmin = 0;
        
    case 32 % Linear function - full rank
        name = 'Linear function - full rank';
        x0 = ones(1,n);
        solu = -ones(1,n);
        fmin = m-n;
        
    case 33 % Linear function - rank 1
        name = 'Linear function - rank 1';
        x0 = ones(1,n);
        warning('Multiple solutions: f=m*(m-1)/(2*(2*m+1)) at any point satisfies sum([1:n].*x)=3/(2*m+1). ');
        solu = [];
        fmin = [];
        
    case 34 % Linear function - rank 1 with 0 columns and rows
        name = 'Linear function - rank 1 with 0 columns and rows';
        x0 = ones(1,n);
        warning('Multiple solutions: f=(m^2+3*m-6)/(2*(2*m-3)) at any point satisfies sum([2:m-1].*x(2:m-1))=3/(2*m-3). ');
        solu = [];
        fmin = [];
        
    case 35 % Chebyquad function (more81-35)
        name = 'Chebyquad';
        x0 = [1:n]/(n+1);
        solu = [];
        if m == n && n <=10
            if n == 8
                fmin = 3.51687e-3;
            elseif n == 10 
                fmin = 6.50395e-3;
            else
                fmin = 0;
            end
        else 
            fmin = [];
        end
end

x0 = x0';
solu = solu';

%% %%%%%% Find minimum if not given %%%%%%
if isempty(solu) && findSolution == 1
    fprintf('Solution unknown. Searching an approximate one...\n');
    
    fh = @(x) vecfcn( x, nprob, m );
    Jh = @(x) vecjac( x, nprob, m );
    pars.fid = fopen('foo','w+'); % dumb way to suppress output
    pars.maxit= 1000;
    pars.tolX = 1e-9;
    pars.tol  = 1e-10; % tolG
    pars.tolFun = 1e-9;
    pars.relFun = 1e-6;
    pars.cg_tol_opt = 1;
    pars.cg_tol = 1e-16;
    pars.walltime = Inf;
    addpath('/home/hjc/Public/Codes/ML_estimation/methods/') % modifiedLM()
    [ solu, info ] = modifiedLM_fh( @(x) funObj(x,fh,Jh ), x0, pars );
    delete foo;

    if strcmpi( info.exit, 'i' )
        disp('Maximal iteration reached.');
    else 
        disp('Solution found.');
    end

    if isempty(fmin)
        fmin = norm( fh(solu) )^2;
        fprintf('Exact minimum is not known. Give an approximate one: \|f\|=%6.4e. (Solution does not guarantee minimum!)\n',fmin);        
    elseif fmin==0
        err = abs( fmin - norm( fh(solu) )^2 );
        if  err > 1e-3
            fprintf('Absolute error: %6.4e\n',err);
            error('Minimum not achieved.');
        end
    else
        relerr = abs( fmin - norm( fh(solu) )^2 ) / abs(fmin);
        if  relerr > 0.05
            fprintf('Relative error: %6.4e\n',relerr);
            error('Minimum not achieved.');
        end
    end

    if ( norm( fh(solu) ) > 1e-5 )
        warning( 'fvec is not close to zero.');
    end
end

function [f,J]=funObj(x, fh, Jh )
f = fh(x);
J = Jh(x);

