function [ problems, pars ] = set_logreg( data, datadir )

pars.fid = 1;
pars.eta = 0.001;
pars.p1 = 0.25;
pars.p2 = 0.75;
pars.c1 = 4;
pars.c2 = 0.25;
pars.mu_0 = 1;
pars.mu_min = 1e-12;
pars.mu_max = 1e5;
pars.maxit = 100;
pars.tolX = 1e-9;
pars.tol = 1e-5;
pars.tolG = pars.tol; % same, used in different context
pars.tolFun = 1e-16;
pars.relFun = 1e-16;
pars.cg_tol_opt = 3; % 2: sqrt, 3:\|.\|
% pars.cg_tol = 1e-4;
pars.cg_maxit = 30;
pars.M = 5;
pars.walltime = Inf;
pars.delta = 1e-1; % ELM_MIX
pars.delta0= 1e-0; % TR

switch ( data )
    case 'a9a'
        problems = process_libsvm_data(sprintf('%s/a9a',datadir));    
    case 'real-sim'
        problems = process_libsvm_data(sprintf('%s/real-sim',datadir));   
    case 'news20'
        problems = process_libsvm_data(sprintf('%s/news20.binary',datadir));
    case 'rcv1'
        problems = process_libsvm_data(sprintf('%s/rcv1_train.binary',datadir));
end



function problems = process_libsvm_data( datapath )

[ b, X ] = libsvmread(datapath); 
   
X = [ X, ones( size(X,1), 1 ) ]; % include the bias term, so that wbar'*xbar = w'*x + b

lambda = 0.25; % 2; % 

% n_fold = 5;
n_fold = 0; % skip cross-validation
n_instance = size(X,1);

problems = cell(n_fold+1,1);
for i = 1:n_fold
    one_fold = fix(size(X,1)/n_fold);
    i_training = setdiff([1:n_instance],one_fold*(i-1)+[1:one_fold]);
    Xi = X(i_training,:);
    bi = b(i_training,:);

    problems{i}.n_feature   = size(Xi,2);
    problems{i}.n_instance  = size(Xi,1);
    problems{i}.lambda      = lambda;    
    problems{i}.X           = Xi;
    problems{i}.b           = bi;
    problems{i}.w0          = rand( size(X,2),1);
end

i = n_fold+1;
problems{i}.n_feature   = size(X,2);
problems{i}.n_instance  = size(X,1);
problems{i}.lambda      = lambda;    
problems{i}.X           = X;
problems{i}.b           = b;
problems{i}.w0          = rand( size(X,2),1);


