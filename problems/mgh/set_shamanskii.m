function [ problems, pars ] = set_shamanskii( Nprob, n )
% select problems from those used in published papers
% 7 (problem 35) should be omitted, no one performs well
% Chol has problem with: 8 (problem 27)

% MGH nonlinear equation system
problemID = [1,13,3,14,7,20,35,27,28,29,26,25,30,31]; % problem 1-14, shamanskii choose problem 8-14

MODEL1 = {'zero','Uniform','Gaussian','Laplace','t','Poisson'}; % 5 in total
MODEL2 = {'L2','L1','L1-L2','Lp','Fair','Huber','Cauchy','Geman-McClure',...
'Welsch','Tukey','Andrew','t'}; % 12 in total

% shamanskii problem set
            
l = 1;
if any(Nprob==8)
    % 8 (problem 27): inf when factor = 10, 100
    problems{l}.pid = problemID( 8 );
    problems{l}.n = n;
    problems{l}.factor = 1;
    problems{l}.name = '8';
    problems{l}.noise_model = MODEL1{4};
    problems{l}.rho_model   = MODEL2{6};
    l = l+1; % next waiting index 
end

Nprob = setdiff(Nprob,8); 
factor = [1,10,100];
for i = 1:length(Nprob)
    for j = 1:length(n)
        for k = 1:length(factor)
            problems{l}.pid = problemID( Nprob(i) );
            problems{l}.n = n(j);
            problems{l}.factor = factor(k);
            problems{l}.name = num2str( Nprob(i) );
            problems{l}.noise_model = MODEL1{4};
            problems{l}.rho_model   = MODEL2{5};
            l = l+1; % next waiting index 
        end
    end
end


pars.fid = 1;
pars.eta = 0.0001;
pars.p1 = 0.25;
pars.p2 = 0.75;
pars.p_accept = 0.75;
pars.c1 = 4;
pars.c2 = 0.25;
pars.mu_0 = 1e-5;
pars.mu_min = 1e-8;
pars.maxit  = 2000;
pars.tolX = 1e-9;
pars.tolG  = 1e-10; % tolG
pars.tolFun = 1e-9;
pars.relFun = 1e-6;
% pars.cg_tol_opt = 1;
% pars.cg_tol = 1e-16;
pars.cg_maxit = 20;
pars.walltime = Inf;


