function [ problems, pars ] = set_cutest( DataDir, problemSet )

pars.fid = 1;
pars.eta = 0.0001;
pars.p1 = 0.25;
pars.p2 = 0.75;
pars.c1 = 4;
pars.c2 = 0.25;
pars.mu_0 = 1e-5;
pars.mu_min = 1e-8;
pars.maxit  = 2000;
pars.tolX = 1e-9;
pars.tolG  = 1e-10; % tolG
pars.tolFun = 1e-9;
pars.relFun = 1e-9;
pars.cg_tol_opt = 3;
% pars.cg_tol = 1e-16;
pars.cg_maxit = 50;
pars.cg_total_it = 100000;
pars.walltime = 3;
pars.delta0 = 1e-3; % newton_cg_tr
pars.delta = 1e-6; % ELM_tcg_wls6
pars.M = 5;

switch ( problemSet )
    case 'noisy_subset'
        MODEL1 = {'zero','Uniform','Gaussian','Laplace','t','Poisson'}; % 6 in total
        MODEL2 = {'L2','L1','L1-L2','Lp','Fair','Huber','Cauchy','Geman-McClure',...
        'Welsch','Tukey','Andrew','t'}; % 12 in total
        
        noise_model = MODEL1{5};
        rho_model = MODEL2{12};
        
        prob_dir = strcat(DataDir,'le1000_prob/');
        prob_names = read_fileContent( strcat( prob_dir,'problemNames') );
        prob_id = 1:100; 
        prob_id = setdiff( prob_id, [32:34,61,70:71,82:84,92,93] );
        
        model_pars.p = 0.5; % Lp
        model_pars.c = 2; % Fair, Cauchy, Welsch, Tukey, Andrew
        model_pars.k = 2; % Huber
        model_pars.nu= 1; % t
        rho = @(x) rhoSum( x, rho_model, model_pars ); % for simplicity, we assume identical parameter
                
        % assume the nonlinear equation is perturbed by a certain noise        
        noise_pars.sigma = 1; % Gaussian
        noise_pars.b = 1; % Laplace
        noise_pars.nu = 1; % t
        noise_pars.lambda = 1; % Poisson

end

problems.dir = prob_dir;
problems.pid = prob_id;
problems.names = prob_names;
problems.rho = rho;
problems.noise_model = noise_model;
problems.noise_pars = noise_pars;



