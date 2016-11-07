HOME = '/home/hjc/Projects/Extended_LM/codes/';
CUTEST = '/home/hjc/.linuxbrew/opt/cutest/libexec/src/matlab/'; % cutest execution required
addpath( CUTEST );
addpath( genpath(HOME) );
DataDir = strcat( HOME,'data/cutest/');

problemSet = {'noisy_subset'};
chosen_set = problemSet{1};

rng(2016);

[ problems, pars ] = set_cutest( DataDir, chosen_set );
% problems.pid = randi(100,1,2);

methods = {'ELM','ELM_WLS','ELM_FIM','ELM_MIX','LBFGS'};
options.methods = methods;
options.savedata = 1;
options.savedir = strcat( HOME, 'temp/cutest/test/' );
mkdir(options.savedir);
test_cutest( problems, pars, options );

