function Script_logreg()
% problems{i}.lambda greatly affects the performance of LBFGS
% 
% clear all;  close all;  clc;

HOME = '/home/hjc/Projects/Extended_LM/codes';
DataDir = strcat(HOME,'/data/logreg');
addpath(genpath( HOME ));

methods = {'ELM_MIX','LBFGS'};
problem_sets = {'a9a','real-sim','news20','rcv1'};
test_set = problem_sets{1};

[ problems, pars ] = set_logreg( test_set, DataDir );

options.savedata = 1;
options.pars = pars;

t0 = tic;
cd( HOME );
mkdir(strcat('temp/logreg/',test_set));
cd(strcat('temp/logreg/',test_set));
test_logreg(problems, methods, options);
toc( t0 )
