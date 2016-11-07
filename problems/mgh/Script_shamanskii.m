% clear all;  close all;  
clc;

PATH = '/home/hjc/Projects/Extended_LM/codes';
addpath(genpath( PATH ));

Nprob = [8:14]; % rd = 0


n = 5000;

[ problems, pars ] = set_shamanskii( Nprob, n );

methods = {'ELM_MIX_exact','ELM_MIX_inexact','ELM_MIX_approxJ_exact','ELM_MIX_approxJ_inexact'};

options.savedata = 1;
options.pars = pars;


cd( PATH );
mkdir temp/shamanskii/test0;
cd temp/shamanskii/test0;
test_mgh( problems, 0, methods, options );
