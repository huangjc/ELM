function results_shamanskii()

close all

HOME = '/home/hjc/Projects/Extended_LM/codes';
addpath( genpath( HOME ));
datapath = strcat( HOME, '/temp/shamanskii/test0/' );
pid = [1:19]; % rec id
methods = {'ELM_MIX_exact','ELM_MIX_inexact','ELM_MIX_approxJ_exact','ELM_MIX_approxJ_inexact'};

show_results_mgh( datapath, pid, methods );