function results_cutest()

close all;

HOME = '/home/hjc/Projects/Extended_LM/codes';
datapath = strcat(HOME,'/temp/cutest/test1/');

pid = [1:100]; % rec id
methods = {'ELM','ELM_WLS','ELM_FIM','ELM_MIX','LBFGS'};

pid = remove_empty_id( datapath, pid );
show_results_cutest( datapath, pid, methods );
