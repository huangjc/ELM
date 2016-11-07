function [ breakflag, exitflag ] = termination_test( fval, gnorm, iter, fdiff, d, t, cg_iter, pars )
% a unified test for all methods

% cutest
%%% a more complicated one
% tolX = pars.tolX;
tolG = pars.tolG;
tolFun = pars.tolFun;
relFun = pars.relFun;
% maxit = pars.maxit;
% cg_total_it = pars.cg_total_it
walltime = pars.walltime;

breakflag = 0;
exitflag = '';
if ~( gnorm>tolG && abs(fval) > tolFun && abs(fdiff)>relFun && t <= walltime )
            % && cg_iter < cg_total_it
            % && iter<maxit && abs(fdiff)>relFun  % 
            % && norm(d)>tolX  % cannot use this, since turnc_cg() may have d = 0 in rejection
    breakflag = 1;
    if gnorm<=tolG
        exitflag = 'g';
%     elseif iter >= maxit
%         exitflag = 'i';
    elseif abs(fval) <= tolFun
        exitflag = 'f';
    elseif t > walltime
        exitflag = 't';
    else % abs(fdiff) <= relFun
        exitflag = 'r';
%     else % norm(d) <= tolX
%         exitflag = 'd';
    end
end

% % logistic
% % simple one
% tol = pars.tol;
% tolX = pars.tolX;
% maxit = pars.maxit;
% 
% breakflag = 0;
% exitflag = '';
% if ~( gnorm>tol && iter<maxit && norm(d)>tolX )
%     breakflag = 1;
%     if gnorm<=tol
%         exitflag = 'g';
%     elseif iter >= maxit
%         exitflag = 'i';
%     else % norm(d) <= tolX
%         exitflag = 'd';
%     end
% end



