function test_cutest( problems, pars, options )

DIR = problems.dir;
pids = problems.pid;
problemNames = problems.names;
rho = problems.rho; % problem model
noise_model = problems.noise_model;
noise_pars = problems.noise_pars;

methods = options.methods;
savedata = options.savedata;

nmthd = length(methods);

rec = cell(1, nmthd);

for kprob = 1 : length( problemNames )
    if any( pids == kprob )
        
        %====== Initialize problem ========================
        cd(strcat( DIR, problemNames{kprob} ));

        prob = cutest_setup();
        name = prob.name;
        n = prob.n;
        m = prob.m;
        x0 = prob.x;  
        % bl = prob.bl;
        % bu = prob.bu;
        
        r0 = res(x0);
        noise = 0.01*max( abs(r0) )*genNoise( size(r0,1), 1, noise_model, noise_pars );
%         noise = genNoise( size(r0,1), 1, noise_model, noise_pars );

        %==================================================
        fprintf('Computing problem %d: %s...\nn = %d, m = %d\n', kprob, name, n, m );
        
        for kmthd = 1:nmthd
            res_n = @(x) res_noisy( x, noise );
            [ x{kmthd}, info{kmthd} ] = algo_batching( methods{kmthd}, res_n, rho, x0, pars );
        end

        for kmthd = 1:nmthd
            rec{kmthd}.pid = kprob;
            rec{kmthd}.name = name;
            rec{kmthd}.n = n;
            rec{kmthd}.m = m;
            rec{kmthd}.x0 = x0;
            rec{kmthd}.r0 = r0;
            rec{kmthd}.mid = methods{kmthd};
            rec{kmthd}.x = x{kmthd};
            rec{kmthd}.info = info{kmthd};
            rec{kmthd}.pars = pars;
        end

        cutest_terminate();
       
    else
        fprintf('Problem skipped: %d...\n', kprob);
        rec = cell(1, nmthd);
    end
    
    % save rec
    % rec is empty if all(Nprob ~= kprob) (thus the stored data is more structured)
    if savedata == 1
        save([options.savedir,num2str(kprob)],'rec');
    end
end


function [ r, J, Jt ] = res( x )
if nargout == 1
    r = cutest_cons(x);
else
    [r,Jm] = cutest_cons(x);
    J = @(x) Jm*x;
    Jt = @(x) Jm.'*x;
end

function [ r, J, Jt ] = res_noisy( x, noise )
% add noise
[ r, J, Jt ] = res( x );
r = r + noise;
