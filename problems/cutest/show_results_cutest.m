function show_results( PATH, pid, methods )

cd(PATH); % direct to PATH to save figures

np = length(pid);
nm = length(methods);

mid = 1:nm;
record = load_data( PATH, pid, mid ); % np x nm cell

% plot_perfprof( record, 'info.time' );
% plot_perfprof( record, 'info.iter' );
plot_perfprof( record, 'info.nF' );
