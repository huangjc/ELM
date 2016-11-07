function results_logreg()

addpath /home/hjc/Public/Codes/MyLib/processing; % lineOption

data_names = {'a9a','real-sim','news20','rcv1'};
dataname = data_names{1};
datapath = ['/home/hjc/Projects/Extended_LM/codes/temp/logreg/',dataname,'/'];

problems = 1; % rec id
methods = {'ELM','LBFGS'};

cd(datapath); % direct to PATH to save figures

np = length(problems);
nm = length(methods);
record = load_data( datapath, problems, 1:nm ); % np x nm cell
nm = length(methods);

Res = getFieldCell( record, [], 'info.res' );
minVal = min( Res, [], 2);
Res = Res - repmat(minVal,1,nm);
fprintf( [ repmat('%6.4e \t', 1, nm ), '\n'], Res' );

%%% performance profile
% convergence test
tau = 1e-3;
PASS = cvg_test( Res, minVal, tau );
       
[ T ] = getFieldCell( record, [], 'info.time' );


%% plot some figures
i_sample = 1;
rec = record(i_sample,:);

t = getFieldCell( record, [], 'info.time' )
t = max(t);
xl = [0:5]'*ceil(t/5);
% xl = [0:10:60]';

% plot function value history
figure(1);
% opt = lineOption(length(method),'color','line');
opt = {'b-','r--','g^-','k*-','cd-'};
for i = 1:length(methods)
    plot( rec{i}.info.t_hist, log10(rec{i}.info.f_hist - minVal(i_sample)), opt{i}, 'linewidth', 2 );
    hold on;
end
hold off;

% yl = get(gca,'YTickLabel');
% set(gca,'YTickLabel',strcat('10^',yl));
ylim([-10,6])
yl =[-10,-5,0,5]';
yl = num2str(yl,'%-3d');
set(gca,'YTickLabel',[]);
text(repmat(-0.25, size(yl,1),1), str2num(yl), strcat('10^{',yl,'}'),...
    'HorizontalAlignment','center','FontSize',12,'FontWeight','Bold')

set(gca,'XTickLabel',[]);
% xl = get(gca,'XTickLabel');
% xl = [0:5]';
set(gca,'XTick',xl);
xlim([0,max(xl)])
text( xl, repmat(-10.5, size(xl,1),1), num2str(xl),...
    'HorizontalAlignment','center','FontSize',12,'FontWeight','Bold')

legend(strrep(methods,'_',' '))
saveas(gcf,[dataname,'_f'],'epsc');

% plot gnorm history
figure(2);
% opt = lineOption(length(method),'color','line');
% opt = {'b^','r*'};
for i = 1:length(methods)
    plot( rec{i}.info.t_hist, log10(rec{i}.info.gnorm_hist), opt{i}, 'linewidth', 2 );
    hold on;
end
hold off;
% yl=get(gca,'YTickLabel');
% set(gca,'YTickLabel',strcat('10^',yl));
ylim([-5,5])
yl =[-5,0,5]';
yl = num2str(yl,'%-3d');
set(gca,'YTickLabel',[]);
text(repmat(-0.25, size(yl,1),1), str2num(yl), strcat('10^{',yl,'}'),...
    'HorizontalAlignment','center','FontSize',12,'FontWeight','Bold')

set(gca,'XTickLabel',[]);
% xl = get(gca,'XTickLabel');
% xl = [0:5]';
set(gca,'XTick',xl);
xlim([0,max(xl)])
text( xl, repmat(-5.3, size(xl,1),1), num2str(xl),...
    'HorizontalAlignment','center','FontSize',12,'FontWeight','Bold')

legend(strrep(methods,'_',' '))
saveas(gcf,[dataname,'_gnorm'],'epsc');


