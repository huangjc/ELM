function plot_perfprof( record, fieldname )

Field = getFieldCell( record, [], fieldname );

PASS = cvg_test( record );
    
% filtered failed problems
q = find( any(PASS,2) ~= 0 );
% q = find( all(PASS,2) ~= 0 );
PASS = PASS(q,:);
Field = Field(q,:);
fprintf(['#(prob)\t #(any passed)', repmat('\t #(pass)', 1, size(PASS,2)), '\n', ...
        '%d\t %d',repmat('\t %d', 1, size(PASS,2)),'\n'], ...
        size(PASS,1), size(q,1), sum(PASS)');

%%% performance profile
Field( PASS==0 ) = NaN;

figure;
perf_profile(Field,1);
% xlabel('t/t_{min}')
% xlabel('iter/iter_{min}')
% xlabel('nfeval/nfeval_{min}')
xl = get(gca,'XTickLabel');
set(gca,'XTickLabel',xl,'FontSize',12,'FontWeight','Bold')
% ylabel('percent solved')
yl = get(gca,'YTickLabel');
set(gca,'YTickLabel',yl,'FontSize',12,'FontWeight','Bold')

% title(fieldname);

methods = getFieldCell( record, 0, 'mid' );
methods = methods(1,:);
% legend('ELM','ELM\_WLS','ELM\_FIM','ELM\_MIX','LBFGS','location','southeast')
legend(methods,'Location','best'); 
saveas(gcf, fieldname, 'epsc');

% printTable(record(q,:));
