function Plot_rhoSum()

close all;

name = {'L1','Tukey','t'};
n = 100;
xx = linspace(-3,3,n)';
% d = zeros(size(xx)); % distribution
f = zeros(size(xx));
g = zeros(size(xx));
h = zeros(size(xx));
w = zeros(size(xx));

addpath('/home/hjc/Public/Codes/MyLib/processing/plot'); % lineOption

mu = 0;     sigma = 0.8;     b = 0.5;       nu = 1.4;  

for i = 1:length(name);
    list{i}.name = name{i};
    list{i}.mu = mu;
    list{i}.sigma = sigma;
    list{i}.nu = nu;  
    list{i}.b = b;
    list{i}.c = 2;
    list{i}.k = 2;
end

% opt = lineOption( length(list), 'color','line' );
opt = {'r-.','g--','b-','k:'};

pars.p = 0.5; % Lp
pars.c = 2.5; % Fair, Cauchy, Welsch, Tukey, Andrew
pars.k = 2; % Huber
pars.nu= 1.4; % t

% hold on 
% axis equal
for i = 1:length(list)
    name = list{i}.name;
%     switch( name )
%         case 'L2'
%             mu = list{i}.mu;    sigma = list{i}.sigma;
%             d = pdf('normal',xx,mu,sigma);
%         case 'L1'
%             mu = list{i}.mu;    b = list{i}.b;
%             d = exp(-abs(xx-mu)./b)./(2*b);
%         case 't'
%             nu = list{i}.nu; 
%             d = pdf('t',xx,nu);
%             pars.nu= nu; % t
%     end

    for j = 1:length(f)
        [ f(j), g(j), Bx ] = rhoSum( xx(j), name, pars );
        h(j) = Bx(1);
        w(j) = g(j)/xx(j);
    end
    
%     figure(1);
%     % subplot(3,4,i)
%     hold on
%     fh=plot(xx,d,opt{i},'LineWidth',2);  
%     hold off
    
    figure(2);
    % subplot(3,4,i)
    hold on
    fh=plot(xx,f,opt{i},'LineWidth',2);  
    hold off
%     % add a Gaussian "background"
%     lim_x = get(gca,'XLim');  lim_y = get(gca,'YLim');
%     hold on; plot(xx,xx.^2./2,'r--'); hold off
%     axis([lim_x(1),lim_x(2),lim_y(1),lim_y(2)]);
%     saveas(fh,['f_',name{i}],'epsc')

    figure(3);
    % subplot(3,4,i)
    hold on
    plot(g,opt{i},'LineWidth',2)
    hold off
    
    figure(4);
    % subplot(3,4,i)
    hold on
    fh=plot(xx,h,opt{i},'LineWidth',2);  
    hold off
%     hold on; plot(xx,zeros(size(xx)),'r--'); hold off
%     saveas(fh,['H_',name{i}],'epsc')

    figure(5);
    % subplot(3,4,i)
    hold on
    plot(w,opt{i},'LineWidth',2)
    hold off
end
% hold off


for i = 1:5
    figure(i)
    lgd = {};
    for j = 1:length(list)
        % lgd{j} = num2str(list{j}.nu);
        lgd{j} = list{j}.name;
    end
%     legend(lgd)
end