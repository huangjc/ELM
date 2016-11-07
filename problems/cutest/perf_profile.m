function hl = perf_profile(T, logplot)
% failed problems should be set as NaN in T
% Input:    T is np x ns, each column forms a plot in the performance
%                         profile

if (nargin < 2)
    logplot = 0; 
end

[np,ns] = size(T);
t_min = min( T, [], 2 );
if min(abs(t_min)) == 0
    error('Error: Divided over 0.');
end
r = T./repmat( t_min, 1, ns );
max_ratio = max(max(r));
    
r( isnan(r) ) = 2*max_ratio;
r = sort(r, 1);

% % print percentage sketch
% nmin(1,:) = sum( r<= 1.01 );
% nmin(2,:) = sum( r<= 4 );
% nmin(3,:) = sum( r<= 16 );
% nmin(4,:) = sum( r<= max_ratio );
% 
% fprintf([ % repmat('\t #(min)', 1, size(PASS,2)), '\n', ...
%         repmat('\t %4.3f', 1, ns),'\n'], ...
%         (nmin/np)' );
    
colors  = ['b' 'r' 'k' 'm' 'c' 'g' 'y'];   lines   = {'-' '-.' '--'};
markers = [ 's' 'o' '^' 'v' 'p' '<' 'x' 'h' '+' 'd' '*' '<' ];

% Plot stair graphs with markers.
hl = zeros(ns,1);
opt = {'b--','g-.','k-.','r-','m-','c--'};
for s = 1:ns
    [xs,ys] = stairs(r(:,s),(1:np)/np);

    % Only plot one marker at the intercept
    if (xs(1)==1)
        vv = find(xs==1,1,'last');
        xs = xs(vv:end);   ys = ys(vv:end);
    end

    sl = mod(s-1,3) + 1; sc = mod(s-1,7) + 1; sm = mod(s-1,12) + 1;
    option1 = opt{s};
    if (logplot)
        hl(s) = semilogx(xs,ys,option1,'linewidth',2);
    else
        hl(s) = plot(xs,ys,option1);
    end
    hold on;
end
hold off

% Axis properties are set so that failures are not shown, but with the
% max_ratio data points shown. This highlights the "flatline" effect.
if (logplot) 
  axis([1 1.1*max_ratio 0 1]);
  twop = max(floor(log2(1.1*max_ratio)),1);
  set(gca,'XTick',2.^[0:twop])
else
  axis([1 1.1*max_ratio 0 1]);
end