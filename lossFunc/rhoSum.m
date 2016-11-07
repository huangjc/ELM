function [ f, g, Bx ] = rhoSum( r, name, pars )
% Bx    function handle for Hessian approximation, convenient for large B,
%       especially when it is expensive to be stored.
% d     diagonal of B

if nargin < 3
    pars.p = 0.5; % Lp
    pars.c = 2; % Fair, Cauchy, Welsch, Tukey, Andrew
    pars.k = 2; % Huber
    pars.nu= 2; % t
end

switch( name )
    case {'gauss','L2'}
        f = r'*r/2;
        if nargout > 1
            g = r;
            if nargout > 2
                Bx = @(x) x;
            end
        end

    case {'laplace','L1'}
        f = sum( abs(r) );
        if nargout > 1
            g = sign(r);
            if nargout > 2
                Bx = @(x) 0*x;
            end
        end
        
    case 'L1-L2'
        root = sqrt(1+r.^2/2);
        f = sum( 2*(root-1) );
        if nargout > 1
            g = r./root;
            if nargout > 2
                Bx = @(x) diag( sparse( 1./(root.^3) ) )*x;
            end
        end
        
    case 'Lp'
        p = pars.p;
        f = sum( abs(r).^p )/p;
        if nargout > 1
            g = sign(r).*abs(r).^(p-1);
            if nargout > 2
                Bx = @(x) diag( sparse( (p-1)*abs(r).^(p-2) ))*x;
            end
        end
        
    case 'Fair'
        c = pars.c;
        roc = abs(r)/c;
        f = c^2*sum( roc - log(1+roc) );
        if nargout > 1
            g = r./(1+roc);
            if nargout > 2
                Bx = @(x) diag( sparse( 1./( (1+roc).^2 ) ))*x;
            end
        end
        
    case 'Huber'
        k = pars.k;
        n = length(r);
        f = zeros(n,1);     g = zeros(n,1);     h = zeros(n,1);
        for i = 1:n
            if abs(r(i)) <= k 
                f(i) = r(i)^2/2;
                if nargout > 1
                    g(i) = r(i);
                    if nargout > 2
                        h(i) = 1;
                    end
                end
            else
                f(i) = k*(abs(r(i))-k/2);
                if nargout > 1
                    g(i) = k*sign(r(i));
                    if nargout > 2
                        h(i) = 0;
                    end
                end
            end
        end
        f = sum(f);
        Bx = @(x) diag( sparse( h ))*x;
    case 'Cauchy'
        c = pars.c;
        roc2 =  (r./c).^2;
        f = c^2/2*sum( log(1 + roc2) );
        if nargout > 1
            g = r./(1+roc2);
            if nargout > 2
                Bx = @(x) diag( sparse( (1-roc2)./( (1+roc2).^2 ) ) )*x;
            end
        end
    case 'Geman-McClure'
        r2 = r.^2;
        f = sum( r2/2./ (1+r2) );
        if nargout > 1
            g = r./( (1+r2).^2 );
            if nargout > 2
                Bx = @(x) diag( sparse( (1-3*r2)./( (1+r2).^3 ) ) )*x;
            end
        end
    case 'Welsch'
        c = pars.c;
        roc2 = (r./c).^2;
        f = c^2/2*sum( 1-exp(-roc2) );
        if nargout > 1
            g = r.*exp( -roc2 );
            if nargout > 2
                Bx = @(x) diag( sparse( (1-2*roc2).*exp( -roc2 ) ))*x;
            end
        end
    case 'Tukey'
        c = pars.c;
        n = length(r);
        f = zeros(n,1);     g = zeros(n,1);     h = zeros(n,1);
        for i = 1:n
            x = r(i);
            if abs(r(i)) <= c 
                f(i) = c^2/6*( 1-(1-(x/c)^2)^3 );
                if nargout > 1
                    g(i) = x*(1-(x/c)^2)^2;
                    if nargout > 2
                        h(i) = (1-(x/c)^2)*((1-5*(x/c)^2));
                    end
                end
            else
                f(i) = c^2/6;
                if nargout > 1
                    g(i) = 0;
                    if nargout > 2
                        h(i) = 0;
                    end
                end
            end
        end
        f = sum(f);
        Bx = @(x) diag( sparse( h ))*x;
    case 'Andrew'
        c = pars.c;
        n = length(r);
        f = zeros(n,1);     g = zeros(n,1);     h = zeros(n,1);
        for i = 1:n
            x = r(i);
            if abs( x ) <= pi*c 
                f(i) = c*(1-cos(x/c));
                if nargout > 1
                    g(i) = sin(x/c);
                    if nargout > 2
                        h(i) = 1/c*cos(x/c);
                    end
                end
            else
                f(i) = 2*c;
                if nargout > 1
                    g(i) = 0;
                    if nargout > 2
                        h(i) = 0;
                    end
                end
            end
        end
        f = sum(f);
        Bx = @(x) diag( sparse( h ))*x;
    case 't'
        nu = pars.nu;
        r2 = r.^2;
        f = sum( log( 1 + r2/nu ) );
        if nargout > 1
            g = 2*r./(nu+r2);
            if nargout > 2
                Bx = @(v) diag( sparse( 2*( nu-r2 )./(( nu+r2 ).^2) )) *v;
            end
        end
end