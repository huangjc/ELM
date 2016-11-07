function [ f, g ] = phi( x, res_h, rho_h )
% [ r, J ] = res_h( x );
% [ f, g ] = rho_h( r );

if nargout == 1
    f = rho_h( res_h( x ) );
else
    [ r, J, Jt ] = res_h( x );
    [ rho, rho_g ] = rho_h( r );
    f = rho;
    g = Jt(rho_g);
end