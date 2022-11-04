function dS = dricbesj( nu, z )
%DRICBESJ calculates first derivative of 
%                      the Riccati-Bessel function of the first kind
% -------------------------------------------------------------------------
%% INPUT:
% -------------------------------------------------------------------------
% nu - the order of the Bessel function
% z  - argument
% -------------------------------------------------------------------------
%% OUTPUT
% -------------------------------------------------------------------------
% dS - first derivative of the Riccati-Bessel function of the first kind
% -------------------------------------------------------------------------
%% CALCULATING dS
% -------------------------------------------------------------------------
dS = 0.5*sqrt( pi/2./z ).* besselj( nu+0.5, z ) ...
       + sqrt( pi*z/2 ).*dbesselj( nu+0.5, z );
% -------------------------------------------------------------------------
end