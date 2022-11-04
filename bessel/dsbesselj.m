function dj = dsbesselj( nu, z )
%DSBESSELJ calculates first derivative of 
%                       the spherical Bessel function of the first kind
% -------------------------------------------------------------------------
%% INPUT:
% -------------------------------------------------------------------------
% nu - the order of the Bessel function
% z  - argument
% -------------------------------------------------------------------------
%% OUTPUT
% -------------------------------------------------------------------------
% dj - first derivative of the spherical Bessel function of the first kind
% -------------------------------------------------------------------------
%% CALCULATING dj
% -------------------------------------------------------------------------
dj = - sbesselj( nu+1, z ) + ( nu./z ).*sbesselj( nu, z );
% -------------------------------------------------------------------------
end