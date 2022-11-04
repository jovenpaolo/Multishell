function xi = ricbesh( nu, k, z )
%RICBESH calculates the Riccati-Bessel function of the third kind
% -------------------------------------------------------------------------
%% INPUT:
% -------------------------------------------------------------------------
% nu - the order of Ricatti-Bessel function
% k  - 1 or 2, default = 1
% z  - argument
% -------------------------------------------------------------------------
%% OUTPUT
% -------------------------------------------------------------------------
% xi - Ricatti-Bessel function of the third kind
% -------------------------------------------------------------------------
%% COPYRIGHT
% -------------------------------------------------------------------------
% Copyright 2020 Ilia Rasskazov, University of Rochester
% -------------------------------------------------------------------------
% Author:        Ilia Rasskazov, irasskaz@ur.rochester.edu
% -------------------------------------------------------------------------
% Organization:  The Institute of Optics, University of Rochester
%                http://www.hajim.rochester.edu/optics/
% -------------------------------------------------------------------------
%% CHECKING INPUT
% -------------------------------------------------------------------------
if nargin == 2
    z = k;
    k = 1;
end
% -------------------------------------------------------------------------
%% CALCULATING xi
% -------------------------------------------------------------------------
xi = sqrt( pi.*z/2 ).*besselh( nu+0.5, k, z );
% -------------------------------------------------------------------------
end