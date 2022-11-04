function [ sc, ab, ex ] = crs_sec( rad, lam, nh, l, T )
%CRS_SEC calculates scattering, absorption and extinction of a general 
%   multilayered spherical particle
% -------------------------------------------------------------------------
%% INPUT
% -------------------------------------------------------------------------
% rad - outer radii for each layer of the sphere
% lam - vacuum wavelength
% nh  - refractive index of the host medium
% l   - numbers of terms in expansion (array or scalar)
% T   - transfer matrix
% -------------------------------------------------------------------------
%% OUTPUT
% -------------------------------------------------------------------------
% sc - scattering
% ab - absorption
% ex - extinction
%       key to structures:  cs - cross-section
%                           ef - efficiency
%                           e  - electric
%                           m  - magnetic
%                           em - e+m
%                           l  - for each l separately 
% -------------------------------------------------------------------------
%% ALLOCATING USEFUL QUANTITIES
% -------------------------------------------------------------------------
cf  = 0.5 * (lam/nh)^2 / pi;                                                % coefficient
gcs = pi * rad(end)^2;                                                      % geometric cross-section
% -------------------------------------------------------------------------
%% CALCULATING EXPANSION COEFFICIENTS
% -------------------------------------------------------------------------
a = -(squeeze(T.te(:,end,2,1))./squeeze(T.te(:,end,1,1))).';
b = -(squeeze(T.tm(:,end,2,1))./squeeze(T.tm(:,end,1,1))).';
% -------------------------------------------------------------------------
%% CONSTRUCTING OUTPUT
% -------------------------------------------------------------------------
sc.csel  = cf.*( 2*l + 1 ).*abs(a).^2;                                      % cross-sections for each mode l
sc.csml  = cf.*( 2*l + 1 ).*abs(b).^2;
sc.cseml = sc.csml + sc.csel;
ex.csel  = cf.*( 2*l + 1 ).*real(a);
ex.csml  = cf.*( 2*l + 1 ).*real(b);
ex.cseml = ex.csel + ex.csml;
ab.csel  = ex.csel - sc.csel;
ab.csml  = ex.csml - sc.csml;
ab.cseml = ab.csel + ab.csml;
% -------------------------------------------------------------------------
sc.cse  = sum( sc.csel );                                                   % total cross-sections
sc.csm  = sum( sc.csml );
sc.csem = sum( sc.cseml );
ex.cse  = sum( ex.csel );
ex.csm  = sum( ex.csml );
ex.csem = sum( ex.cseml );
ab.cse  = ex.cse - sc.cse;
ab.csm  = ex.csm - sc.csm;
ab.csem = ab.cse + ab.csm;
% -------------------------------------------------------------------------
sc.efel  = sc.csel / gcs;                                                   % efficiencies for each mode
sc.efml  = sc.csml / gcs;
sc.efeml = sc.cseml / gcs;
ex.efel  = ex.csel / gcs;
ex.efml  = ex.csml / gcs;
ex.efeml = ex.cseml / gcs;
ab.efel  = ab.csel / gcs;
ab.efml  = ab.csml / gcs;
ab.efeml = ab.cseml / gcs;
% -------------------------------------------------------------------------
sc.efe  = sc.cse / gcs;                                                     % total efficiencies
sc.efm  = sc.csm / gcs;
sc.efem = sc.csem / gcs;
ex.efe  = ex.cse / gcs;
ex.efm  = ex.csm / gcs;
ex.efem = ex.csem / gcs;
ab.efe  = ab.cse / gcs;
ab.efm  = ab.csm / gcs;
ab.efem = ab.csem / gcs;
% -------------------------------------------------------------------------
end