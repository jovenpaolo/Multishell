% this script calculates scattering, absorption and extinction efficiencies
% of Au@Fe core-shell sphere and shows the contribution of multipoles
% -------------------------------------------------------------------------
function [sct,abt,ext] = PeakQ( radFe, radAu, lam)

l   = 1:3;                                                                  % up to 3-rd multipole
% -------------------------------------------------------------------------
% setting up refractive indices
load( 'Ag_JC.mat' );
load( 'Fe304_Q.mat');
% loading refractive indices
nfe = complex( interp1( Fe3O4Q(:,1), Fe3O4Q(:,2), lam, 'spline' ),...       % interpolating tabulated data for Au
                 interp1( Fe3O4Q(:,1), Fe3O4Q(:,3), lam, 'spline') );          
nau   = complex( interp1( Ag_JC(:,1), Ag_JC(:,2), lam, 'spline' ),...       % interpolating tabulated data for Au
                 interp1( Ag_JC(:,1), Ag_JC(:,3), lam, 'spline') );          
nh    = 1.33;                                                               % water host
% -------------------------------------------------------------------------
% setting up particle
rad = [radFe,radFe+radAu]*1e-9;                                             % radii of shells, meters
mu  = ones(1,3);                                                            % permeabilities, including the host medium     
% -------------------------------------------------------------------------
% setting up output
abel = zeros( numel( lam ), numel( l ) ); 
abml = zeros( numel( lam ), numel( l ) );
abt  = zeros( numel( lam ), 1 );
exel = zeros( numel( lam ), numel( l ) ); 
exml = zeros( numel( lam ), numel( l ) );
ext  = zeros( numel( lam ), 1 );
scel = zeros( numel( lam ), numel( l ) ); 
scml = zeros( numel( lam ), numel( l ) );
sct  = zeros( numel( lam ), 1 );
% -------------------------------------------------------------------------
%% CALCULATING SCATTERING, ABSORPTION AND EXTINCTION
% -------------------------------------------------------------------------
parfor (i = 1 : numel(lam),6)
    nb = [nfe(i),nau(i), nh];                                               
    T = t_mat( rad, nb, mu, lam(i), l );                                    % transfer matrices
    [ sc, ab, ex ] = crs_sec( rad, lam(i), nh, l, T );                      % getting sc, ab, ex
    abt(i)    = ab.efem;  sct(i)    = sc.efem;  ext(i)    = ex.efem;
end
lam1 = lam*1e9;
end
