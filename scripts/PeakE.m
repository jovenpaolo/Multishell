function [ Emax] = PeakE( radFe, radAu, lam)
%% INPUT
% -------------------------------------------------------------------------
% setting up incident illumination
% -------------------------------------------------------------------------
% setting up particles
rad = [radFe,radFe+radAu]*1e-9;                                            
% -------------------------------------------------------------------------
% setting up refractive indices
load( 'Au_JC.mat' );                                                       % loading Au Johnson and Christy
load( 'Fe304_Q.mat');
% loading refractive indices
%nfe = 1.45;                                                               % silica shell
nfe = complex( interp1( Fe3O4Q(:,1), Fe3O4Q(:,2), lam, 'spline' ),...      % interpolating tabulated data for Au
                 interp1( Fe3O4Q(:,1), Fe3O4Q(:,3), lam, 'spline') ); 
nau = complex( interp1( Au_JC(:,1), Au_JC(:,2), lam, 'spline' ),...        % interpolating data for Au
               interp1( Au_JC(:,1), Au_JC(:,3), lam, 'spline') ).';   
nh = 1.33;                                                                 % host medium
% -------------------------------------------------------------------------
% setting up mesh for near-field
xl = linspace(-70,70,700)*1e-9;                                            % 140x140nm^2 area in XOY plane
yl = linspace(-70,70,700)*1e-9;
zl = 0;
[X,Y,Z] = meshgrid(xl,yl,zl);
% -------------------------------------------------------------------------
%% CALCULATING NEAR-FIELDS AND PLOTTING RESULTS
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% calculating    
    nb = [repmat([nau,nfe],1,1),nh];                                       % refractive indices
    mu = ones(1,2*1+1);                                                    % permeabilities        
    l = 1 : l_conv( rad(end), nh, lam, 'near' );                           % defining l required for convergence
    T = t_mat( rad, nb, mu, lam, l );                                      % transfer-matrices
    [~,~,I] = near_fld(rad, nb, mu, lam, l, T, X,Y,Z);                     % returning field intensities
% -------------------------------------------------------------------------
Emax= max(max(I.E));
end