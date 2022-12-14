function F = F_lossless( A, B, x, sm, jl, hl, l_max )
%F_LOSSLESS calculates function F, which enters final expressions
%           for the total energy within lossless shells
% -------------------------------------------------------------------------
%% INPUT
% -------------------------------------------------------------------------
% A, B  - expansion coefficients
% x     - size parameter
% sm    - indices of the lossless shells
% jl    - spherical Bessel functions of the 1st kind
% hl    - spherical Hankel functions of the 1st kind
% l_max - truncation for the l
% -------------------------------------------------------------------------
%% OUTPUT
% -------------------------------------------------------------------------
% F - structure with functions F.xyyyz with the following notation:
%   x: 
%       - 'e' for electric
%       - 'm' for magnetic
%   yyy:
%       - 'min' for l-1 order
%       - 'sam' for l   order
%       - 'pls' for l+1 order
%   z:
%       - 'u' for upper boundary of definitive integral
%       - 'b' for lower boundary of definitive integral
% -------------------------------------------------------------------------
%% CALCULATING USEFUL VALUES
% -------------------------------------------------------------------------
xu     = repmat( x.u(sm), l_max, 1 );
xb     = repmat( x.b(sm), l_max, 1 );
l_2    = repmat( 2*(0:l_max+1)' + 1, 1, numel( x.u(sm) ) );
l_2min = l_2(1:end-2,:);
l_2sam = l_2(2:end-1,:);
l_2pls = l_2(3:end,  :);
% -------------------------------------------------------------------------
%% CONSTRUCTING OUTPUT
% -------------------------------------------------------------------------
F.eminu = xu.*( abs( A.e(:,sm).*jl.u(1:end-3,sm) + B.e(:,sm).*hl.u(1:end-3,sm) ).^2 + ...
                abs( A.e(:,sm).*jl.u(2:end-2,sm) + B.e(:,sm).*hl.u(2:end-2,sm) ).^2 ) - ...
 l_2min.*real( conj( A.e(:,sm).*jl.u(1:end-3,sm) + B.e(:,sm).*hl.u(1:end-3,sm) ).* ...
                   ( A.e(:,sm).*jl.u(2:end-2,sm) + B.e(:,sm).*hl.u(2:end-2,sm) ) );
F.mminu = xu.*( abs( A.m(:,sm).*jl.u(1:end-3,sm) + B.m(:,sm).*hl.u(1:end-3,sm) ).^2 + ...
                abs( A.m(:,sm).*jl.u(2:end-2,sm) + B.m(:,sm).*hl.u(2:end-2,sm) ).^2 ) - ...
 l_2min.*real( conj( A.m(:,sm).*jl.u(1:end-3,sm) + B.m(:,sm).*hl.u(1:end-3,sm) ).* ...
                   ( A.m(:,sm).*jl.u(2:end-2,sm) + B.m(:,sm).*hl.u(2:end-2,sm) ) );
% -------------------------------------------------------------------------
F.esamu = xu.*( abs( A.e(:,sm).*jl.u(2:end-2,sm) + B.e(:,sm).*hl.u(2:end-2,sm) ).^2 + ...
                abs( A.e(:,sm).*jl.u(3:end-1,sm) + B.e(:,sm).*hl.u(3:end-1,sm) ).^2 ) - ...
 l_2sam.*real( conj( A.e(:,sm).*jl.u(2:end-2,sm) + B.e(:,sm).*hl.u(2:end-2,sm) ).* ...
                   ( A.e(:,sm).*jl.u(3:end-1,sm) + B.e(:,sm).*hl.u(3:end-1,sm) ) );
F.msamu = xu.*( abs( A.m(:,sm).*jl.u(2:end-2,sm) + B.m(:,sm).*hl.u(2:end-2,sm) ).^2 + ...
                abs( A.m(:,sm).*jl.u(3:end-1,sm) + B.m(:,sm).*hl.u(3:end-1,sm) ).^2 ) - ...
 l_2sam.*real( conj( A.m(:,sm).*jl.u(2:end-2,sm) + B.m(:,sm).*hl.u(2:end-2,sm) ).* ...
                   ( A.m(:,sm).*jl.u(3:end-1,sm) + B.m(:,sm).*hl.u(3:end-1,sm) ) );
% -------------------------------------------------------------------------
F.eplsu = xu.*( abs( A.e(:,sm).*jl.u(3:end-1,sm) + B.e(:,sm).*hl.u(3:end-1,sm) ).^2 + ...
                abs( A.e(:,sm).*jl.u(4:end  ,sm) + B.e(:,sm).*hl.u(4:end  ,sm) ).^2 ) - ...
 l_2pls.*real( conj( A.e(:,sm).*jl.u(3:end-1,sm) + B.e(:,sm).*hl.u(3:end-1,sm) ).* ...
                   ( A.e(:,sm).*jl.u(4:end  ,sm) + B.e(:,sm).*hl.u(4:end  ,sm) ) );
F.mplsu = xu.*( abs( A.m(:,sm).*jl.u(3:end-1,sm) + B.m(:,sm).*hl.u(3:end-1,sm) ).^2 + ...
                abs( A.m(:,sm).*jl.u(4:end  ,sm) + B.m(:,sm).*hl.u(4:end  ,sm) ).^2 ) - ...
 l_2pls.*real( conj( A.m(:,sm).*jl.u(3:end-1,sm) + B.m(:,sm).*hl.u(3:end-1,sm) ).* ...
                   ( A.m(:,sm).*jl.u(4:end  ,sm) + B.m(:,sm).*hl.u(4:end  ,sm) ) );                            
% -------------------------------------------------------------------------
F.eminb = xb.*( abs( A.e(:,sm).*jl.b(1:end-3,sm) + B.e(:,sm).*hl.b(1:end-3,sm) ).^2 + ...
                abs( A.e(:,sm).*jl.b(2:end-2,sm) + B.e(:,sm).*hl.b(2:end-2,sm) ).^2 ) - ...
 l_2min.*real( conj( A.e(:,sm).*jl.b(1:end-3,sm) + B.e(:,sm).*hl.b(1:end-3,sm) ).* ...
                   ( A.e(:,sm).*jl.b(2:end-2,sm) + B.e(:,sm).*hl.b(2:end-2,sm) ) );
F.mminb = xb.*( abs( A.m(:,sm).*jl.b(1:end-3,sm) + B.m(:,sm).*hl.b(1:end-3,sm) ).^2 + ...
                abs( A.m(:,sm).*jl.b(2:end-2,sm) + B.m(:,sm).*hl.b(2:end-2,sm) ).^2 ) - ...
 l_2min.*real( conj( A.m(:,sm).*jl.b(1:end-3,sm) + B.m(:,sm).*hl.b(1:end-3,sm) ).* ...
                   ( A.m(:,sm).*jl.b(2:end-2,sm) + B.m(:,sm).*hl.b(2:end-2,sm) ) );
% -------------------------------------------------------------------------
F.esamb = xb.*( abs( A.e(:,sm).*jl.b(2:end-2,sm) + B.e(:,sm).*hl.b(2:end-2,sm) ).^2 + ...
                abs( A.e(:,sm).*jl.b(3:end-1,sm) + B.e(:,sm).*hl.b(3:end-1,sm) ).^2 ) - ...
 l_2sam.*real( conj( A.e(:,sm).*jl.b(2:end-2,sm) + B.e(:,sm).*hl.b(2:end-2,sm) ).* ...
                   ( A.e(:,sm).*jl.b(3:end-1,sm) + B.e(:,sm).*hl.b(3:end-1,sm) ) );
F.msamb = xb.*( abs( A.m(:,sm).*jl.b(2:end-2,sm) + B.m(:,sm).*hl.b(2:end-2,sm) ).^2 + ...
                abs( A.m(:,sm).*jl.b(3:end-1,sm) + B.m(:,sm).*hl.b(3:end-1,sm) ).^2 ) - ...
 l_2sam.*real( conj( A.m(:,sm).*jl.b(2:end-2,sm) + B.m(:,sm).*hl.b(2:end-2,sm) ).* ...
                   ( A.m(:,sm).*jl.b(3:end-1,sm) + B.m(:,sm).*hl.b(3:end-1,sm) ) );
% -------------------------------------------------------------------------
F.eplsb = xb.*( abs( A.e(:,sm).*jl.b(3:end-1,sm) + B.e(:,sm).*hl.b(3:end-1,sm) ).^2 + ...
                abs( A.e(:,sm).*jl.b(4:end  ,sm) + B.e(:,sm).*hl.b(4:end  ,sm) ).^2 ) - ...
 l_2pls.*real( conj( A.e(:,sm).*jl.b(3:end-1,sm) + B.e(:,sm).*hl.b(3:end-1,sm) ).* ...
                   ( A.e(:,sm).*jl.b(4:end  ,sm) + B.e(:,sm).*hl.b(4:end  ,sm) ) );
F.mplsb = xb.*( abs( A.m(:,sm).*jl.b(3:end-1,sm) + B.m(:,sm).*hl.b(3:end-1,sm) ).^2 + ...
                abs( A.m(:,sm).*jl.b(4:end  ,sm) + B.m(:,sm).*hl.b(4:end  ,sm) ).^2 ) - ...
 l_2pls.*real( conj( A.m(:,sm).*jl.b(3:end-1,sm) + B.m(:,sm).*hl.b(3:end-1,sm) ).* ...
                   ( A.m(:,sm).*jl.b(4:end  ,sm) + B.m(:,sm).*hl.b(4:end  ,sm) ) );                            
% ------------------------------------------------------------------------- 
end