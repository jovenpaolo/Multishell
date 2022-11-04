function l_max = l_conv( rad, ref, lam, fld )
%l_conv calculates l_max required for truncation
% -------------------------------------------------------------------------
%% INPUT
% -------------------------------------------------------------------------
% rad - outer radii of the sphere
% ref - refractive index of the host medium
% lam - vacuum wavelength (can be array with all wavelengths)
% fld - regime: 'far' or 'near' field
% -------------------------------------------------------------------------
%% OUTPUT
% -------------------------------------------------------------------------
% l_max - truncation l calculated via Wiscombe criterion [1] for far field
%           or with amended criterion for a near field [2]
% -------------------------------------------------------------------------
%% CALCULATING SIZE PARAMETER
% -------------------------------------------------------------------------
x = 2*pi*rad*ref/min(lam,[],'all');
% -------------------------------------------------------------------------
%% CONSTRUCTING OUTPUT
% -------------------------------------------------------------------------
switch fld
    case 'near'
        l_max = x + 11*x^(1/3) + 1;
    case 'far'
        if x <= 8
            l_max = x + 4*x^(1/3) + 1;
        elseif x < 4200
            l_max = x + 4.05*x^(1/3) + 2;
        else
            error('too large sphere. exiting');
        end
end