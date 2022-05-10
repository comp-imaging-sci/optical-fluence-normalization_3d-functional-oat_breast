function mask = BreastMask(z, rho, f_ellp)
% MASK = BREASTMASK(Z, RHO, F_ELLP) creates a mask MASK by assigning the 
% value "1" to voxels inside the breast boundary, that is defined by a fit 
% object F_ELLP (elliptical curve-fitting result), and "0" outside.
%
% Input:
%   Z:          z-coordinates over a grid, a 3D array [voxel]
%   RHO:        polar angle over a grid, a 3D array [voxel]
%   F_ELLP:     fit object, i.e., elliptical curve-fitting result
%
% Output:
%   MASK:       breast mask, a 3D array
%
% -------------------------------------------------------------------------
% Author:   Seonyeong Park (sp33@illinois.edu)
% Date:     Jan. 24, 2022
%
% This function is part of optical-fluence-normalization_3d-functional-oat_
% breast.
%
% Copyright (C) 2022 Seonyeong Park
% License:  GNU General Public License version 3, Please see 'LICENSE' for
%           details.
%

mask = ones(size(rho));
bndr_rho = reshape(sqrt(f_ellp(z)), size(z));
bndr_rho = real(bndr_rho); % Exclude imaginary part
mask(rho > bndr_rho) = 0;
mask = logical(mask);
