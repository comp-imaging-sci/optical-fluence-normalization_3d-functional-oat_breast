function [r, theta, phi, rho] = SphericalCoord(dim_size, origin)
% [R, THETA, PHI, RHO] = SPHERICALCOORD(DIM_SIZE, ORIGIN)
% computes spherical coordinates of each voxel in a domain defined by its
% dimension size DIM_SIZE, based on the given origin coordinates ORIGIN. 
%
% Input:
%   DIM_SIZE:   dimension size. The dimension size of outputs R, THETA, 
%               PHI, and RHO is same as DIM_SIZE.
%   ORIGIN:     origin coordinates of the domain
%
% Output:
%   R:          radial distance from ORIGIN over a grid, a 3D array [voxel]
%   THETA:      polar angle over a grid, a 3D array [degree]
%   PHI:        azimuthal angle over a grid, a 3D array [degree]
%   RHO:        radial distance from z-axis over a grid, a 3D array [voxel]
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

[x, y, z] = meshgrid(1:dim_size(1), 1:dim_size(2), 1:dim_size(3));
x0 = x - origin(1); % x-coordinate [voxel]
y0 = y - origin(2); % y-coordinate [voxel]
z0 = z - origin(3); % z-coordinate [voxel]

r     = sqrt(x0.^2 + y0.^2 + z0.^2); % Distance from origin [voxel]
theta = acosd(z0./ r);               % Polar angle [degree]
phi   = angle(x0 + 1i.*y0).*180./pi; % Azimuthal angle [degree]
rho   = sqrt(x0.^2 + y0.^2);         % Distance from z-axis [voxel]
