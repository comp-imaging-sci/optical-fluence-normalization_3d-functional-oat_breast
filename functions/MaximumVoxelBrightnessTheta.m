function [theta_range, maxval, maxidx_mask] = MaximumVoxelBrightnessTheta(img, theta)
% [THETA_RANGE, MAXVAL, MAXINX_MASK] = MAXIMUMVOXELBRIGHTNESSTHETA(IMG,
% THETA) extracts maximum voxel brightness MAXVAL at each polar angle THETA 
% in the given three-dimensional (3D) image IMG and its index mask 
% MAXIDX_MASK.
%
% Input:
%   IMG:            3D image, a 3D array
%   THETA:          polar angle over a grid, a 3D array [degree]
%
% Output:
%   THETA_RANGE:    polar angle range for extraction of maximum voxel 
%                   brightness [degree]
%   MAXVAL:         maximum voxel brightness at each polar angle 
%   MAXIDX_MASK:    index mask of MAXVAL, a 3D array
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

inc = 1; % Discrete polar angle step [degree]

theta_range = (floor(min(theta(:))):inc:ceil(max(theta(:))))';
maxval = zeros(length(theta_range), 1);
maxidx_mask = zeros(size(img));

for theta_i = 1:length(theta_range)
    % Set region of interest
    img_at_theta = img;
    img_at_theta(theta >  theta_range(theta_i) + inc/2) = 0;
    img_at_theta(theta <= theta_range(theta_i) - inc/2) = 0;
    
    maxval(theta_i) = max(img_at_theta(:));
    maxidx_mask(img_at_theta == maxval(theta_i)) = 1;
end
