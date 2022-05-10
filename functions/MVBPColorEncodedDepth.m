function [mvbpx, mvbpy, mvbpz] = MVBPColorEncodedDepth(img, voxel_size, depth, max_depth)
% [MVBPX, MVBPY, MVBPZ] = MVBPCOLORENCODEDDEPTH(IMG, VOXEL_SIZE, DEPTH,
% MAX_DEPTH) plots maximum voxel brightness projections of a
% three-dimensional(3D) image IMG along x, y, and z-axis that are 
% color-encoded by depth DEPTH.
% 
% Input:
%   IMG:        3D image
%   VOXEL_SIZE: voxel size
%   DEPTH:      depth map [voxel]
%   MAX_DEPTH:  maximum depth to visualize [voxel]
%
% Output:
%   MVBPX:      maximum voxel brightness projection along x-axis that is
%               color-encoded by depth
%   MVBPY:      maximum voxel brightness projection along y-axis that is
%               color-encoded by depth
%   MVBPZ:      maximum voxel brightness projection along z-axis that is
%               color-encoded by depth
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

% Get dimension size [voxel]
[Nx, Ny, Nz] = size(img);

% Normalize brightness
img = img./max(img(:));

% Colormap
cmap = jet(round(max_depth*1.5));
cmap = cmap(end - max_depth:end, :);

% Depth colormap
depth(depth > max_depth) = max_depth;
d_cmap = cmap(round(depth) + 1,:);
d_cmap = reshape(d_cmap, [Nx, Ny, Nz, 3]);

% 4D image color-encoded by depth
imgc = cat(4, img, img, img);
imgc = imgc.*d_cmap;

% Maximum voxel brightness projection along x-axis, color-encoded by depth
[~, mvbx_idx] = max(img, [], 1);
mvbx_idx = squeeze(mvbx_idx);
mvbpx = zeros(Ny, Nz, 3);
for i = 1:Ny
    for j = 1:Nz
        mvbpx(i, j, :) = imgc(mvbx_idx(i, j), i, j, :);
    end
end
figure; imshow(mvbpx);
    
% Maximum voxel brightness projection along y-axis, color-encoded by depth
[~, mvby_idx] = max(img, [], 2);
mvby_idx = squeeze(mvby_idx);
mvbpy = zeros(Nx, Nz, 3);
for i = 1:Nx
    for j = 1:Nz
        mvbpy(i, j, :) = imgc(i, mvby_idx(i, j), j, :);
    end
end
figure; imshow(mvbpy);

% Maximum voxel brightness projection along z-axis, color-encoded by depth
[~, mvbz_idx] = max(img, [], 3);
mvbpz = zeros(Nx, Ny, 3);
for i = 1:Nx
    for j = 1:Ny
        mvbpz(i, j, :) = imgc(i, j, mvbz_idx(i, j), :);
    end
end
figure; imshow(mvbpz);

% Color bar
NShade = 30;
cbar = zeros(NShade, max_depth + 1, 3);
for i = 1:25
    cbar(i, :, :) = cmap.*(i/25);
end
for i = 1:5
    cbar(25 + i, :, :) = cmap.*(1 + i*2/25);
end
cbar = rot90(cbar);

% Plot color bar
figure; imshow(cbar); axis on
xlabel('Normalized brightness', 'Interpreter', 'latex');
ylabel('Depth from breast surface [mm]', 'Interpreter', 'latex');
xticks([1, NShade]); xticklabels([0, 1]);
yticks(1:max_depth/5:max_depth + 1);
yticklabels({ ...
    num2str(max_depth*voxel_size), ...
    num2str(max_depth*voxel_size*4/5), ...
    num2str(max_depth*voxel_size*3/5), ...
    num2str(max_depth*voxel_size*2/5), ...
    num2str(max_depth*voxel_size/5), ...
    '0'});
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 14, 'FontName', 'Times');
