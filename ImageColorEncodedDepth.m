% -------------------------------------------------------------------------
%                       Images Color-Encoded by Depth
% -------------------------------------------------------------------------
% Author:   Seonyeong Park (sp33@illinois.edu)
% Date:     Jan. 24, 2022
%
% This project is to plot and save maximum voxel brightness projections of
% three-dimensional (3D) images along x, y, and z-axis. Before running this
% code, the detected blood vessel masks should be saved by running
% 'ArteryVeinDetectionClassification.m'. The method of the software is
% detailed in [Park2022].
%
% Reference:
%   [Park2022] Seonyeong Park, Frank J. Brooks, Umberto Villa, Richard Su,
%           Mark A. Anastasio, Alexander A. Oraevsky, "Normalization of
%           optical fluence distribution for three-dimensional functional
%           optoacoustic tomography of the breast," J. Biomed. Opt. 27(3)
%           036001 (16 March 2022)
%           https://doi.org/10.1117/1.JBO.27.3.036001
%
% Copyright (C) 2022 Seonyeong Park
% License:  GNU General Public License version 3, Please see 'LICENSE' for
%           details.
%

addpath('functions');   % Add a path of functions
addpath('frangi_filter_version2a'); % Add a path of vessel enhancement filter

voxel_size = 0.25; % Voxel size [mm]
depth_vis = 30;   % Maximum depth for visualization [mm]

% File name
filename_p0_w757 = fullfile('data', 'RECON_NOISY_w757_FBP.mat');
filename_p0_w800 = fullfile('data', 'RECON_NOISY_w800_FBP.mat');
filename_p0_w850 = fullfile('data', 'RECON_NOISY_w850_FBP.mat');
filename_phi0_est_w757 = fullfile('data', 'phi0_est_w757.mat');
filename_phi0_est_w800 = fullfile('data', 'phi0_est_w800.mat');
filename_phi0_est_w850 = fullfile('data', 'phi0_est_w850.mat');
filename_phia_est_w757 = fullfile('data', 'phia_est_w757.mat');
filename_phia_est_w800 = fullfile('data', 'phia_est_w800.mat');
filename_phia_est_w850 = fullfile('data', 'phia_est_w850.mat');
filename_vessel_filter_non_w757_w800_w850 = ...
    fullfile('data', 'vessel_filter_non_w757_w800_w850.mat');
filename_vessel_filter_on_w757_w800_w850 = ...
    fullfile('data', 'vessel_filter_on_w757_w800_w850.mat');
filename_depth_w757 = fullfile('data', 'depth_w757.mat');
filename_depth_w800 = fullfile('data', 'depth_w800.mat');
filename_depth_w850 = fullfile('data', 'depth_w850.mat');
filename_bv_non_mvbpx_w757 = ...
    fullfile('data', 'vessels_color-encoded_depth_non_mvbpx_w757.png');
filename_bv_non_mvbpy_w757 = ...
    fullfile('data', 'vessels_color-encoded_depth_non_mvbpy_w757.png');
filename_bv_non_mvbpz_w757 = ...
    fullfile('data', 'vessels_color-encoded_depth_non_mvbpz_w757.png');
filename_bv_non_mvbpx_w800 = ...
    fullfile('data', 'vessels_color-encoded_depth_non_mvbpx_w800.png');
filename_bv_non_mvbpy_w800 = ...
    fullfile('data', 'vessels_color-encoded_depth_non_mvbpy_w800.png');
filename_bv_non_mvbpz_w800 = ...
    fullfile('data', 'vessels_color-encoded_depth_non_mvbpz_w800.png');
filename_bv_non_mvbpx_w850 = ...
    fullfile('data', 'vessels_color-encoded_depth_non_mvbpx_w850.png');
filename_bv_non_mvbpy_w850 = ...
    fullfile('data', 'vessels_color-encoded_depth_non_mvbpy_w850.png');
filename_bv_non_mvbpz_w850 = ...
    fullfile('data', 'vessels_color-encoded_depth_non_mvbpz_w850.png');
filename_bv_on_mvbpx_w757 = ...
    fullfile('data', 'vessels_color-encoded_depth_on_mvbpx_w757.png');
filename_bv_on_mvbpy_w757 = ...
    fullfile('data', 'vessels_color-encoded_depth_on_mvbpy_w757.png');
filename_bv_on_mvbpz_w757 = ...
    fullfile('data', 'vessels_color-encoded_depth_on_mvbpz_w757.png');
filename_bv_on_mvbpx_w800 = ...
    fullfile('data', 'vessels_color-encoded_depth_on_mvbpx_w800.png');
filename_bv_on_mvbpy_w800 = ...
    fullfile('data', 'vessels_color-encoded_depth_on_mvbpy_w800.png');
filename_bv_on_mvbpz_w800 = ...
    fullfile('data', 'vessels_color-encoded_depth_on_mvbpz_w800.png');
filename_bv_on_mvbpx_w850 = ...
    fullfile('data', 'vessels_color-encoded_depth_on_mvbpx_w850.png');
filename_bv_on_mvbpy_w850 = ...
    fullfile('data', 'vessels_color-encoded_depth_on_mvbpy_w850.png');
filename_bv_on_mvbpz_w850 = ...
    fullfile('data', 'vessels_color-encoded_depth_on_mvbpz_w850.png');


% Load estimated intial pressure, reconstructed image
fprintf('Loading reconstructed images...\n');
load(filename_p0_w757, 'recon'); p0_est_w757 = recon;
load(filename_p0_w800, 'recon'); p0_est_w800 = recon;
load(filename_p0_w850, 'recon'); p0_est_w850 = recon;

clear recon;

% Load depth map
fprintf('Loading depth maps...\n');
load(filename_depth_w757); depth_w757 = depth;
doi_w757 = depth_w757; % Depth of interest
doi_w757(doi_w757 > depth_vis/voxel_size)= depth_vis/voxel_size;
load(filename_depth_w800); depth_w800 = depth;
doi_w800 = depth_w800; % Depth of interest
doi_w800(doi_w800 > depth_vis/voxel_size)= depth_vis/voxel_size;
load(filename_depth_w850); depth_w850 = depth;
doi_w850 = depth_w850; % Depth of interest
doi_w850(doi_w850 > depth_vis/voxel_size)= depth_vis/voxel_size;

clear depth;

% Load estimated incident optical fluence and optical attenuation
fprintf('Loading estimated incident optical fluence...\n');
load(filename_phi0_est_w757); phi0_est_w757 = phi0_est;
load(filename_phi0_est_w800); phi0_est_w800 = phi0_est;
load(filename_phi0_est_w850); phi0_est_w850 = phi0_est;

fprintf('Loading estimated optical attenuation...\n');
load(filename_phia_est_w757, 'f_phia');
phia_est_w757 = exp(-f_phia.b.*doi_w757); % Only for visualization purpose
load(filename_phia_est_w800, 'f_phia');
phia_est_w800 = exp(-f_phia.b.*doi_w800); % Only for visualization purpose
load(filename_phia_est_w850, 'f_phia'); 
phia_est_w850 = exp(-f_phia.b.*doi_w850); % Only for visualization purpose

clear phi0_est f_phia;

% Optical fluence normalization
fprintf('Applying optical fluence normalization...\n');
mu_a_est_w757 = p0_est_w757./phi0_est_w757./phia_est_w757;
mu_a_est_w800 = p0_est_w800./phi0_est_w800./phia_est_w800;
mu_a_est_w850 = p0_est_w850./phi0_est_w850./phia_est_w850;

% Load detected blood vessel mask
fprintf('Loading detected blood vessel masks...\n');
load(filename_vessel_filter_non_w757_w800_w850, 'mask_non_w757_w800_w850');
load(filename_vessel_filter_on_w757_w800_w850, 'mask_on_w757_w800_w850');

% Apply detected blood vessel mask
fprintf('Applying detected blood vessel masks...\n');
p0_est_w757 = p0_est_w757.*mask_non_w757_w800_w850;
p0_est_w800 = p0_est_w800.*mask_non_w757_w800_w850;
p0_est_w850 = p0_est_w850.*mask_non_w757_w800_w850;
mu_a_est_w757 = mu_a_est_w757.*mask_on_w757_w800_w850;
mu_a_est_w800 = mu_a_est_w800.*mask_on_w757_w800_w850;
mu_a_est_w850 = mu_a_est_w850.*mask_on_w757_w800_w850;

% Normalize brightness
fprintf('Normalizing voxel brightness...\n');
p0_est_w757 = p0_est_w757./max(p0_est_w757(:));
p0_est_w800 = p0_est_w800./max(p0_est_w800(:));
p0_est_w850 = p0_est_w850./max(p0_est_w850(:));
mu_a_est_w757 = mu_a_est_w757./max(mu_a_est_w757(:));
mu_a_est_w800 = mu_a_est_w800./max(mu_a_est_w800(:));
mu_a_est_w850 = mu_a_est_w850./max(mu_a_est_w850(:));


% Plot maximum voxel brightness projection (MVBP), color encoded by depth
fprintf('Calculating maximum voxel brightness projections color-encoded by depth...\n');
% No optical fluence normalization
p0_est_w757 = p0_est_w757/0.6; % Clip brightness for visualization
p0_est_w757(p0_est_w757 > 1) = 1;
[mvbpx_non_w757, mvbpy_non_w757, mvbpz_non_w757] = ...
    MVBPColorEncodedDepth(p0_est_w757, voxel_size, depth_w757, ...
    depth_vis/voxel_size);

p0_est_w800 = p0_est_w800/0.6; % Clip brightness for visualization
p0_est_w800(p0_est_w800 > 1) = 1;
[mvbpx_non_w800, mvbpy_non_w800, mvbpz_non_w800] = ...
    MVBPColorEncodedDepth(p0_est_w800, voxel_size, depth_w800, ...
    depth_vis/voxel_size);

p0_est_w850 = p0_est_w850/0.6; % Clip brightness for visualization
p0_est_w850(p0_est_w850 > 1) = 1;
[mvbpx_non_w850, mvbpy_non_w850, mvbpz_non_w850] = ...
    MVBPColorEncodedDepth(p0_est_w850, voxel_size, depth_w850, ...
    depth_vis/voxel_size);

% Optical fluence normalization
mu_a_est_w757 = mu_a_est_w757/0.6; % Clip brightness for visualization
mu_a_est_w757(mu_a_est_w757 > 1) = 1;
[mvbpx_on_w757, mvbpy_on_w757, mvbpz_on_w757] = ...
    MVBPColorEncodedDepth(mu_a_est_w757, voxel_size, depth_w757, ...
    depth_vis/voxel_size);

mu_a_est_w800 = mu_a_est_w800/0.6; % Clip brightness for visualization
mu_a_est_w800(mu_a_est_w800 > 1) = 1;
[mvbpx_on_w800, mvbpy_on_w800, mvbpz_on_w800] = ...
    MVBPColorEncodedDepth(mu_a_est_w800, voxel_size, depth_w800, ...
    depth_vis/voxel_size);

mu_a_est_w850 = mu_a_est_w850/0.6; % Clip brightness for visualization
mu_a_est_w850(mu_a_est_w850 > 1) = 1;
[mvbpx_on_w850, mvbpy_on_w850, mvbpz_on_w850] = ...
    MVBPColorEncodedDepth(mu_a_est_w850, voxel_size, depth_w850, ...
    depth_vis/voxel_size);


% Save MVBP, color encoded by depth
fprintf('Saving maximum voxel brightness projections color-encoded by depth...\n');
% No optical fluence normalization
imwrite(mvbpx_non_w757, filename_bv_non_mvbpx_w757);
imwrite(mvbpy_non_w757, filename_bv_non_mvbpy_w757);
imwrite(mvbpz_non_w757, filename_bv_non_mvbpz_w757);
imwrite(mvbpx_non_w800, filename_bv_non_mvbpx_w800);
imwrite(mvbpy_non_w800, filename_bv_non_mvbpy_w800);
imwrite(mvbpz_non_w800, filename_bv_non_mvbpz_w800);
imwrite(mvbpx_non_w850, filename_bv_non_mvbpx_w850);
imwrite(mvbpy_non_w850, filename_bv_non_mvbpy_w850);
imwrite(mvbpz_non_w850, filename_bv_non_mvbpz_w850);

% Optical fluence normalization
imwrite(mvbpx_on_w757, filename_bv_on_mvbpx_w757);
imwrite(mvbpy_on_w757, filename_bv_on_mvbpy_w757);
imwrite(mvbpz_on_w757, filename_bv_on_mvbpz_w757);
imwrite(mvbpx_on_w800, filename_bv_on_mvbpx_w800);
imwrite(mvbpy_on_w800, filename_bv_on_mvbpy_w800);
imwrite(mvbpz_on_w800, filename_bv_on_mvbpz_w800);
imwrite(mvbpx_on_w850, filename_bv_on_mvbpx_w850);
imwrite(mvbpy_on_w850, filename_bv_on_mvbpy_w850);
imwrite(mvbpz_on_w850, filename_bv_on_mvbpz_w850);
