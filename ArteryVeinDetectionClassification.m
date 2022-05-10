% -------------------------------------------------------------------------
%                 Artery Vein Detection and Classification
% -------------------------------------------------------------------------
% Author:   Seonyeong Park (sp33@illinois.edu)
% Date:     Jan. 24, 2022
%
% This project is to detect blood vessels via estimation of total
% hemoglobin concentration and to classify them into arteries and veins via
% estimation of oxygen saturation, from three-dimensional (3D) optoacoustic
% tomography (OAT) images of the breast at two and three illumination
% wavelengths. The total hemoglobin concentration and oxygen saturation are
% estimated using spectral linear unmixing. The molar extinction
% coefficients of deoxy- and oxyhemoglobin used in the spectral linear
% unmixing are from the reference [OMLC]. Multiscale vessel enhancement
% filtering, also known as Frangi filtering [Frangi], and Otsu thresholding
% are applied to the estimate of total hemoglobin concentration to detect
% the blood vessels. Before running this code, the function package
% of the multiscale vessel enhancement filter [Fraingi] should be
% downloaded, and the 'eig3volume.c' file should be compiled ('mex
% eig3volume.c'). The method of the software is detailed in [Park2022].
%
% Reference:
%   [Park2022] Seonyeong Park, Frank J. Brooks, Umberto Villa, Richard Su,
%           Mark A. Anastasio, Alexander A. Oraevsky, "Normalization of
%           optical fluence distribution for three-dimensional functional
%           optoacoustic tomography of the breast," J. Biomed. Opt. 27(3)
%           036001 (16 March 2022)
%           https://doi.org/10.1117/1.JBO.27.3.036001
%
%   [OMLC]  https://omlc.org/spectra/hemoglobin/summary.html
%
%   [Frangi] Dirk-Jan Kroon (2022). Hessian based Frangi Vesselness filter
%           (https://www.mathworks.com/matlabcentral/fileexchange/24409-
%           hessian-based-frangi-vesselness-filter), MATLAB Central File
%           Exchange. Retrieved April 15, 2022.
%
% Copyright (C) 2022 Seonyeong Park
% License:  GNU General Public License version 3, Please see 'LICENSE' for
%           details.
%

addpath('functions');   % Add a path of functions
addpath('frangi_filter_version2a'); % Add a path of vessel enhancement filter

flag_fig  = true;       % Flag to plot
% flag_fig  = false;      % Flag to plot

% Molar extinction coffient matix
eps = [156.048,	56.8; ...   % 757 nm
    76.172, 81.6; ...       % 800 nm
    69.132, 105.8];         % 850 nm

% True oxygen saturation
s_artery    = 0.97;     % Artery
s_vein      = 0.7;      % Vein
threshold = mean([s_artery, s_vein]);

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
% No optical fluence normalization
filename_vessel_filter_non_w757_w800 = ...
    fullfile('data', 'vessel_filter_non_w757_w800.mat');
filename_vessel_filter_non_w757_w850 = ...
    fullfile('data', 'vessel_filter_non_w757_w850.mat');
filename_vessel_filter_non_w800_w850 = ...
    fullfile('data', 'vessel_filter_non_w800_w850.mat');
filename_vessel_filter_non_w757_w800_w850 = ...
    fullfile('data', 'vessel_filter_non_w757_w800_w850.mat');
% Optical fluence normalization
filename_vessel_filter_on_w757_w800 = ...
    fullfile('data', 'vessel_filter_on_w757_w800.mat');
filename_vessel_filter_on_w757_w850 = ...
    fullfile('data', 'vessel_filter_on_w757_w850.mat');
filename_vessel_filter_on_w800_w850 = ...
    fullfile('data', 'vessel_filter_on_w800_w850.mat');
filename_vessel_filter_on_w757_w800_w850 = ...
    fullfile('data', 'vessel_filter_on_w757_w800_w850.mat');


% Load estimated intial pressure, reconstructed image
fprintf('Loading reconstructed images...\n');
load(filename_p0_w757, 'recon'); p0_est_w757 = recon;
load(filename_p0_w800, 'recon'); p0_est_w800 = recon;
load(filename_p0_w850, 'recon'); p0_est_w850 = recon;
clear recon;

% Load estimated incident optical fluence and optical attenuation
fprintf('Loading estimated incident optical fluence...\n');
load(filename_phi0_est_w757); phi0_est_w757 = phi0_est;
load(filename_phi0_est_w800); phi0_est_w800 = phi0_est;
load(filename_phi0_est_w850); phi0_est_w850 = phi0_est;
clear phi0_est;

fprintf('Loading estimated optical attenuation...\n');
load(filename_phia_est_w757); phia_est_w757 = phia_est;
load(filename_phia_est_w800); phia_est_w800 = phia_est;
load(filename_phia_est_w850); phia_est_w850 = phia_est;
clear phia_est;

% Estimate optical absorption coefficient
fprintf('Estimating optical absorption coefficient...\n');
mu_a_est_w757 = smooth3(p0_est_w757./phi0_est_w757./phia_est_w757);
mu_a_est_w800 = smooth3(p0_est_w800./phi0_est_w800./phia_est_w800);
mu_a_est_w850 = smooth3(p0_est_w850./phi0_est_w850./phia_est_w850);

%% Spectral linear unmixing
fprintf('Computing Spectral linear unmixing...\n');

% No optical normalization
fprintf('    With no optical normaliation...\n');
[s_est_non_w757_w800, c_thb_est_non_w757_w800] = ... % 757, 800 nm
    LinearUnmixing2([eps(1, :); eps(2, :)], p0_est_w757, p0_est_w800);
[s_est_non_w757_w850, c_thb_est_non_w757_w850] = ... % 757, 850 nm
    LinearUnmixing2([eps(1, :); eps(3, :)], p0_est_w757, p0_est_w850);
[s_est_non_w800_w850, c_thb_est_non_w800_w850] = ... % 800, 850 nm
    LinearUnmixing2([eps(2, :); eps(3, :)], p0_est_w800, p0_est_w850);
[s_est_non_w757_w800_w850, c_thb_est_non_w757_w800_w850] = ... % 757, 800, 850 nm
    LinearUnmixing3(eps, p0_est_w757, p0_est_w800, p0_est_w850);

% Optical normalization
fprintf('    With optical normaliation...\n');
[s_est_on_w757_w800, c_thb_est_on_w757_w800] = ... % 757, 800 nm
    LinearUnmixing2([eps(1, :); eps(2, :)], mu_a_est_w757, mu_a_est_w800);
[s_est_on_w757_w850, c_thb_est_on_w757_w850] = ... % 757, 850 nm
    LinearUnmixing2([eps(1, :); eps(3, :)], mu_a_est_w757, mu_a_est_w850);
[s_est_on_w800_w850, c_thb_est_on_w800_w850] = ... % 800, 850 nm
    LinearUnmixing2([eps(2, :); eps(3, :)], mu_a_est_w800, mu_a_est_w850);
[s_est_on_w757_w800_w850, c_thb_est_on_w757_w800_w850] = ... % 757, 800, 850 nm
    LinearUnmixing3(eps, mu_a_est_w757, mu_a_est_w800, mu_a_est_w850);

%% Multiscale vessel enhancement filtering
fprintf('Applying multiscale vessel filtering...\n');
options.FrangiScaleRange = [1, 5];
options.FrangiScaleRatio = 1;
options.BlackWhite = false;

% No optical normalization
fprintf('    With no optical normaliation...\n');
[filter_non_w757_w800, ~, ~, ~, ~] = ...
    FrangiFilter3D(c_thb_est_non_w757_w800, options); % 757, 800 nm
[filter_non_w757_w850, ~, ~, ~, ~] = ...
    FrangiFilter3D(c_thb_est_non_w757_w850, options); % 757, 850 nm
[filter_non_w800_w850, ~, ~, ~, ~] = ...
    FrangiFilter3D(c_thb_est_non_w800_w850, options); % 800, 850 nm
[filter_non_w757_w800_w850, ~, ~, ~, ~] = ...
    FrangiFilter3D(c_thb_est_non_w757_w800_w850, options); % 757, 800, 850 nm

% Optical normalization
fprintf('    With optical normaliation...\n');
[filter_on_w757_w800, ~, ~, ~, ~] = ...
    FrangiFilter3D(c_thb_est_on_w757_w800, options); % 757, 800 nm
[filter_on_w757_w850, ~, ~, ~, ~] = ...
    FrangiFilter3D(c_thb_est_on_w757_w850, options); % 757, 850 nm
[filter_on_w800_w850, ~, ~, ~, ~] = ...
    FrangiFilter3D(c_thb_est_on_w800_w850, options); % 800, 850 nm
[filter_on_w757_w800_w850, ~, ~, ~, ~] = ...
    FrangiFilter3D(c_thb_est_on_w757_w800_w850, options); % 757, 800, 850 nm

%% Otsu thresholding
fprintf('Applying Otsu thresholding...\n');

% No optical normalization
fprintf('    With no optical normaliation...\n');
mask_non_w757_w800 = imbinarize( ...
    filter_non_w757_w800./max(filter_non_w757_w800(:)));
mask_non_w757_w850 = imbinarize( ...
    filter_non_w757_w850./max(filter_non_w757_w850(:)));
mask_non_w800_w850 = imbinarize( ...
    filter_non_w800_w850./max(filter_non_w800_w850(:)));
mask_non_w757_w800_w850 = imbinarize(...
    filter_non_w757_w800_w850./max(filter_non_w757_w800_w850(:)));

% Optical normalization
fprintf('    With optical normaliation...\n');
mask_on_w757_w800 = imbinarize( ...
    filter_on_w757_w800./max(filter_on_w757_w800(:)));
mask_on_w757_w850 = imbinarize( ...
    filter_on_w757_w850./max(filter_on_w757_w850(:)));
mask_on_w800_w850 = imbinarize( ...
    filter_on_w800_w850./max(filter_on_w757_w850(:)));
mask_on_w757_w800_w850 = imbinarize( ...
    filter_on_w757_w800_w850./max(filter_on_w757_w800_w850(:)));


% Save vessel enhancement filtering result and vessel mask
fprintf('Saving vessel enhancement filtering result and vessel mask...\n');

% No optical normalization
fprintf('    With no optical normaliation...\n');
save(filename_vessel_filter_non_w757_w800, ...
    'filter_non_w757_w800', 'mask_non_w757_w800');
save(filename_vessel_filter_non_w757_w850, ...
    'filter_non_w757_w850', 'mask_non_w757_w850');
save(filename_vessel_filter_non_w800_w850, ...
    'filter_non_w800_w850', 'mask_non_w800_w850');
save(filename_vessel_filter_non_w757_w800_w850, ...
    'filter_non_w757_w800_w850', 'mask_non_w757_w800_w850');

% Optical normalization
fprintf('    With optical normaliation...\n');
save(filename_vessel_filter_on_w757_w800, ...
    'filter_on_w757_w800', 'mask_on_w757_w800');
save(filename_vessel_filter_on_w757_w850, ...
    'filter_on_w757_w850', 'mask_on_w757_w850');
save(filename_vessel_filter_on_w800_w850, ...
    'filter_on_w800_w850', 'mask_on_w800_w850');
save(filename_vessel_filter_on_w757_w800_w850, ...
    'filter_on_w757_w800_w850', 'mask_on_w757_w800_w850');

if flag_fig == 1
    % Plot detected blooe vessels
    % No optical fluence normalization
    figure; % 757, 800 nm
    subplot(1, 3, 1); imagesc(squeeze(max(mask_non_w757_w800,[],1)));
    xlabel('$z$ [voxel]', 'Interpreter', 'latex');
    ylabel('$y$ [voxel]', 'Interpreter', 'latex');
    colorbar('TickLabelInterpreter', 'latex'); colormap hot; axis image;
    set(gca, 'FontSize', 14, 'FontName', 'Times', ...
        'TickLabelInterpreter', 'latex');
    subplot(1, 3, 2); imagesc(squeeze(max(mask_non_w757_w800,[],2)));
    xlabel('$z$ [voxel]', 'Interpreter', 'latex');
    ylabel('$x$ [voxel]', 'Interpreter', 'latex');
    colorbar('TickLabelInterpreter', 'latex'); colormap hot; axis image;
    set(gca, 'FontSize', 14, 'FontName', 'Times', ...
        'TickLabelInterpreter', 'latex');
    subplot(1, 3, 3); imagesc(squeeze(max(mask_non_w757_w800,[],3)));
    xlabel('$y$ [voxel]', 'Interpreter', 'latex');
    ylabel('$x$ [voxel]', 'Interpreter', 'latex');
    colorbar('TickLabelInterpreter', 'latex'); colormap hot; axis image;
    set(gca, 'FontSize', 14, 'FontName', 'Times', ...
        'TickLabelInterpreter', 'latex');
    sgtitle({'Detected blood vessels, No optical fluence normalization', ...
        '757 nm and 800 nm'}, 'Interpreter', 'latex');
    
    figure; % 757, 850 nm
    subplot(1, 3, 1); imagesc(squeeze(max(mask_non_w757_w850,[],1)));
    xlabel('$z$ [voxel]', 'Interpreter', 'latex');
    ylabel('$y$ [voxel]', 'Interpreter', 'latex');
    colorbar('TickLabelInterpreter', 'latex'); colormap hot; axis image;
    set(gca, 'FontSize', 14, 'FontName', 'Times', ...
        'TickLabelInterpreter', 'latex');
    subplot(1, 3, 2); imagesc(squeeze(max(mask_non_w757_w850,[],2)));
    xlabel('$z$ [voxel]', 'Interpreter', 'latex');
    ylabel('$x$ [voxel]', 'Interpreter', 'latex');
    colorbar('TickLabelInterpreter', 'latex'); colormap hot; axis image;
    set(gca, 'FontSize', 14, 'FontName', 'Times', ...
        'TickLabelInterpreter', 'latex');
    subplot(1, 3, 3); imagesc(squeeze(max(mask_non_w757_w850,[],3)));
    xlabel('$y$ [voxel]', 'Interpreter', 'latex');
    ylabel('$x$ [voxel]', 'Interpreter', 'latex');
    colorbar('TickLabelInterpreter', 'latex'); colormap hot; axis image;
    set(gca, 'FontSize', 14, 'FontName', 'Times', ...
        'TickLabelInterpreter', 'latex');
    sgtitle({'Detected blood vessels, No optical fluence normalization', ...
        '757 nm and 850 nm'}, 'Interpreter', 'latex');
    
    figure; % 800, 850 nm
    subplot(1, 3, 1); imagesc(squeeze(max(mask_non_w800_w850,[],1)));
    xlabel('$z$ [voxel]', 'Interpreter', 'latex');
    ylabel('$y$ [voxel]', 'Interpreter', 'latex');
    colorbar('TickLabelInterpreter', 'latex'); colormap hot; axis image;
    set(gca, 'FontSize', 14, 'FontName', 'Times', ...
        'TickLabelInterpreter', 'latex');
    subplot(1, 3, 2); imagesc(squeeze(max(mask_non_w800_w850,[],2)));
    xlabel('$z$ [voxel]', 'Interpreter', 'latex');
    ylabel('$x$ [voxel]', 'Interpreter', 'latex');
    colorbar('TickLabelInterpreter', 'latex'); colormap hot; axis image;
    set(gca, 'FontSize', 14, 'FontName', 'Times', ...
        'TickLabelInterpreter', 'latex');
    subplot(1, 3, 3); imagesc(squeeze(max(mask_non_w800_w850,[],3)));
    xlabel('$y$ [voxel]', 'Interpreter', 'latex');
    ylabel('$x$ [voxel]', 'Interpreter', 'latex');
    colorbar('TickLabelInterpreter', 'latex'); colormap hot; axis image;
    set(gca, 'FontSize', 14, 'FontName', 'Times', ...
        'TickLabelInterpreter', 'latex');
    sgtitle({'Detected blood vessels, No optical fluence normalization', ...
        '800 nm and 850 nm'}, 'Interpreter', 'latex');
    
    figure; % 757, 800, 850 nm
    subplot(1, 3, 1); imagesc(squeeze(max(mask_non_w757_w800_w850,[],1)));
    xlabel('$z$ [voxel]', 'Interpreter', 'latex');
    ylabel('$y$ [voxel]', 'Interpreter', 'latex');
    colorbar('TickLabelInterpreter', 'latex'); colormap hot; axis image;
    set(gca, 'FontSize', 14, 'FontName', 'Times', ...
        'TickLabelInterpreter', 'latex');
    subplot(1, 3, 2); imagesc(squeeze(max(mask_non_w757_w800_w850,[],2)));
    xlabel('$z$ [voxel]', 'Interpreter', 'latex');
    ylabel('$x$ [voxel]', 'Interpreter', 'latex');
    colorbar('TickLabelInterpreter', 'latex'); colormap hot; axis image;
    set(gca, 'FontSize', 14, 'FontName', 'Times', ...
        'TickLabelInterpreter', 'latex');
    subplot(1, 3, 3); imagesc(squeeze(max(mask_non_w757_w800_w850,[],3)));
    xlabel('$y$ [voxel]', 'Interpreter', 'latex');
    ylabel('$x$ [voxel]', 'Interpreter', 'latex');
    colorbar('TickLabelInterpreter', 'latex'); colormap hot; axis image;
    set(gca, 'FontSize', 14, 'FontName', 'Times', ...
        'TickLabelInterpreter', 'latex');
    sgtitle({'Detected blood vessels, No optical fluence normalization', ...
        '757 nm, 800 nm, and 850 nm'}, 'Interpreter', 'latex');
    
    
    % Optical fluence normalization
    figure; % 757, 800 nm
    subplot(1, 3, 1); imagesc(squeeze(max(mask_on_w757_w800,[],1)));
    xlabel('$z$ [voxel]', 'Interpreter', 'latex');
    ylabel('$y$ [voxel]', 'Interpreter', 'latex');
    colorbar('TickLabelInterpreter', 'latex'); colormap hot; axis image;
    set(gca, 'FontSize', 14, 'FontName', 'Times', ...
        'TickLabelInterpreter', 'latex');
    subplot(1, 3, 2); imagesc(squeeze(max(mask_on_w757_w800,[],2)));
    xlabel('$z$ [voxel]', 'Interpreter', 'latex');
    ylabel('$x$ [voxel]', 'Interpreter', 'latex');
    colorbar('TickLabelInterpreter', 'latex'); colormap hot; axis image;
    set(gca, 'FontSize', 14, 'FontName', 'Times', ...
        'TickLabelInterpreter', 'latex');
    subplot(1, 3, 3); imagesc(squeeze(max(mask_on_w757_w800,[],3)));
    xlabel('$y$ [voxel]', 'Interpreter', 'latex');
    ylabel('$x$ [voxel]', 'Interpreter', 'latex');
    colorbar('TickLabelInterpreter', 'latex'); colormap hot; axis image;
    set(gca, 'FontSize', 14, 'FontName', 'Times', ...
        'TickLabelInterpreter', 'latex');
    sgtitle({'Detected blood vessels, Optical fluence normalization', ...
        '757 nm and 800 nm'}, 'Interpreter', 'latex');
    
    figure; % 757, 850 nm
    subplot(1, 3, 1); imagesc(squeeze(max(mask_on_w757_w850,[],1)));
    xlabel('$z$ [voxel]', 'Interpreter', 'latex');
    ylabel('$y$ [voxel]', 'Interpreter', 'latex');
    colorbar('TickLabelInterpreter', 'latex'); colormap hot; axis image;
    set(gca, 'FontSize', 14, 'FontName', 'Times', ...
        'TickLabelInterpreter', 'latex');
    subplot(1, 3, 2); imagesc(squeeze(max(mask_on_w757_w850,[],2)));
    xlabel('$z$ [voxel]', 'Interpreter', 'latex');
    ylabel('$x$ [voxel]', 'Interpreter', 'latex');
    colorbar('TickLabelInterpreter', 'latex'); colormap hot; axis image;
    set(gca, 'FontSize', 14, 'FontName', 'Times', ...
        'TickLabelInterpreter', 'latex');
    subplot(1, 3, 3); imagesc(squeeze(max(mask_on_w757_w850,[],3)));
    xlabel('$y$ [voxel]', 'Interpreter', 'latex');
    ylabel('$x$ [voxel]', 'Interpreter', 'latex');
    colorbar('TickLabelInterpreter', 'latex'); colormap hot; axis image;
    set(gca, 'FontSize', 14, 'FontName', 'Times', ...
        'TickLabelInterpreter', 'latex');
    sgtitle({'Detected blood vessels, Optical fluence normalization', ...
        '757 nm and 850 nm'}, 'Interpreter', 'latex');
    
    figure; % 800, 850 nm
    subplot(1, 3, 1); imagesc(squeeze(max(mask_on_w800_w850,[],1)));
    xlabel('$z$ [voxel]', 'Interpreter', 'latex');
    ylabel('$y$ [voxel]', 'Interpreter', 'latex');
    colorbar('TickLabelInterpreter', 'latex'); colormap hot; axis image;
    set(gca, 'FontSize', 14, 'FontName', 'Times', ...
        'TickLabelInterpreter', 'latex');
    subplot(1, 3, 2); imagesc(squeeze(max(mask_on_w800_w850,[],2)));
    xlabel('$z$ [voxel]', 'Interpreter', 'latex');
    ylabel('$x$ [voxel]', 'Interpreter', 'latex');
    colorbar('TickLabelInterpreter', 'latex'); colormap hot; axis image;
    set(gca, 'FontSize', 14, 'FontName', 'Times', ...
        'TickLabelInterpreter', 'latex');
    subplot(1, 3, 3); imagesc(squeeze(max(mask_on_w800_w850,[],3)));
    xlabel('$y$ [voxel]', 'Interpreter', 'latex');
    ylabel('$x$ [voxel]', 'Interpreter', 'latex');
    colorbar('TickLabelInterpreter', 'latex'); colormap hot; axis image;
    set(gca, 'FontSize', 14, 'FontName', 'Times', ...
        'TickLabelInterpreter', 'latex');
    sgtitle({'Detected blood vessels, Optical fluence normalization', ...
        '800 nm and 850 nm'}, 'Interpreter', 'latex');
    
    figure; % 757, 800, 850 nm
    subplot(1, 3, 1); imagesc(squeeze(max(mask_on_w757_w800_w850,[],1)));
    xlabel('$z$ [voxel]', 'Interpreter', 'latex');
    ylabel('$y$ [voxel]', 'Interpreter', 'latex');
    colorbar('TickLabelInterpreter', 'latex'); colormap hot; axis image;
    set(gca, 'FontSize', 14, 'FontName', 'Times', ...
        'TickLabelInterpreter', 'latex');
    subplot(1, 3, 2); imagesc(squeeze(max(mask_on_w757_w800_w850,[],2)));
    xlabel('$z$ [voxel]', 'Interpreter', 'latex');
    ylabel('$x$ [voxel]', 'Interpreter', 'latex');
    colorbar('TickLabelInterpreter', 'latex'); colormap hot; axis image;
    set(gca, 'FontSize', 14, 'FontName', 'Times', ...
        'TickLabelInterpreter', 'latex');
    subplot(1, 3, 3); imagesc(squeeze(max(mask_on_w757_w800_w850,[],3)));
    xlabel('$y$ [voxel]', 'Interpreter', 'latex');
    ylabel('$x$ [voxel]', 'Interpreter', 'latex');
    colorbar('TickLabelInterpreter', 'latex'); colormap hot; axis image;
    set(gca, 'FontSize', 14, 'FontName', 'Times', ...
        'TickLabelInterpreter', 'latex');
    sgtitle({'Detected blood vessels, Optical fluence normalization', ...
        '757 nm, 800 nm, and 850 nm'}, 'Interpreter', 'latex');
end

% Classify arteries and veins
fprintf('Classifying arteries and veins...\n');
% No optical normalization
fprintf('    With no optical normaliation...\n');
% 757, 800 nm
s_est_non_w757_w800_av = s_est_non_w757_w800;
s_est_non_w757_w800_av(s_est_non_w757_w800 > threshold)  = s_artery;
s_est_non_w757_w800_av(s_est_non_w757_w800 <= threshold) = s_vein;
s_est_non_w757_w800_av = s_est_non_w757_w800_av.*mask_non_w757_w800;
% 757, 850 nm
s_est_non_w757_w850_av = s_est_non_w757_w850;
s_est_non_w757_w850_av(s_est_non_w757_w850 > threshold)  = s_artery;
s_est_non_w757_w850_av(s_est_non_w757_w850 <= threshold) = s_vein;
s_est_non_w757_w850_av = s_est_non_w757_w850_av.*mask_non_w757_w850;
% 800, 850 nm
s_est_non_w800_w850_av = s_est_non_w800_w850;
s_est_non_w800_w850_av(s_est_non_w800_w850 > threshold)  = s_artery;
s_est_non_w800_w850_av(s_est_non_w800_w850 <= threshold) = s_vein;
s_est_non_w800_w850_av = s_est_non_w800_w850_av.*mask_non_w800_w850;
% 757, 800, 850 nm
s_est_non_w757_w800_w850_av = s_est_non_w757_w800_w850;
s_est_non_w757_w800_w850_av(s_est_non_w757_w800_w850 > threshold)  = s_artery;
s_est_non_w757_w800_w850_av(s_est_non_w757_w800_w850 <= threshold) = s_vein;
s_est_non_w757_w800_w850_av = s_est_non_w757_w800_w850_av.*mask_non_w757_w800_w850;

% Optical normalization
fprintf('    With optical normaliation...\n');
% 757, 800 nm
s_est_on_w757_w800_av = s_est_on_w757_w800;
s_est_on_w757_w800_av(s_est_on_w757_w800 > threshold)  = s_artery;
s_est_on_w757_w800_av(s_est_on_w757_w800 <= threshold) = s_vein;
s_est_on_w757_w800_av = s_est_on_w757_w800_av.*mask_on_w757_w800;
% 757, 850 nm
s_est_on_w757_w850_av = s_est_on_w757_w850;
s_est_on_w757_w850_av(s_est_on_w757_w850 > threshold)  = s_artery;
s_est_on_w757_w850_av(s_est_on_w757_w850 <= threshold) = s_vein;
s_est_on_w757_w850_av = s_est_on_w757_w850_av.*mask_on_w757_w850;
% 800, 850 nm
s_est_on_w800_w850_av = s_est_on_w800_w850;
s_est_on_w800_w850_av(s_est_on_w800_w850 > threshold)  = s_artery;
s_est_on_w800_w850_av(s_est_on_w800_w850 <= threshold) = s_vein;
s_est_on_w800_w850_av = s_est_on_w800_w850_av.*mask_on_w800_w850;
% 757, 800, 850 nm
s_est_on_w757_w800_w850_av = s_est_on_w757_w800_w850;
s_est_on_w757_w800_w850_av(s_est_on_w757_w800_w850 > threshold)  = s_artery;
s_est_on_w757_w800_w850_av(s_est_on_w757_w800_w850 <= threshold) = s_vein;
s_est_on_w757_w800_w850_av = s_est_on_w757_w800_w850_av.*mask_on_w757_w800_w850;

if flag_fig == 1
    % Plot classified arteries and veins
    % No optical fluence normalization
    figure;
    subplot(1, 3, 1);
    imagesc(squeeze(max(s_est_non_w757_w800_av,[],1)), [0.5, 1]);
    xlabel('$z$ [voxel]', 'Interpreter', 'latex');
    ylabel('$y$ [voxel]', 'Interpreter', 'latex');
    colorbar('TickLabelInterpreter', 'latex'); colormap jet; colorbar; axis image;
    set(gca, 'FontSize', 14, 'FontName', 'Times', ...
        'TickLabelInterpreter', 'latex');
    subplot(1, 3, 2);
    imagesc(squeeze(max(s_est_non_w757_w800_av,[],2)), [0.5, 1]);
    xlabel('$z$ [voxel]', 'Interpreter', 'latex');
    ylabel('$x$ [voxel]', 'Interpreter', 'latex');
    colorbar('TickLabelInterpreter', 'latex'); colormap jet; colorbar; axis image;
    set(gca, 'FontSize', 14, 'FontName', 'Times', ...
        'TickLabelInterpreter', 'latex');
    subplot(1, 3, 3);
    imagesc(squeeze(max(s_est_non_w757_w800_av,[],3)), [0.5, 1]);
    xlabel('$y$ [voxel]', 'Interpreter', 'latex');
    ylabel('$x$ [voxel]', 'Interpreter', 'latex');
    colorbar('TickLabelInterpreter', 'latex'); colormap jet; colorbar; axis image;
    set(gca, 'FontSize', 14, 'FontName', 'Times', ...
        'TickLabelInterpreter', 'latex');
    sgtitle({'Classified arteries and veins, No optical fluence normalization', ...
        '757 nm and 800 nm'}, 'Interpreter', 'latex');
    
    figure;
    subplot(1, 3, 1);
    imagesc(squeeze(max(s_est_non_w757_w850_av,[],1)), [0.5, 1]);
    xlabel('$z$ [voxel]', 'Interpreter', 'latex');
    ylabel('$y$ [voxel]', 'Interpreter', 'latex');
    colorbar('TickLabelInterpreter', 'latex'); colormap jet; colorbar; axis image;
    set(gca, 'FontSize', 14, 'FontName', 'Times', ...
        'TickLabelInterpreter', 'latex');
    subplot(1, 3, 2);
    imagesc(squeeze(max(s_est_non_w757_w850_av,[],2)), [0.5, 1]);
    xlabel('$z$ [voxel]', 'Interpreter', 'latex');
    ylabel('$x$ [voxel]', 'Interpreter', 'latex');
    colorbar('TickLabelInterpreter', 'latex'); colormap jet; colorbar; axis image;
    set(gca, 'FontSize', 14, 'FontName', 'Times', ...
        'TickLabelInterpreter', 'latex');
    subplot(1, 3, 3);
    imagesc(squeeze(max(s_est_non_w757_w850_av,[],3)), [0.5, 1]);
    xlabel('$y$ [voxel]', 'Interpreter', 'latex');
    ylabel('$x$ [voxel]', 'Interpreter', 'latex');
    colorbar('TickLabelInterpreter', 'latex'); colormap jet; colorbar; axis image;
    set(gca, 'FontSize', 14, 'FontName', 'Times', ...
        'TickLabelInterpreter', 'latex');
    sgtitle({'Classified arteries and veins, No optical fluence normalization', ...
        '757 nm and 850 nm'}, 'Interpreter', 'latex');
    
    figure;
    subplot(1, 3, 1);
    imagesc(squeeze(max(s_est_non_w800_w850_av,[],1)), [0.5, 1]);
    xlabel('$z$ [voxel]', 'Interpreter', 'latex');
    ylabel('$y$ [voxel]', 'Interpreter', 'latex');
    colorbar('TickLabelInterpreter', 'latex'); colormap jet; colorbar; axis image;
    set(gca, 'FontSize', 14, 'FontName', 'Times', ...
        'TickLabelInterpreter', 'latex');
    subplot(1, 3, 2);
    imagesc(squeeze(max(s_est_non_w800_w850_av,[],2)), [0.5, 1]);
    xlabel('$z$ [voxel]', 'Interpreter', 'latex');
    ylabel('$x$ [voxel]', 'Interpreter', 'latex');
    colorbar('TickLabelInterpreter', 'latex'); colormap jet; colorbar; axis image;
    set(gca, 'FontSize', 14, 'FontName', 'Times', ...
        'TickLabelInterpreter', 'latex');
    subplot(1, 3, 3);
    imagesc(squeeze(max(s_est_non_w800_w850_av,[],3)), [0.5, 1]);
    xlabel('$y$ [voxel]', 'Interpreter', 'latex');
    ylabel('$x$ [voxel]', 'Interpreter', 'latex');
    colorbar('TickLabelInterpreter', 'latex'); colormap jet; colorbar; axis image;
    set(gca, 'FontSize', 14, 'FontName', 'Times', ...
        'TickLabelInterpreter', 'latex');
    sgtitle({'Classified arteries and veins, No optical fluence normalization', ...
        '800 nm and 850 nm'}, 'Interpreter', 'latex');
    
    figure;
    subplot(1, 3, 1);
    imagesc(squeeze(max(s_est_non_w757_w800_w850_av,[],1)), [0.5, 1]);
    xlabel('$z$ [voxel]', 'Interpreter', 'latex');
    ylabel('$y$ [voxel]', 'Interpreter', 'latex');
    colorbar('TickLabelInterpreter', 'latex'); colormap jet; colorbar; axis image;
    set(gca, 'FontSize', 14, 'FontName', 'Times', ...
        'TickLabelInterpreter', 'latex');
    subplot(1, 3, 2);
    imagesc(squeeze(max(s_est_non_w757_w800_w850_av,[],2)), [0.5, 1]);
    xlabel('$z$ [voxel]', 'Interpreter', 'latex');
    ylabel('$x$ [voxel]', 'Interpreter', 'latex');
    colorbar('TickLabelInterpreter', 'latex'); colormap jet; colorbar; axis image;
    set(gca, 'FontSize', 14, 'FontName', 'Times', ...
        'TickLabelInterpreter', 'latex');
    subplot(1, 3, 3);
    imagesc(squeeze(max(s_est_non_w757_w800_w850_av,[],3)), [0.5, 1]);
    xlabel('$y$ [voxel]', 'Interpreter', 'latex');
    ylabel('$x$ [voxel]', 'Interpreter', 'latex');
    colorbar('TickLabelInterpreter', 'latex'); colormap jet; colorbar; axis image;
    set(gca, 'FontSize', 14, 'FontName', 'Times', ...
        'TickLabelInterpreter', 'latex');
    sgtitle({'Classified arteries and veins, No optical fluence normalization', ...
        '757 nm, 800 nm, and 850 nm'}, 'Interpreter', 'latex');
    
    
    % Optical fluence normalization
    figure;
    subplot(1, 3, 1);
    imagesc(squeeze(max(s_est_on_w757_w800_av,[],1)), [0.5, 1]);
    xlabel('$z$ [voxel]', 'Interpreter', 'latex');
    ylabel('$y$ [voxel]', 'Interpreter', 'latex');
    colorbar('TickLabelInterpreter', 'latex'); colormap jet; colorbar; axis image;
    set(gca, 'FontSize', 14, 'FontName', 'Times', ...
        'TickLabelInterpreter', 'latex');
    subplot(1, 3, 2);
    imagesc(squeeze(max(s_est_on_w757_w800_av,[],2)), [0.5, 1]);
    xlabel('$z$ [voxel]', 'Interpreter', 'latex');
    ylabel('$x$ [voxel]', 'Interpreter', 'latex');
    colorbar('TickLabelInterpreter', 'latex'); colormap jet; colorbar; axis image;
    set(gca, 'FontSize', 14, 'FontName', 'Times', ...
        'TickLabelInterpreter', 'latex');
    subplot(1, 3, 3);
    imagesc(squeeze(max(s_est_on_w757_w800_av,[],3)), [0.5, 1]);
    xlabel('$y$ [voxel]', 'Interpreter', 'latex');
    ylabel('$x$ [voxel]', 'Interpreter', 'latex');
    colorbar('TickLabelInterpreter', 'latex'); colormap jet; colorbar; axis image;
    set(gca, 'FontSize', 14, 'FontName', 'Times', ...
        'TickLabelInterpreter', 'latex');
    sgtitle({'Classified arteries and veins, Optical fluence normalization', ...
        '757 nm and 800 nm'}, 'Interpreter', 'latex');
    
    figure;
    subplot(1, 3, 1);
    imagesc(squeeze(max(s_est_on_w757_w850_av,[],1)), [0.5, 1]);
    xlabel('$z$ [voxel]', 'Interpreter', 'latex');
    ylabel('$y$ [voxel]', 'Interpreter', 'latex');
    colorbar('TickLabelInterpreter', 'latex'); colormap jet; colorbar; axis image;
    set(gca, 'FontSize', 14, 'FontName', 'Times', ...
        'TickLabelInterpreter', 'latex');
    subplot(1, 3, 2);
    imagesc(squeeze(max(s_est_on_w757_w850_av,[],2)), [0.5, 1]);
    xlabel('$z$ [voxel]', 'Interpreter', 'latex');
    ylabel('$x$ [voxel]', 'Interpreter', 'latex');
    colorbar('TickLabelInterpreter', 'latex'); colormap jet; colorbar; axis image;
    set(gca, 'FontSize', 14, 'FontName', 'Times', ...
        'TickLabelInterpreter', 'latex');
    subplot(1, 3, 3);
    imagesc(squeeze(max(s_est_on_w757_w850_av,[],3)), [0.5, 1]);
    xlabel('$y$ [voxel]', 'Interpreter', 'latex');
    ylabel('$x$ [voxel]', 'Interpreter', 'latex');
    colorbar('TickLabelInterpreter', 'latex'); colormap jet; colorbar; axis image;
    set(gca, 'FontSize', 14, 'FontName', 'Times', ...
        'TickLabelInterpreter', 'latex');
    sgtitle({'Classified arteries and veins, Optical fluence normalization', ...
        '757 nm and 850 nm'}, 'Interpreter', 'latex');
    
    figure;
    subplot(1, 3, 1);
    imagesc(squeeze(max(s_est_on_w800_w850_av,[],1)), [0.5, 1]);
    xlabel('$z$ [voxel]', 'Interpreter', 'latex');
    ylabel('$y$ [voxel]', 'Interpreter', 'latex');
    colorbar('TickLabelInterpreter', 'latex'); colormap jet; colorbar; axis image;
    set(gca, 'FontSize', 14, 'FontName', 'Times', ...
        'TickLabelInterpreter', 'latex');
    subplot(1, 3, 2);
    imagesc(squeeze(max(s_est_on_w800_w850_av,[],2)), [0.5, 1]);
    xlabel('$z$ [voxel]', 'Interpreter', 'latex');
    ylabel('$x$ [voxel]', 'Interpreter', 'latex');
    colorbar('TickLabelInterpreter', 'latex'); colormap jet; colorbar; axis image;
    set(gca, 'FontSize', 14, 'FontName', 'Times', ...
        'TickLabelInterpreter', 'latex');
    subplot(1, 3, 3);
    imagesc(squeeze(max(s_est_on_w800_w850_av,[],3)), [0.5, 1]);
    xlabel('$y$ [voxel]', 'Interpreter', 'latex');
    ylabel('$x$ [voxel]', 'Interpreter', 'latex');
    colorbar('TickLabelInterpreter', 'latex'); colormap jet; colorbar; axis image;
    set(gca, 'FontSize', 14, 'FontName', 'Times', ...
        'TickLabelInterpreter', 'latex');
    sgtitle({'Classified arteries and veins, Optical fluence normalization', ...
        '800 nm and 850 nm'}, 'Interpreter', 'latex');
    
    figure;
    subplot(1, 3, 1);
    imagesc(squeeze(max(s_est_on_w757_w800_w850_av,[],1)), [0.5, 1]);
    xlabel('$z$ [voxel]', 'Interpreter', 'latex');
    ylabel('$y$ [voxel]', 'Interpreter', 'latex');
    colorbar('TickLabelInterpreter', 'latex'); colormap jet; colorbar; axis image;
    set(gca, 'FontSize', 14, 'FontName', 'Times', ...
        'TickLabelInterpreter', 'latex');
    subplot(1, 3, 2);
    imagesc(squeeze(max(s_est_on_w757_w800_w850_av,[],2)), [0.5, 1]);
    xlabel('$z$ [voxel]', 'Interpreter', 'latex');
    ylabel('$x$ [voxel]', 'Interpreter', 'latex');
    colorbar('TickLabelInterpreter', 'latex'); colormap jet; colorbar; axis image;
    set(gca, 'FontSize', 14, 'FontName', 'Times', ...
        'TickLabelInterpreter', 'latex');
    subplot(1, 3, 3);
    imagesc(squeeze(max(s_est_on_w757_w800_w850_av,[],3)), [0.5, 1]);
    xlabel('$y$ [voxel]', 'Interpreter', 'latex');
    ylabel('$x$ [voxel]', 'Interpreter', 'latex');
    colorbar('TickLabelInterpreter', 'latex'); colormap jet; colorbar; axis image;
    set(gca, 'FontSize', 14, 'FontName', 'Times', ...
        'TickLabelInterpreter', 'latex');
    sgtitle({'Classified arteries and veins, Optical fluence normalization', ...
        '757 nm, 800 nm, and 850 nm'}, 'Interpreter', 'latex');
end