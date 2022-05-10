% -------------------------------------------------------------------------
%                       Optical Fluence Normalization
% -------------------------------------------------------------------------
% Author:   Seonyeong Park (sp33@illinois.edu)
% Date:     Jan. 24, 2022
%
% This project is to estimate and compensate for the non-uniform optical 
% fluence distribution in three-dimensional (3D) optoacoustic tomography 
% (OAT) of the breast, based on reasonable assumptions regarding breast 
% anatomy and optical properties. The non-uniform incident optical fluence
% is estimated based on the illumination geometry in the OAT system, and 
% the depth-dependent optical attenuation is approximated using the
% Beer-Lambert law.
% The current version of the software is based on the reference imaging
% system, LOUISA-3D of TomoWave Laboratories (Houston, Texas). The method
% of the software is detailed in [Park2022].
%
% Reference:
%   [Park2022] Seonyeong Park, Frank J. Brooks, Umberto Villa, Richard Su,
%           Mark A. Anastasio, Alexander A. Oraevsky, "Normalization of
%           optical fluence distribution for three-dimensional functional
%           optoacoustic tomography of the breast," J. Biomed. Opt. 27(3)
%           036001 (16 March 2022)
%           https://doi.org/10.1117/1.JBO.27.3.036001
%
% Copyright (C) 2022 Seonyeong Park, 
% License:  GNU General Public License version 3, Please see 'LICENSE' for
%           details.
%

addpath('functions');   % Add a path of functions

lambs = [757, 800, 850]; % Wavelangths [nm]

voxel_size = 0.25; % Voxel size [mm]

flag_fig  = true;  % Flag to plot
% flag_fig  = false; % Flag to plot

depth_vis = 30;   % Maximum depth for visualization [mm]
L = 2; % Polynomial degree for curve-fitting, Eq. (2) in Reference


for lamb_i = 1:length(lambs)
    lamb = lambs(lamb_i); % Wavelength
    fprintf('-------------------------------------------------------------------------\n');
    fprintf('Wavelength: %d nm\n', lamb);
    fprintf('-------------------------------------------------------------------------\n');

    % File name
    % Reconstructed image
    filename_recon  = fullfile('data', ['RECON_NOISY_w', ...
        num2str(lamb), '_FBP.mat']);
    
    % Breast mask
    filename_breast_mask = fullfile('data', ['mask_breast_w', ...
        num2str(lamb), '.mat']);

    % Depth map
    filename_depth = fullfile('data', ['depth_w', num2str(lamb), '.mat']);
    
    % Estimated incident optical fluence
    filename_phi0_est = fullfile('data', ['phi0_est_w', num2str(lamb), ...
        '.mat']);

    % Estimated optical attenuation
    filename_phia_est = fullfile('data', ['phia_est_w', num2str(lamb), ...
        '.mat']);
    
    
    % Load reconstructed image
    fprintf('Loading reconstructed image...\n');
    load(filename_recon, 'recon');

    % Normalize brightness
    recon(recon < 0) = 0; recon = recon./max(recon(:));

    % Get dimension size [voxel]
    [Nx, Ny, Nz] = size(recon); 

    if flag_fig == 1
        % Plot maximum voxel brightness projections of reconstructed image
        figure; imagesc(squeeze(max(recon, [], 2))); colormap hot;
        colorbar('TickLabelInterpreter', 'latex'); truesize;
        xlabel('$z$ [voxel]', 'Interpreter', 'latex');
        ylabel('$x$ [voxel]', 'Interpreter', 'latex');
        title('$\bf{\hat{\alpha}}$, MVBP along $y$-axis', ...
        'Interpreter', 'latex');
        set(gca, 'FontSize', 14, 'FontName', 'Times', ...
            'TickLabelInterpreter', 'latex');
    end


    %% Compute spherical coordinates of each voxel
    fprintf('Computing spherical coordinates of each voxel...\n');

    O = [Nx + 1, Ny + 1, 2*Nz + 1]./2; % Origin coordinates [voxel]

    % Compute polar angle (THETA) and radial distance from z-axis (RHO)
    [~, THETA, ~, RHO] = SphericalCoord([Nx, Ny, Nz], O); % [voxel]

    %% COMPENSATION FOR NON-UNIFORM INCIDENT OPTICAL FLUENCE
    fprintf('\nCOMPENSATION FOR NON-UNIFORM INCIDENT OPTICAL FLUENCE\n');
    fprintf('[Step 1] Maximum voxel brightness extraction at each polar angle\n');

    % Extracts maximum voxel brightness at each polar angle in the
    % reconstructed image
    [theta, max_vb_theta, ~] = MaximumVoxelBrightnessTheta(recon, THETA);


    fprintf('[Step 2] Non-uniform illumination estimation using polynomial curve-fitting\n');

    theta_arla = 160; % Polar angle of areola [degree]

    % Fit L-degree polynomial curve to maximum voxel brightness according to 
    % the discretized polar angles (the region containing the nipple and areola
    % is excluded.), Eq. (2) in Reference
    f_phi0 = polyfit( ...
        theta(theta <= theta_arla), max_vb_theta(theta <= theta_arla), L);

    if flag_fig == 1
        phi0_est_theta  = polyval(f_phi0, theta); % This is just for a plot

        % Plot maximum voxel brightness according to the polar angles
        figure;
        plot(theta, max_vb_theta, 'k-', 'LineWidth', 1.5); hold on;
        plot(theta, phi0_est_theta, 'r-', 'LineWidth', 2); hold on;
        plot([theta_arla, theta_arla], [0, 1], 'b--', 'LineWidth', 2); hold on;
        xlim([90, 180]); ylim([0, 1]);
        xlabel('Polar angle $\theta$ [degree]', 'Interpreter', 'latex');
        ylabel('Normalized brightness $\bf{\hat{\alpha}}$', ...
            'Interpreter', 'latex');
        lgnd1 = 'Maximum voxel brightness $\bf{\hat{\alpha}_{max}}$';
        lgnd2 = 'Estimated incident optical flunece $\bf{\hat{\phi}_0}$';
        lgnd3 = 'Maximum $\theta$ for curve-fitting';
        legend({lgnd1, lgnd2, lgnd3}, ...
            'Interpreter', 'latex', 'Location', 'northoutside');
        set(gca, 'FontSize', 14, 'FontName', 'Times', ...
            'TickLabelInterpreter', 'latex');

        clear phi0_est_theta lgnd1 lgnd2 lgnd3
    end

    % Estimated incident optical fluence, Eq. (2) in Reference
    phi0_est = polyval(f_phi0, THETA);
    phi0_est = phi0_est./max(phi0_est(:));

    % Save estimated incident optical fluence
    fprintf('    Saving estimated incident optical fluence...\n');
    save(filename_phi0_est, 'phi0_est', 'f_phi0');

    clear theta max_vb_theta


    fprintf('[Step 3] Compensation for non-uniform illumination\n');

    % Compensate for the non-uniform incident optical fluence
    reconN0 = recon./phi0_est;              % Eq. (3) in Reference
    reconN0 = reconN0./max(reconN0(:));     % Normalize brightness

    if flag_fig == 1
        % Plot maximum voxel brightness projections of reconstructed image
        % after compensate for the non-uniform incident optical fluence
        figure; imagesc(squeeze(max(reconN0, [], 2))); colormap hot;
        colorbar('TickLabelInterpreter', 'latex'); truesize;
        xlabel('$z$ [voxel]', 'Interpreter', 'latex');
        ylabel('$x$ [voxel]', 'Interpreter', 'latex');
        title('$\bf{\hat{\alpha_{N0}}}$, MVBP along $y$-axis', ...
            'Interpreter', 'latex');
        set(gca, 'FontSize', 14, 'FontName', 'Times', ...
            'TickLabelInterpreter', 'latex');
    end


    %% ESTIMATION OF BREAST SURFACE AND DEPTH
    fprintf('\nESTIMATION OF BREAST SURFACE AND DEPTH\n');
    fprintf('[Step 1] Extraction of voxels near breast surface\n');

    % 3D median filter to reduce noise, Eq. (4) in Reference
    reconN0_p = medfilt3(reconN0);

    % Square operation to enhance contrast, Eq. (4) in Reference
    reconN0_p = reconN0_p.^2;

    % Binarize using Otsu's threshold
    reconN0_p = imbinarize(reconN0_p);


    fprintf('[Step 2] Estimation of breast surface using elliptical curve-fitting\n');

    % Estimate radius on each z-slice, i.e., x-y plane, Eq. (5) in Reference
    rho = zeros(Nz, 1);
    for zi = 1:Nz
        if sum(reconN0_p(:, :, zi), 'all')
            [xn, yn] = find(reconN0_p(:, :, zi));
            rho(zi) = max(sqrt((xn - O(1)).^2 + (yn - O(2)).^2));
        end
    end

    clear reconN0_p xn yn

    % Squared elliptical equation for curve fitting to prevent complex values
    % (z - z_C)^2/a^2 + rho^2/b^2 = 1
    % rho^2 = (b/a)^2*(a^2 - (z - z_C)^2)
    zn  = (1:Nz)';
    zn  = zn(rho > 0);  % Exclude z-slices where the breast does not appear
    rho = rho(rho > 0); % Exclude z-slices where the breast does not appear

    ellp = fittype('(b/a)^2*(a^2 - (x - c)^2)', ...
        'independent', 'x', 'dependent', 'y');  % Ellipse equation 
    opts = fitoptions(ellp);
    opts.Lower      = [Nz - find(rho > 0, 1, 'first'), max(rho), Nz];
    opts.StartPoint = [Nz - find(rho > 0, 1, 'first'), max(rho), Nz];
    opts.Robust = 'LAR';
    f_ellp = fit(zn, rho.^2, ellp, opts);       % Elliptical curve fitting

    clear opts

    % Estimated radius on each z-slice, i.e., x-y plane after the elliptical
    % curve fitting
    z = (1:Nz)';
    rho_est = sqrt(f_ellp(z));
    rho_est = real(rho_est); % Exclude imaginary part

    if flag_fig == 1
        % Plot elliptical curve fitting result
        figure;
        plot(zn, rho, 'bo', 'LineWidth', 0.5, 'MarkerSize', 6); hold on;
        plot([z(find(rho_est > 0, 1, 'first')); z(rho_est > 0)], ...
            [-rho_est(find(rho_est > 0, 1, 'first')); rho_est(rho_est > 0)], ...
            '-', 'Color', 'r', 'LineWidth', 2.5); hold on;
        xlabel('$z$ [voxel]', 'Interpreter', 'latex');
        ylabel('$\rho$ [voxel]', 'Interpreter', 'latex');
        xlim([1, Nz]); ylim([0, Nx/2]);
        legend({'$\rho_z$', '$\hat{\rho_z}$'}, 'Interpreter', 'latex', ...
            'Location', 'northwest');
        set(gca, 'FontSize', 14, 'FontName', 'Times', ...
            'TickLabelInterpreter', 'latex', 'DataAspectRatio', [1, 1, 1]);
    end

    clear zn rho

    % Breat mask
    [~, ~, z] = meshgrid(1:Nx, 1:Ny, 1:Nz);
    mask = BreastMask(z, RHO, f_ellp);

    % Save breast mask
    fprintf('    Saving breast mask...\n');
    save(filename_breast_mask, 'mask');
    
    clear z


    fprintf('[Step 3] Estimation of depth\n');

    % Exclude unneccessary computation
    RHO_slc = squeeze(RHO(:,:,1));
    rho_roi = RHO_slc.*squeeze(mask(:,:,end));

    % Consider all discrete rho values in z-rho grid
    rho_uniq = unique(rho_roi(rho_roi > 0));
    [Z, RHO_uniq] = meshgrid(1:Nz, rho_uniq);

    clear rho_roi

    % Mask in the z-rho grid
    mask_zrho = zeros(size(RHO_uniq));
    for zi = 1:Nz, mask_zrho(RHO_uniq(:,zi) <= rho_est(zi), zi) = 1; end

    % Minimum distance from voxels to the elliptical curve
    dist_zrho = DistancePointEllipse(Z, RHO_uniq, f_ellp); % [voxel]
    dist_zrho = dist_zrho.*mask_zrho; % [voxel]

    clear Z RHO_uniq rho_est mask_zrho

    % Map z-rho grid to Cartesian grid for a depth map
    rho_idx = zeros(size(RHO_slc));
    for i = 1:length(rho_uniq), rho_idx(RHO_slc == rho_uniq(i)) = i; end

    depth = zeros(Nx, Ny, Nz); % Depth map [voxel]
    for zi = 1:Nz
        for xi = 1:Nx
            for yi = 1:Ny
                if rho_idx(xi, yi) > 0
                    depth(xi, yi, zi) = dist_zrho(rho_idx(xi, yi), zi);
                end
            end
        end
    end
    
    % Save depth map
    fprintf('    Saving voxel depth map...\n');
    save(filename_depth, 'depth');
    
    clear RHO_slc rho_uniq dist_zrho rho_idx xi yi zi


    %% COMPENSATION FOR THE EFFECTIVE OPTICAL ATTENUATION
    fprintf('\nCOMPENSATION FOR THE EFFECTIVE OPTICAL ATTENUATION\n');
    fprintf('[Step 1] Maximum voxel brightness extraction of depth\n');

    max_depth = round(max(depth(:))); % Maximum depth [voxel]

    % Maximum voxel brightness at each discretized depth in the reconstructed
    % image, Eq. (6) in Reference
    max_vb_depth = zeros(max_depth + 1, 1);  

    delta_d = 1; % Depth increment [voxel]

    for d = 0:max_depth
        recon_bin = recon(depth <= d + delta_d/2 & depth > d - delta_d/2);
        max_vb_depth(d + 1) = max(recon_bin);
    end
    clear delta_d recon_bin


    fprintf('[Step 2] Estimation of optical attenuation using Beer-Lambert law\n');

    % Beer-Lambert law, Eq. (7) in Reference
    beer_lambert_law = fittype('a*exp(-b*x)'); 
    d = double(0:max_depth)'; % [voxel]

    % Local maxima for exponential curve-fitting (the function 'islocalmax' is
    % available in MATLAB R2017b or later.)
    lm_dix = islocalmax(max_vb_depth, 'MinSeparation', 10); 
    f_phia = fit(d(lm_dix), max_vb_depth(lm_dix), beer_lambert_law, ...
        'Lower', [0, 0], 'Upper', [1, 1], 'StartPoint', [0, 0]);

    % Estimated  effective optical attenuation [cm^(-1)]
    mu_eff = f_phia.b / voxel_size * 10;
    fprintf('	Estimated mu_eff: %.4f [1/cm]\n', mu_eff);

    if flag_fig == 1
        phia_est_d = f_phia(d); % This is just for a plot

        % Plot maximum voxel brightness according to the depth
        figure;
        plot(d.*voxel_size, max_vb_depth, 'k-', 'LineWidth', 1.5); hold on;
        plot(d.*voxel_size, phia_est_d, 'b-', 'LineWidth', 2);
        xlim([0, max_depth] .* voxel_size); ylim([0, 1]);
        xlabel('Depth [mm]','Interpreter','latex');
        ylabel('Normalized brightness $\bf{\hat{\alpha}}$','Interpreter','latex');
        lgnd1 = 'Maximum voxel brightness $\bf{\hat{\alpha}_{BV}}$';
        lgnd2 = 'Estimated optical attenuation $\bf{\hat{\phi}_a}$';
        legend({lgnd1, lgnd2}, 'Interpreter','latex', 'Location', 'northeast');
        set(gca, 'FontSize', 14, 'FontName', 'Times', 'TickLabelInterpreter', 'latex');
        txt = text(10, 0.6, ['Estimated $\mu_{eff}$: ', num2str(mu_eff), ' cm$^{-1}$'], 'interpreter','latex');
        set(txt, 'FontSize', 14);

        clear phia_est_d lgnd1 lgnd2 txt
    end

    % Estimated optical attenuation
    phia_est = exp(-f_phia.b.*depth);
    
    % Save estimated optical attenuation
    fprintf('    Saving estimated optical attenuation...\n');
    save(filename_phia_est, 'phia_est', 'f_phia');
    
    doi = depth; % Depth of interest for visualization [voxel]
    doi(doi > depth_vis/voxel_size)= depth_vis/voxel_size; % [voxel]
    phia_est = exp(-f_phia.b.*doi); % Only for visualization purpose

    clear d max_vb_depth lm_dix


    fprintf('[Step 3] Compensation for optical attenuation\n\n');

    % Compensate for the optical attenuation
    reconN = reconN0./phia_est;             % Eq. (8) in Reference
    reconN = reconN./max(reconN(:));        % Normalize brightness

    if flag_fig == 1
        % Plot maximum voxel brightness projections of reconstructed image
        % after optical fluence normalization
        figure; imagesc(squeeze(max(reconN, [], 2))); colormap hot;
        colorbar('TickLabelInterpreter', 'latex'); truesize;
        xlabel('$z$ [voxel]', 'Interpreter', 'latex');
        ylabel('$x$ [voxel]', 'Interpreter', 'latex');
        title('$\bf{\hat{\alpha_{N}}}$, MVBP along $y$-axis', ...
            'Interpreter', 'latex');
        set(gca, 'FontSize', 14, 'FontName', 'Times', ...
            'TickLabelInterpreter', 'latex');
    end
end
