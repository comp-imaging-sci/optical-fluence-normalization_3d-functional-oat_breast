function [s, c_thb] = LinearUnmixing2(epsilon, img1, img2)
% [S, C_THB] = LINEARUNMIXING2(EPSILON, IMG1, IMG2) computes oxygen
% saturation S and total hemoglobin concentration C_THB from two images
% IMG1 and IMG2 at different wavelengths using spectral linear unmixing.
%
% Input:
%   EPSILON:    2 x 2 molar extinction coefficient matrix of deoxy
%               [EPSILON(:, 1)] and oxy [EPSILON(:, 2)] -hemoglobin at 
%               wavelength 1 [EPSILON(1, :)] and wavelength 2 
%               [EPSILON(2, :)]
%   IMG1:       image at wavelength 1
%   IMG2:       image at wavelength 2
%
% Output:
%   S:          oxygen saturation
%   C_THB:      total hemoglobin concentration
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

    coeff = inv(epsilon);
    c_hb    = coeff(1,1).*img1 + coeff(1,2).*img2;
    c_hbo2  = coeff(2,1).*img1 + coeff(2,2).*img2;

    c_thb   = c_hb + c_hbo2;
    s       = c_hbo2./c_thb;
    s(s < 0) = 0;
    s(s > 1) = 1;
end