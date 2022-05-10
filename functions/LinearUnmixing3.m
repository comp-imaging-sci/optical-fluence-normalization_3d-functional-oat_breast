function [s, c_thb] = LinearUnmixing3(epsilon, img1, img2, img3)
% [S, C_THB] = LINEARUNMIXING3(EPSILON, IMG1, IMG2, IMG3) computes oxygen
% saturation S and total hemoglobin concentration C_THB from three images
% IMG1, IMG2, and IMG3 at different wavelengths using spectral linear
% unmixing.
%
% Input:
%   EPSILON:    3 x 2  molar extinction coefficient matrix of deoxy
%               [EPSILON(:, 1)] and oxy [EPSILON(:, 2)] -hemoglobin at 
%               wavelength 1 [EPSILON(1, :)], wavelength 2 [EPSILON(2, :)],
%               and [EPSILON(3, :)]
%   IMG1:       image at wavelength 1
%   IMG2:       image at wavelength 2
%   IMG3:       image at wavelength 3
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

    coeff = pinv(epsilon);
    c_hb    = coeff(1,1).*img1 + coeff(1,2).*img2 + coeff(1,3).*img3;
    c_hbo2  = coeff(2,1).*img1 + coeff(2,2).*img2 + coeff(2,3).*img3;
        
    c_thb   = c_hb + c_hbo2;
    s       = c_hbo2./c_thb;
    s(s < 0) = 0;
    s(s > 1) = 1;
end