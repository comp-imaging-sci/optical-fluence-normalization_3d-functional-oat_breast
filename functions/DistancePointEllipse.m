function dist = DistancePointEllipse(x, y, f_ellp)
% DIST = DISTANCEPOINTELLIPSE(X, Y, F_ELLP) calculates a distnace DIST from
% a point at X, Y to an elliptical curve, that is defined by a fit object
% F_ELLP (elliptical curve-fitting result).
%
% Input:
%   X:          x-coordinates, a 3D array [voxel]
%   Y:          y-coordinates, a 3D array [voxel]
%   F_ELLP:     fit object, i.e., elliptical curve-fitting result
%
% Output: 
%   DIST:       distance from a point at (X, Y) to an elliptical curve
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

a = f_ellp.a;
b = f_ellp.b;
c = f_ellp.c;

[Nx, Ny] = size(x);
x = x(:) - c;
y = y(:);

n          = length(x);    % Number of points
x_proj     = zeros(n, 1);  % x-coordinate of projection points onto ellipse
y_proj     = zeros(n, 1);  % y-coordinate of projection points onto ellipse
tolerance  = 1e-9;         % Tolerance

% Circle
if abs((a - b)/a) < tolerance
    theta = angle(x + 1i*y);
    x_proj = a*cos(theta);
    y_proj = b*sin(theta);

% Ellipse
else
    aa = a^2;
    bb = b^2;
    tolerance_a  = tolerance*a;
    tolerance_b  = tolerance*b;
    tolerance_aa = tolerance*aa;

    x_abs = abs(x);
    y_abs = abs(y);
    
    t_init = max(a*(x_abs - a), b*(y_abs - b));

    for i = 1:n
        u = x_abs(i);
        v = y_abs(i);
        if u == 0
            z1 = 1;
        else
            z1 = sign(x(i));
        end
        if v == 0
            z2 = 1;
        else
            z2 = sign(y(i));
        end

        % Point on the minor axis
        if u < tolerance_a
            if y(i) < 0
                x_proj(i) = 0;
                y_proj(i) = -b;
            else
                x_proj(i) = 0;
                y_proj(i) = b;
            end
            continue
        end

        % Point on the major axis
        if v < tolerance_b
            if u < a - bb/a
                x_proj(i) = z1*aa*u/(aa - bb);
                y_proj(i) = z2*b*sqrt(1 - (aa*u/(aa - bb)/a)^2);
            else
                x_proj(i) = z1*a;
                y_proj(i) = 0;
            end
            continue
        end

        % Point not on the axis
        % Newton's method
        t = t_init(i);
        for iter = 1:100
            t_aa = t + aa;
            t_bb = t + bb;
            pp1 = (u*a/t_aa)^2;
            pp2 = (v*b/t_bb)^2;
            F  = pp1 + pp2 - 1;
            if F < 0
                break
            end
            F_p = 2*(pp1/t_aa + pp2/t_bb);
            ratio = F/F_p;
            if ratio < tolerance_aa
                break
            end
            t = t + ratio;
        end
        x_proj(i) = x(i)*aa/t_aa;
        y_proj(i) = sign(y(i))*b*sqrt(1 - (x_proj(i)/a)^2);
    end
end
dist = reshape(sqrt((x_proj - x).^2 + (y_proj - y).^2), [Nx, Ny]);
