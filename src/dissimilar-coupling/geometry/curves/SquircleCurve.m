classdef SquircleCurve < SmoothCurveClosed
%SQUIRCLECURVE Rounded square: r(theta) = R*(1 - eps*cos(4*theta)).
%
%   A closed, convex, smooth curve with four-fold symmetry. r is minimal on
%   the axes, so the flats face the axes and the high-curvature corners sit
%   on the diagonals.
%
%   eps concentrates curvature; it does not approach a true square. At a
%   flat's midpoint the polar curvature numerator r^2 + 2r'^2 - r r'' is
%   R^2 (1 - eps)(1 - 17 eps), so the curve loses convexity at eps = 1/17.
%   find_closest_point is unreliable on non-convex closed curves, so eps
%   must stay below that, which caps the radius ratio at about 1.12
%   (sqrt(2) for a true square).
%
%   center, radius and eps are stored, as in CircularCurve, so that
%   FluidStructureCoupling.get_manufactured_solution can build the normal
%   as a function of (x, y); a normal in the curve parameter would collide
%   with the time symbol.
%
%   See also: CircularCurve, SmoothCurveClosed

    properties
        center
        radius
        eccentricity   % the eps of r = R(1 - eps*cos(4 theta))
    end

    methods
        function obj = SquircleCurve(x0y0, R, eps_squircle, options)
            arguments
                x0y0 double = [0.5; 0.5]
                R double = 0.4
                eps_squircle double = 0.04
                options.mapping = "closest";
            end
            assert(eps_squircle < 1/17, "SquircleCurve:nonConvex", ...
                "eps = %.4f is at or past the convexity wall 1/17 = %.4f.", ...
                eps_squircle, 1/17);
            r_t = sprintf("%.15g*(1 - %.15g*cos(8*pi*t))", R, eps_squircle);
            obj = obj@SmoothCurveClosed( ...
                sprintf("%.15g + %s*cos(2*pi*t)", x0y0(1), r_t), ...
                sprintf("%.15g + %s*sin(2*pi*t)", x0y0(2), r_t), ...
                [0 1], mapping = options.mapping);
            obj.center = x0y0;
            obj.radius = R;
            obj.eccentricity = eps_squircle;
        end

        function A = enclosed_area(obj)
            % Exact: A = 1/2 int r^2 dtheta = pi R^2 (1 + eps^2/2).
            A = pi * obj.radius^2 * (1 + obj.eccentricity^2 / 2);
        end
    end
end
