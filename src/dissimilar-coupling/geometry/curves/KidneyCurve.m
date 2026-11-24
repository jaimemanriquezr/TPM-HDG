classdef KidneyCurve < LevelSetCurve
    % Kidney-shaped zero contour of the quartic level-set function
    % phi = (2*r^2 - X)^2 - r^2 + offset, with X = x + center and
    % r^2 = X^2 + y^2.
    %
    %   curve = KidneyCurve(offset, center)
    %
    % offset (default 0.1) is the constant term of phi and center
    % (default 0.5) the shift of x in phi; both are kept as properties.
    % The bounding box is [-0.7 0.7 -0.8 0.8], which contains the default
    % curve; the curve must be discretized (discretize) before it is used
    % as an interface.
    %
    % See also: LevelSetCurve
    properties
        offset
        center
    end

    methods
        function obj = KidneyCurve(offset, center)
            arguments
                offset double = 0.1;
                center double = 0.5;
            end
            a = center;
            b = offset;

            phi = @(x, y) (2*((x+a).^2 + y.^2) - (x+a)).^2 - ((x+a).^2 + y.^2) + b;

            obj = obj@LevelSetCurve(@(X) phi(X(1, :), X(2, :)), ...
                box = [-0.7 0.7 -0.8 0.8]);   % contains the default kidney
            obj.offset = b;
            obj.center = a;
        end
    end
end
