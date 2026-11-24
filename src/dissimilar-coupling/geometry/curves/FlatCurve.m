classdef FlatCurve < ParametrizedCurve
    % Straight segment gamma(t) = origin + t*direction, t in t_range, with
    % a unit direction. The closest-point map is the orthogonal projection,
    % the tangent is constant, and the unit normal is the tangent rotated
    % by +90 degrees (upward for a left-to-right horizontal line).
    %
    % See also: FlatCurveHorizontal, FlatCurveVertical
    properties (SetAccess = immutable)
        origin (2,1) double
        direction (2,1) double
    end

    methods
        function obj = FlatCurve(origin, direction, t_range)
            arguments
                origin (2,1) double = [0; 0.5]
                direction (2,1) double = [1; 0]
                t_range (1,2) double = [0 1]
            end
            d = direction / norm(direction);
            obj = obj@ParametrizedCurve(@(t) origin(1) + d(1)*t, ...
                                        @(t) origin(2) + d(2)*t, ...
                                        @(t1, t2) abs(t2 - t1), t_range);
            obj.origin = origin;
            obj.direction = d;
            obj.tangent = @(t) d .* t.^0;
        end

        function v_t = map_points_to_interface(obj, v_x, options)
            arguments
                obj FlatCurve
                v_x (2, :, :) double
                options.coordinates = "local";
                options.normals = [];
            end
            % orthogonal projection: t = d . (x - origin), clamped to the
            % segment (beyond an end the closest point is that end)
            v_t = sum(obj.direction .* (v_x - obj.origin), 1);
            v_t = min(max(v_t, obj.t_range(1)), obj.t_range(2));

            switch options.coordinates
                case "local"
                    return
                case "physical"
                    v_t = obj.get_2d_points(v_t);
            end
        end

        function n = unit_normal(obj, t)
            % the direction rotated by +90 degrees, page-tolerant
            n = ([0 -1; 1 0] * obj.direction) .* t.^0;
        end
    end
end
