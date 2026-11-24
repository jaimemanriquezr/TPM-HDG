classdef WallEndedCurve < DeformedCurve
    % A DeformedCurve of an open reference whose ends slide along the
    % vertical walls x = walls(1) and x = walls(2) of the box: points on a
    % wall map to the end parameter on that wall, all others to their
    % closest point. A fluid mesh carved from the box then ends its
    % interface chain at the curve's ends even where Sigma_h meets the wall
    % obliquely, which the closest-point map alone does not give (a wall
    % vertex at height delta above an end of slope m maps to t ~ m delta).
    %
    %   curve = WallEndedCurve(reference, d, walls = [0 1], DeformedCurve options...)
    %
    % See also: DeformedCurve, lagrangian_carved_march
    properties (SetAccess = immutable)
        walls (1, 2) double = [0 1]
    end

    methods
        function obj = WallEndedCurve(reference, displacement, walls, varargin)
            obj = obj@DeformedCurve(reference, displacement, varargin{:});
            obj.walls = walls;
        end

        function v_t = map_points_to_interface(obj, v_x, varargin)
            v_t = map_points_to_interface@DeformedCurve(obj, v_x, varargin{:});
            on_wall = abs(v_x(1, :, :) - reshape(obj.walls, 1, 1, 1, 2)) < 1e-12;
            v_t(any(on_wall(:, :, :, 1), 1)) = obj.t_range(1);
            v_t(any(on_wall(:, :, :, 2), 1)) = obj.t_range(2);
        end
    end
end
