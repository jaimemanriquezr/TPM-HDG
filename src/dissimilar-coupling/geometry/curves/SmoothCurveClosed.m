classdef SmoothCurveClosed < SmoothCurve
    % SmoothCurve whose parameter wraps around (is_periodic is true), so
    % a mesh face may straddle the seam at the ends of t_range.
    %
    %   curve = SmoothCurveClosed(x_t, y_t, t_range, mapping="closest")
    %
    % Arguments are those of SmoothCurve; the defaults describe a closed
    % curve about (0.5, 0.5) on t in [0, 1].
    %
    % See also: SmoothCurve, CircularCurve, AirfoilCurve
    methods
        function obj = SmoothCurveClosed(x_t, y_t, t_range, options)
            arguments
                x_t string = "0.5 + .2*(1 + .5*sin(2*pi*t))*cos(2*pi*t)";
                y_t string = "0.5 + .4*(1 + .5*cos(2*pi*t))*sin(2*pi*t)";
                t_range double = [0 1]
                options.mapping = "closest";
            end
            obj = obj@SmoothCurve(x_t, y_t, t_range, mapping=options.mapping)
        end

        function tf = is_periodic(~)
            % The parameter wraps: a face may straddle the seam, and the
            % generic map_mesh_to_interface splits it there
            tf = true;
        end

        % find_closest_point is inherited from ParametrizedCurve, which
        % brackets globally against a sampling of the whole period and so
        % also covers the seam.
    end
end
