classdef CircularCurve < SmoothCurveClosed
    % Circle gamma(t) = center + radius*(cos 2*pi*t, sin 2*pi*t), t in
    % [0, 1], traversed counterclockwise from (center + radius*e_1).
    %
    %   curve = CircularCurve(center, radius, mapping="closest")
    %
    % center is a 2-vector (default [0.5; 0.5]) and radius a scalar
    % (default 0.1); both are kept as properties. mapping is the
    % SmoothCurve point-to-interface map ("closest", "vertical",
    % "normal").
    %
    % See also: SmoothCurveClosed
    properties
        center
        radius
    end

    methods
        function obj = CircularCurve(x0y0, r, options)
            arguments
                x0y0 double = [0.5; 0.5]
                r double = 0.1;
                options.mapping = "closest";
            end
            x_t = sprintf("%.17g + %.17g*cos(2*pi*t)", x0y0(1), r);
            y_t = sprintf("%.17g + %.17g*sin(2*pi*t)", x0y0(2), r);
            obj = obj@SmoothCurveClosed(x_t, y_t, [0 1], mapping=options.mapping)
            obj.center = x0y0;
            obj.radius = r;
        end
    end
end
