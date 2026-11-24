classdef SmoothCurve < ParametrizedCurve
    % Interface given by symbolic coordinate expressions x(t), y(t),
    % with the tangent, sym_normal (the tangent rotated by -90 degrees)
    % and arc length derived from them.
    %
    %   curve = SmoothCurve(x_t, y_t, t_range, mapping="closest")
    %
    % x_t, y_t   strings in the symbol t (defaults "t" and
    %            "0.5 + .2*sin(2*pi*t)").
    % t_range    parameter interval (default [0 1]).
    % mapping    rule of map_points_to_interface: "closest" (closest
    %            point of the curve), "vertical" (t = x), or "normal"
    %            (along the given normals, one per point, to the curve
    %            point with the smallest |r| <= reach; map_mesh_to_interface
    %            passes the averaged face normals of the interface vertices).
    %
    % See also: ParametrizedCurve, SmoothCurveClosed
    properties
        sym_expr
        mapping
    end

    properties (Hidden)
        sym_normal
        dist_grad
    end
    methods
        function obj = SmoothCurve(x_t, y_t, t_range, options)
            arguments
                x_t string = "t";
                y_t string = "0.5 + .2*sin(2*pi*t)";
                t_range double = [0 1]
                options.mapping = "closest";
            end
            x_sym_expr = str2sym(x_t);
            x_fun = matlabFunction(x_sym_expr, "Vars", sym("t"));
            x_t = @(t) (t.^0) .* feval(x_fun, t);

            y_sym_expr = str2sym(y_t);
            y_fun = matlabFunction(y_sym_expr, "Vars", sym("t"));
            y_t = @(t) (t.^0) .* feval(y_fun, t);

            dx_fun = matlabFunction(diff(x_sym_expr), "Vars", sym("t"));
            dy_fun = matlabFunction(diff(y_sym_expr), "Vars", sym("t"));

            dt = @(t1, t2) integral(@(t) sqrt(dx_fun(t).^2 + dy_fun(t).^2), t1, t2);
            obj = obj@ParametrizedCurve(x_t, y_t, dt, t_range)
            obj.tangent = @(t) [t.^0 .* dx_fun(t); t.^0 .* dy_fun(t)];
            obj.sym_expr = [x_sym_expr; y_sym_expr];
            obj.sym_normal = -[-diff(y_sym_expr); diff(x_sym_expr)] / sqrt(diff(x_sym_expr)^2 + diff(y_sym_expr)^2);

            obj.dist_grad = @(x0, y0) matlabFunction((x_sym_expr - x0)*diff(x_sym_expr) + (y_sym_expr - y0)*diff(y_sym_expr));
            obj.mapping = options.mapping;
        end

        function v_t = map_points_to_interface(obj, v_x, options)
            arguments
                obj SmoothCurve
                v_x (2, :, :) double
                options.normals (2, :) double = double.empty(2, 0);
                options.reach (1, :) double = Inf;
                options.coordinates = "local";
            end
            switch obj.mapping
                case "closest"
                    v_t = zeros([1 size(v_x, [2 3])]);
                    x0 = v_x(1, :);
                    y0 = v_x(2, :);
                    for i = 1:numel(v_t)
                        v_tx = find_closest_point(obj, x0(i), y0(i));
                        v_t(i) = max(obj.t_range(1), min(v_tx, obj.t_range(2)));
                    end
                case "vertical"
                    v_t = v_x(1, :, :);
                case "normal"
                    n_points = numel(v_x)/2;
                    if size(options.normals, 2) ~= n_points
                        error("SmoothCurve:normalsRequired", ...
                            "mapping = ""normal"" needs one normal per point (%d), got %d.", ...
                            n_points, size(options.normals, 2));
                    end
                    v_t = find_normal_intersection(obj, reshape(v_x, 2, []), ...
                        options.normals, options.reach);
                    v_t = reshape(v_t, [1 size(v_x, [2 3])]);
            end
            assert(obj.t_range(1) <= min(v_t(:)) && max(v_t(:)) <= obj.t_range(2));

            switch options.coordinates
                case "local"
                    return
                case "physical"
                    v_t = obj.get_2d_points(reshape(v_t, [1 size(v_x, [2 3])]));
            end
        end
        % find_closest_point is inherited from ParametrizedCurve, which
        % searches globally: the distance to a curve is multimodal wherever
        % the curve passes near itself.
    end
end
