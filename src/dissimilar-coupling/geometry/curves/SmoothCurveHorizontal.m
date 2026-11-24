classdef SmoothCurveHorizontal < SmoothCurve
    % Graph y = y(t) of a symbolic expression over x = t, with the upward
    % unit normal.
    %
    %   curve = SmoothCurveHorizontal(y_t, t_range, mapping="vertical")
    %
    % y_t        string in the symbol t (default "0.5 + -.1*sin(2*pi*t)").
    % t_range    parameter interval (default [0 1]).
    % mapping    rule of map_points_to_interface: "vertical" (t = x, the
    %            default), "closest" (closest point of the curve) or
    %            "normal" (along the given normals, as SmoothCurve).
    %
    % See also: SmoothCurve, FlatCurveHorizontal
    methods
        function obj = SmoothCurveHorizontal(y_t, t_range, options)
            arguments
                y_t string = "0.5 + -.1*sin(2*pi*t)";
                t_range double = [0 1]
                options.mapping = "vertical";
            end
            obj = obj@SmoothCurve("t", y_t, t_range)

            sym_expr = str2sym(y_t);

            dy_fun = matlabFunction(diff(sym_expr), "Vars", sym("t"));
            obj.arc_length = @(t1, t2) integral(@(t) sqrt(1 + dy_fun(t).^2), t1, t2);
            obj.tangent = @(t) [t.^0; t.^0 .* dy_fun(t)];
            obj.sym_expr = sym_expr;
            obj.sym_normal = [-diff(sym_expr); 1] / sqrt(1 + diff(sym_expr)^2);

            syms t
            obj.dist_grad = @(x0, y0) matlabFunction((t - x0) + (sym_expr - y0)*diff(sym_expr));
            obj.mapping = options.mapping;
        end

        function v_t = map_points_to_interface(obj, v_x, options)
            arguments
                obj SmoothCurveHorizontal
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
                    % Bracket from the curve's own parameter range, padded by
                    % 10% of its length.
                    L = obj.t_range(2) - obj.t_range(1);
                    bracket = [obj.t_range(1) - 0.1*L, obj.t_range(2) + 0.1*L];
                    for i = 1:numel(v_t)
                        v_tx = fzero(obj.dist_grad(x0(i), y0(i)), bracket);
                        v_t(i) = max(obj.t_range(1), min(v_tx, obj.t_range(2)));
                    end
                case "vertical"
                    v_t = v_x(1, :, :);
                case "normal"
                    v_t = map_points_to_interface@SmoothCurve(obj, v_x, ...
                        normals = options.normals, reach = options.reach);
            end

            assert(obj.t_range(1) <= min(v_t(:)) && max(v_t(:)) <= obj.t_range(2));

            switch options.coordinates
                case "local"
                    return
                case "physical"
                    v_t = obj.get_2d_points(reshape(v_t, [1 size(v_x, [2 3])]));
            end
        end

        function n = unit_normal(obj, t)
            % Upward normal, the tangent rotated by +90 degrees:
            % [-y'; 1]/sqrt(1 + y'^2). Matches sym_normal,
            % FlatCurveHorizontal.unit_normal and the n_S of the
            % manufactured solutions; the inherited
            % ParametrizedCurve.unit_normal (outward convention of closed
            % curves) would point downward.
            arguments
                obj SmoothCurveHorizontal
                t (1, :, :) double
            end
            T = obj.tangent(reshape(t, 1, []));
            N = [-T(2, :); T(1, :)] ./ vecnorm(T);
            sz = size(t);
            n = reshape(N, [2, sz(2:end)]);
        end

    end
end
