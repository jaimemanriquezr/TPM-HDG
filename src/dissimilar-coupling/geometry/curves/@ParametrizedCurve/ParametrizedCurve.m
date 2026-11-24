classdef ParametrizedCurve < AbstractCurve
    % Curve gamma(t) = (x_coord(t), y_coord(t)), t in t_range, with the
    % arc-length function arc_length(t1, t2) and the tangent gamma'(t).
    % Provides find_closest_point (global search on a cached sampling of
    % the curve, refined by fminbnd), find_normal_intersection (the same
    % sampling, along given lines), unit_normal, and the mesh-to-curve
    % machinery of AbstractCurve.
    %
    %   curve = ParametrizedCurve(x_coord, y_coord, arc_len, t_range)
    %
    % x_coord, y_coord  function handles of t (defaults t and 0.5).
    % arc_len           handle @(t1, t2) with the arc length between two
    %                   parameters (default abs(t2 - t1)).
    % t_range           parameter interval (default [0 1]).
    %
    % is_periodic is false; closed subclasses override it.
    %
    % See also: FlatCurve, SmoothCurve, LevelSetCurve
    properties
        t_range double
        x_coord function_handle
        y_coord function_handle
        arc_length function_handle
        tangent function_handle = @(t) [t.^0; 0*t]   % gamma'(t); default matches (t, const)
    end

    properties (Hidden)
        % Dense sampling of the curve, cached by the constructor and used
        % by find_closest_point to pick the right basin before refining.
        sample_t double = double.empty(1, 0)
        sample_points double = double.empty(2, 0)
    end

    methods
        function obj = ParametrizedCurve(x_coord, y_coord, arc_len, t_range)
            arguments
                x_coord function_handle = @(t) t
                y_coord function_handle = @(t) 0.5 * t.^0
                arc_len function_handle = @(t1, t2) abs(t2 - t1)
                t_range double = [0 1]
            end
            obj.t_range = t_range;
            obj.x_coord = x_coord;
            obj.y_coord = y_coord;
            obj.arc_length = arc_len;

            % Sampling for find_closest_point's global step, cached since
            % the query runs point by point. The closing point of a
            % periodic curve is detected by coincidence rather than by
            % is_periodic, to avoid subclass dispatch in the constructor.
            t = linspace(t_range(1), t_range(2), 8001);
            P = [x_coord(t); y_coord(t)];
            if norm(P(:, end) - P(:, 1)) <= 1e-12 * max(1, norm(P(:, 1)))
                t(end) = [];
                P(:, end) = [];
            end
            obj.sample_t = t;
            obj.sample_points = P;
        end

        function t = find_closest_point(obj, x0, y0)
            % Parameter of the globally closest point of the curve. The
            % distance is multimodal wherever the curve passes near
            % itself, so the cached sampling picks the basin and fminbnd
            % refines inside that bracket only.
            d2 = (obj.sample_points(1, :) - x0).^2 ...
               + (obj.sample_points(2, :) - y0).^2;
            [~, nearest] = min(d2);
            dt = obj.sample_t(2) - obj.sample_t(1);
            lo = obj.sample_t(nearest) - dt;
            hi = obj.sample_t(nearest) + dt;
            period = diff(obj.t_range);
            if ~is_periodic(obj)
                % An open curve has nothing beyond its ends to search.
                lo = max(lo, obj.t_range(1));
                hi = min(hi, obj.t_range(2));
            end
            distance = @(t) (x0 - obj.x_coord(t)).^2 + (y0 - obj.y_coord(t)).^2;
            t = fminbnd(distance, lo, hi, ...
                optimset('TolX', 1e-12, 'MaxIter', 2e3, 'MaxFunEvals', 2e3));
            if is_periodic(obj)
                t = obj.t_range(1) + mod(t - obj.t_range(1), period);
            else
                t = min(max(t, obj.t_range(1)), obj.t_range(2));
            end
        end

        function t = find_normal_intersection(obj, X, N, reach)
            % Parameter where each line X(:, i) + r*N(:, i) meets the
            % curve, taking the intersection with the smallest |r|.
            % Crossings are located as sign changes of the cross product
            % (gamma(t) - x) x n over the cached sampling (on an open
            % curve extended by 10% of t_range past either end, where the
            % ray of an end vertex may land; t is then clamped) and refined
            % by fzero. Errors when no intersection has |r| <= reach(i).
            arguments
                obj ParametrizedCurve
                X (2, :) double
                N (2, :) double
                reach (1, :) double
            end
            T = obj.sample_t;
            P = obj.sample_points;
            t0 = obj.t_range(1);
            period = diff(obj.t_range);
            if is_periodic(obj)
                T = [T, t0 + period];
                P = [P, P(:, 1)];
            else
                pad_lo = linspace(t0 - 0.1*period, t0, 801);
                pad_hi = linspace(t0 + period, t0 + 1.1*period, 801);
                T = [pad_lo(1:end-1), T, pad_hi(2:end)];
                P = [obj.get_2d_points(pad_lo(1:end-1)), P, obj.get_2d_points(pad_hi(2:end))];
            end
            n_points = size(X, 2);
            reach = reach .* ones(1, n_points);
            t = zeros(1, n_points);
            for i = 1:n_points
                x = X(:, i);
                n = N(:, i);
                cross_n = @(s) (obj.x_coord(s) - x(1))*n(2) - (obj.y_coord(s) - x(2))*n(1);
                g = (P(1, :) - x(1))*n(2) - (P(2, :) - x(2))*n(1);
                crossings = find(g(1:end-1).*g(2:end) <= 0);
                [r_best, t_best] = deal(Inf, NaN);
                for j = crossings
                    % re-evaluated: a sample value at roundoff level may
                    % differ in sign from the pointwise one
                    g_ends = [cross_n(T(j)), cross_n(T(j + 1))];
                    if prod(sign(g_ends)) < 0
                        t_j = fzero(cross_n, T(j:j + 1));
                    else
                        [~, at] = min(abs(g_ends));
                        t_j = T(j - 1 + at);
                    end
                    r = (obj.get_2d_points(t_j) - x).' * n;
                    if abs(r) < abs(r_best)
                        [r_best, t_best] = deal(r, t_j);
                    end
                end
                if ~(abs(r_best) <= reach(i))
                    error("ParametrizedCurve:noNormalIntersection", ...
                        "The normal line through (%g, %g) meets the curve " + ...
                        "at no |r| <= %g (closest at |r| = %g).", ...
                        x(1), x(2), reach(i), abs(r_best));
                end
                t(i) = t_best;
            end
            if is_periodic(obj)
                t = t0 + mod(t - t0, period);
            else
                t = min(max(t, obj.t_range(1)), obj.t_range(2));
            end
        end

        plot(obj, ax, options)

        function tf = is_periodic(~)
            % Closed interfaces override this (see SmoothCurveClosed):
            % their parameter wraps, so a face may straddle the seam.
            tf = false;
        end

        function n = unit_normal(obj, t)
            % unit normal at parameter values, page-tolerant: the tangent
            % rotated by -90 degrees, matching the SmoothCurve sym_normal
            % convention (outward on the closed curves)
            arguments
                obj ParametrizedCurve
                t (1, :, :) double
            end
            T = obj.tangent(reshape(t, 1, []));
            N = [T(2, :); -T(1, :)] ./ vecnorm(T);
            sz = size(t);
            n = reshape(N, [2, sz(2:end)]);
        end

        function xy = local2physical(obj, local_coords)
            arguments
                obj ParametrizedCurve
                local_coords (1, :)
            end
            xy = [obj.x_coord(local_coords); obj.y_coord(local_coords)];
        end

        function area = get_partition_area(obj, partition)
            arguments
                obj ParametrizedCurve
                partition struct
            end
            % arc length on Sigma of each patch of an induced partition
            % (get_partitions), from its interface-parameter endpoints
            points = reshape(partition.vertices_interface, 2, []);
            n_patches = size(points, 2);
            area = zeros(1, n_patches);

            get_length = @(p) obj.arc_length(p(1), p(2));
            for e = 1:n_patches
                area(e) = get_length(points(:, e));
            end
        end

        function p = get_2d_points(obj, t)
            arguments
                obj ParametrizedCurve
                t (1, :, :) double
            end
            p = cat(1, obj.x_coord(t), obj.y_coord(t));
        end

    end
end
