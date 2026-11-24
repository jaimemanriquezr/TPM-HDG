classdef LevelSetCurve < AbstractCurve
    % Implicitly defined curve: the zero contour of a level-set function.
    % Canonical calling convention is columnwise points, phi(X) with X
    % (2, n) -> (1, n); the constructor also accepts the cell form
    % phi({x; y}) and wraps it. Sign convention: phi < 0 is inside.
    %
    % projection_damping (in (0, 1], default 1) is the damping of every
    % Newton projection the curve makes when no damping is passed to
    % project explicitly: the TPM transfer paths (map_points_to_interface),
    % the canonical interface points (get_2d_points), the vertices of the
    % induced parametrization (discretize) and the single-mesh transfer
    % users of project. 1 is plain Newton; see project.
    properties
        level_set_function
        gradient            % (2, n) -> (2, n); empty = central differences
        box double = [-1 1 -1 1]   % bounding box [x1 x2 y1 y2] for contour extraction and plotting
        projection_damping (1,1) double {mustBePositive, ...
            mustBeLessThanOrEqual(projection_damping, 1)} = 1
    end

    properties (Hidden)
        % induced parametrization (set by discretize): ordered polyline of
        % the zero contour and its cumulative arc-length coordinate
        polyline double = double.empty(2, 0)
        arc double = double.empty(1, 0)
    end

    properties (Dependent)
        t_range
    end

    methods
        function obj = LevelSetCurve(phi, options)
            arguments
                phi function_handle = @(X) X(1, :).^2 + X(2, :).^2 - 0.5^2;
                options.gradient function_handle = function_handle.empty
                options.box (1, 4) double = [-1 1 -1 1]
                options.projection_damping (1,1) double = 1
            end
            % accept the cell form phi({x; y}) by probing, wrap to array form
            try
                phi([0; 0]);
            catch
                phi = @(X) phi({X(1, :); X(2, :)});
            end
            obj.level_set_function = phi;
            obj.gradient = options.gradient;
            obj.box = options.box;
            obj.projection_damping = options.projection_damping;
        end

        function z = evaluate(obj, X)
            % level-set values at columnwise points: (2, ...) -> (1, ...)
            % (trailing page dims preserved)
            sz = size(X);
            z = obj.level_set_function(reshape(X, 2, []));
            z = reshape(z, [1, sz(2:end)]);
        end

        function G = grad(obj, X)
            % gradient at columnwise points (2, n) -> (2, n); analytic when
            % declared, central differences otherwise
            if ~isempty(obj.gradient)
                G = obj.gradient(X);
            else
                h = 1e-6 * max(1, vecnorm(X));
                e1 = [1; 0]; e2 = [0; 1];
                G = [obj.evaluate(X + h.*e1) - obj.evaluate(X - h.*e1);
                     obj.evaluate(X + h.*e2) - obj.evaluate(X - h.*e2)] ./ (2*h);
            end
        end

        function N = normal_at(obj, X)
            % outward unit normal (phi < 0 inside) at columnwise points
            G = obj.grad(X);
            N = G ./ vecnorm(G);
        end

        function n = unit_normal(obj, t)
            % outward unit normal at parameter values (the query-surface
            % contract), page-tolerant; grad-based, so independent of the
            % contour orientation
            arguments
                obj LevelSetCurve
                t (1, :, :) double
            end
            sz = size(t);
            P = obj.get_2d_points(reshape(t, 1, []));
            n = reshape(obj.normal_at(P), [2, sz(2:end)]);
        end

        function [P, converged] = project(obj, X, options)
            % Newton projection onto {phi = 0} along grad phi (closest
            % point to first order): the iteration
            % x <- x - damping * phi(x) grad/|grad|^2, vectorized.
            % damping = 1 is plain Newton. For damping < 1 the iterates
            % follow the field line of grad phi through x (explicit Euler
            % with step `damping` on dx/ds = -phi grad/|grad|^2), landing
            % at its end on {phi = 0} up to O(damping); convergence is
            % linear with rate 1 - damping, so the default iteration cap
            % grows like 1/damping. The damping defaults to the curve's
            % projection_damping property.
            arguments
                obj LevelSetCurve
                X double
                options.tolerance (1,1) double = 1e-12
                options.max_iterations double {mustBeScalarOrEmpty} = []
                options.damping (1,1) double {mustBePositive, mustBeLessThanOrEqual(options.damping, 1)} = obj.projection_damping
            end
            if isempty(options.max_iterations)
                options.max_iterations = ceil(50/options.damping);
            end
            sz = size(X);
            X = reshape(X, 2, []);
            P = X;
            for it = 1:options.max_iterations
                z = obj.evaluate(P);
                if all(abs(z) <= options.tolerance), break; end
                G = obj.grad(P);
                P = P - options.damping * z .* G ./ sum(G.^2, 1);
            end
            converged = abs(obj.evaluate(P)) <= max(options.tolerance, 1e3*eps(vecnorm(P)));
            assert(all(converged), "LevelSetCurve:projection", ...
                "%d point(s) failed to project onto the zero contour", nnz(~converged));
            P = reshape(P, sz);
        end

        plot(obj, ax, options)
        obj = discretize(obj, options)

        % ------------------------------------------------------------ %
        % induced parametrization: the ParametrizedCurve query surface
        % on the discretized zero contour (arc-length parameter)
        % ------------------------------------------------------------ %

        function tr = get.t_range(obj)
            assert_discretized(obj);
            tr = [0, obj.arc(end)];
        end

        function tf = is_periodic(~)
            tf = true;   % scope: closed zero contours only
        end

        function p = get_2d_points(obj, t)
            % point on {phi = 0} at arc-length t: polyline interpolation
            % refined by projection, so canonical points lie on the curve
            arguments
                obj LevelSetCurve
                t (1, :, :) double
            end
            assert_discretized(obj);
            sz = size(t);
            s = mod(t(:).', obj.arc(end));
            V = [obj.polyline, obj.polyline(:, 1)];
            p = obj.project(interp1(obj.arc, V.', s).');
            p = reshape(p, [2, sz(2:end)]);
        end

        function T = tangent(obj, t)
            % unit tangent gamma'(t) at arc-length t (|gamma'| = 1 by
            % construction): grad-phi normal rotated onto the direction
            % of increasing arc length
            arguments
                obj LevelSetCurve
                t (1, :) double
            end
            P = obj.get_2d_points(t);
            N = obj.normal_at(P);
            T = [N(2, :); -N(1, :)];
            % orient along increasing s via the local polyline chord
            V = [obj.polyline, obj.polyline(:, 1)];
            D = interp1(obj.arc(1:end-1) + diff(obj.arc)/2, ...
                diff(V, 1, 2).', mod(t(:).', obj.arc(end)), ...
                "nearest", "extrap").';
            T = T .* sign(sum(T .* D, 1));
        end

        function t = map_points_to_interface(obj, X)
            % arc-length coordinate of the Newton projection of X onto the
            % curve (closest point to first order), page-tolerant:
            % (2, ...) -> (1, ...)
            assert_discretized(obj);
            sz = size(X);
            P = obj.project(reshape(X, 2, []));
            % nearest polyline segment, then its interior parameter
            V = obj.polyline;
            W = [V(:, 2:end), V(:, 1)] - V;                 % (2, m) segments
            seg_len2 = sum(W.^2, 1);
            dx = permute(P(1, :), [2 1]) - V(1, :);         % (n, m)
            dy = permute(P(2, :), [2 1]) - V(2, :);
            alpha = min(max((dx.*W(1, :) + dy.*W(2, :)) ./ seg_len2, 0), 1);
            d2 = (dx - alpha.*W(1, :)).^2 + (dy - alpha.*W(2, :)).^2;
            [~, j] = min(d2, [], 2);
            n = size(P, 2);
            lin = sub2ind(size(alpha), (1:n).', j);
            t = obj.arc(j).' + alpha(lin) .* sqrt(seg_len2(j).');
            t = reshape(t.', [1, sz(2:end)]);
        end
    end

    methods (Access = private)
        function assert_discretized(obj)
            assert(~isempty(obj.polyline), "LevelSetCurve:undiscretized", ...
                ["the induced parametrization is not built: call " ...
                 "curve = discretize(curve) (check the box property) " ...
                 "before using this curve as an interface"]);
        end
    end
end
