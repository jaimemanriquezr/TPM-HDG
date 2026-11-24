classdef DeformedCurve < ParametrizedCurve
    % Image of a reference curve under a displacement field, on the
    % reference parameter:
    %
    %   x(t) = X(t) + d(X(t)),   X = reference.get_2d_points(t).
    %
    % It is the current interface Sigma(t) of a Lagrangian solid whose
    % reference interface Sigma_0 is X, paired with it through t by the
    % per-side curve pair {DeformedCurve, reference} (interface_curves,
    % InterfaceProblem.fsi).
    %
    %   curve = DeformedCurve(reference, displacement, mapping = "closest", step = [])
    %   curve = DeformedCurve(reference, d, parametric = true, derivative = d_prime)
    %
    % reference     ParametrizedCurve X(t); t_range and periodicity are its.
    % displacement  handle d(X), 2 x n points -> 2 x n displacements; with
    %               parametric = true, d(t), 1 x n parameters -> 2 x n
    %               displacements: x(t) = X(t) + d(t), e.g. a discrete
    %               displacement evaluated piecewise along Sigma_0.
    % derivative    handle d'(t) (parametric only): the tangent is then
    %               X'(t) + d'(t), with no finite differences, which a
    %               piecewise d needs (a stencil straddling a break would
    %               mix two pieces).
    % mapping       rule of map_points_to_interface: "closest" (closest
    %               point of x) or "normal" (along given normals, as
    %               SmoothCurve).
    % step          finite-difference step in t; [] is 2e-4*|t_range|.
    %
    % tangent: x'(t) by the fourth-order central difference
    %   (x(t-2h) - 8 x(t-h) + 8 x(t+h) - x(t+2h)) / (12 h),
    % evaluated through the reference's parametrization beyond the ends of
    % an open t_range; error O(h^4) plus round-off O(eps/h).
    % unit_normal: the tangent rotated in the same sense in which the
    % reference rotates its tangent into its normal, so a pair keeps one
    % orientation (upward / outward when the reference's is).
    % arc_length: the integral of |x'(t)|.
    %
    % See also: ParametrizedCurve, interface_curves, InterfaceProblem.fsi
    properties (SetAccess = immutable)
        reference     % ParametrizedCurve
        displacement  % function_handle
        parametric    % logical: displacement reads t instead of X(t)
    end

    properties
        mapping (1,1) string {mustBeMember(mapping, ["closest", "normal"])} = "closest"
    end

    methods
        function obj = DeformedCurve(reference, displacement, options)
            arguments
                reference (1,1) ParametrizedCurve
                displacement (1,1) function_handle
                options.mapping (1,1) string = "closest"
                options.step double {mustBeScalarOrEmpty} = []
                options.parametric (1,1) logical = false
                options.derivative function_handle {mustBeScalarOrEmpty} = function_handle.empty
            end
            assert(options.parametric || isempty(options.derivative), ...
                "DeformedCurve:derivative", "derivative needs parametric = true.");
            position = @(t) DeformedCurve.deform(reference, displacement, t, options.parametric);
            h = options.step;
            if isempty(h), h = 2e-4*diff(reference.t_range); end
            % x'(t) on a row of parameters, 2 x n
            x_prime = @(t) (position(t - 2*h) - 8*position(t - h) ...
                + 8*position(t + h) - position(t + 2*h)) / (12*h);
            if ~isempty(options.derivative)
                x_prime = @(t) reference.tangent(t) + 0*t + options.derivative(t);
            end
            speed = @(t) vecnorm(x_prime(reshape(t, 1, [])), 2, 1);
            obj = obj@ParametrizedCurve(@(t) DeformedCurve.component(position, t, 1), ...
                @(t) DeformedCurve.component(position, t, 2), ...
                @(t1, t2) integral(speed, t1, t2), reference.t_range);
            obj.reference = reference;
            obj.displacement = displacement;
            obj.parametric = options.parametric;
            obj.mapping = options.mapping;
            obj.tangent = @(t) x_prime(reshape(t, 1, []));
        end

        function p = get_2d_points(obj, t)
            arguments
                obj DeformedCurve
                t (1, :, :) double
            end
            p = DeformedCurve.deform(obj.reference, obj.displacement, t, obj.parametric);
        end

        function tf = is_periodic(obj)
            tf = is_periodic(obj.reference);
        end

        function t = find_closest_point(obj, x0, y0)
            % ParametrizedCurve's global search, whose fminbnd stops at a
            % relative accuracy of about sqrt(eps), polished by
            % Gauss-Newton steps on |x(t) - x0|^2 (quadratic for points
            % on the curve, where the TPM vertices of a mesh with
            % vertices on Sigma lie); a step that clearly moves the
            % point away ends the polish (at a corner of the curve the
            % minimum sits on the kink, and a step along either one-sided
            % tangent overshoots and can walk far from it)
            t = find_closest_point@ParametrizedCurve(obj, x0, y0);
            t0 = obj.t_range(1);
            period = diff(obj.t_range);
            distance = @(t) sum((obj.get_2d_points(t) - [x0; y0]).^2);
            d = distance(t);
            for it = 1:3
                g = obj.tangent(t);
                t_new = t - sum((obj.get_2d_points(t) - [x0; y0]) .* g) / sum(g.^2);
                if is_periodic(obj)
                    t_new = t0 + mod(t_new - t0, period);
                else
                    t_new = min(max(t_new, obj.t_range(1)), obj.t_range(2));
                end
                d_new = distance(t_new);
                % reject only a clear increase: off the curve the last
                % steps change d below its round-off
                if ~(d_new <= d*(1 + 1e-8) + eps), break, end
                [t, d] = deal(t_new, d_new);
            end
        end

        function n = unit_normal(obj, t)
            % page-tolerant; the reference's rotation sense, pointwise
            arguments
                obj DeformedCurve
                t (1, :, :) double
            end
            t_row = reshape(t, 1, []);
            T = obj.tangent(t_row);
            T0 = obj.reference.tangent(t_row) + 0*t_row;
            N0 = reshape(obj.reference.unit_normal(t_row), 2, []) + 0*t_row;
            sense = sign(T0(1, :).*N0(2, :) - T0(2, :).*N0(1, :));
            N = sense .* [-T(2, :); T(1, :)] ./ vecnorm(T, 2, 1);
            sz = size(t);
            n = reshape(N, [2, sz(2:end)]);
        end

        function v_t = map_points_to_interface(obj, v_x, options)
            arguments
                obj DeformedCurve
                v_x (2, :, :) double
                options.normals (2, :) double = double.empty(2, 0);
                options.reach (1, :) double = Inf;
                options.coordinates = "local";
            end
            switch obj.mapping
                case "closest"
                    v_t = zeros([1 size(v_x, [2 3])]);
                    for i = 1:numel(v_t)
                        v_t(i) = find_closest_point(obj, v_x(1, i), v_x(2, i));
                    end
                case "normal"
                    n_points = numel(v_x)/2;
                    if size(options.normals, 2) ~= n_points
                        error("DeformedCurve:normalsRequired", ...
                            "mapping = ""normal"" needs one normal per point (%d), got %d.", ...
                            n_points, size(options.normals, 2));
                    end
                    v_t = find_normal_intersection(obj, reshape(v_x, 2, []), ...
                        options.normals, options.reach);
                    v_t = reshape(v_t, [1 size(v_x, [2 3])]);
            end
            switch options.coordinates
                case "local"
                    return
                case "physical"
                    v_t = obj.get_2d_points(v_t);
            end
        end
    end

    methods (Static, Access = private)
        function p = deform(reference, displacement, t, parametric)
            % x(t) = X(t) + d(X(t)), or X(t) + d(t), 2 x size(t)(2:end)
            X = reference.get_2d_points(t);
            if parametric
                p = X + reshape(displacement(reshape(t, 1, [])), size(X));
            else
                p = X + reshape(displacement(reshape(X, 2, [])), size(X));
            end
        end

        function c = component(position, t, i)
            p = position(t);
            c = reshape(p(i, :), size(t));
        end
    end
end
