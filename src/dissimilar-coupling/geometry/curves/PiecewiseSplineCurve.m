classdef PiecewiseSplineCurve < ParametrizedCurve
%PIECEWISESPLINECURVE A curve made of cubic-spline pieces joined at corners.
%   Piece i is the not-a-knot cubic spline (MATLAB spline) through the
%   points P_i (2 x m_i) at the knots s_i (1 x m_i, monotone) of its own
%   coordinate s. Consecutive pieces must share their end points, so the
%   curve is continuous; its tangent may jump at the joins (corners). With
%   closed = true the last piece must also end where the first begins, and
%   the curve is periodic.
%
%   Parameter: one unit per piece, t in [0, n_pieces]; on piece i, t in
%   [i-1, i] maps affinely to s from s_i(1) to s_i(end). On a closed curve
%   the parameter wraps (t = n_pieces is t = 0). gamma'(t) is the
%   exact derivative of the spline (so |gamma'| is the arc-length element,
%   as get_sliver_jacobian expects). The normal is the tangent rotated by
%   -90 degrees (ParametrizedCurve), scaled by normal_sign.
%
%Usage:
%   pieces = struct("s", {s1, s2}, "P", {P1, P2});
%   curve = PiecewiseSplineCurve(pieces);
%   loop = PiecewiseSplineCurve(pieces, closed = true);
%
%See also: ParametrizedCurve, PolygonalCurve
    properties
        pieces struct                   % s, P per piece, as given
        normal_sign (1,1) double = 1
        closed (1,1) logical = false
    end

    properties (Hidden)
        position_pp cell                % pp-form of P(s), one per piece
        derivative_pp cell              % pp-form of dP/ds
        s_start double                  % s at the start / end of each piece
        s_end double
    end

    methods
        function obj = PiecewiseSplineCurve(pieces, options)
            arguments
                pieces
                options.normal_sign (1,1) double {mustBeMember(options.normal_sign, [-1, 1])} = 1
                options.closed (1,1) logical = false
            end
            if iscell(pieces), pieces = [pieces{:}]; end
            n = numel(pieces);
            [pos, der] = deal(cell(1, n));
            [s0, s1] = deal(zeros(1, n));
            for i = 1:n
                s = reshape(double(pieces(i).s), 1, []);
                P = double(pieces(i).P);
                assert(size(P, 1) == 2 && size(P, 2) == numel(s), ...
                    "PiecewiseSplineCurve:size", "piece %d: P must be 2 x numel(s)", i);
                assert(all(diff(s) > 0) || all(diff(s) < 0), ...
                    "PiecewiseSplineCurve:knots", "piece %d: knots must be strictly monotone", i);
                if i > 1
                    gap = norm(P(:, 1) - pieces(i-1).P(:, end));
                    assert(gap <= 1e-12*max(1, norm(P(:, 1))), "PiecewiseSplineCurve:gap", ...
                        "pieces %d and %d do not meet (gap %.3e)", i-1, i, gap);
                end
                [ss, o] = sort(s);
                pos{i} = spline(ss, P(:, o));
                [breaks, coefs, ~, order, d] = unmkpp(pos{i});
                der{i} = mkpp(breaks, coefs(:, 1:order-1) .* (order-1:-1:1), d);
                s0(i) = s(1);  s1(i) = s(end);
            end
            if options.closed
                gap = norm(pieces(1).P(:, 1) - pieces(n).P(:, end));
                assert(gap <= 1e-12*max(1, norm(pieces(1).P(:, 1))), "PiecewiseSplineCurve:gap", ...
                    "a closed curve must end where it begins (gap %.3e)", gap);
            end

            % piece index and local coordinate s of each parameter value
            locate = @(t) PiecewiseSplineCurve.locate(t, n, s0, s1, options.closed);
            x_of = @(t) PiecewiseSplineCurve.component(pos, locate, t, 1);
            y_of = @(t) PiecewiseSplineCurve.component(pos, locate, t, 2);
            tangent = @(t) PiecewiseSplineCurve.derivative(der, locate, t, s1 - s0);
            % arc length: |gamma'| integrated piece by piece
            arc_of = @(t1, t2) PiecewiseSplineCurve.arc(tangent, t1, t2);

            obj = obj@ParametrizedCurve(x_of, y_of, arc_of, [0 n]);
            obj.pieces = pieces;
            obj.normal_sign = options.normal_sign;
            obj.closed = options.closed;
            obj.position_pp = pos;
            obj.derivative_pp = der;
            obj.s_start = s0;
            obj.s_end = s1;
            obj.tangent = tangent;
        end

        function tf = is_periodic(obj)
            tf = obj.closed;
        end

        function n = unit_normal(obj, t)
            n = obj.normal_sign * unit_normal@ParametrizedCurve(obj, t);
        end

        function v_t = map_points_to_interface(obj, v_x, options)
            % closest-point projection (base-class find_closest_point)
            arguments
                obj PiecewiseSplineCurve
                v_x (2, :, :) double
                options.normals (2, :) double = double.empty(2, 0) %#ok<INUSA>
                options.coordinates = "local"
            end
            v_t = zeros([1 size(v_x, [2 3])]);
            x0 = v_x(1, :); y0 = v_x(2, :);
            for i = 1:numel(v_t)
                v_t(i) = find_closest_point(obj, x0(i), y0(i));
            end
        end
    end

    methods (Static, Access = private)
        function [piece, s] = locate(t, n, s0, s1, closed)
            % a t on a join takes the piece starting there, except at the
            % end of an open curve; a closed curve wraps t into [0, n)
            if closed
                t = mod(t, n);
            else
                t = min(max(t, 0), n);
            end
            piece = min(floor(t) + 1, n);
            tau = t - (piece - 1);
            s = s0(piece) + tau .* (s1(piece) - s0(piece));
        end

        function v = component(pp, locate, t, c)
            [piece, s] = locate(t);
            v = zeros(size(t));
            for i = unique(piece(:)).'
                mask = piece == i;
                P = ppval(pp{i}, s(mask));
                v(mask) = P(c, :);
            end
        end

        function T = derivative(der, locate, t, ds_dt)
            [piece, s] = locate(reshape(t, 1, []));
            T = zeros(2, numel(s));
            for i = unique(piece(:)).'
                mask = piece == i;
                T(:, mask) = ppval(der{i}, s(mask)) * ds_dt(i);
            end
        end

        function L = arc(tangent, t1, t2)
            speed = @(t) reshape(vecnorm(tangent(t), 2, 1), size(t));
            a = min(t1, t2);  b = max(t1, t2);
            joins = ceil(a):floor(b);
            L = integral(speed, a, b, Waypoints = joins(joins > a & joins < b), ...
                AbsTol = 1e-14, RelTol = 1e-12);
        end
    end
end
