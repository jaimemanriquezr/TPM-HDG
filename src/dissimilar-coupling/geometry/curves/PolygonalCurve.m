classdef PolygonalCurve < ParametrizedCurve
%POLYGONALCURVE The linear interpolant of an interface, used as the interface.
%   A piecewise-linear curve through a given vertex list. Meshing both
%   subdomains to the interpolant of Sigma and passing the same polygon to
%   the solver gives zero-length transfer paths and interface integrals on
%   the chords: the non-TPM baseline.
%
%   A chord departs from a curve of curvature kappa by kappa*h^2/8, so this
%   baseline caps the error at O(h^2) regardless of the polynomial degree.
%
%   Parameter: one unit per segment, t in [0, n_segments]. On a closed
%   polygon the parameter wraps (is_periodic is true), so the seam is
%   crossed like any vertex. The tangent, and hence the normal, is
%   piecewise constant and jumps at every vertex.
%
%Usage:
%   V = centre + radius*[cos(th); sin(th)];      % th = linspace(0,2*pi,N+1)
%   poly = PolygonalCurve(V);
%
%See also: ParametrizedCurve, fitted_split, CircularCurve
    properties
        vertices double        % 2 x (n_segments+1), first point repeated at the end
        chords   double        % 1 x n_segments, segment lengths
        closed   (1,1) logical
        normal_sign (1,1) double = 1   % +1: the tangent rotated by -90 degrees (ParametrizedCurve); -1 flips it
    end

    methods
        function obj = PolygonalCurve(vertices, options)
            arguments
                vertices double {mustBeNonempty}
                % A closed polygon accepts the vertex list either with or
                % without the repeated closing point.
                options.closed (1,1) logical = true
                % Orientation of unit_normal; interface_polygon sets it to
                % agree with the curve being interpolated.
                options.normal_sign (1,1) double {mustBeMember(options.normal_sign, [-1, 1])} = 1
            end
            V = vertices;
            if options.closed
                if norm(V(:, end) - V(:, 1)) <= 1e-14 * max(1, norm(V(:, 1)))
                    V(:, end) = [];
                end
                V = [V, V(:, 1)];
            end
            n_seg = size(V, 2) - 1;
            assert(n_seg >= 2, "PolygonalCurve:tooFewVertices", ...
                "need at least two segments, got %d", n_seg);
            s = 0:n_seg;
            chords = vecnorm(diff(V, 1, 2));
            cum = [0, cumsum(chords)];

            if options.closed
                wrap = @(t) mod(t, n_seg);      % t = n_seg lands on vertex 1
            else
                wrap = @(t) min(max(t, 0), n_seg);
            end
            x_of = @(t) interp1(s, V(1, :), wrap(t), "linear");
            y_of = @(t) interp1(s, V(2, :), wrap(t), "linear");
            % Arc length is exact on a polygon: interpolate the cumulative
            % chord length. It clamps rather than wraps, so arc_length(0,
            % n_seg) is the whole perimeter.
            len_of = @(t) interp1(s, cum, min(max(t, 0), n_seg), "linear");
            arc_of = @(t1, t2) abs(len_of(t2) - len_of(t1));

            obj = obj@ParametrizedCurve(x_of, y_of, arc_of, [0 n_seg]);
            obj.vertices = V;
            obj.chords = chords;
            obj.closed = options.closed;
            obj.normal_sign = options.normal_sign;
            % gamma'(t) is the segment vector: constant per segment, and
            % its norm is the chord length, so arc length and parameter
            % agree the way the base class expects.
            obj.tangent = @(t) segment_tangent(V, n_seg, wrap(t));
        end

        function tf = is_periodic(obj)
            tf = obj.closed;
        end

        function n = unit_normal(obj, t)
            n = obj.normal_sign * unit_normal@ParametrizedCurve(obj, t);
        end

        function v_t = map_points_to_interface(obj, v_x, options)
            % Closest-point projection (base-class find_closest_point):
            % the cached sampling brackets the segment and fminbnd refines
            % inside it, exactly, since the distance is quadratic along a
            % chord.
            arguments
                obj PolygonalCurve
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
end

%=========================================================================%
function T = segment_tangent(V, n_seg, t)
% Segment index for each parameter value. A t sitting exactly on a vertex
% takes the segment starting there, except at the very end, which has no
% following segment and takes the one arriving.
idx = floor(t) + 1;
idx = min(max(idx, 1), n_seg);
T = V(:, idx + 1) - V(:, idx);
T = reshape(T, 2, []);
end
