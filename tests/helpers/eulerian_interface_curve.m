function [curve, info] = eulerian_interface_curve(reference, mesh, k, U, options)
%EULERIAN_INTERFACE_CURVE  The current interface of an Eulerian solid from its displacement.
%
%   [curve, info] = eulerian_interface_curve(reference, mesh, k, U, ...
%       breaks = [], spacing = [], velocity = [], time_step = 0)
%
%   Sigma = {x : x - u(x) in Sigma_0}, Sigma_0 = reference (a
%   ParametrizedCurve, the solid's reference boundary or a part of it),
%   as the DeformedCurve x(t) = X(t) + d(t) on the reference parameter:
%   at sample parameters t_j the point x_j solves x = X(t_j) + u(x) by
%   fixed-point iteration from x = X(t_j) (a contraction while |grad u| <
%   1), and d is the cubic spline through d_j = x_j - X(t_j), one spline
%   per piece between consecutive breaks (corners of Sigma_0, where d is
%   continuous but its derivative is not), the samples wrapped around on a
%   periodic reference. u is the Dubiner field U (2*d2 x n_el, x over y)
%   on mesh; a point outside the mesh takes the polynomial of the element
%   of the nearest vertex whose closure it is nearest (an extension of at
%   most the gap). U may instead be a handle u(x) (2 x n -> 2 x n), the
%   exact displacement (mesh and k unused).
%   spacing: the sample spacing in arc length of Sigma_0 ([] = a quarter
%   of the mean face length of mesh).
%   velocity, time_step: a velocity field V (same layout as U, or a
%   handle) moves the samples to x_j + time_step V(x_j), the one-step
%   forward prediction of the material points (default: none).
%   zero_on: rows [a b] of reference parameter intervals (closed) on which
%   d = 0 (a clamped part of Sigma_0, or a single point a = b).
%
%   info: t (samples), x (2 x n points), residual (max |x - X - u(x)|),
%   iterations, outside (samples evaluated by extension), max_step (the
%   largest move x_j - X(t_j)).
%
%See also: DeformedCurve, eulerian_hyperelastic_fsi_march
arguments
    reference ParametrizedCurve
    mesh
    k
    U
    options.breaks double = []
    options.spacing double = []
    options.velocity = []
    options.time_step (1,1) double = 0
    options.tolerance (1,1) double = 1e-13
    options.max_iterations (1,1) double = 50
    options.zero_on (:,2) double = zeros(0, 2)
end
tr = reference.t_range;
periodic = is_periodic(reference);
breaks = unique([tr(1), options.breaks(:).', tr(2)]);
if isa(U, "function_handle")
    evaluate_u = @(x) handle_values(U, x);
else
    locator = point_locator(mesh, k);
    evaluate_u = @(x) locator(U, x);
end
spacing = options.spacing;
if isempty(spacing)
    spacing = mean(mesh.area)/4;
end
% samples per piece, ends included (shared ends counted once)
t = zeros(1, 0);
for p = 1:numel(breaks) - 1
    n = max(8, ceil(reference.arc_length(breaks(p), breaks(p + 1))/spacing));
    tp = linspace(breaks(p), breaks(p + 1), n + 1);
    if p > 1, tp(1) = []; end
    t = [t, tp]; %#ok<AGROW>
end
if periodic, t(end) = []; end
X = reshape(reference.get_2d_points(t), 2, []);
fixed = false(1, numel(t));
for i = 1:size(options.zero_on, 1)
    fixed = fixed | (t >= options.zero_on(i, 1) - 1e-12 & t <= options.zero_on(i, 2) + 1e-12);
end
free = ~fixed;
x = X;
it = 0;
outside = false(1, numel(t));
residual = 0;
if any(free)
    for it = 1:options.max_iterations
        u = evaluate_u(x(:, free));
        x_new = X(:, free) + u;
        step = max(vecnorm(x_new - x(:, free)), [], "all");
        x(:, free) = x_new;
        if step <= options.tolerance*max(1, max(abs(X), [], "all")), break, end
    end
    [u, outside(free)] = evaluate_u(x(:, free));
    residual = max(vecnorm(x(:, free) - X(:, free) - u), [], "all");
    if ~isempty(options.velocity) && options.time_step ~= 0
        if isa(options.velocity, "function_handle")
            v = options.velocity(x(:, free));
        else
            v = locator(options.velocity, x(:, free));
        end
        x(:, free) = x(:, free) + options.time_step*v;
    end
end
d = x - X;
% one cubic spline per piece; on a periodic reference the samples wrap
% around the seam (three on each side) so the spline is smooth there
pieces = cell(1, numel(breaks) - 1);
for p = 1:numel(breaks) - 1
    if periodic && numel(breaks) == 2
        period = diff(tr);
        m = 3;
        ts = [t(end - m + 1:end) - period, t, t(1:m) + period];
        ds = [d(:, end - m + 1:end), d, d(:, 1:m)];
    else
        in = t >= breaks(p) - 1e-14 & t <= breaks(p + 1) + 1e-14;
        if periodic && p == numel(breaks) - 1
            ts = [t(in), tr(2)];
            ds = [d(:, in), d(:, 1)];
        else
            ts = t(in);
            ds = d(:, in);
        end
    end
    pieces{p} = spline(ts, ds);
end
displacement = @(s) piecewise(pieces, breaks, s, periodic, tr, false);
derivative = @(s) piecewise(pieces, breaks, s, periodic, tr, true);
curve = DeformedCurve(reference, displacement, parametric = true, derivative = derivative);
info = struct("t", t, "x", x, "residual", residual, "iterations", it, "outside", nnz(outside), ...
    "max_step", max(vecnorm(d), [], "all"));
end

function [u, outside] = handle_values(U, x)
u = U(x);
outside = false(1, size(x, 2));
end

function v = piecewise(pieces, breaks, s, periodic, tr, derivative)
% the spline of the piece containing each parameter (2 x numel(s)), or
% its derivative
s = reshape(s, 1, []);
if periodic
    s = tr(1) + mod(s - tr(1), diff(tr));
else
    s = min(max(s, tr(1)), tr(2));
end
v = zeros(2, numel(s));
p = discretize(s, breaks);
p(isnan(p)) = numel(pieces);
for i = unique(p)
    sel = p == i;
    pp = pieces{i};
    if derivative
        pp = spline_derivative(pp);
    end
    v(:, sel) = reshape(ppval(pp, s(sel)), 2, []);
end
end

function dp = spline_derivative(pp)
[b, c, l, o, d] = unmkpp(pp);
dc = c(:, 1:o - 1).*(o - 1:-1:1);
dp = mkpp(b, dc, d);
end

function locate = point_locator(mesh, k)
% a handle (U, x) -> [u(x) (2 x n), outside (1 x n)] for Dubiner fields
% U (2*d2 x n_el) on mesh
T = triangulation(mesh.elements(1:3, :).', mesh.coordinates.');
attachments = vertexAttachments(T);
d2 = nchoosek(k + 2, 2);
    function [u, outside] = evaluate(U, x)
        K = pointLocation(T, x.');
        outside = isnan(K).';
        for i = find(outside)
            % the element of the nearest vertex whose barycentric
            % coordinates at x are least negative
            candidates = attachments{nearestNeighbor(T, x(:, i).')};
            B = cartesianToBarycentric(T, candidates(:), repmat(x(:, i).', numel(candidates), 1));
            [~, j] = max(min(B, [], 2));
            K(i) = candidates(j);
        end
        rs = mesh.get_2d_local(reshape(x, 2, 1, []), K(:));
        P = dubiner_basis(k, reshape(rs(1, 1, :), 1, []), reshape(rs(2, 1, :), 1, []));   % n x d2
        u = [sum(P.*U(1:d2, K).', 2), sum(P.*U(d2 + (1:d2), K).', 2)].';
    end
locate = @evaluate;
end
