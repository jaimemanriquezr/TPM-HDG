function rule = foot_face_quadrature(old_mesh, new_mesh, faces, foot, n_points, integration)
%FOOT_FACE_QUADRATURE  A quadrature on new faces of piecewise polynomials of an old mesh evaluated at the feet.
%
%   rule = foot_face_quadrature(old_mesh, new_mesh, faces, foot, n_points, integration)
%
%   The face analogue of foot_quadrature: a rule for the integrals over the
%   new faces e (faces, a column of indices of new_mesh) of products of
%   polynomials of e with x -> u(foot(x)), u a piecewise polynomial of
%   old_mesh. A face is parametrized by t in [0, 1] from its first vertex
%   to its second. rule has one entry per point:
%     face  the position in faces of the point's face;
%     t     the point's parameter on it (a column);
%     L, r  the old element whose polynomial is evaluated at foot(x), and
%           the reference coordinates of foot(x) in it;
%     w     the weight, relative to the length of the face (a column).
%   integration:
%     "quadrature"  n_points Gauss points on each face, each foot
%       evaluated with the old element that holds it (locate_feet);
%     "clipped"  exact integration of the piecewise polynomial integrand
%       when the foot is affine on each face: the segment foot(e) (from the
%       feet of its end points) is clipped against the old elements it
%       meets (a parameter interval per element), with n_points Gauss
%       points on each interval. The part of foot(e) outside the old mesh
%       takes the polynomial of the face's source, extrapolated: the old
%       element that holds (or, outside the old mesh, has the nearest
%       centroid to) the foot of the midpoint of e. The rule is the
%       source's on all of e plus, on each interval of an old element L
%       other than the source, L's points with weight w and the source's
%       with weight -w. Exact for integrands of degree 2 n_points - 1 on
%       each interval. A foot that is not affine on some face is an error
%       ("foot_face_quadrature:notAffine").
%
%See also: foot_quadrature, transport_solid_state
arguments
    old_mesh SimplicialGrid
    new_mesh SimplicialGrid
    faces (:,1) double
    foot function_handle
    n_points (1,1) double
    integration (1,1) string {mustBeMember(integration, ["clipped", "quadrature"])}
end
n = numel(faces);
[tg, wg] = lgwt(n_points, 0, 1);
tg = tg(:);  wg = wg(:);
xa = new_mesh.coordinates(1:2, new_mesh.faces(1, faces));
xb = new_mesh.coordinates(1:2, new_mesh.faces(2, faces));
on_face = @(p, a, b, t) a(:, p) + (b(:, p) - a(:, p)).*t(:).';
rule.face = repelem(1:n, n_points);
rule.t = repmat(tg, n, 1);
rule.w = repmat(wg, n, 1);
if integration == "quadrature"
    [rule.L, rule.r] = locate_feet(old_mesh, foot(on_face(rule.face, xa, xb, rule.t)));
    return
end

% the segments foot(e) and the source of each face
ya = foot(xa);  yb = foot(xb);
y_m = foot((xa + xb)/2);
scale = max(abs([ya, yb]), [], "all") + max(abs([xa, xb]), [], "all");
assert(max(abs(y_m - (ya + yb)/2), [], "all") <= 1e-12*scale, "foot_face_quadrature:notAffine", ...
    "The foot is not affine on every face.");
source = locate_feet(old_mesh, y_m);

% the parameter interval of foot(e) in each old element L it may meet
N_old = size(old_mesh.elements, 2);
T_old = reshape(old_mesh.coordinates(1:2, old_mesh.elements), 2, 3, N_old);
lo = min(ya, yb);  hi = max(ya, yb);
lo_old = reshape(min(T_old, [], 2), 2, N_old);  hi_old = reshape(max(T_old, [], 2), 2, N_old);
slack = 1e-12*scale;
meets = lo(1, :).' <= hi_old(1, :) + slack & lo_old(1, :) <= hi(1, :).' + slack ...
      & lo(2, :).' <= hi_old(2, :) + slack & lo_old(2, :) <= hi(2, :).' + slack;   % n x N_old
meets(sub2ind([n, N_old], 1:n, source)) = false;
[F, L] = find(meets);
F = F.';  L = L.';
t0 = zeros(1, numel(F));  t1 = ones(1, numel(F));
signed = (T_old(1, 2, :) - T_old(1, 1, :)).*(T_old(2, 3, :) - T_old(2, 1, :)) ...
       - (T_old(2, 2, :) - T_old(2, 1, :)).*(T_old(1, 3, :) - T_old(1, 1, :));
orientation = sign(reshape(signed(L), 1, []));
for e = 1:3
    % side(p) >= 0 inside: side(ya + t (yb - ya)) = s_a + t (s_b - s_a)
    p = reshape(T_old(:, e, L), 2, []);
    d = reshape(T_old(:, mod(e, 3) + 1, L), 2, []) - p;
    side = @(y) orientation.*(d(1, :).*(y(2, F) - p(2, :)) - d(2, :).*(y(1, F) - p(1, :)));
    s_a = side(ya);  s_b = side(yb);
    slope = s_b - s_a;
    bound = -s_a./slope;
    t0(slope > 0) = max(t0(slope > 0), bound(slope > 0));
    t1(slope < 0) = min(t1(slope < 0), bound(slope < 0));
    t1(slope == 0 & s_a < 0) = -Inf;
end
% intervals of zero length (a vertex, or none) contribute nothing
keep = t1 - t0 > 1e-12;
F = F(keep);  L = L(keep);  t0 = t0(keep);  t1 = t1(keep);
t_I = t0 + tg.*(t1 - t0);                                              % n_points x n_I
w_I = wg.*(t1 - t0);
f_I = repelem(F, n_points);

% the source on all of e; on each interval, L in and the source out
rule.L = [source(rule.face), repelem(L, n_points), source(f_I)];
rule.face = [rule.face, f_I, f_I];
rule.t = [rule.t; t_I(:); t_I(:)];
rule.w = [rule.w; w_I(:); -w_I(:)];
y = on_face(rule.face, ya, yb, rule.t);
rule.r = reshape(old_mesh.get_2d_local(reshape(y, 2, 1, []), rule.L(:)), 2, []);
end
