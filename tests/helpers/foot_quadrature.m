function rule = foot_quadrature(old_mesh, new_mesh, elements, foot, q, integration)
%FOOT_QUADRATURE  A quadrature on new elements of piecewise polynomials of an old mesh evaluated at the feet.
%
%   rule = foot_quadrature(old_mesh, new_mesh, elements, foot, q, integration)
%
%   A rule for the integrals over the new elements K (elements, a column of
%   indices of new_mesh) of products of polynomials of K with x -> u(foot(x)),
%   u a piecewise polynomial of old_mesh. foot maps points x (2 x n) to
%   their feet (2 x n). rule has one column (or entry) per point:
%     element  the position in elements of the point's element K;
%     x, rx    the point and its reference coordinates in K;
%     L, r     the old element whose polynomial is evaluated at foot(x),
%              and the reference coordinates of foot(x) in it;
%     w        the weight, relative to the area of K (a column).
%   integration:
%     "quadrature"  the positive-weight rule exact to degree q on each K
%       (collapsed_gauss_rule), each foot evaluated with the old element
%       that holds it (locate_feet: outside the old mesh, the element with
%       the nearest centroid, extrapolated): a quadrature of an integrand
%       that is discontinuous where the feet cross old element faces;
%     "clipped"  exact integration of the piecewise polynomial integrand
%       when the foot is affine on each K. With foot_K(x) = M_K x + b_K the
%       affine interpolant of foot at the vertices of K,
%         int_K u(foot(x)) v(x) dx = sum_L int_{K ∩ foot_K^-1(L)} u_L(foot_K(x)) v(x) dx,
%       the overlaps foot_K(K) ∩ L clipped as convex polygons (Sutherland-
%       Hodgman), fanned into triangles, mapped back by foot_K^-1, and
%       integrated with the rule exact to degree q. The part of foot_K(K)
%       outside the old mesh takes the polynomial of the source of K,
%       extrapolated: the old element with the same vertices when K is
%       present in both meshes, otherwise the old element that holds (or,
%       outside the old mesh, has the nearest centroid to) the foot of the
%       centroid of K. The rule is the source's on all of K plus, on each
%       overlap with an old element L other than the source, L's points
%       with weight w and the source's with weight -w. Exact for integrands
%       of degree q on each overlap, e.g. q = deg u + deg v. A foot that
%       is not affine on some K is an error ("foot_quadrature:notAffine").
%
%See also: foot_face_quadrature, foot_projection, transport_solid_state, collapsed_gauss_rule
arguments
    old_mesh SimplicialGrid
    new_mesh SimplicialGrid
    elements (:,1) double
    foot function_handle
    q (1,1) double
    integration (1,1) string {mustBeMember(integration, ["clipped", "quadrature"])}
end
n = numel(elements);
[rs, w] = collapsed_gauss_rule(q);
nq = numel(w);
x = reshape(new_mesh.get_2d_physical(repmat(rs, 1, 1, n), elements), 2, []);
rule.element = repelem(1:n, nq);
rule.x = x;
rule.rx = repmat(rs, 1, n);
rule.w = repmat(w(:), n, 1);
if integration == "quadrature"
    [rule.L, rule.r] = locate_feet(old_mesh, foot(x));
    return
end

% the affine foot of each element, from its vertices
X = reshape(new_mesh.coordinates(1:2, new_mesh.elements(:, elements)), 2, 3, n);
Y = reshape(foot(reshape(X, 2, [])), 2, 3, n);
E = [X(:, 2, :) - X(:, 1, :), X(:, 3, :) - X(:, 1, :)];               % 2 x 2 x n
M = pagemrdivide([Y(:, 2, :) - Y(:, 1, :), Y(:, 3, :) - Y(:, 1, :)], E);
b = Y(:, 1, :) - pagemtimes(M, X(:, 1, :));
c = mean(X, 2);
y_c = foot(reshape(c, 2, []));
scale = max(abs(Y), [], "all") + max(abs(X), [], "all");
assert(max(abs(y_c - reshape(pagemtimes(M, c) + b, 2, [])), [], "all") <= 1e-12*scale, ...
    "foot_quadrature:notAffine", "The foot is not affine on every element.");
area_K = abs(E(1, 1, :).*E(2, 2, :) - E(1, 2, :).*E(2, 1, :))/2;      % 1 x 1 x n
det_M = abs(M(1, 1, :).*M(2, 2, :) - M(1, 2, :).*M(2, 1, :));

% the source of each element
N_old = size(old_mesh.elements, 2);
T_old = reshape(old_mesh.coordinates(1:2, old_mesh.elements), 2, 3, N_old);
c_old = reshape(mean(T_old, 2), 2, N_old);
[same, source] = ismembertol(reshape(c, 2, []).', c_old.', 1e-14, "ByRows", true, "DataScale", 1);
source = source(:).';
if any(~same)
    source(~same) = locate_feet(old_mesh, y_c(:, ~same));
end

% the overlaps foot_K(K) ∩ L, L other than the source
lo = reshape(min(Y, [], 2), 2, n);  hi = reshape(max(Y, [], 2), 2, n);
lo_old = reshape(min(T_old, [], 2), 2, N_old);  hi_old = reshape(max(T_old, [], 2), 2, N_old);
slack = 1e-12*scale;
meets = lo(1, :).' <= hi_old(1, :) + slack & lo_old(1, :) <= hi(1, :).' + slack ...
      & lo(2, :).' <= hi_old(2, :) + slack & lo_old(2, :) <= hi(2, :).' + slack;   % n x N_old
meets(sub2ind([n, N_old], 1:n, source)) = false;
[K, L] = find(meets);
K = K.';  L = L.';
[P, n_v] = clip_polygons(counterclockwise(Y(:, :, K)), counterclockwise(T_old(:, :, L)));

% fan triangulation of each overlap: (P_1, P_j, P_j+1), j = 2..n_v-1
pair = zeros(1, 0);  a = zeros(2, 0);  bb = a;  cc = a;
for j = 2:size(P, 2) - 1
    has = find(n_v > j);
    pair = [pair, has]; %#ok<AGROW>
    a = [a, reshape(P(:, 1, has), 2, [])]; %#ok<AGROW>
    bb = [bb, reshape(P(:, j, has), 2, [])]; %#ok<AGROW>
    cc = [cc, reshape(P(:, j + 1, has), 2, [])]; %#ok<AGROW>
end
area_T = abs((bb(1, :) - a(1, :)).*(cc(2, :) - a(2, :)) - (bb(2, :) - a(2, :)).*(cc(1, :) - a(1, :)))/2;
% overlaps of zero area (a shared edge or vertex) contribute nothing
owner = K(pair);
keep = area_T > 1e-14*reshape(det_M(owner).*area_K(owner), 1, []);
pair = pair(keep);  owner = owner(keep);
a = a(:, keep);  bb = bb(:, keep);  cc = cc(:, keep);  area_T = area_T(keep);
y_T = reshape(a, 2, 1, []) + reshape(bb - a, 2, 1, []).*rs(1, :) + reshape(cc - a, 2, 1, []).*rs(2, :);
x_T = pagemldivide(M(:, :, owner), y_T - b(:, :, owner));               % 2 x nq x n_tri
w_T = w(:) .* reshape(area_T ./ reshape(det_M(owner).*area_K(owner), 1, []), 1, []);   % nq x n_tri
rx_T = zeros(2, 0);
if ~isempty(owner)
    rx_T = reshape(new_mesh.get_2d_local(x_T, elements(owner)), 2, []);
end
x_T = reshape(x_T, 2, []);
e_T = repelem(owner, nq);

% the source on all of K; on each overlap, L in and the source out
rule.L = [source(rule.element), repelem(L(pair), nq), source(e_T)];
y = [foot(x), reshape(y_T, 2, []), reshape(y_T, 2, [])];
rule.element = [rule.element, e_T, e_T];
rule.x = [rule.x, x_T, x_T];
rule.rx = [rule.rx, rx_T, rx_T];
rule.w = [rule.w; w_T(:); -w_T(:)];
rule.r = reshape(old_mesh.get_2d_local(reshape(y, 2, 1, []), rule.L(:)), 2, []);
end

function T = counterclockwise(T)
% triangles 2 x 3 x n with their vertices in counterclockwise order
signed = (T(1, 2, :) - T(1, 1, :)).*(T(2, 3, :) - T(2, 1, :)) ...
       - (T(2, 2, :) - T(2, 1, :)).*(T(1, 3, :) - T(1, 1, :));
flip = signed < 0;
T(:, [2 3], flip) = T(:, [3 2], flip);
end

function [P, n] = clip_polygons(subject, clip)
% Sutherland-Hodgman, vectorized over pages: the convex polygons
% subject ∩ clip, both counterclockwise triangles (2 x 3 x n). P is
% 2 x 6 x n (the intersection of two triangles has at most six vertices),
% n the vertex count of each page (below 3 when they do not overlap).
% (LinearElasticityProblem's shifted_projection clips the same way.)
n_pages = size(subject, 3);
MAX = 6;
P = zeros(2, MAX, n_pages);
P(:, 1:3, :) = subject;
n = 3*ones(1, n_pages);
pages = 1:n_pages;
for e = 1:3
    a = reshape(clip(:, e, :), 2, []);
    b = reshape(clip(:, mod(e, 3) + 1, :), 2, []);
    side = @(p) (b(1, :) - a(1, :)).*(p(2, :) - a(2, :)) ...
              - (b(2, :) - a(2, :)).*(p(1, :) - a(1, :));
    Q = zeros(2, MAX, n_pages);
    n_new = zeros(1, n_pages);
    for i = 1:MAX
        active = i <= n;
        % the previous vertex, cyclically (vertex n before vertex 1)
        j_prev = (i - 1)*ones(1, n_pages);
        if i == 1, j_prev = max(n, 1); end
        current = reshape(P(:, i, :), 2, []);
        prior = P(:, (pages - 1)*MAX + j_prev);
        d_cur = side(current);
        d_pri = side(prior);
        in_cur = d_cur >= 0;
        in_pri = d_pri >= 0;
        crossing = active & (in_cur ~= in_pri);
        lambda = d_pri ./ (d_pri - d_cur + ~crossing);
        X = prior + (current - prior).*lambda;
        % the crossing point first, then the current vertex when inside
        n_new = n_new + crossing;
        Q = put(Q, X, n_new, crossing);
        keep = active & in_cur;
        n_new = n_new + keep;
        Q = put(Q, current, n_new, keep);
    end
    P = Q;
    n = n_new;
end
end

function Q = put(Q, X, slot, mask)
% Q(:, slot(p), p) = X(:, p) on the pages where mask holds
[~, MAX, n_pages] = size(Q);
p = find(mask);
Q = reshape(Q, 2, MAX*n_pages);
Q(:, (p - 1)*MAX + slot(p)) = X(:, p);
Q = reshape(Q, 2, MAX, n_pages);
end
