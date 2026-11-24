function V = foot_projection(old_mesh, new_mesh, U, deg_from, deg_to, foot, options)
%FOOT_PROJECTION  L2 projection on a new mesh of piecewise polynomials of an old mesh evaluated at the feet.
%
%   V = foot_projection(old_mesh, new_mesh, U, deg_from, deg_to, foot, ...
%       elements = all, quadrature_degree = deg_from + deg_to + 2, ...
%       foot_integration = "clipped")
%
%   U holds Dubiner fields of degree deg_from on old_mesh stacked by rows
%   (component blocks of nchoosek(deg_from + 2, 2) rows), one column per
%   old element; foot maps points x of the new configuration (2 x n) to
%   their feet y (2 x n). V holds, for each new element K in elements (one
%   column each), the dofs of the L2 projection onto Dubiner P_deg_to(K) of
%   x -> U(foot(x)), component by component. The integrals are computed by
%   foot_quadrature with degree quadrature_degree and foot_integration:
%   "clipped" (default; an affine foot on each K) integrates the piecewise
%   polynomial U(foot(.)) exactly on the overlaps of foot(K) with the old
%   elements, the part of foot(K) outside the old mesh taking K's source
%   polynomial, extrapolated; "quadrature" is a positive-weight rule on K
%   (collapsed_gauss_rule) with each foot evaluated in the old element
%   that holds it (locate_feet; outside the old mesh, the element with the
%   nearest centroid, extrapolated), exact when U(foot(.)) is a polynomial
%   of degree at most quadrature_degree - deg_to on K (an affine foot and
%   a foot image inside one old element), and otherwise a quadrature of a
%   discontinuous integrand.
%
%See also: foot_quadrature, locate_feet, transport_solid_state
arguments
    old_mesh SimplicialGrid
    new_mesh SimplicialGrid
    U double
    deg_from (1,1) double
    deg_to (1,1) double
    foot function_handle
    options.elements (:,1) double = (1:size(new_mesh.elements, 2)).'
    options.quadrature_degree double {mustBeScalarOrEmpty} = []
    options.foot_integration (1,1) string {mustBeMember(options.foot_integration, ...
        ["clipped", "quadrature"])} = "clipped"
end
q = options.quadrature_degree;
if isempty(q), q = deg_from + deg_to + 2; end
elements = options.elements;
d_from = nchoosek(deg_from + 2, 2);
d_to = nchoosek(deg_to + 2, 2);
n_c = size(U, 1)/d_from;
V = zeros(n_c*d_to, numel(elements));
if isempty(elements), return, end
rule = foot_quadrature(old_mesh, new_mesh, elements, foot, q, options.foot_integration);
[rs, w] = collapsed_gauss_rule(q);
P = dubiner_basis(deg_to, rs);
mass = P.'*(w(:).*P);                                                    % d_to x d_to
B = dubiner_basis(deg_from, rule.r);                                     % n_points x d_from
test = rule.w .* dubiner_basis(deg_to, rule.rx);                         % n_points x d_to
S = sparse(rule.element, 1:numel(rule.w), 1, numel(elements), numel(rule.w));
for c = 0:n_c - 1
    values = sum(B .* U(c*d_from + (1:d_from), rule.L).', 2);
    V(c*d_to + (1:d_to), :) = mass \ (S*(test.*values)).';
end
end
