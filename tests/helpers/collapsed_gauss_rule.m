function [rs, w] = collapsed_gauss_rule(q)
%COLLAPSED_GAUSS_RULE  A positive-weight quadrature rule on the reference triangle.
%
%   [rs, w] = collapsed_gauss_rule(q)
%
%   The tensor Gauss-Legendre rule on [0, 1]^2 mapped to the reference
%   triangle r, s >= 0, r + s <= 1 by (u, v) -> (u, v (1 - u)) (Duffy),
%   with ceil((q + 2)/2) points per direction: exact to degree q, every
%   weight positive, sum(w) = 1 (weights relative to the triangle's area,
%   as get_2d_reference_quadrature's). rs is 2 x n, w a column. Unlike the
%   Grundmann-Moeller rules, whose weights change sign, the rule is a
%   stable quadrature of a discontinuous integrand (a field evaluated at
%   the feet of a step, which cross element faces).
%
%See also: get_2d_reference_quadrature, foot_projection, transport_solid_state
arguments
    q (1,1) double {mustBeNonnegative}
end
n = max(ceil((q + 2)/2), 1);
[u, wu] = lgwt(n, 0, 1);
[U, V] = meshgrid(u, u);
[WU, WV] = meshgrid(wu, wu);
rs = [U(:).'; V(:).'.*(1 - U(:).')];
w = WU(:).*WV(:).*(1 - U(:));
w = w/sum(w);
end
