function grid = graded_cavity_grid(N, top_ratio)
%GRADED_CAVITY_GRID  Crossed N x N grid of the unit square, graded toward the lid y = 1.
%
%   grid = graded_cavity_grid(N, top_ratio)
%
%   The vertices of StructuredSimplicialGrid(N, N, "crossed") are mapped by
%   y = 1 - (1 - eta) (r + (1 - r)(1 - eta)), r = top_ratio in (0, 1]: a
%   monotone quadratic with dy/deta = r at the lid and 2 - r at the bottom,
%   so the rows next to y = 1 are r times the uniform height 1/N (r = 1 is
%   the uniform grid). x is unchanged. The grid is rebuilt from the mapped
%   vertices and the same connectivity, so tags and geometry are those of
%   a freshly constructed grid. Used where a body approaches the lid
%   closer than the uniform element height (the carved fluid must keep
%   elements in the gap).
%
%See also: StructuredSimplicialGrid, cavity_problems
arguments
    N (1,1) double {mustBeInteger, mustBePositive}
    top_ratio (1,1) double {mustBeInRange(top_ratio, 0, 1, "exclude-lower")} = 0.25
end
uniform = StructuredSimplicialGrid(N, N, "crossed");
X = uniform.coordinates;
s = 1 - X(2, :);
X(2, :) = 1 - s.*(top_ratio + (1 - top_ratio)*s);
grid = SimplicialGrid(X, uniform.elements, mesh_generator = "array");
end
