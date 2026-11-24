function mesh = fill_band(obj, interface, options)
%-------------------------------------------------------------------------
%Description: Closes the cut band of a carved mesh: the strip between its
%             staircase boundary and Sigma is meshed (mesh2d, staircase
%             nodes kept, so the strip meets the carved elements node to
%             node) and merged in. The result is fitted to Sigma, and the
%             grid structure survives away from it.
%
%             Carved elements with a vertex within drop*h of Sigma are
%             removed first; otherwise the strip holds slivers where Sigma
%             grazes a staircase corner.
%
%             Assumes the stenosis-type layout: Sigma a graph y(x) over the
%             mesh's x-range (x_coord(t) = t), and the carved boundary off
%             the bounding box a single staircase path reaching both
%             vertical walls. Only "Interface" is tagged (on Sigma); the
%             caller tags the walls, as it does for any "array" mesh.
%
%Usage:       [~, carved] = bg.split(interface, "cut_elements", "discard");
%             fluid = carved.fill_band(interface, h = 0.05);
%
%See also: split, submesh
%-------------------------------------------------------------------------
arguments
    obj SimplicialGrid
    interface ParametrizedCurve   % a graph over x: x_coord(t) = t
    % Grid spacing: sets Sigma's sampling and the drop distance. Empty
    % takes the median face length, which on a structured grid is h.
    options.h double {mustBeScalarOrEmpty} = []
    % Drop distance in multiples of h (vertical distance to Sigma).
    options.drop (1,1) double {mustBeNonnegative} = 0.3
end
h = options.h;
if isempty(h), h = median(obj.area); end
V = obj.coordinates;
E = obj.elements;
tol = 1e-9*max(1, max(abs(V(:))));
x_lo = min(V(1, :));  x_hi = max(V(1, :));

t_check = linspace(x_lo, x_hi, 101);
assert(max(abs(interface.x_coord(t_check) - t_check)) < tol, ...
    "fill_band:notGraph", "Sigma must be a graph over x (x_coord(t) = t).");
sigma = @(x) interface.y_coord(x);

if options.drop > 0
    near = abs(sigma(V(1, :)) - V(2, :)) < options.drop*h;
    E = E(:, ~any(near(E), 1));
end
y_lo = min(V(2, :));  y_hi = max(V(2, :));
used = unique(E(:));
renumber = zeros(1, size(V, 2));
renumber(used) = 1:numel(used);
V = V(:, used);
E = renumber(E);

% Staircase: boundary edges off the bounding box, ordered x_lo -> x_hi.
edges = sort([E([1 2], :), E([2 3], :), E([3 1], :)], 1).';
[unique_edges, ~, ic] = unique(edges, "rows");
boundary = unique_edges(accumarray(ic, 1) == 1, :);
on_line = @(c, value) all(abs(reshape(V(c, boundary.'), 2, []) - value) < tol, 1).';
% The far wall is the horizontal side away from Sigma; the near one would
% hold the staircase's own flat runs.
if mean(sigma(V(1, :)) - V(2, :)) > 0, y_far = y_lo; else, y_far = y_hi; end
on_box = on_line(1, x_lo) | on_line(1, x_hi) | on_line(2, y_far);
stair = boundary(~on_box, :);
G = graph(stair(:, 1), stair(:, 2));
deg = degree(G);
ends = find(deg == 1);
assert(numel(ends) == 2 && all(deg <= 2), "fill_band:staircase", ...
    "The carved boundary off the bounding box is not a single path.");
[~, first] = min(V(1, ends));
path = shortestpath(G, ends(first), ends(3 - first));
assert(numel(path) == nnz(deg), "fill_band:staircase", ...
    "The staircase path misses some of its own vertices.");
P = V(:, path);
assert(abs(P(1, 1) - x_lo) < tol && abs(P(1, end) - x_hi) < tol, ...
    "fill_band:staircase", "The staircase does not reach both walls.");

% Strip polygon: staircase (x_lo -> x_hi), then Sigma back (x_hi -> x_lo);
% the wall segments close it. "preserve" keeps the polygon edges unsplit.
xs = linspace(x_hi, x_lo, round((x_hi - x_lo)/h) + 1);
polygon = [P, [xs; sigma(xs)]];
n = size(polygon, 2);
[node, ~, tri] = refine2(polygon.', [1:n; 2:n, 1].', {}, h, ...
    ref1 = "preserve", disp = Inf);
assert(all(ismembertol(polygon.', node, tol, "ByRows", true, "DataScale", 1)), ...
    "fill_band:mesher", "mesh2d moved a polygon node.");

% Merge, identify shared nodes, orient counterclockwise.
X = [V, node.'];
T = [E, size(V, 2) + tri.'];
[~, keep, map] = uniquetol(X.', tol, "ByRows", true, "DataScale", 1);
X = X(:, keep);
T = reshape(map(T), 3, []);
signed = (X(1, T(2, :)) - X(1, T(1, :))).*(X(2, T(3, :)) - X(2, T(1, :))) ...
       - (X(1, T(3, :)) - X(1, T(1, :))).*(X(2, T(2, :)) - X(2, T(1, :)));
T([2 3], signed < 0) = T([3 2], signed < 0);

mesh = SimplicialGrid(X, T, mesh_generator = "array");
mesh = apply_tag(mesh, "Interface", @(Y) abs(Y(2, :) - sigma(Y(1, :))) < tol);
end
