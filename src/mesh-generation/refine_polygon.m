function [new_vertices, new_connectivity] = refine_polygon(vertices, connectivity, boundary_markers, sub_edge_length)
%REFINE_POLYGON  Subdivide the edges of a polygon into segments of a given length.
%
%   [new_vertices, new_connectivity] = refine_polygon(vertices, ...
%       connectivity, boundary_markers, sub_edge_length)
%
%   vertices         2-by-n_vertices.
%   connectivity     2-by-n_edges, start and end vertex of each edge.
%   boundary_markers 1-by-n_edges, the marker of each edge.
%   sub_edge_length  1-by-n_edges target length per edge; NaN leaves the
%                    edge unrefined. An edge of length L is split into
%                    ceil(L/sub_edge_length) equal segments.
%
%   new_vertices is the vertex list followed by the inserted points;
%   new_connectivity is 3-by-n_segments (start, end, marker), the
%   sub-segments inheriting the marker of their edge.
arguments
    vertices (2, :) double
    connectivity (2, :) double
    boundary_markers (1, :) double
    sub_edge_length (1, :) double
end
N_VERTICES = size(vertices, 2);
N_EDGES = size(connectivity, 2);

x_start  = vertices(:, connectivity(1, :));
x_end = vertices(:, connectivity(2, :));

dx = vecnorm(x_end - x_start);
no_refinement = isnan(sub_edge_length);
sub_edge_length(no_refinement) = dx(no_refinement);
insert_num = ceil(dx ./ sub_edge_length) - 1;

new_vertices = zeros(2, N_VERTICES + sum(insert_num));
new_connectivity = zeros(3, N_EDGES + sum(insert_num));

new_vertices(:, 1:N_VERTICES) = vertices;
new_connectivity(:, 1:N_EDGES) = [connectivity; boundary_markers];

insert_cum = cumsum([0, insert_num]);
for i = 1:N_EDGES
    n_points = insert_num(i);
    refined_points = zeros(2, n_points + 2);
    refined_points(1, :) = linspace(x_start(1, i), x_end(1, i), n_points + 2);
    refined_points(2, :) = linspace(x_start(2, i), x_end(2, i), n_points + 2);
    refined_points = refined_points(:, 2:end-1);

    insert_range = N_VERTICES + (insert_cum(i)+1:insert_cum(i+1));
    new_vertices(:, insert_range) = refined_points;
    new_connectivity(1, [i insert_range]) = [connectivity(1, i) insert_range];
    new_connectivity(2, [i insert_range]) = [insert_range connectivity(2, i)];
    new_connectivity(3, [i insert_range]) = boundary_markers(i);
end
end
