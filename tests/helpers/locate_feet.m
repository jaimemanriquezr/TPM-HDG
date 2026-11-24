function [elements, rs] = locate_feet(mesh, y)
%LOCATE_FEET  The elements of a mesh that hold given points, and the points' reference coordinates.
%
%   [elements, rs] = locate_feet(mesh, y)
%
%   y is 2 x n. elements (1 x n) is, for each point, an element of mesh
%   that contains it (triangulation.pointLocation: a point on a shared face
%   or vertex gets one of its elements); a point outside the mesh gets the
%   element with the nearest centroid, whose polynomial is then
%   extrapolated (the source rule of remap_solid_state). rs (2 x n) holds
%   the reference coordinates of each point in its element.
%
%See also: foot_projection, transport_solid_state
arguments
    mesh SimplicialGrid
    y (2,:) double
end
T = triangulation(mesh.elements.', mesh.coordinates(1:2, :).');
elements = pointLocation(T, y.').';
outside = find(isnan(elements));
if ~isempty(outside)
    c = squeeze(mean(reshape(mesh.coordinates(1:2, mesh.elements), 2, 3, []), 2));   % 2 x N
    for j = outside
        [~, elements(j)] = min((c(1, :) - y(1, j)).^2 + (c(2, :) - y(2, j)).^2);
    end
end
rs = reshape(mesh.get_2d_local(reshape(y, 2, 1, []), elements(:)), 2, []);
end
