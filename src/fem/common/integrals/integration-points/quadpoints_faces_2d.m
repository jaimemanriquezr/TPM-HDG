function [xy, t, w, F] = quadpoints_faces_2d(mesh, order, face_indices, options)
%QUADPOINTS_FACES_2D  Quadrature points on the faces of a mesh.
%
%   [xy, t, w, F] = quadpoints_faces_2d(mesh, order, face_indices, ...
%                                       output = "cell")
%
%   Maps the reference rule get_1d_reference_quadrature(order) to the
%   chosen faces (default all); F is the vector from the first to the
%   second vertex of each face. With output "cell" xy is split into
%   {x; y}.
arguments
    mesh AbstractGrid
    order double = 0
    face_indices (:, 1) double = 1:size(mesh.faces, 2)
    options.output = "cell"
end
[t, w] = get_1d_reference_quadrature(order);

xy_1 = permute(mesh.coordinates(:, mesh.faces(1, face_indices)), [1 3 2]);
xy_2 = permute(mesh.coordinates(:, mesh.faces(2, face_indices)), [1 3 2]);
F = xy_2 - xy_1;

xy = pagemtimes(F, t) + xy_1;
if options.output == "cell"
    xy = {reshape(xy(1, :), [], length(face_indices));
        reshape(xy(2, :), [], length(face_indices))};
end
end
