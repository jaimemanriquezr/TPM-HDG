function xy_points = get_1d_physical(obj, t_points, face_list)
%GET_1D_PHYSICAL  Physical points at parameters along mesh faces.
%
%   xy_points = get_1d_physical(mesh, t_points, face_list)
%
%   t_points  parameters in [0, 1], n_t-by-n_faces (or reshaped to it).
%   face_list indices of the faces, one per column of t_points.
%   xy_points 2-by-n_t-by-n_faces, x_1 + t*(x_2 - x_1) along each face
%             from its first to its second vertex.
arguments
    obj SimplicialGrid
    t_points double
    face_list (1, :) double
end
assert(all(0 <= t_points(:)) && all(t_points(:) <= 1), ...
    "Local coordinates outside of [0, 1] range!")
t_points = reshape(t_points, [], length(face_list));
assert(numel(face_list) == size(t_points, 2), ...
    "Number of input points and faces do not match!")
    
n_patches = size(t_points, 2);
xy_points = zeros(2, size(t_points, 1), n_patches);
for e = 1:n_patches
    F = face_list(e);
    xy_end = obj.coordinates(1:2, obj.faces(1:2, F));
    t = reshape(t_points(:, e), 1, []);
    xy_points(:, :, e) = xy_end(:, 1) + t.*(xy_end(:, 2) - xy_end(:, 1));
end

end