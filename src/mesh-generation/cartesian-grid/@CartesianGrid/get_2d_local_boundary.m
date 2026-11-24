function rs_points = get_2d_local_boundary(obj, t_points, element_list, side_number)
%GET_2D_LOCAL_BOUNDARY  Reference coordinates of points on element sides.
%
%   rs_points = get_2d_local_boundary(mesh, t_points, element_list, side_number)
%
%   t_points     parameters in [0, 1], n_t-by-n, one column per pair.
%   element_list elements, n-by-1.
%   side_number  sides (1 to 4; others give NaN), n-by-1.
%   rs_points    2-by-n_t-by-n reference coordinates on the given side of
%                each element. Every side runs in the +r or +s direction
%                (the face parametrization), so the element is immaterial:
%                side 1 (r, 0), 2 (1, s), 3 (r, 1), 4 (0, s).
%
%See also: SimplicialGrid/get_2d_local_boundary, get_physical_points
arguments
    obj CartesianGrid %#ok<INUSA>
    t_points double
    element_list (:,1) double
    side_number (:,1) double
end
assert(all(0 <= t_points(:)) && all(t_points(:) <= 1), ...
    "Local coordinates outside of [0, 1] range!")
assert(numel(element_list) == size(t_points, 2), ...
    "Local points array does not match elements in quantity!")
assert(numel(element_list) == numel(side_number), ...
    "Element and side numbers do not match in quantity!")

t = permute(t_points, [3 1 2]);          % 1 x n_t x n
side = reshape(side_number, 1, 1, []);   % 1 x 1 x n

r = (side == 1 | side == 3).*t + (side == 2);
s = (side == 2 | side == 4).*t + (side == 3);
bad = ~ismember(side, 1:4);
r(:, :, bad) = nan;
s(:, :, bad) = nan;
rs_points = [r; s];                      % 2 x n_t x n
end
