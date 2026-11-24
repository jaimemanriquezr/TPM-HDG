function rs_points = get_2d_local_boundary(obj, t_points, element_list, side_number)
%GET_2D_LOCAL_BOUNDARY  Reference coordinates of points on element sides.
%
%   rs_points = get_2d_local_boundary(mesh, t_points, element_list, side_number)
%
%   t_points     parameters in [0, 1], n_t-by-n, one column per pair.
%   element_list elements, n-by-1.
%   side_number  sides (1, 2 or 3; others give NaN), n-by-1.
%   rs_points    2-by-n_t-by-n reference coordinates on the given side of
%                each element; the permutation flag of (side, element)
%                fixes the orientation of the parameter along the side.
arguments
    obj SimplicialGrid
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

% one permutation flag per (side, element) pair, 1 x n
mu = reshape(obj.permutation(sub2ind(size(obj.permutation), side_number, element_list)), 1, []);
t = permute(t_points, [3 1 2]);          % 1 x n_t x n
mu = reshape(mu, 1, 1, []);              % 1 x 1 x n
side = reshape(side_number, 1, 1, []);   % 1 x 1 x n

% along the side: a + b*t with (a, b) set by the permutation flag
along = (2 - mu) + (2*mu - 3).*t;
r = (side ~= 1).*along;
s = (side == 1).*along + (side == 3).*((mu - 1) + (3 - 2*mu).*t);
bad = ~ismember(side, [1 2 3]);
r(:, :, bad) = nan;
s(:, :, bad) = nan;
rs_points = [r; s];                      % 2 x n_t x n
end
