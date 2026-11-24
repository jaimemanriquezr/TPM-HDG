function [rs_points, element_list] = get_2d_local(obj, xy_points, element_list)
%GET_2D_LOCAL  Reference coordinates of physical points.
%
%   [rs, elements] = get_2d_local(mesh, xy_points)
%   rs = get_2d_local(mesh, xy_points, element_list)
%
%   Inverts the element maps of get_affine_map. With an element list,
%   page e of xy_points (2-by-n-by-n_pages) is mapped with the element
%   element_list(e). Without one, xy_points is a single page (2-by-n) and
%   each point is located by its breakpoint intervals in x and y, returning
%   the element indices. Unlike SimplicialGrid, points on faces are found:
%   a point on an interior face goes to the element on its east/north side
%   (the last element in each direction keeps its closing face). Points
%   outside the grid get NaN entries.
%
%See also: get_affine_map, SimplicialGrid/get_2d_local
arguments
    obj CartesianGrid
    xy_points (2,:,:) double
    element_list (:,1) double = []
end
if isempty(element_list)
    assert(size(xy_points, 3) == 1)
    i = discretize(xy_points(1, :), obj.x);
    j = discretize(xy_points(2, :), obj.y);
    element_list = i + obj.Nx*(j - 1);          % NaN outside the grid
    found = ~isnan(element_list);
    rs_points = nan(size(xy_points));
    [F, x0] = get_affine_map(obj, element_list(found), inverse = true);
    rs_points(:, found) = reshape(pagemtimes(F, ...
        reshape(xy_points(:, found) - x0, 2, 1, [])), 2, []);
else
    assert(size(xy_points, 3) == length(element_list))
    [F, x0] = get_affine_map(obj, element_list, inverse = true);
    rs_points = pagemtimes(F, xy_points - reshape(x0, 2, 1, []));
end
end
