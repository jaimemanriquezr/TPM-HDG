function [rs_points, element_list] = get_2d_local(obj, xy_points, element_list)
%GET_2D_LOCAL  Reference coordinates of physical points.
%
%   [rs, elements] = get_2d_local(mesh, xy_points)
%   rs = get_2d_local(mesh, xy_points, element_list)
%
%   Inverts the element affine maps. With an element list, page e of
%   xy_points (2-by-n-by-n_pages) is mapped with the element
%   element_list(e). Without one, xy_points is a single page (2-by-n) and
%   each point is located in the element whose interior contains it
%   (strictly: points on faces and vertices are not found), returning the
%   element indices; the entries of unlocated points are NaN.
%
%See also: get_affine_map, get_physical_points
arguments
    obj SimplicialGrid
    xy_points (2,:,:) double
    element_list(:,1) double = []
end
if isempty(element_list)
    assert(size(xy_points, 3) == 1)
    [B, x3] = obj.get_affine_map(output="array");
    rs_points = nan(size(xy_points));
    element_list = nan(1, size(xy_points, 2));
    [local_points, valid_element, is_inside] = deal(cell(size(B, 3), 1));
    parfor K = 1:size(B, 3)
        l_coords = B(:,:,K) \ (xy_points - x3(:,K));

        r = l_coords(1, :);
        s = l_coords(2, :);
        is_interior = (r > 0) & (s > 0) & (r+s < 1);
        if any(is_interior)
            valid_point = find(is_interior);
            local_points{K} = l_coords(:, valid_point);
            is_inside{K} = valid_point;
            valid_element{K} = K * (valid_point.^0);
        else
            is_inside{K} = 0;
            valid_element{K} = 0;
        end
    end
    is_inside = horzcat(is_inside{:});
    valid_element = horzcat(valid_element{:});
    valid_element = valid_element(is_inside > 0);
    is_inside = is_inside(is_inside > 0);

    local_points = horzcat(local_points{:});
    rs_points(:, is_inside) = local_points;
    element_list(:, is_inside) = valid_element;
else
    assert(size(xy_points, 3) == length(element_list))
    [B, x3] = obj.get_affine_map(element_list, output="array");
    % one 2x2 solve per page, all elements at once
    rs_points = pagemldivide(B, xy_points - reshape(x3, 2, 1, []));
end
end
