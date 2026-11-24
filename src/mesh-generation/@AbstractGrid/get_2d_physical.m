function xy_points = get_2d_physical(obj, rs_points, element_list)
    % Map local coordinates to physical coordinates. Two modes, mirroring
    % get_2d_local:
    %   no element_list: rs_points (2 x N), broadcast through the affine
    %       map of every element -> (2 x N x n_elements);
    %   element_list:    rs_points (2 x N x n_pages), page k mapped by the
    %       affine map of element_list(k) -> (2 x N x n_pages).
    arguments
        obj AbstractGrid
        rs_points (2, :, :) double
        element_list (:, 1) double = [];
    end
    if isempty(obj)
        xy_points = rs_points;
        return;
    end
    if isempty(element_list)
        assert(size(rs_points, 3) == 1)
        [F, xy_0] = get_affine_map(obj, output="array");
    else
        assert(size(rs_points, 3) == length(element_list))
        [F, xy_0] = get_affine_map(obj, element_list, output="array");
    end
    xy_points = pagemtimes(F, rs_points) + reshape(xy_0, 2, 1, []);
end
