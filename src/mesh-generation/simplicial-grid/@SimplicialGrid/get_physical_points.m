function x_physical = get_physical_points(obj, x_local)
    %GET_PHYSICAL_POINTS  Physical images of reference points on all elements or faces.
    %
    %   x_physical = get_physical_points(mesh, x_local)
    %
    %   x_local  2-by-m reference points of the element, or 1-by-m
    %            parameters in [0, 1] along a face.
    %   x_physical  2-by-m-by-n_elements (2 rows) or 2-by-m-by-n_faces (1
    %            row), the images under the element affine maps or along
    %            the faces, from their first to their second vertex.
    %
    %See also: get_affine_map, get_1d_physical
    arguments
        obj SimplicialGrid
        x_local (:, :) double
    end
    
    input_dim = size(x_local, 1);
    switch input_dim
        case 2
            [B, x_0] = obj.get_affine_map(output="array");
            x_physical = pagemtimes(B, x_local) + reshape(x_0, 2, 1, []);
        case 1
            x = reshape(obj.coordinates(:, obj.faces(1:2, :)), 2, 2, []);
            B = diff(x, 1, 2);
            x_0 = x(:, 1, :);
            x_physical = B .* x_local + x_0;
    end
end