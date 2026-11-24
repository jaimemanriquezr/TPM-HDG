function p = evaluate_boundary(obj, points, mesh, options)
%EVALUATE_BOUNDARY  Evaluate a function on the faces of every element.
%
%   p = evaluate_boundary(u, points, mesh, output = "array")
%
%   points  1-by-n parameters on the reference face (one row fewer than
%           the geometric dimension of the mesh).
%   A function defined on elements is evaluated at the sides of each
%   element, parametrized as the reference element of mesh.ElementType
%   ("Simplex" or "Quad"); a function defined on faces is evaluated on
%   each face, with mesh.permutation selecting the orientation. The output option is
%   forwarded to evaluate.
    arguments
        obj         AbstractFEFunction
        points      double
        mesh        AbstractGrid = SimplicialGrid.empty
        options.output string = "array"
    end
    assert(size(points, 1) == (mesh.geometric_dimension - 1))
    if isa(obj, 'FEField')
        b_dim = obj.GeometricInputDimension;
    else
        [~, b_dim] = read_signature(obj);
    end
    op_args = namedargs2cell(options);
    switch b_dim
        case 1
            points_boundary = {points; 
                               1 - points};
            index_matrix = mesh.permutation.';
        case 2
            switch mesh.ElementType
                case "Simplex"
                    points_boundary = {[0*points; 1 - points]; 
                                       [1 - points; 0*points]; 
                                       [1 - points; points]};
                case "Quad"
                    points_boundary = {[points; 0*points]; 
                                       [1 + 0*points; points]; 
                                       [points; 1 + 0*points]; 
                                       [0*points; points]};
            end
            index_matrix = [];
    end
    s_boundary = cat(2, points_boundary{:});
    p = split(obj.evaluate(s_boundary, mesh, op_args{:}), ...
              length(points_boundary));
    if size(p, 6) == size(mesh.faces, 2)
        n_sides = size(mesh.boundaries, 1);
        p_1 = permute(reshape(p(1, :, :, :, :, mesh.boundaries), 1, size(p, 2), size(p, 3), 1, 1, n_sides, []), [6 7 2 3 1 4 5]);
        p_2 = permute(reshape(p(2, :, :, :, :, mesh.boundaries), 1, size(p, 2), size(p, 3), 1, 1, n_sides, []), [6 7 2 3 1 4 5]);
        p = (mesh.permutation == 1).*p_1 + (mesh.permutation==2).*p_2;
        p = permute(p, [3 4 5 6 2 1]);
    else
        p = index_into(p, index_matrix);
    end
end

function x = split(x, n)
    size_x = size(x);
    size_x = [size_x(1)/n n size_x(2:end)];
    x = reshape(x, size_x);
    
    x_dims = 1:length(size(x));
    x = permute(x, x_dims([2 1 3:end]));
end

function x = index_into(x, A)
    if isempty(A)
        % move the boundary-piece dim to the end; dimension-general so
        % tensor-valued evaluations (7-D) pass through
        nd = max(ndims(x), 6);
        x = permute(x, [2:nd 1]);
    else
        x_dims = 1:length(size(x));
        x = reshape(x(A, :), [size(A, 1) size(A, 2) size(x, x_dims(2:end))]);
        x = permute(x, [3 4 5 6 1 2]);
    end
end