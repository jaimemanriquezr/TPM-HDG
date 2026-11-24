function p = evaluate(obj, points, mesh, options)
%EVALUATE  Values of a discrete field at points.
%
%   p = evaluate(obj, points, mesh, element_list=, coordinate_system=, output=)
%
%   The options are passed to obj.basis.evaluate. If obj has no dof, p is
%   the basis evaluation itself. Otherwise p contracts it with the dof:
%   n_points-by-1-by-(value size)-by-n_el, with one page per column of
%   obj.dof, or per entry of element_list when that is given.
    arguments
        obj FEField
        points double
        mesh AbstractGrid = SimplicialGrid.empty;
        options.element_list = [];
        options.coordinate_system = "local";
        options.output = "array";
    end
    op_args = namedargs2cell(options);
    p = obj.basis.evaluate(points, mesh, op_args{:});
    if isempty(obj.dof)
        return;
    end
    dim_dof = obj.basis.Dimension;

    dofs = reshape(obj.dof, dim_dof, 1, 1, 1, []);
    out_dim = obj.basis.GeometricOutputDimension;
    if isscalar(out_dim) && out_dim == 1
        out_dim = {1 1};
    else
        out_dim = num2cell(reshape(out_dim, 1, []));
    end
    n_points = size(points, 2);
    if isempty(options.element_list)
        n_elements = size(obj.dof, 2);
        p = reshape(pagemtimes(p, dofs), n_points, 1, out_dim{:}, n_elements, []);
    else
        dofs = dofs(:, :, :, :, options.element_list);
        dofs = reshape(dofs, size(dofs, 1), 1, []);
        p = reshape(pagemtimes(p, dofs), n_points, 1, out_dim{:}, []);
    end
end 