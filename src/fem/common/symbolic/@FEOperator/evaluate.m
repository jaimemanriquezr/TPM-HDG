function F_u = evaluate(obj, u, points, mesh, options)
%EVALUATE  Apply a pointwise operator to an FE function at given points.
%
%   F_u = evaluate(obj, u, points, mesh, element_list = [], ...
%                  coordinate_system = "local", output = "array")
%
%   u is evaluated at the points and obj.apply_action is applied to the
%   result; a spatially dependent action also receives the physical
%   points. output = "cell" splits the value dimensions into a cell.
arguments
    obj FEOperator
    u {mustBeA(u, ["AbstractFEFunction", "AbstractFEBasis"])}
    points double
    mesh AbstractGrid = SimplicialGrid.empty;
    options.element_list = [];
    options.coordinate_system = "local";
    options.output = "array";
end
u_ = u.evaluate(points, mesh, ...
    element_list=options.element_list, ...
    coordinate_system=options.coordinate_system);

if isa(obj, 'IdentityOperator')
    F_u = u_;
elseif obj.is_spatial_dependent
    if options.coordinate_system == "physical"
        x_physical = points;
    else
        x_physical = mesh.get_2d_physical(points, options.element_list);
    end
    F_u = obj.apply_action(u_, x_physical);
else
    F_u = obj.apply_action(u_);
end

switch options.output
    case "array"
        return
    case "cell"
        out_sz = size(F_u, 3:4);
        F_u = mat2cell(F_u, size(F_u, 1), size(F_u, 2), ...
            ones(out_sz(1), 1), ones(out_sz(2), 1), size(F_u, 5));
end
end
