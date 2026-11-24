function f = evaluate(obj, points, mesh, options)
%EVALUATE  Values of a Datum at points.
%
%   f = evaluate(obj, points, mesh, coordinate_system=, output=)
%
%   points are reference points mapped to physical points through mesh
%   (coordinate_system = "local", the default) or physical points
%   (coordinate_system = "physical", in which case mesh is not needed).
%
%   With output = "array" (default) f is n_points-by-1-by-m-by-n-by-n_el,
%   where m-by-n is the value size of the Datum; with output = "cell" it
%   is the cell array returned by the Datum's get_value.
    arguments
        obj         Datum
        points      {mustBeNumericOrCell}
        mesh        AbstractGrid = SimplicialGrid.empty;
        options.coordinate_system string = "local"
        options.output string = "array"
        options.element_list = [];      % accepted for uniformity; unused
    end
    switch options.coordinate_system
        case "physical"
            x_physical = points;
        case "local"
            x_physical = get_physical_points(mesh, points);
    end
    f = obj.get_value(x_physical);
    switch options.output
        case "array"
            f = cat(2, f{:});
            f = reshape(f, size(f, 1), 1, value_size(obj, 1), value_size(obj, 2), []);
        case "cell"
            return
    end
end