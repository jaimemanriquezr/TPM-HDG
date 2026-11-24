function varargout = evaluate(obj, points, mesh, options)
%EVALUATE  Values (and derivatives) of a tensor basis at points.
%
%   [P, ...] = evaluate(obj, points, mesh, element_list=, coordinate_system=, output=)
%
%   Evaluates obj.ScalarBasis (same points, mesh and options) and
%   multiplies every output by each of the linear maps, concatenating
%   the results along the second dimension: component c occupies the
%   columns (c-1)*d+1 : c*d, d being the scalar Dimension. With
%   output = "array" (default) each output is
%   N-by-Dimension-by-(map size); with "cell" it is split into a cell
%   array by map entry. For cell array points the result is a cell array
%   of the same size.
    arguments
        obj TensorFEBasis
        points {mustBeNumericOrCell}
        mesh AbstractGrid = SimplicialGrid.empty
        options.element_list = [];
        options.coordinate_system = "local";
        options.output = "array";
    end
    basis_name = class(obj.ScalarBasis);
    basis_fcn = sprintf('%s>%s.evaluate_basis', basis_name, basis_name);
    p = cell(max(nargout(basis_fcn), 1), 1);
    [p{:}] = obj.ScalarBasis.evaluate(points, mesh, ...
                coordinate_system=options.coordinate_system, ...
                element_list=options.element_list);
    
    varargout = cell(1, numel(p));
    if isnumeric(points)
        for i = 1:nargout
            P = cellfun(@(A) p{i} .* reshape(A, [1 1 size(A)]), ...
                        obj.LinearMaps, ...
                        "UniformOutput", false);
            P = cat(2, P{:});
            switch options.output
                case "array"
                    varargout{i} = P;
                case "cell"
                    out_sz = num2cell(obj.GeometricOutputDimension);
                    out_sz = cellfun(@(s) ones(s, 1), out_sz, 'UniformOutput', false);
                    P = mat2cell(P, size(P, 1), size(P, 2), out_sz{:});
                    varargout{i} = squeeze(P);
            end
        end
    else
        for i = 1:nargout
            P = cell(numel(points), 1);
            for j = 1:numel(points)
                P{j} = cellfun(@(A) p{i}{j} .* reshape(A, [1 1 size(A)]), ...
                            obj.LinearMaps, ...
                            "UniformOutput", false);
                P{j} = cat(2, P{j}{:});
            end
            switch options.output
                case "array"
                    varargout{i} = P;
                case "cell"
                    out_sz = num2cell(obj.GeometricOutputDimension);
                    out_sz = cellfun(@(s) ones(s, 1), out_sz, 'UniformOutput', false);
                    for j = 1:numel(points)
                        P{j} = squeeze(mat2cell(P{j}, size(P{j}, 1), size(P{j}, 2), out_sz{:}));
                    end
                    varargout{i} = P;
            end
        end
    end
end