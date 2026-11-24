function varargout = evaluate(obj, points, mesh, options)
%EVALUATE  Values (and derivatives) of a scalar basis at points.
%
%   [P, ...] = evaluate(obj, points, mesh, coordinate_system=, output=, element_list=)
%
%   Calls the static evaluate_basis of the concrete class, so the number of
%   outputs is that of evaluate_basis (values, then derivatives).
%   points is a GeometricInputDimension-by-N array, or a cell array of
%   such arrays, in which case every output is a cell array of the same
%   size. coordinate_system = "local" (default) takes reference points;
%   "physical" maps them to reference points through mesh, which is then
%   required. A basis with GeometricInputDimension 0 is evaluated at the
%   point count only. For Piola-transformed bases the inverse affine
%   map of mesh is passed to evaluate_basis.
    arguments
        obj         ScalarFEBasis
        points      {mustBeNumericOrCell}
        mesh        AbstractGrid = SimplicialGrid.empty
        options.coordinate_system string = "local"
        options.output string = "array"
        options.element_list = [];
    end
    basis_name = class(obj);
    basis_fcn = sprintf('%s>%s.evaluate_basis', basis_name, basis_name);
    varargout = cell(1, max(nargout(basis_fcn), 1));
    
    if isnumeric(points)
        if obj.GeometricInputDimension == 0
            % scalar-supported basis: constant in space, evaluate at the
            % point count only (any point-array shape is accepted)
            [varargout{:}] = obj.evaluate_basis(obj.Order, zeros(size(points, 2), 1));
            return
        end
        assert(size(points, 1) == obj.GeometricInputDimension)
        switch options.coordinate_system
            case "physical"
                assert(~isempty(mesh), "No mesh specified!");
                switch obj.GeometricInputDimension
                    case 2
                        local_points = mesh.get_2d_local(points, options.element_list);
                end
            otherwise % "local" or "vertex"
                local_points = points;
        end
        args = mat2cell(local_points, ones(obj.GeometricInputDimension, 1));
        if obj.IsPiolaTransformed
            args{end+1} = get_affine_map(mesh, inverse=true);
        end
        [varargout{:}] = obj.evaluate_basis(obj.Order, args{:});
    else
        op_args = namedargs2cell(options);
        var_out = cell(size(points));
        for i = 1:numel(points)
            var_out{i} = varargout;
            [var_out{i}{:}] = obj.evaluate(points{i}, mesh, op_args{:});
        end
        
        for j = 1:length(varargout)
            varargout{j} = cellfun(@(C) C{j}, var_out, 'UniformOutput', false);
        end
    end
end