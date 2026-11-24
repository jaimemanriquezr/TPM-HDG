function dp = gradient(obj)
%GRADIENT  Gradient of a scalar- or vector-valued field.
%
%   dp = gradient(obj)
%
%   dp is an FEFunction whose evaluation at (x, mesh) gives, as for
%   FEField.evaluate, the gradient of obj's basis functions (one column
%   each) when obj has no dof (trial and test functions), and the
%   gradient of the field itself (a single column) when it has dof; in
%   the latter case mesh defaults to obj.mesh. Matrix-valued fields are
%   not supported: evaluating dp then prints a message and errors (its
%   output is not assigned).
    dp = FEFunction(@(varargin) evaluate_gradient(obj, varargin{:}), ...
                    "Gradient-" + obj.signature);
end

function dP = evaluate_gradient(u, x, mesh)
    arguments
        u FEField
        x double
        mesh AbstractGrid = SimplicialGrid.empty;
    end
    if isempty(mesh) && ~isempty(u.dof)
        mesh = u.mesh;
    end
    if isa(u.basis, 'ScalarFEBasis')
        dP = u.basis.evaluate_gradient(x, mesh);
    else
        if u.basis.GeometricOutputDimension(2) > 1
            disp("No gradient for tensors (yet)!")
            return
        end
        dp = p_transpose(u.basis.ScalarBasis.evaluate_gradient(x, mesh));
        A_times_dp = @(A) pagemtimes(A, 'none', dp, 'transpose');
        dP = apply_fun(A_times_dp, u.basis.LinearMaps);
    end
    dP = apply_dof(u, dP);
end

function dP = apply_dof(u, dP)
    % combine the per-basis values with the field's dof, as FEField.evaluate
    if ~isempty(u.dof)
        dP = sum(dP .* reshape(u.dof, 1, size(u.dof, 1), 1, 1, []), 2);
    end
end

function x = p_transpose(x)
    x_dim = 1:length(size(x));
    x = permute(x, [3 4 1 2 x_dim(5:end)]);
end

function fun_array = apply_fun(fun, cell_array)
    fun_cell = cellfun(@(x) p_transpose(fun(x)), cell_array, ...
                        'UniformOutput', false);
    fun_array = cat(2, fun_cell{:});
end