function dp = divergence(obj)
%DIVERGENCE  Divergence of a vector- or matrix-valued field.
%
%   dp = divergence(obj)
%
%   dp is an FEFunction whose evaluation at (x, mesh) gives, as for
%   FEField.evaluate, the divergence of obj's basis functions (one column
%   each) when obj has no dof (trial and test functions), and the
%   divergence of the field itself (a single column) when it has dof; in
%   the latter case mesh defaults to obj.mesh. Matrix-valued fields are
%   differentiated row-wise. For a PEERSBasis the bubble part has zero
%   divergence. Evaluating it for a scalar basis raises an error.
    dp = FEFunction(@(varargin) evaluate_divergence(obj, varargin{:}), ...
                    "Divergence-" + obj.signature);
end

function dP = evaluate_divergence(u, x, mesh)
    arguments
        u FEField
        x double
        mesh AbstractGrid = SimplicialGrid.empty;
    end
    if isempty(mesh) && ~isempty(u.dof)
        mesh = u.mesh;
    end
    if isa(u.basis, 'ScalarFEBasis')
        error('Divergence operator undefined for scalar functions.')
    elseif isa(u.basis, 'PEERSBasis')
        sub_basis = DubinerBasis(u.basis.Order).to_tensor();
        dp = p_transpose(sub_basis.ScalarBasis.evaluate_gradient(x, mesh));
        dp_row_A = @(A) pagemtimes(A, dp);
        dP = apply_fun(dp_row_A, sub_basis.LinearMaps); 
       
        dim_bubble = u.basis.Dimension - sub_basis.Dimension;
        div_curl_bubble = zeros(size(x, 2), dim_bubble, 2, 1, size(dP, 5));
        dP = [dP div_curl_bubble];
    else
        dp = p_transpose(u.basis.ScalarBasis.evaluate_gradient(x, mesh));
        if u.GeometricOutputDimension(2) == 1
            dp_dot_A = @(A) pagemtimes(dp, 'transpose', A, 'none');
            dP = apply_fun(dp_dot_A, u.basis.LinearMaps); 
        else
            dp_row_A = @(A) pagemtimes(A, dp);
            dP = apply_fun(dp_row_A, u.basis.LinearMaps);
        end
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