classdef LagrangeBasis < ScalarFEBasis
    % Continuous Lagrange basis on the triangle, evaluated with
    % lagrange_eval (only the linear nodal functions are implemented).
    %
    %   basis = LagrangeBasis(k)
    %
    % k is the polynomial order (default empty, giving Dimension 0).
    % evaluate_basis(k, r, s) returns [P, dPx, dPy], the values and first
    % derivatives at the natural coordinates (r, s).
    %
    % See also: ScalarFEBasis, DubinerBasis, lagrange_eval
    properties (Constant, Hidden)
        ElementType = "Simplex"
        GeometricInputDimension = 2
        GeometricOutputDimension = 1
        IsDiscontinuous = false
        IsPiolaTransformed = false
    end
    
    methods
        function obj = LagrangeBasis(k)
            arguments
                k double = [];
            end
            obj.Order = k;
        end

        function E = compute_l2_error(obj, dof, mesh, options)
            % Element-based L2 error of a field with nodal values dof (one
            % column per element), as DubinerBasis.compute_l2_error.
            arguments
                obj LagrangeBasis
                dof double
                mesh SimplicialGrid
                options.exact_solution function_handle = @(x,y) 0*(x.*y);
            end
            [xy, rs, w] = quadpoints_elements_2d(mesh, obj.Order + 1, output="array");
            P = obj.evaluate(rs) * dof;
            P_ex = reshape(options.exact_solution(xy(1,:,:), xy(2,:,:)), size(P));
            E = sqrt(sum(mesh.volume .* (w'*(P - P_ex).^2)));
        end
    end

    methods (Static)
        function [P, dPx, dPy] = evaluate_basis(k, r, s)
            arguments
                k double
                r double
                s double = r(2, :);
            end
            if nargin < 3
                r = r(1, :);
            end

            if nargout > 1
                [P, dPx, dPy] = lagrange_eval(k, r, s);
            else
                P = lagrange_eval(k, r, s);
            end
        end
    end
end