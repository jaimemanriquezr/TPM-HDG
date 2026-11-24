classdef ConstantBasis < ScalarFEBasis
    % The one-dof constant space on the whole domain: the scalar-supported
    % global unknown of the support taxonomy (Lagrange multipliers, global
    % means). Evaluates to 1 at every point; one degree of freedom total.
    properties (Constant, Hidden)
        ElementType = "Simplex"
        GeometricInputDimension = 0
        GeometricOutputDimension = 1
        IsDiscontinuous = false
        IsPiolaTransformed = false
    end

    methods
        function obj = ConstantBasis()
            obj.Order = 0;   % Dimension = nchoosek(2, 2) = 1
        end

        function E = compute_l2_error(obj, dof, mesh, options) %#ok<INUSD,STOUT>
            error('ConstantBasis:NotImplemented', ...
                'A scalar-supported global unknown has no L2 error field.');
        end

        % scalar support overrides the geometric enumeration: one global
        % dof, seen by every element
        function dof_list = get_global_dof(~, ~)
            dof_list = 1;
        end

        function [rows, cols, rows_rhs] = get_assembly_dof(~, mesh)
            [rows, cols, rows_rhs] = deal(ones(1, 1, size(mesh.elements, 2)));
        end
    end

    methods (Static)
        function [P, dPx, dPy] = evaluate_basis(k, r, s)
            arguments
                k double %#ok<INUSA>
                r double
                s double = r %#ok<INUSA>
            end
            P = ones(numel(r), 1);
            [dPx, dPy] = deal(zeros(numel(r), 1));
        end
    end
end
