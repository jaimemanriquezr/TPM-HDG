classdef DubinerBasis < ScalarFEBasis
    % Discontinuous Dubiner polynomial basis of total degree up to k on
    % the triangle, Dimension = nchoosek(k+2, 2).
    %
    %   basis = DubinerBasis(k)
    %
    % k is the polynomial order (default empty, giving Dimension 0).
    % [P, dPx, dPy, dPxx, dPyy, dPxy] = DubinerBasis.evaluate_basis(k, r, s)
    % returns the values and first and second derivatives at the natural
    % coordinates (r, s); a scalar k selects all Dimension functions, a
    % vector k the functions with those indices (dubiner_eval).
    %
    % See also: ScalarFEBasis, LegendreBasis, dubiner_eval
    properties (Constant, Hidden)
        ElementType = "Simplex"
        GeometricInputDimension = 2
        GeometricOutputDimension = 1
        IsDiscontinuous = true
        IsPiolaTransformed = false
    end

    methods
        function obj = DubinerBasis(k)
            arguments
                k double = [];
            end
            obj.Order = k;
        end
        E = compute_l2_error(obj, dof, mesh, options)
    end

    methods (Static)
        function [P, dPx, dPy, dPxx, dPyy, dPxy] = evaluate_basis(k, r, s)
            arguments
                k double
                r double
                s double = r(2, :);
            end
            if nargin < 3
                r = r(1, :);
            end

            if isscalar(k)
                d2 = nchoosek(k+2, 2);
                k_list = 1:d2;
            else
                k_list = k;
            end

            if nargout > 3
                [P, dPx, dPy, dPxx, dPyy, dPxy] = dubiner_eval(k_list, r, s);
            elseif nargout > 1
                [P, dPx, dPy] = dubiner_eval(k_list, r, s);
            else
                P = dubiner_eval(k_list, r, s);
            end
        end
    end
end
