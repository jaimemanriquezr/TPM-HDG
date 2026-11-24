classdef LegendreBasis < ScalarFEBasis
    % Discontinuous Legendre polynomial basis of degree up to k on the
    % faces (one-dimensional support), Dimension = k + 1.
    %
    %   basis = LegendreBasis(k)
    %
    % k is the polynomial order (default empty, giving Dimension 0).
    % evaluate_basis(k, s) returns [P, dP], the values and derivatives
    % at the points s of the k + 1 basis functions (legendre_eval).
    %
    % See also: ScalarFEBasis, DubinerBasis
    properties (Constant, Hidden)
        ElementType = "Simplex"
        GeometricInputDimension = 1
        GeometricOutputDimension = 1
        IsDiscontinuous = true
        IsPiolaTransformed = false
    end

    methods
        function obj = LegendreBasis(k)
            arguments
                k double = [];
            end
            obj.Order = k;
        end
        E = compute_l2_error(obj, dof, mesh, options)

    end

    methods (Static)

        function [P, dP] = evaluate_basis(k, s)
            arguments
                k double
                s (:, 1) double
            end
            d1 = k + 1;
            [P, dP] = deal(zeros(length(s), d1));
            for j = 1:d1
                [P(:,j), dP(:,j)] = legendre_eval(j, s);
            end
        end
    end
end
