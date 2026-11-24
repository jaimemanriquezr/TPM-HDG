classdef (Abstract) ScalarFEBasis < AbstractFEBasis
    % Abstract base of the scalar-valued bases (DubinerBasis,
    % LegendreBasis, LagrangeBasis, ConstantBasis). Subclasses implement
    % the static evaluate_basis and compute_l2_error.
    %
    % Dimension = nchoosek(Order + GeometricInputDimension,
    % GeometricInputDimension), or 0 when Order is empty.
    %
    %   A_basis = basis * A    TensorFEBasis(basis, {A}) for a double A
    %                          (also A * basis, which errors for a
    %                          non-double A)
    %   to_vector(basis, d)    TensorFEBasis with the d unit-vector maps
    %                          (d defaults to GeometricInputDimension)
    %   to_tensor(basis, maps) TensorFEBasis with the maps given as a
    %                          double (one map) or a cell; by default the
    %                          d^2 matrices I(:,j)*I(:,i)'
    %
    % See also: AbstractFEBasis, TensorFEBasis
    properties (Abstract, Constant, Hidden)
        ElementType string
        GeometricInputDimension double
        GeometricOutputDimension double
        IsPiolaTransformed logical
        IsDiscontinuous logical
    end

    %--------------------------------------------------------------------------%
    methods (Abstract)
        error = compute_l2_error(obj, dof, mesh, options);
    end

    methods (Abstract, Static)
        varargout = evaluate_basis(varargin)
    end
    %--------------------------------------------------------------------------%

    properties
        Order double
    end

    properties (Dependent)
        Dimension
    end

    methods
        function dim_space = get.Dimension(obj)
            dim_space = return_Dimension(obj);
        end

        function A_basis = mtimes(basis, A)
            arguments
                basis
                A
            end
            if ~isa(basis, "ScalarFEBasis")
                [A, basis] = deal(basis, A);
                if ~isa(A, "double")
                    error("Invalid argument. Must be double.");
                end
            end
            A_basis = TensorFEBasis(basis, {A});
        end

        function vector_basis = to_vector(obj, d)
            arguments
                obj ScalarFEBasis
                d = obj.GeometricInputDimension;
            end
            A = mat2cell(eye(d), d, ones(1, d)).';
            vector_basis = TensorFEBasis(obj, A);
        end

        function tensor_basis = to_tensor(obj, linear_maps)
            arguments
                obj ScalarFEBasis
                linear_maps = []
            end
            if isempty(linear_maps)
                d = obj.GeometricInputDimension;
                I = eye(d);
                [R, C] = meshgrid(1:d);
                A = arrayfun(@(i, j) I(:,j)*I(:,i)', R(:), C(:), 'UniformOutput', 0);
            elseif isa(linear_maps, "double")
                A = {linear_maps};
            else
                A = linear_maps;
            end
            tensor_basis = TensorFEBasis(obj, A);
        end
    end

    methods (Access = protected)
        function dim_space = return_Dimension(obj)
            if isempty(obj.Order)
                dim_space = 0;
            else
                dim_space = nchoosek(obj.Order + obj.GeometricInputDimension, obj.GeometricInputDimension);
            end
        end
    end
end
