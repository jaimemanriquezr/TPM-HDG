classdef TensorFEBasis < AbstractFEBasis
    % Vector- or tensor-valued basis built from a scalar basis and a list
    % of constant linear maps: the functions are phi_i * A_c for every
    % scalar function phi_i and every map A_c.
    %
    %   basis = TensorFEBasis(scalar_basis, linear_maps)
    %
    % scalar_basis is a ScalarFEBasis; linear_maps is a cell array of
    % matrices of equal size (default {eye(2)}), otherwise an assertion
    % fails. Dimension = numel(linear_maps) * scalar_basis.Dimension and
    % GeometricOutputDimension is the size of the maps; Order,
    % ElementType, GeometricInputDimension, IsDiscontinuous and
    % IsPiolaTransformed are those of the scalar basis. Usually built by
    % ScalarFEBasis/to_vector, to_tensor or mtimes.
    %
    % See also: ScalarFEBasis, AbstractFEBasis
    properties (Hidden)
        ScalarBasis
        LinearMaps
    end

    properties (Dependent)
        Order
        Dimension
        GeometricOutputDimension
    end

    properties (Dependent, Hidden)
        ElementType
        GeometricInputDimension
        IsDiscontinuous
        IsPiolaTransformed
    end

    methods
        function obj = TensorFEBasis(basis, linear_maps)
            arguments
                basis ScalarFEBasis
                linear_maps cell = {eye(2)};
            end
            obj.ScalarBasis = basis;
            if numel(linear_maps) > 1
                linear_map_sizes = cellfun(@size, linear_maps, 'UniformOutput', 0);
                assert(isequal(linear_map_sizes{:}));
            end
            obj.LinearMaps = linear_maps;
        end
        varargout = evaluate(obj, points, mesh, options)

        function k = get.Order(obj)
            k = obj.ScalarBasis.Order;
        end

        function dim = get.GeometricOutputDimension(obj)
            dim = size(obj.LinearMaps{1});
        end

        function dim = get.Dimension(obj)
            dim = numel(obj.LinearMaps) * obj.ScalarBasis.Dimension;
        end

        function el_type = get.ElementType(obj)
            el_type = obj.ScalarBasis.ElementType;
        end

        function dim = get.GeometricInputDimension(obj)
            dim = obj.ScalarBasis.GeometricInputDimension;
        end

        function el_type = get.IsDiscontinuous(obj)
            el_type = obj.ScalarBasis.IsDiscontinuous;
        end

        function el_type = get.IsPiolaTransformed(obj)
            el_type = obj.ScalarBasis.IsPiolaTransformed;
        end
    end

    methods (Access = private)
        function [ids, stride, cols_ids] = scalar_enumeration(obj, mesh, kind)
            % the wrapped element enumerates itself; the component stride
            % is its total dof count on the mesh (entity count from the
            % support codimension)
            sdim = obj.ScalarBasis.Dimension;
            switch mesh.geometric_dimension - obj.ScalarBasis.GeometricInputDimension
                case 1
                    stride = sdim * size(mesh.faces, 2);
                otherwise
                    stride = sdim * size(mesh.elements, 2);
            end
            if kind == "global"
                ids = obj.ScalarBasis.get_global_dof(mesh);
                cols_ids = [];
            else
                [ids, cols_ids, ~] = obj.ScalarBasis.get_assembly_dof(mesh);
            end
        end
    end
end
