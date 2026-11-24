classdef FEField < AbstractFEFunction
    %FEFIELD  A discrete field: an FE basis and its degrees of freedom.
    %
    %   u = FEField(basis, dof) builds the field with the coefficients dof
    %   in the AbstractFEBasis basis. basis may also be an FESpace, in which
    %   case its Basis and Mesh are taken. Both arguments are optional
    %   (defaults: an empty DubinerBasis and empty dof), and FEField(basis)
    %   without dof is the placeholder used for trial and test functions
    %   (TrialFunction(u), TestFunction(u)).
    %
    %   Stored dof are one column per mesh entity (per face or element, all
    %   components contiguous). A global dof vector of a tensor basis is
    %   taken in component-major numbering and converted on assignment.
    %
    %   Methods:
    %     evaluate(u, points, mesh, ...)     values at points
    %     plot(u, mesh, ...)                 plot the field on a mesh
    %     compute_l2_error(u, mesh, exact)   L2 error against a Datum
    %     gradient(u), divergence(u)         FEFunction of the derivative
    %     norm(u)                            L2 norm
    %     u + v, u - v, -u, +u               arithmetic on fields
    %
    %See also: Datum, FESolution, AbstractFEBasis
    properties
        basis AbstractFEBasis
    end

    properties (Dependent)
        dof         % double, Dimension-by-n_entities
    end

    properties (Access = protected)
        % the stored coefficients, one column per mesh entity (set.dof
        % converts global vectors into this layout)
        StoredDof
    end

    properties (Hidden)
        name string = string.empty;
        mesh AbstractGrid = SimplicialGrid.empty;
    end

    properties (Dependent, Hidden)
        Dimension
        GeometricInputDimension
        GeometricOutputDimension
        order
        id
        local_dof
    end

    methods
        function obj = FEField(basis, dof)
            arguments
                basis = DubinerBasis.empty;
                dof = [];
            end
            if isa(basis, "FESpace")
                obj.basis = basis.Basis;
                obj.mesh = basis.Mesh;
            else
                obj.basis = basis;
            end
            obj.dof = dof;
            obj.signature = "FEField";

            if isempty(obj.basis)
                return
            else
                if ~isempty(obj.basis.id)
                    obj.signature = obj.signature + "-" + obj.basis.id;
                end
                obj.signature = obj.signature ...
                    + ":" + obj.basis.GeometricInputDimension ...
                    + ":" + obj.basis.Order;
            end
        end

        function dof = get.dof(obj)
            dof = obj.StoredDof;
        end

        function obj = set.dof(obj, new_dof)
            space = obj.basis;
            % Stored layout is per-entity columns (per face/element, all
            % components contiguous). Global dof vectors of tensor bases
            % arrive in component-major numbering ([all comp_1; all
            % comp_2; ...], stride = the scalar count over the mesh --
            % the enumeration convention), so convert on the way in.
            if isempty(new_dof)
                obj.StoredDof = new_dof;
            elseif isa(space, 'TensorFEBasis') ...
                    && isvector(new_dof) ...
                    && numel(new_dof) > space.Dimension
                pdim = space.ScalarBasis.Dimension;
                n_comp = numel(space.LinearMaps);
                v = reshape(new_dof, pdim, [], n_comp);
                obj.StoredDof = reshape(permute(v, [1 3 2]), pdim*n_comp, []);
            else
                obj.StoredDof = reshape(new_dof, space.Dimension, []);
            end
        end

        function dof = get.local_dof(obj)
            n_elems = size(obj.mesh.elements, 2);
            dof = reshape(obj.dof(:, obj.mesh.boundaries), [], 1, n_elems);
        end

        function id = get.id(obj)
            id = obj.basis.id;
        end

        function dim = get.Dimension(obj)
            dim = obj.basis.Dimension;
        end

        function dim = get.order(obj)
            dim = obj.basis.Order;
        end

        function dim = get.GeometricInputDimension(obj)
            dim = obj.basis.GeometricInputDimension;
        end

        function dim = get.GeometricOutputDimension(obj)
            dim = obj.basis.GeometricOutputDimension;
        end

        minus_u = uminus(u, x);
        plus_u = uplus(u, x);
        u_plus_x = plus(u, x);
        u_minus_x = minus(u, x);

        function norm_u = norm(u, p)
            arguments
                u FEField
                p string = "L2";
            end
            switch p
                case "L2"
                    norm_u = compute_l2_error(u, u.mesh);
            end
        end

        function u = TrialFunction(obj)
            arguments
                obj FEField
            end
            u = FEField(obj.basis);
        end

        function v = TestFunction(obj)
            arguments
                obj FEField
            end
            v = FEField(obj.basis);
        end
    end
end
