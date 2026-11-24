classdef FESpace
    % Finite element space: mesh x element (x continuity), in the Ciarlet
    % sense. The Basis is the reference triple (domain implicit in its
    % GeometricInputDimension, P and the functionals fixed by the basis
    % choice); the space adds the mesh, the gluing rule (Continuity;
    % currently only "discontinuous") and the global dof enumeration.
    %
    % Global enumeration (get_assembly_dof, get_global_dof,
    % global_dof_layout) is a space-level operation; the basis-tier
    % get_global_dof / get_assembly_dof supply each entity's dof ordering.
    %
    % The Mesh may be empty at construction (a space template); solve
    % binds it before assembly, and enumeration requires a bound space.
    properties
        Basis     AbstractFEBasis
        Mesh        AbstractGrid = SimplicialGrid.empty
        Continuity  string {mustBeMember(Continuity, "discontinuous")} = "discontinuous"
        % Statically condense this space's dofs out of the global system
        % (declared by the variational problem, never inferred)
        Condense    logical = false
        % Local coefficient directions of this (condensed) space that stay global
        % unknowns instead of being eliminated: a cell of index vectors into the
        % space's own dofs, each meaning the 0/1 combination v = sum_i e_i (the
        % first index is the pivot). The coefficient along v becomes one
        % element-supported global unknown, laid out after the problem's global
        % spaces (condense_global_system). Only meaningful with Condense = true.
        UncondensedModes cell = {}
    end

    properties (Dependent)
        Order
        Dimension
        GeometricInputDimension
        IsPiolaTransformed
        Support
    end

    methods
        function obj = FESpace(basis, mesh, options)
            % FESpace(basis, mesh?) - a basis array yields an FESpace
            % array (one space per element, shared mesh).
            arguments
                basis   AbstractFEBasis  = AbstractFEBasis.empty
                mesh    AbstractGrid = SimplicialGrid.empty
                options.Continuity string = "discontinuous"
                options.Condense logical = false
                options.UncondensedModes cell = {}
            end
            if numel(basis) > 1
                assert(isempty(options.UncondensedModes), "FESpace:uncondensedModes", ...
                    "set UncondensedModes on the element of the array it belongs to.");
                for i = numel(basis):-1:1
                    obj(i) = FESpace(basis(i), mesh, ...
                        Continuity=options.Continuity, Condense=options.Condense);
                end
                obj = reshape(obj, size(basis));
                return
            end
            obj.Basis = basis;
            obj.Mesh = mesh;
            obj.Continuity = options.Continuity;
            obj.Condense = options.Condense;
            obj.UncondensedModes = options.UncondensedModes;
        end

        function L = global_dof_layout(spaces, mesh)
            % Block layout of the global unknown vector plus the
            % per-element scatter/gather map. Footprint row order is
            % d-fastest within a face slot, face slots stacked.
            arguments
                spaces FESpace
                mesh   AbstractGrid
            end
            n = numel(spaces);
            n_el = size(mesh.elements, 2);
            L.support = [spaces.Support];
            [L.count, L.offset] = deal(zeros(1, n));
            L.footprint = cell(1, n);
            off = 0;
            for i = 1:n
                el = spaces(i).Basis;
                d = el.Dimension;
                % footprints derive from each element's own enumeration
                % (per-face contiguous for scalar traces, component-major
                % for tensor elements)
                dofA = el.get_global_dof(mesh);  % per-entity columns
                switch L.support(i)
                    case "face"
                        L.footprint{i} = reshape(dofA(:, mesh.boundaries), [], n_el) + off;
                        L.count(i) = d * size(mesh.faces, 2);
                    case "element"
                        L.footprint{i} = dofA + off;
                        L.count(i) = d * n_el;
                    case "scalar"
                        L.footprint{i} = ones(1, n_el) + off;
                        L.count(i) = d;
                end
                L.offset(i) = off;
                off = off + L.count(i);
            end
            L.n_dof = off;
        end

        % --- DSL function factories (delegate to the elements) ---

        function varargout = get_trial_functions(spaces)
            [varargout{1:max(nargout, 1)}] = get_trial_functions([spaces.Basis]);
        end

        function varargout = get_test_functions(spaces)
            [varargout{1:max(nargout, 1)}] = get_test_functions([spaces.Basis]);
        end

        function varargout = get_field_functions(spaces, id)
            arguments
                spaces FESpace
            end
            arguments (Repeating)
                id double
            end
            [varargout{1:max(nargout, 1)}] = get_field_functions([spaces.Basis], id{:});
        end

        function k = get.Order(obj)
            k = obj.Basis.Order;
        end
        function d = get.Dimension(obj)
            d = obj.Basis.Dimension;
        end
        function d = get.GeometricInputDimension(obj)
            d = obj.Basis.GeometricInputDimension;
        end
        function b = get.IsPiolaTransformed(obj)
            b = obj.Basis.IsPiolaTransformed;
        end
        function s = get.Support(obj)
            s = obj.Basis.Support;
        end
    end
end
