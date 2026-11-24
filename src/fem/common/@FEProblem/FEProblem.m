classdef FEProblem
    % Variational problem: an array of FESpaces (each declaring whether
    % its dofs are statically condensed), an AffineOperator holding the
    % equations, and BoundaryConditions. Two ways to provide the algebra
    % share one solve:
    %   - declarative: set Operator (symbolic forms and/or pre-assembled
    %     blocks); solve assembles, condenses and recovers.
    %   - operator overrides (low-level path): subclass and implement
    %     get_matrices / get_assemble_indices / get_boundary_conditions /
    %     get_local_solutions + static get_discrete_spaces, and assign
    %     the FieldsLocal/FieldsHybrid metadata.
    properties
        Label   string          = "FEProblem"
        Order   double          = []
        Data    struct          = struct()

        Mesh    AbstractGrid    = SimplicialGrid.empty
        Spaces  FESpace         = FESpace.empty
        Operator AffineOperator = AffineOperator.empty
        BoundaryConditions      = []

        FieldsLocal  struct     = struct.empty
        FieldsHybrid struct     = struct.empty
    end

    properties (Hidden)
        AdditionalConstraints   function_handle = @no_constraints
    end

    properties (Hidden, Dependent)
        % views of Spaces by condensation, and the id bookkeeping the
        % assembly partitions blocks with
        LocalSpaces
        GlobalSpaces
        Space
        Id
    end

    methods
        function obj = FEProblem(label, options)
            arguments
                label                          string  = "FEProblem"
                options.Order                  double  = []
                options.Mesh                   AbstractGrid = SimplicialGrid.empty
                options.Data                   struct  = struct()
                options.AdditionalConstraints  function_handle = @no_constraints
                options.Spaces                 FESpace = FESpace.empty
                options.Operator               AffineOperator = AffineOperator.empty
                options.BoundaryConditions             = []
            end
            obj.Label = label;
            for p = string(fieldnames(options)).'
                obj.(p) = options.(p);
            end
        end

        % --- operator overrides (low-level path; implement in subclasses) ---

        function varargout = get_matrices(obj, varargin) %#ok<STOUT>
            error('FEProblem:NotImplemented', ...
                '%s must implement get_matrices for the low-level path.', class(obj));
        end

        function varargout = get_assemble_indices(obj, varargin) %#ok<STOUT>
            error('FEProblem:NotImplemented', ...
                '%s must implement get_assemble_indices for the low-level path.', class(obj));
        end

        function [x_hat, dir_dof, g_hat] = get_boundary_conditions(obj, mesh, g_D, q_N) %#ok<INUSD,STOUT>
            error('FEProblem:NotImplemented', ...
                '%s must implement get_boundary_conditions for the low-level path.', class(obj));
        end

        function [volume_dof, hybrid_dof] = get_local_solutions(obj, x_hat, local_operators, mesh) %#ok<INUSD,STOUT>
            error('FEProblem:NotImplemented', ...
                '%s must implement get_local_solutions for the low-level path.', class(obj));
        end

        function n = get_ndof(pde, mesh, order)
            arguments (Input)
                pde   (1,1)   FEProblem
                mesh          AbstractGrid
                order (1,1)   double = pde.Order
            end
            [lhs_rows, ~, ~] = get_assemble_indices(pde, mesh, order);
            n = max(lhs_rows(:));
        end

        % --- dependent views ---

        function spaces = get.LocalSpaces(obj)
            spaces = obj.Spaces([obj.Spaces.Condense]);
        end

        function spaces = get.GlobalSpaces(obj)
            spaces = obj.Spaces(~[obj.Spaces.Condense]);
        end

        function space = get.Space(obj)
            % assembly order: condensed spaces first, then retained -
            % trial/test ids follow this ordering
            space = [obj.LocalSpaces(:); obj.GlobalSpaces(:)]';
        end

        function id = get.Id(obj)
            % The struct keeps the geometric names field/hybrid: they index
            % the solution side, and `global` is a MATLAB keyword.
            n_local  = length(obj.LocalSpaces);
            n_global = length(obj.GlobalSpaces);
            id = struct("field", 1:n_local, "hybrid", n_local + (1:n_global));
        end
    end
end

function [H, RHS] = no_constraints(~, ~, H, RHS); end
