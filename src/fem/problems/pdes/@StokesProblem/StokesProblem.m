classdef StokesProblem < FEProblem
%-------------------------------------------------------------------------
%Description: HDG discretization of the Stokes problem:
%             Find (L,u,p) in V^(2x2) x V^2 x V such that
%               L − ∇u = 0,   −∇·(νL − pI) = f,   ∇·u = 0   in Ω
%             Hybrid unknowns: velocity trace û and per-element pressure mean ρ.
%             The zero-mean pressure constraint (∫p = 0) is imposed as an
%             additional Lagrange multiplier via ZeroMeanPressure.
%-------------------------------------------------------------------------
    methods
        function obj = StokesProblem(k, label, options)
            arguments (Input)
                k       double          {mustBeScalarOrEmpty} = []                         % polynomial degree
                label   string          {mustBeScalarOrEmpty} = "Stokes (Velocity-Pressure)"  % display label
                options.viscosity       (1,1)   double        = 1                          % kinematic viscosity ν
                options.domain_measure  (1,1)   double        = 1                          % |Ω|, for normalising pressure mean
                options.source_term     (2,1)   cell          = {@(x,y) 0*(x.*y); @(x,y) 0*(x.*y)}  % body force {f_x, f_y}
                options.datum_dirichlet (2,1)   cell          = {@(x,y) 0*(x.*y); @(x,y) 0*(x.*y)}  % Dirichlet data {g_x, g_y}
                options.datum_neumann   (2,1)   cell          = {@(x,y) 0*(x.*y); @(x,y) 0*(x.*y)}  % Neumann data {q_x, q_y}
                % "constrained": the zero-mean pressure multiplier (all-Dirichlet
                % velocity); "pinned": a traction boundary fixes the level, no constraint
                options.pressure_mean   (1,1)   string {mustBeMember(options.pressure_mean, ["constrained", "pinned"])} = "constrained"
                % TPM: x on Gamma_h -> x_bar on Gamma (handle or LevelSetCurve); the
                % Dirichlet datum is then transferred along the paths (apply_transfer_conditions)
                options.transfer_mapping {mustBeA(options.transfer_mapping, ["function_handle", "LevelSetCurve"])} = function_handle.empty
            end
            obj = obj@FEProblem(label, ...
                Order=k, ...
                AdditionalConstraints = @ZeroMeanPressure, ...
                Data = options);
            obj.FieldsLocal = [struct(name="L_h",       description="velocity gradient");
                                 struct(name="u_h",       description="velocity");
                                 struct(name="p_h",       description="pressure")];
            obj.FieldsHybrid = [struct(name="\hat{u}_h", description="velocity trace");
                                 struct(name="\rho_h",    description="pressure mean")];
        end
    end

    methods (Static)
        spaces   = get_discrete_spaces(k)
        data = get_manufactured_solution(U, p, kinematic_viscosity, filename)
    end
end
