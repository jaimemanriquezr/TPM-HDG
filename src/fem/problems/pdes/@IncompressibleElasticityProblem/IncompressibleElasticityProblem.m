classdef IncompressibleElasticityProblem < FEProblem
%-------------------------------------------------------------------------
%Description: HDG discretization of the incompressible elasticity problem
%             in stress-velocity form (PEERS element):
%             Find (σ,u,γ) in S^(2x2) x V^2 x A^(2x2) such that
%               (1/2μ)dev(σ) − ∇u + γ = 0,   −∇·σ = f,   ∇·u = 0   in Ω
%             Hybrid unknowns: velocity trace û and per-element pressure mean ρ.
%             The zero-mean pressure constraint is imposed via ZeroMeanPressure.
%             Supports steady (density=0) and transient (density>0) problems.
%
%Time stepping: the θ-scheme on the momentum equation (time_scheme
%             "euler" -> θ = 1, "cranknicolson" -> θ = 1/2, or an explicit
%             theta; see get_theta), written as backward Euler with
%             Δt -> θΔt plus the history ((1-θ)/θ) ρ M a⁰. The history of a
%             step is Data.state, the struct built at t = 0 by initial_state
%             and advanced by get_kinematic_state:
%                 problem.Data.source_term = f(0);
%                 state = initial_state(problem, mesh, k);
%                 for n = 1:n_steps
%                     problem.Data.source_term = f(n*dt);   % and the BC data
%                     problem.Data.state = state;
%                     solution = solve(problem, mesh, k);
%                     state = get_kinematic_state(problem, solution, state);
%                 end
%             state.fluid holds the volume dof arrays (2*d2 x N) velocity
%             and acceleration, the skeleton array velocity_skeleton
%             (2*d1 x F), pressure (d2 x N) and stress (PEERS dofs). The
%             acceleration is u_t (Eulerian, on the fixed fluid domain).
%             [] (or the scalar 0) means zero; an empty state means all
%             histories are zero: build the state with initial_state to
%             start from the initial data. initial_velocity (handles) and
%             datum_dirichlet_rate are read only by initial_state, which
%             solves the acceleration problem for (u_t, p, σ) at t = 0 only
%             when θ < 1 or with_pressure = true (otherwise acceleration,
%             pressure and stress are []).
%-------------------------------------------------------------------------
    methods
        function obj = IncompressibleElasticityProblem(k, label, options)
            arguments (Input)
                k       double          {mustBeScalarOrEmpty} = []                           % polynomial degree
                label   string          {mustBeScalarOrEmpty} = "Stokes (Stress-Velocity)"  % display label

                options.density         (1,1)   double = 0.0   % ρ (0 → steady; >0 → transient)
                options.time            (1,1)   double = 0.0   % initial time t0: seeds state.time (initial_state, get_kinematic_state); not read by the solve
                options.time_step       (1,1)   double = 1e-1  % time step Δt
                options.initial_velocity (2,1)  cell   = {@(x,y) 0*(x.*y); @(x,y) 0*(x.*y)}  % u(0) {u_x, u_y} (initial_state)
                % theta-method for the momentum equation: time_scheme is a
                % preset ("euler" -> θ=1, "cranknicolson" -> θ=1/2); an
                % explicit finite theta overrides it (see get_theta).
                options.time_scheme     (1,1)   string {mustBeMember(options.time_scheme, ["euler", "cranknicolson"])} = "euler"
                options.theta           (1,1)   double = NaN   % θ ∈ (0,1]; NaN -> use time_scheme preset
                options.state                   struct = struct.empty  % step history (initial_state / get_kinematic_state)
                options.viscosity       (1,1)   double = 1e-3  % shear viscosity μ
                options.domain_measure  (1,1)   double = 1.0   % |Ω|, for normalising pressure mean
                % Pressure-level closure. "constrained": the zero-mean
                % pressure constraint is enforced (required under
                % all-Dirichlet velocity BCs). "pinned": the boundary
                % conditions determine the level and no constraint is
                % applied.
                options.pressure_mean   (1,1)   string {mustBeMember(options.pressure_mean, ["constrained", "pinned"])} = "constrained"
                options.source_term     (2,1)   cell   = {@(x,y) 0*(x.*y); @(x,y) 0*(x.*y)}  % body force {f_x, f_y}
                options.datum_dirichlet (2,1)   cell   = {@(x,y) 0*(x.*y); @(x,y) 0*(x.*y)}  % Dirichlet data {g_x, g_y}
                options.datum_neumann   (2,1)   cell   = {@(x,y) 0*(x.*y); @(x,y) 0*(x.*y)}  % Neumann data {q_x, q_y}
                % ∂_t g of the Dirichlet data at Data.time: the essential
                % datum of the acceleration problem (initial_state)
                options.datum_dirichlet_rate (2,1) cell = {@(x,y) 0*(x.*y); @(x,y) 0*(x.*y)}
                % TPM: x on Gamma_h -> x_bar on Gamma (handle or LevelSetCurve); the
                % Dirichlet datum is then transferred along the paths (apply_transfer_conditions)
                options.transfer_mapping {mustBeA(options.transfer_mapping, ["function_handle", "LevelSetCurve"])} = function_handle.empty
                % false: keep S = mu*(mesh tau) under convection; the caller sets tau(beta) on the mesh
                options.convective_stabilization (1,1) logical = true
            end
            obj = obj@FEProblem(label, ...
                Order=k, ...
                AdditionalConstraints = @ZeroMeanPressure, ...
                Data = options);
            obj.FieldsLocal = [struct(name="\sigma_h", description="shear stress");
                                 struct(name="u_h",      description="velocity");
                                 struct(name="\gamma_h", description="rotation")];
            obj.FieldsHybrid = [struct(name="\hat{u}_h", description="velocity trace");
                                 struct(name="\rho_h",    description="pressure mean")];
        end
    end

    methods (Static)
        spaces   = get_discrete_spaces(k);
        data = get_manufactured_solution(U, p, shear_viscosity, rho, time, filename);
    end
end
