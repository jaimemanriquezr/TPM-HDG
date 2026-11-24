classdef LinearElasticityProblem < FEProblem
    %-------------------------------------------------------------------------
    %Description: HDG discretization of the linear elasticity problem (PEERS element):
    %             Find (σ,u,γ) in S^(2x2) x V^2 x A^(2x2) such that
    %               ζ⁻¹σ − ε(u) + γ = 0,   −∇·σ = f   in Ω
    %             where ζ is the compliance tensor (inverse of Hooke's law).
    %             Hybrid unknown: displacement trace û on mesh skeleton.
    %             Supports static (density=0) and transient (density>0)
    %             problems, with ρ ü added to the momentum equation.
    %
    %Time stepping: time_scheme = "newmark" (Newmark-β, parameters
    %             newmark_beta, newmark_gamma) or "onestep" (central
    %             difference; its first step, selected by a state whose
    %             displacement_previous is unset, is the Taylor start
    %             from (u₀, u̇₀)). The history of a step is
    %             Data.state, the struct built at t = 0 by initial_state and
    %             advanced by get_kinematic_state:
    %                 problem.Data.source_term = f(0);
    %                 state = initial_state(problem, mesh, k);
    %                 for n = 1:n_steps
    %                     problem.Data.source_term = f(n*dt);   % and the BC data
    %                     problem.Data.state = state;
    %                     solution = solve(problem, mesh, k);
    %                     state = get_kinematic_state(problem, solution, state);
    %                 end
    %             state.solid holds the
    %             volume dof arrays (2*d2 x N) displacement, velocity,
    %             acceleration, displacement_previous, the skeleton arrays
    %             (2*d1 x F) displacement_skeleton, velocity_skeleton,
    %             acceleration_skeleton, and the stress/rotation dofs of the
    %             displacement, velocity and acceleration (stress, rotation,
    %             stress_velocity, ..., rotation_acceleration). [] (or the
    %             scalar 0) means zero, except that displacement_previous = []
    %             selects the onestep Taylor start (march from a zero u⁻¹
    %             with 0 or a zeros array). With density > 0 and no state, every
    %             history is zero and the solve is the steady problem
    %             with the reaction ρ/Δt² u. The initial_displacement and
    %             initial_velocity handles are read only by initial_state: a
    %             time-dependent solve without state starts from rest
    %             whatever they hold.
    %
    %Constitutive source: an optional field Data.constitutive_source
    %             (PEERS dofs per element, (4 d2 + d2*) x N) holds the
    %             moments (r, τ)_K of a right-hand side r of the
    %             constitutive row, ζ⁻¹σ − ε(u) + γ = r; absent or empty,
    %             r = 0.
    %-------------------------------------------------------------------------
    methods
        function obj = LinearElasticityProblem(k, label, options)
            arguments (Input)
                k       double          {mustBeScalarOrEmpty} = []                  % polynomial degree
                label   string          {mustBeScalarOrEmpty} = "Linear Elasticity" % display label

                options.density                         double = 0.0    % ρ (0 → static; >0 → transient)
                options.time                (1,1)       double = 0.0    % current time t
                options.time_step           (1,1)       double = 1e-1   % Newmark time step Δt
                options.initial_displacement (2,1)      cell   = {@(x,y) 0*(x.*y); @(x,y) 0*(x.*y)}  % u(0) {u_x, u_y} (initial_state)
                options.initial_velocity     (2,1)      cell   = {@(x,y) 0*(x.*y); @(x,y) 0*(x.*y)}  % u̇(0) {v_x, v_y} (initial_state)
                options.time_scheme         (1,1)       string {mustBeMember(options.time_scheme, ["onestep", "newmark"])} = "onestep"
                options.newmark_beta        (1,1)       double = 0.25   % β (Newmark)
                options.newmark_gamma       (1,1)       double = 0.5    % γ (Newmark)
                options.state                           struct = struct.empty  % step history (initial_state / get_kinematic_state)
                options.hooke_modulus       (1,1)       double = 0.0    % b₀ (Winkler / Robin stiffness)
                % Condense the PEERS solid with its constant isotropic stress
                % coefficient kept global (FESpace.UncondensedModes on the stress
                % space): the same discrete solution, without the lambda-sized
                % element penalties in the condensed system (read by
                % FluidStructureCoupling.interface_problem; solve(problem, mesh, k)
                % ignores it)
                options.factored_condensation (1,1)     logical = false
                options.young_modulus       (1,1)       double = 1e4    % E
                options.poisson_ratio       (1,1)       double = .1     % ν
                options.source_term         (2,1)       cell   = {@(x,y) 0*(x.*y); @(x,y) 0*(x.*y)}  % body force {f_x, f_y}
                options.datum_dirichlet     (2,1)       cell   = {@(x,y) 0*(x.*y); @(x,y) 0*(x.*y)}  % Dirichlet data {g_x, g_y}
                options.datum_neumann       (2,1)       cell   = {@(x,y) 0*(x.*y); @(x,y) 0*(x.*y)}  % Neumann data {q_x, q_y}
            end
            obj = obj@FEProblem(label, Order=k, Data=options);
            obj.FieldsLocal = [struct(name="\sigma_h", description="stress");
                struct(name="u_h",      description="displacement");
                struct(name="\gamma_h", description="rotation")];
            obj.FieldsHybrid = [struct(name="\hat{u}_h", description="displacement trace")];
        end
    end

    methods (Static)
        spaces   = get_discrete_spaces(k);
        data = get_manufactured_solution(u, young_modulus, poisson_ratio, hooke_modulus, density, time, filename)

        function [lambda, mu] = get_lame_parameters(E, nu)
            % Converts Young's modulus E and Poisson ratio ν to Lamé parameters λ and μ.
            lambda = E*nu/((1 + nu)*(1 - 2*nu));
            mu = E/2/(1 + nu);
        end

        function ix = isotropic_stress_dofs(k)
            % PEERS local dofs of the constant isotropic stress sigma = phi_1 I:
            % the constant sigma_11 and sigma_22 Dubiner functions, [1, 3 d2 + 1]
            d2 = nchoosek(k + 2, 2);
            ix = [1, 3*d2 + 1];
        end

    end
end
