classdef BiotProblem < FEProblem
    %-------------------------------------------------------------------------
    %Description: HDG discretization of Biot poroelasticity (PEERS total
    %             stress, Darcy flux), one implicit time step:
    %             find (σ, η, γ, q, p) with
    %               A(σ + α p I) − ∇η + γ = 0,   ρ_P ∂ₜₜη − ∇·σ + b₀ η = f,
    %               κ⁻¹ q + ∇p = 0,   s₀ ∂ₜp + α ∂ₜ∇·η + ∇·q = g   in Ω,
    %             σ = 2μ ε(η) + λ (∇·η) I − α p I the total stress, A the
    %             compliance tensor; quasi-static when density ρ_P = 0.
    %             Hybrid unknown: one three-component trace (η̂_x, η̂_y, p̂)
    %             on the mesh skeleton.
    %
    %             The algebra is BiotMatrices on the declarative path:
    %             get_declarative_problem builds the FEProblem, solve
    %             solves it. ∂ₜw ≈ (w − w*)/cdt, with (cdt, w*) from
    %             biot_step_history(time_scheme, time_step, state.biot);
    %             an empty state is the zero history (cdt = time_step).
    %             The inertia applies it twice: v = (η − η*)/cdt and
    %             ∂ₜₜη ≈ (v − v*)/cdt, with the volume displacement and
    %             velocity of state.biot (initial_velocity at t = 0).
    %
    %             Boundary conditions, first match: Data.boundary_conditions
    %             (componentwise specs on η̂, and on p̂ when a spec sets
    %             boundary_condition's p; otherwise p̂ is natural, zero flux), a
    %             transfer_mapping (transferred Dirichlet data), otherwise
    %             strong Dirichlet data datum_dirichlet = {η_x; η_y; p} on
    %             the faces tagged "Dirichlet".
    %-------------------------------------------------------------------------
    methods
        function obj = BiotProblem(k, label, options)
            arguments (Input)
                k     double {mustBeScalarOrEmpty} = []
                label string {mustBeScalarOrEmpty} = "Biot (quasi-static)"
                options.young_modulus       (1,1) double = 1.0
                options.poisson_ratio       (1,1) double = 0.25
                options.biot_alpha          (1,1) double {mustBeInRange(options.biot_alpha, 0, 1)} = 1.0
                options.storage             (1,1) double {mustBeNonnegative} = 0.0
                options.diffusion           (2,2) double = eye(2)          % kappa
                options.hooke_modulus       (1,1) double {mustBeNonnegative} = 0.0
                options.density             (1,1) double {mustBeNonnegative} = 0.0   % rho_P
                options.stabilization_darcy (1,1) double {mustBePositive} = 1.0
                options.stabilization_solid (1,1) double {mustBePositive} = 1.0   % tau_s / tau
                options.time                (1,1) double = 0.0
                options.time_step           (1,1) double {mustBePositive, mustBeFinite} = 1.0
                options.time_scheme         (1,1) string {mustBeMember(options.time_scheme, ["euler", "bdf2"])} = "euler"
                options.state                     struct = struct.empty   % the coupled state (reads state.biot)
                options.source_term          (2,1) cell = {@(x,y) 0*(x.*y); @(x,y) 0*(x.*y)}          % f_p
                options.source_mass          (1,1) function_handle = @(x,y) 0*(x.*y)                   % g_p
                options.initial_displacement (2,1) cell = {@(x,y) 0*(x.*y); @(x,y) 0*(x.*y)}
                options.initial_velocity     (2,1) cell = {@(x,y) 0*(x.*y); @(x,y) 0*(x.*y)}
                options.initial_pressure     (1,1) function_handle = @(x,y) 0*(x.*y)
                options.datum_dirichlet      (3,1) cell = {@(x,y) 0*(x.*y); @(x,y) 0*(x.*y); @(x,y) 0*(x.*y)}   % {eta_x; eta_y; p}
                options.transfer_mapping {mustBeA(options.transfer_mapping, ["function_handle", "LevelSetCurve"])} = function_handle.empty
            end
            obj = obj@FEProblem(label, Order=k, Data=options);
            obj.FieldsLocal = [struct(name="\sigma_h", description="total stress");
                struct(name="\eta_h",   description="displacement");
                struct(name="\gamma_h", description="rotation");
                struct(name="q_h",      description="Darcy flux");
                struct(name="p_h",      description="pore pressure")];
            obj.FieldsHybrid = struct(name="(\hat{\eta}_h, \hat{p}_h)", ...
                description="displacement and pressure trace");
        end
    end

    methods (Static)
        spaces = get_discrete_spaces(k)
        mms = get_manufactured_solution(displacement, pressure, options)

        function factor = darcy_tau_factor(stabilization_darcy, kappa)
            % tau_p / tau: the Darcy stabilization scales with the
            % permeability, tau_p = stabilization_darcy * ||kappa|| * tau
            % (spectral norm), in the element (BiotMatrices) and in the
            % interface rows (InterfaceProblem.stokes_biot)
            factor = stabilization_darcy * norm(kappa);
        end

        function sigma_eff = effective_stress(stress, pressure, alpha, k)
            % stress + alpha*p*I on PEERS dofs: the diagonal Dubiner
            % blocks (rows 1:d2 and 3*d2 + (1:d2)) get alpha*pressure;
            % pressure 0 or [] leaves the stress unchanged
            sigma_eff = stress;
            if isempty(pressure) || isequal(pressure, 0) || alpha == 0
                return
            end
            d2 = nchoosek(k + 2, 2);
            rows = [1:d2, 3*d2 + (1:d2)];
            sigma_eff(rows, :) = sigma_eff(rows, :) + alpha*[pressure; pressure];
        end
    end
end
