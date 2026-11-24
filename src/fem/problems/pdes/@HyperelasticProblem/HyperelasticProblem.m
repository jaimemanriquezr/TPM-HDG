classdef HyperelasticProblem < LinearElasticityProblem
    %-------------------------------------------------------------------------
    %Description: HDG discretization of the geometrically nonlinear
    %             (hyperelastic) solid in the total Lagrangian form, on the
    %             PEERS unknowns (σ, u, γ; û) of LinearElasticityProblem,
    %             with σ the first Piola–Kirchhoff stress:
    %               ζ⁻¹(σ − N(∇u)) − ∇u + γ = 0,   (σ − N(∇u)) symmetric,
    %               ρ ü − ∇·σ + b₀ u = f                in Ω (reference),
    %             where P(I + G) = Cε(G) + N(G) splits the stress law into
    %             Hooke's law and a remainder N, quadratic at G = 0, and
    %             γ = skew ∇u. The remainder is evaluated at the discrete
    %             gradient L_h of (u_h, û_h), the PEERS field with
    %               (L_h, τ) = −(u_h, ∇·τ) + <û_h, τ n>_∂K   for τ in S_h,
    %             and enters only the constitutive and weak-symmetry rows,
    %               (ζ⁻¹σ, τ) + (u, ∇·τ) + (γ, τ) − (ζ⁻¹N(L_h), τ) = <û, τn>
    %               (σ, η) = (N(L_h), η);
    %             the momentum row, the numerical flux
    %             σ̂n = σ n − τ(u − û) and the transmission (trace) rows are
    %             those of the linear class. With N ≡ 0 the scheme is
    %             LinearElasticityProblem, and it differs from it by
    %             O(load²) at small loads.
    %             Materials (plane strain, Lamé parameters from E, ν):
    %               "saint_venant_kirchhoff"  P = F (λ tr E I + 2μ E),
    %                                         E = (FᵀF − I)/2
    %               "neo_hookean"             P = μ(F − F⁻ᵀ) + λ log(J) F⁻ᵀ
    %               "neo_hookean_eulerian"    the Cauchy stress of the same
    %                                         law, J⁻¹ P Fᵀ, as a function of
    %                                         the Eulerian displacement
    %                                         gradient G = ∇ₓu, F = (I − G)⁻¹,
    %                                         for the Eulerian solid; its
    %                                         remainder N is symmetric, so the
    %                                         weak-symmetry row is the linear one
    %               "saint_venant_kirchhoff_eulerian"  the same for the
    %                                         Saint Venant–Kirchhoff law,
    %                                         σ = J⁻¹ F S Fᵀ
    %             Neumann data are reference tractions P n.
    %
    %Solve:       solve(problem, mesh, k) runs Newton's method. Each Newton
    %             step is one solve of the low-level path whose local
    %             operators (get_matrices) are the linearization at the
    %             iterate Data.iterate (struct with displacement and
    %             displacement_skeleton, the layout of state.solid), written
    %             for the new iterate itself, so the hybrid system, its
    %             static condensation and the boundary conditions are those
    %             of the linear class. Without an iterate the first step is
    %             the linear solve. Stops when the residual of the full
    %             (local and global) nonlinear system (newton_residual),
    %             relative to the linear residual of the zero iterate, is
    %             below newton_tolerance; the Newton report is the third
    %             output of solve.
    %
    %Stabilization: stabilization = "constant" uses mesh.tau (the linear
    %             class's τ). "adaptive" sets, per face, the rule of Terrana,
    %             Nguyen, Bonet and Peraire (CMAME 352, 2019, eq. 18),
    %               τ_F = ρ_max(∂P/∂F) / (2 L_c),
    %             ρ_max the largest eigenvalue of the tangent over the face's
    %             quadrature points at the iterate, L_c = characteristic_length
    %             a fixed length (the size of a coarse mesh resolving the
    %             solution), linear value (μ + λ)/L_c in plane strain. With
    %             stabilization_update = "lagged" τ is recomputed after every
    %             Newton step and held fixed in the Jacobian; "frozen" keeps
    %             the value of the initial iterate for the whole solve.
    %
    %Time stepping: time_scheme = "newmark" only, with the state struct,
    %             initial_state and get_kinematic_state of
    %             LinearElasticityProblem (same fields, same loop). The
    %             inertial terms are linear and enter as in the linear class;
    %             the Newton iteration of a step starts from the state.
    %
    %Viscosity:   solid_viscosity = μ_v > 0 (transient only) adds the
    %             Kelvin–Voigt stress of the Green–Lagrange strain rate,
    %               σ = P(F) + F S_v,   S_v = 2μ_v Ė = μ_v(ḞᵀF + FᵀḞ),
    %             objective (it vanishes on rigid rotations) and, at small
    %             strain, 2μ_v ε(u̇). It joins the stress remainder N
    %             (nonlinear_terms), evaluated at the discrete gradient L_h
    %             and at its rate, the Newmark velocity gradient
    %               L̇_h = γ/(βΔt) (L_h − L_h⁰) + (1 − γ/β) L̇_h⁰
    %                     + Δt (1 − γ/(2β)) L̈_h⁰,
    %             the history from Data.state (velocity, acceleration and
    %             their skeletons); the Newton step linearizes both. The
    %             initial state is the elastic one: its consistent
    %             acceleration omits ∇·(F S_v) at t = 0, exact for u̇(0) = 0.
    %             μ_v = 0 is the scheme without viscosity.
    %             viscous_model selects the viscous stress (viscous_stress):
    %               "material"  (default) the stress above, F S_v; its
    %                           Cauchy stress is (2μ_v/J) B D B;
    %               "spatial"   the Newtonian Cauchy stress 2μ_v D,
    %                           D = sym(Ḟ F⁻¹) the full (not deviatoric)
    %                           rate of deformation of the compressible
    %                           solid, as first Piola stress
    %                           σ_v = 2μ_v J D F⁻ᵀ. In plane strain
    %                           2μ_v D = 2μ_v dev D + μ_v tr(D) I: the
    %                           full D adds a bulk viscous part
    %                           μ_v tr(D) I, absent for an
    %                           incompressible solid (tr D = 0).
    %             Both are objective and agree to O(|∇u| |∇u̇|).
    %-------------------------------------------------------------------------
    methods
        function obj = HyperelasticProblem(k, label, options)
            arguments (Input)
                k       double          {mustBeScalarOrEmpty} = []             % polynomial degree
                label   string          {mustBeScalarOrEmpty} = "Hyperelasticity" % display label

                options.material            (1,1)       string {mustBeMember(options.material, ...
                    ["saint_venant_kirchhoff", "neo_hookean", "neo_hookean_eulerian", ...
                    "saint_venant_kirchhoff_eulerian"])} = "saint_venant_kirchhoff"
                options.newton_tolerance    (1,1)       double = 1e-10  % relative residual tolerance
                options.newton_max_iterations (1,1)     double = 25     % Newton iteration cap
                options.iterate                         struct = struct.empty  % Newton iterate: displacement, displacement_skeleton
                % "constant": τ = mesh.tau, as in the linear class; "adaptive":
                % τ_F = ρ_max(∂P/∂F)/(2 characteristic_length) per face
                % (adaptive_tau), at the Newton iterate ("lagged": recomputed
                % after every step) or at the initial iterate ("frozen")
                options.stabilization       (1,1)       string {mustBeMember(options.stabilization, ["constant", "adaptive"])} = "constant"
                options.characteristic_length (1,1)     double = 0.5    % L_c of the adaptive τ (fixed under refinement)
                options.stabilization_update (1,1)      string {mustBeMember(options.stabilization_update, ["lagged", "frozen"])} = "lagged"
                options.solid_viscosity     (1,1)       double {mustBeNonnegative} = 0  % μ_v (Kelvin–Voigt, transient only)
                options.viscous_model       (1,1)       string {mustBeMember(options.viscous_model, ["material", "spatial"])} = "material"  % F S_v or 2μ_v J D F⁻ᵀ
                options.density                         double = 0.0    % ρ (0 → static; >0 → transient)
                options.time                (1,1)       double = 0.0    % current time t
                options.time_step           (1,1)       double = 1e-1   % Newmark time step Δt
                options.initial_displacement (2,1)      cell   = {@(x,y) 0*(x.*y); @(x,y) 0*(x.*y)}  % u(0) (initial_state)
                options.initial_velocity     (2,1)      cell   = {@(x,y) 0*(x.*y); @(x,y) 0*(x.*y)}  % u̇(0) (initial_state)
                options.time_scheme         (1,1)       string {mustBeMember(options.time_scheme, "newmark")} = "newmark"
                options.newmark_beta        (1,1)       double = 0.25   % β (Newmark)
                options.newmark_gamma       (1,1)       double = 0.5    % γ (Newmark)
                options.state                           struct = struct.empty  % step history (initial_state / get_kinematic_state)
                options.hooke_modulus       (1,1)       double = 0.0    % b₀ (Winkler / Robin stiffness)
                options.young_modulus       (1,1)       double = 1e4    % E
                options.poisson_ratio       (1,1)       double = .1     % ν
                options.source_term         (2,1)       cell   = {@(x,y) 0*(x.*y); @(x,y) 0*(x.*y)}  % body force {f_x, f_y}
                options.datum_dirichlet     (2,1)       cell   = {@(x,y) 0*(x.*y); @(x,y) 0*(x.*y)}  % Dirichlet data {g_x, g_y}
                options.datum_neumann       (2,1)       cell   = {@(x,y) 0*(x.*y); @(x,y) 0*(x.*y)}  % reference traction {q_x, q_y}
            end
            nonlinear = ["material", "newton_tolerance", "newton_max_iterations", "iterate", ...
                "stabilization", "characteristic_length", "stabilization_update", "solid_viscosity", "viscous_model"];
            linear = rmfield(options, nonlinear);
            arguments_linear = namedargs2cell(linear);
            obj = obj@LinearElasticityProblem(k, label, arguments_linear{:});
            for name = nonlinear
                obj.Data.(name) = options.(name);
            end
        end
    end

    methods (Static)
        data = get_manufactured_solution(u, young_modulus, poisson_ratio, hooke_modulus, density, time, material, solid_viscosity)
        [P, dP] = first_piola(G, lambda, mu, material)
        P = first_piola_symbolic(G, lambda, mu, material)
        [P, dP] = viscous_stress(G, G_dot, viscosity, rate, model)
        [R_u, R_hat] = discrete_gradient_operators(mesh, k)
    end
end
