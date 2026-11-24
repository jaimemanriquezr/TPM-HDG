classdef HyperelasticGradientProblem < FEProblem
    %-------------------------------------------------------------------------
    %Description: HDG discretization of the geometrically nonlinear
    %             (hyperelastic) solid in the total Lagrangian form with the
    %             displacement gradient as the mixed variable (the
    %             formulation of Nguyen and Peraire, J. Comput. Phys. 231
    %             (2012), eqs. 64a–64e, compressible case, written for the
    %             displacement u = φ − x, so F = I + H):
    %             Find (H_h, u_h, û_h), H_h a full (nonsymmetric) P_k tensor,
    %               (H_h, E) + (u_h, ∇·E) − <û_h, E n> = 0
    %               (P(I + H_h), ∇w) − <P̂n, w> + (ρü_h + b₀u_h, w) = (f, w)
    %               Σ_K <P̂n, μ>_∂K = <g_N, μ>_Γ_N
    %             with the numerical (Piola) traction
    %               P̂n = P(I + H_h) n − τ (u_h − û_h).
    %             Materials (plane strain, Lamé parameters from E, ν):
    %               "saint_venant_kirchhoff", "neo_hookean" (the laws of
    %               HyperelasticProblem.first_piola) and "linear"
    %               (P = Cε(H), the linear HDG elasticity of the same form).
    %
    %Stabilization: τ constant on each face, from Data.stabilization:
    %               "constant"  τ = stabilization_constant μ / characteristic_length
    %                           (Nguyen–Peraire eq. 65);
    %               "mesh_size" τ = stabilization_constant μ / h_F, h_F the face
    %                           length (the O(1/h) scaling);
    %               "adaptive"  τ = ρ_max(∂P/∂F) / (2 characteristic_length),
    %                           ρ_max the largest tangent eigenvalue over the
    %                           face's quadrature points of both neighbours at
    %                           the Newton iterate (Terrana, Nguyen, Bonet,
    %                           Peraire, CMAME 352 (2019), eq. 18), lagged:
    %                           recomputed after every Newton step;
    %               "cockburn_shen" τ_KF = C_KF,o / h_F on the elements of
    %                           T_h^(-,nb) and τ = 1 elsewhere, both times
    %                           stabilization_scale (Cockburn and Shen,
    %                           Results in Applied Mathematics 1 (2019),
    %                           Theorem 2.4 and Sec. 2.5.5 without the
    %                           global constant C_KF,θ; T_h^(-,nb): the
    %                           elements where ∂P/∂F at the iterate has a
    %                           negative eigenvalue, and their vertex
    %                           neighbours), "lagged" or "frozen"
    %                           (stabilization_update, see solve). The
    %                           degree is uniform: their raise to k + 1 on
    %                           T_h^(-,nb) is not implemented.
    %             The default, "mesh_size" with C = 40, is the robust choice:
    %             with τ = O(1) ("constant", "adaptive") Newton finds
    %             spurious solutions at moderate strains, and the constant
    %             needed grows with the strain (C = 10 suffices for
    %             |∇u| ≲ 0.16, not for 0.3).
    %
    %Solve:       solve(problem, mesh, k) runs Newton's method; each step
    %             is one low-level solve (static condensation of (H, u) onto
    %             û) of the linearization at Data.iterate (struct with
    %             fields gradient, 4*d2 x N, displacement, 2*d2 x N, and
    %             displacement_skeleton, 2*d1 x F), written for the new
    %             iterate. Without an iterate the first step is the linear
    %             solve. load_steps > 1 ramps the source, Dirichlet and
    %             Neumann data in equal fractions (each step from the
    %             previous one), halving a fraction that fails.
    %
    %Time stepping: time_scheme = "newmark" (β, γ), with Data.state built by
    %             initial_state and advanced by get_kinematic_state:
    %                 state = initial_state(problem, mesh, k);
    %                 for n = 1:n_steps
    %                     problem.Data.source_term = f(n*dt);  % and the BC data
    %                     problem.Data.state = state;
    %                     solution = solve(problem, mesh, k);
    %                     state = get_kinematic_state(problem, solution, state);
    %                 end
    %             state.solid holds the volume arrays (2*d2 x N)
    %             displacement, velocity, acceleration, the skeleton arrays
    %             (2*d1 x F) displacement_skeleton, velocity_skeleton,
    %             acceleration_skeleton, and gradient (4*d2 x N, H_h).
    %             traction(problem, mesh, solution) returns the Piola
    %             traction P̂n on the element faces.
    %-------------------------------------------------------------------------
    methods
        function obj = HyperelasticGradientProblem(k, label, options)
            arguments (Input)
                k       double          {mustBeScalarOrEmpty} = []
                label   string          {mustBeScalarOrEmpty} = "Hyperelasticity (gradient)"

                options.material            (1,1)       string {mustBeMember(options.material, ...
                    ["saint_venant_kirchhoff", "neo_hookean", "linear"])} = "saint_venant_kirchhoff"
                options.stabilization       (1,1)       string {mustBeMember(options.stabilization, ...
                    ["constant", "mesh_size", "adaptive", "cockburn_shen"])} = "mesh_size"
                options.stabilization_constant (1,1)    double = 40     % C of the constant and mesh_size rules
                options.characteristic_length (1,1)     double = 1      % ℓ = L_c (fixed under refinement)
                options.stabilization_scale (1,1)       double = 1      % stress scale of the cockburn_shen rule
                options.stabilization_update (1,1)      string {mustBeMember(options.stabilization_update, ["lagged", "frozen"])} = "lagged"
                options.newton_tolerance    (1,1)       double = 1e-10  % relative residual tolerance
                options.newton_max_iterations (1,1)     double = 25     % Newton iteration cap
                options.load_steps          (1,1)       double = 1      % load fractions of a steady solve
                options.iterate                         struct = struct.empty  % Newton iterate
                options.density                         double = 0.0    % ρ (0 → static; >0 → transient)
                options.time                (1,1)       double = 0.0    % current time t
                options.time_step           (1,1)       double = 1e-1   % Newmark time step Δt
                options.initial_displacement (2,1)      cell   = {@(x,y) 0*(x.*y); @(x,y) 0*(x.*y)}  % u(0) (initial_state)
                options.initial_velocity     (2,1)      cell   = {@(x,y) 0*(x.*y); @(x,y) 0*(x.*y)}  % u̇(0) (initial_state)
                options.time_scheme         (1,1)       string {mustBeMember(options.time_scheme, "newmark")} = "newmark"
                options.newmark_beta        (1,1)       double = 0.25   % β (Newmark)
                options.newmark_gamma       (1,1)       double = 0.5    % γ (Newmark)
                options.state                           struct = struct.empty  % step history
                options.hooke_modulus       (1,1)       double = 0.0    % b₀ (Winkler / Robin stiffness)
                options.young_modulus       (1,1)       double = 1e4    % E
                options.poisson_ratio       (1,1)       double = .1     % ν
                options.source_term         (2,1)       cell   = {@(x,y) 0*(x.*y); @(x,y) 0*(x.*y)}  % body force {f_x, f_y}
                options.datum_dirichlet     (2,1)       cell   = {@(x,y) 0*(x.*y); @(x,y) 0*(x.*y)}  % Dirichlet data {g_x, g_y}
                options.datum_neumann       (2,1)       cell   = {@(x,y) 0*(x.*y); @(x,y) 0*(x.*y)}  % reference traction {q_x, q_y}
            end
            obj = obj@FEProblem(label, Order=k, Data=options);
            obj.FieldsLocal = [struct(name="H_h", description="displacement gradient");
                struct(name="u_h", description="displacement")];
            obj.FieldsHybrid = struct(name="\hat{u}_h", description="displacement trace");
        end
    end

    methods (Static)
        function spaces = get_discrete_spaces(k)
            % Volume: full tensor Dubiner P_k^(2x2) for H (components
            % 11, 21, 12, 22), vector Dubiner P_k^2 for u. Hybrid: vector
            % Legendre P_k^2 for û.
            if isempty(k); spaces = struct.empty; return; end
            spaces = struct( ...
                "volume", [DubinerBasis(k).to_tensor(), DubinerBasis(k).to_vector(2)], ...
                "hybrid", LegendreBasis(k).to_vector(2));
        end

        function data = get_manufactured_solution(u, young_modulus, poisson_ratio, hooke_modulus, density, time, material)
            % Manufactured data from an exact displacement u(x, y, t)
            % (strings): f = ρ ü − ∇·P(I + ∇u) + b₀ u. Returns handles at
            % t = time: gradient (2x2, ∇u, field 1), u (2x1, field 2), u_t,
            % u_tt, f (2x1) and piola (2x2, P).
            arguments
                u             = ["x^2*y"; "x*y^2"]
                young_modulus = 1
                poisson_ratio = 0.3
                hooke_modulus = 0.0
                density       = 0.0
                time          = 0.0
                material      (1,1) string = "saint_venant_kirchhoff"
            end
            [lambda, mu] = LinearElasticityProblem.get_lame_parameters(young_modulus, poisson_ratio);
            syms x y t real
            u = str2sym(u(:));
            G = jacobian(u, [x y]);
            if material == "linear"
                P = lambda*trace(G)*eye(2) + mu*(G + G.');
            else
                P = HyperelasticProblem.first_piola_symbolic(G, lambda, mu, material);
            end
            f = density*diff(u, t, 2) - [divergence(P(1, :), [x y]); divergence(P(2, :), [x y])] ...
                + hooke_modulus*u;
            at = @(e) subs(e, t, time);
            data = struct("young_modulus", young_modulus, "poisson_ratio", poisson_ratio, ...
                "gradient", {handles(at(G))}, "u", {handles(at(u))}, ...
                "u_t", {handles(at(diff(u, t)))}, "u_tt", {handles(at(diff(u, t, 2)))}, ...
                "f", {handles(at(f))}, "piola", {handles(at(P))});
        end

        function [P, dP] = first_piola(G, lambda, mu, material)
            % P(I + G) and ∂P/∂G pointwise (n x 4 and n x 4 x 4, components
            % 11, 21, 12, 22); "linear" is Hooke's law Cε(G).
            if material == "linear"
                trG = G(:, 1) + G(:, 4);
                P = [lambda*trG + 2*mu*G(:, 1), mu*(G(:, 2) + G(:, 3)), ...
                     mu*(G(:, 2) + G(:, 3)), lambda*trG + 2*mu*G(:, 4)];
                D = [lambda + 2*mu, 0, 0, lambda; 0, mu, mu, 0; 0, mu, mu, 0; lambda, 0, 0, lambda + 2*mu];
                dP = repmat(reshape(D, 1, 4, 4), size(G, 1), 1, 1);
            else
                [P, dP] = HyperelasticProblem.first_piola(G, lambda, mu, material);
            end
        end
    end
end

function h = handles(e)
% vectorized @(x, y) handles of a symbolic array
syms x y real
h = cell(size(e));
for i = 1:numel(e)
    g = matlabFunction(e(i), Vars = [x, y]);
    h{i} = @(x, y) g(x, y) + 0*x;
end
end
