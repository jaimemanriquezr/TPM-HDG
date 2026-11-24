%-------------------------------------------------------------------------
%Description: Consistent initial acceleration of the Newmark solid
%             (LinearElasticityProblem.initial_state).
%
%             Manufactured displacements u(x, y, t) on fitted unit-square
%             meshes, k = 1, 2, 3, with a Winkler term (hooke_modulus) and
%             density ρ:
%               (A) separable, sin(t + π/4) w(x, y)   (u_tt = -u)
%               (B) non-separable in (x, y, t)
%               (C) polynomial of degree k in (x, y), quadratic in t
%             initial_state under newmark sets (u⁰, û⁰, σ⁰, γ⁰) to the HDG
%             elliptic projection of u(0): the static solution with source
%             g = -∇·Cε(u(0)) + b₀ u(0) and Dirichlet datum u(0). Its local
%             momentum row turns the consistent acceleration into
%             ρ a⁰ = Π (f(0) - g) = ρ Π u_tt(0), so a⁰ is the projection
%             itself, not an O(h^k) approximation of it (the L² projections
%             of u(0) leave an O(h^k) face term <(σ⁰ - σ(0)) n, v>_∂K in
%             that row).
%             Checks (EOC = rate of the last two meshes):
%               (1) the volume acceleration a⁰ against Π u_tt(0): round-off
%                   for A, B, C;
%               (2) the skeleton acceleration â⁰ (face mean of the traces
%                   of a⁰) against P_M u_tt(0): rate k + 1 for A, B
%                   (EOC ≥ k + 0.7), round-off for C (this also pins the
%                   face orientation);
%               (3) one Newmark step from initial_state against the same
%                   step from the state with a⁰ = Π u_tt(0),
%                   â⁰ = P_M u_tt(0): round-off for A and C;
%               (4) for C, the stress σ⁰ is the exact σ(0);
%               (5) for C, the elliptic projection reproduces the
%                   polynomial data: u⁰ = Π u(0), û⁰ = P_M u(0),
%                   v⁰ = Π u_t(0), v̂⁰ = P_M u_t(0) (round-off);
%               (6) σ⁰, γ⁰ equal lift_stress(u⁰, û⁰) (round-off);
%               (7) zero initial data give a zero state, and a handle that
%                   does not evaluate symbolically raises
%                   LinearElasticityProblem:initialDataNotSymbolic.
%             "Round-off" is relative, tolerance 1e-10.
%             The face errors use the discrete norm (Σ_F |F| |e_F|²)^½,
%             equivalent to the L²(skeleton) norm.
%-------------------------------------------------------------------------
n_failed = 0;

E = 2.0; nu = 0.3; K = 0.5; rho = 1.3;
dt = 1e-2;
get_solutions = @LinearElasticityProblem.get_manufactured_solution;
families = struct( ...
    "name", {"A", "B", "C"}, ...
    "u", {["sin(t + pi/4)*sin(pi*x)*cos(pi*y)"; "sin(t + pi/4)*(x^2 + 1)*exp(y)"], ...
          ["sin(pi*x)*sin(pi*y + t)"; "cos(pi*x*t)*sin(pi*y)"], ...
          []});

for k = 1:3
    if k < 3, Ns = [4, 8, 16, 32]; else, Ns = [4, 8, 16]; end
    for fam = families
        u_str = fam.u;
        if fam.name == "C"
            u_str = ["(1 + t + t^2/2)*(x^" + k + " + 2*x*y^" + (k-1) + " + y)";
                     "(1 + t + t^2/2)*(y^" + k + " - x^" + (k-1) + "*y + 3*x)"];
        end
        data0 = get_solutions(u_str, E, nu, K, rho, 0);
        data1 = get_solutions(u_str, E, nu, K, rho, dt);
        u_tt = second_time_derivative_at_zero(u_str);

        [e_a, e_ah, e_u1, e_s, e_p, e_l] = deal(nan(size(Ns)));
        for i = 1:numel(Ns)
            mesh = StructuredSimplicialGrid(Ns(i), Ns(i)).add_stabilization(1);
            pde = LinearElasticityProblem(k, young_modulus = E, poisson_ratio = nu, ...
                hooke_modulus = K, density = rho, time_scheme = "newmark", time_step = dt, ...
                initial_displacement = data0.u, initial_velocity = data0.u_t, ...
                source_term = data0.f);
            state = initial_state(pde, mesh, k);

            a_ref  = L2proj2D(mesh, u_tt, k);
            ah_ref = L2projSkeleton1d(mesh, u_tt, k);
            e_a(i)  = volume_norm(mesh, k, state.solid.acceleration - a_ref) / volume_norm(mesh, k, a_ref);
            e_ah(i) = face_norm(mesh, state.solid.acceleration_skeleton - ah_ref) / face_norm(mesh, ah_ref);

            % one Newmark step from the computed and from the projected a⁰
            pde.Data.source_term = data1.f;
            pde.Data.datum_dirichlet = data1.u;
            pde.Data.state = state;
            u1 = solve(pde, mesh, k).Unknowns.field(2).dof;
            reference = state;
            reference.solid.acceleration = a_ref;
            reference.solid.acceleration_skeleton = ah_ref;
            pde.Data.state = reference;
            u1_ref = solve(pde, mesh, k).Unknowns.field(2).dof;
            e_u1(i) = volume_norm(mesh, k, u1 - u1_ref) / volume_norm(mesh, k, u1_ref);

            % the stress of the projection is the lifted one
            [sigma_l, gamma_l] = pde.lift_stress(mesh, k, state.solid.displacement, ...
                state.solid.displacement_skeleton);
            e_l(i) = max(max(abs(state.solid.stress - sigma_l), [], 'all') ...
                / max(abs(sigma_l), [], 'all'), max(abs(state.solid.rotation - gamma_l), [], 'all') ...
                / max(abs(state.solid.stress), [], 'all'));

            if fam.name == "C"
                rel = @(a, b) max(abs(a - b), [], 'all') / max(abs(b), [], 'all');
                e_p(i) = max([rel(state.solid.displacement, L2proj2D(mesh, data0.u, k)), ...
                    rel(state.solid.displacement_skeleton, L2projSkeleton1d(mesh, data0.u, k)), ...
                    rel(state.solid.velocity, L2proj2D(mesh, data0.u_t, k)), ...
                    rel(state.solid.velocity_skeleton, L2projSkeleton1d(mesh, data0.u_t, k))]);
                d2 = nchoosek(k + 2, 2);
                sigma_ref = [L2proj2D(mesh, data0.sigma{1, 1}, k); L2proj2D(mesh, data0.sigma{2, 1}, k);
                             L2proj2D(mesh, data0.sigma{1, 2}, k); L2proj2D(mesh, data0.sigma{2, 2}, k);
                             zeros(size(state.solid.stress, 1) - 4*d2, size(mesh.elements, 2))];
                e_s(i) = max(abs(state.solid.stress - sigma_ref), [], 'all') ...
                    / max(abs(sigma_ref), [], 'all');
            end
        end

        fprintf("\nk = %d, family %s\n", k, fam.name);
        fprintf("   N   |a0 - Pi u_tt|  EOC   |ah0 - P_M u_tt|  EOC   |u1 - u1(Pi u_tt)|  EOC\n");
        eoc = @(e) [nan, compute_eoc(e(:), 1./Ns).'];
        [r_a, r_ah, r_u1] = deal(eoc(e_a), eoc(e_ah), eoc(e_u1));
        for i = 1:numel(Ns)
            fprintf("  %3d   %.3e  %5.2f    %.3e  %5.2f      %.3e  %5.2f\n", ...
                Ns(i), e_a(i), r_a(i), e_ah(i), r_ah(i), e_u1(i), r_u1(i));
        end
        tag = sprintf("k = %d, family %s: ", k, fam.name);
        fprintf("  max relative |(sigma0, gamma0) - lift_stress(u0, uh0)| = %.3e\n", max(e_l));
        n_failed = n_failed + check(tag + "a0 = Pi u_tt (round-off)", max(e_a) <= 1e-10);
        n_failed = n_failed + check(tag + "stress = lift_stress(u0, uh0)", max(e_l) <= 1e-10);
        if fam.name == "C"
            fprintf("  max relative |sigma0 - sigma(0)| = %.3e\n", max(e_s));
            fprintf("  max relative |(u0, uh0, v0, vh0) - projections| = %.3e\n", max(e_p));
            n_failed = n_failed + check(tag + "ah0 = P_M u_tt (round-off)", max(e_ah) <= 1e-10);
            n_failed = n_failed + check(tag + "first step from a0 = from Pi u_tt", max(e_u1) <= 1e-10);
            n_failed = n_failed + check(tag + "stress = sigma(0)", max(e_s) <= 1e-10);
            n_failed = n_failed + check(tag + "elliptic projection reproduces the data", max(e_p) <= 1e-10);
        else
            n_failed = n_failed + check(tag + "ah0 EOC >= k + 0.7", r_ah(end) >= k + 0.7);
            if fam.name == "A"
                n_failed = n_failed + check(tag + "first step from a0 = from Pi u_tt", max(e_u1) <= 1e-10);
            end
        end
    end
end

% zero data; data that is not symbolic
mesh = StructuredSimplicialGrid(4, 4).add_stabilization(1);
pde = LinearElasticityProblem(1, density = rho, time_scheme = "newmark", time_step = dt);
state = initial_state(pde, mesh, 1);
zero_state = all(cellfun(@(v) all(v == 0, 'all'), struct2cell(state.solid)));
n_failed = n_failed + check("zero initial data: zero state", zero_state);
pde.Data.initial_displacement = {@(x, y) ones(size(x)); @(x, y) 0*x};
try
    initial_state(pde, mesh, 1);
    error_id = "";
catch err
    error_id = string(err.identifier);
end
n_failed = n_failed + check("non-symbolic data raise initialDataNotSymbolic", ...
    error_id == "LinearElasticityProblem:initialDataNotSymbolic");

fprintf("\n%d check(s) failed.\n", n_failed);

%=========================================================================%
function u_tt = second_time_derivative_at_zero(u_str)
% {∂_tt u_x; ∂_tt u_y} at t = 0 as handles of (x, y)
syms x y t real
u = str2sym(u_str);
u_tt = cell(2, 1);
for c = 1:2
    h = matlabFunction(subs(diff(u(c), t, 2), t, 0), Vars = [x, y]);
    u_tt{c} = @(x, y) h(x, y) + 0*x;
end
end

function n = volume_norm(mesh, k, e)
% L² norm of the 2*d2 x N volume dofs e (x block over y block)
M = mass_matrix(mesh, k);
d2 = size(M, 1);
n = 0;
for c = 1:2
    e_c = reshape(e((c-1)*d2 + (1:d2), :), d2, 1, []);
    n = n + sum(pagemtimes(pagemtimes(e_c, 'transpose', M, 'none'), e_c), 'all');
end
n = sqrt(n);
end

function n = face_norm(mesh, e)
% (Σ_F |F| |e_F|²)^½ of the 2*d1 x F skeleton dofs e
n = sqrt(sum(mesh.area .* sum(e.^2, 1)));
end

function failed = check(name, condition)
failed = ~condition;
if failed
    fprintf(2, "  FAIL  %s\n", name);
else
    fprintf("  ok    %s\n", name);
end
end
