%-------------------------------------------------------------------------
%Description: Flux-balanced lift of a volume displacement
%             (LinearElasticityProblem.balanced_lift, balanced_acceleration).
%
%             Meshes: a structured unit square (8 x 8; one side tagged
%             Neumann) and the solid of a carved and snapped disc (centre
%             (0.5, 0.5), radius 0.2, background 16 x 16 "crossed", every
%             boundary face Dirichlet); k = 1, 2; E = 10, 1000; ν = 0.3,
%             b₀ = 0.5, ρ = 1. Data from a manufactured displacement.
%             Checks ("round-off" = relative, max-norm over dofs, 1e-10;
%             the acceleration checks relative L², 1e-9, see (b)):
%               (a) for a static solve, balanced_lift(u_h) returns the
%                   solve's (σ, γ, û) to round-off, with the boundary data
%                   from the problem's handles and with boundary_trace =
%                   the solve's own û;
%               (b) after a Newmark step (from initial_state, and from a
%                   random history), balanced_acceleration(u¹) with the
%                   new-time source and data returns the step's
%                   acceleration a¹ (both boundary routes), relative to
%                   the largest term of the momentum balance
%                   M a¹ = f - K u¹, max(|a¹|, |M⁻¹ K u¹|): a¹ is a
%                   difference of large terms (|M⁻¹ K u¹| reaches
%                   ~3e3 |a¹| at E = 1000), which the error relative to
%                   |a¹| alone (printed in brackets) inherits;
%               (c) with zero data, b₀ = 0 and random u,
%                   -(u, ρ a)_L2 = c(σ, σ) + |u - û|²_S > 0, with
%                   c(σ, σ) = c1 (σ, σ) - c2 (tr σ, tr σ) the
%                   complementary energy of the lifted stress, and
%                   (u₁, ρ a(u₂)) = (u₂, ρ a(u₁));
%               (d) agreement with the Newmark pseudo-step route (a step of
%                   length τ = 1e-6 from (u, 0, 0), then
%                   consistent_acceleration of its (σ, γ, û)), homogeneous
%                   and with data, to the pseudo-step's accuracy (1e-4).
%-------------------------------------------------------------------------
n_failed = 0;

nu = 0.3; K = 0.5; rho = 1; dt = 1e-2;
u_str = ["0.1*sin(2*x + y)*cos(t)"; "0.1*cos(x - y)*(1 + t^2)"];
get_solutions = @LinearElasticityProblem.get_manufactured_solution;
curve = CircularCurve([0.5; 0.5], 0.2);
meshes = {square_mesh(8), disc_mesh(16, curve)};
mesh_names = ["square", "disc"];
rng(7);

fprintf("(b) columns: |a - a1| / max(|a1|, |M^-1 K u1|); in brackets |a - a1| / |a1|\n");
fprintf("mesh    k  E      (a) handles  (a) traces   (b) handles             (b) traces              (c) energy   (c) symm.    (d) homog.   (d) data\n");
for im = 1:numel(meshes)
    mesh = meshes{im};
    for k = 1:2
        for E = [10, 1000]
            tag = sprintf("%s k = %d E = %g: ", mesh_names(im), k, E);
            data0 = get_solutions(u_str, E, nu, K, rho, 0);
            data1 = get_solutions(u_str, E, nu, K, rho, dt);
            neumann = {data0.sigma{1, 1}; data0.sigma{2, 1}};

            %% (a) static solve
            static = LinearElasticityProblem(k, young_modulus = E, poisson_ratio = nu, ...
                hooke_modulus = K, source_term = data0.f, datum_dirichlet = data0.u, ...
                datum_neumann = neumann);
            [sigma_s, u_s, gamma_s, uh_s] = unknowns(solve(static, mesh, k));
            [sigma, gamma, u_hat] = static.balanced_lift(mesh, k, u_s);
            e_a1 = max([rel(sigma, sigma_s), rel(gamma, gamma_s), rel(u_hat, uh_s)]);
            [sigma, gamma, u_hat] = static.balanced_lift(mesh, k, u_s, boundary_trace = uh_s);
            e_a2 = max([rel(sigma, sigma_s), rel(gamma, gamma_s), rel(u_hat, uh_s)]);

            %% (b) Newmark steps
            pde = LinearElasticityProblem(k, young_modulus = E, poisson_ratio = nu, ...
                hooke_modulus = K, density = rho, time_scheme = "newmark", time_step = dt, ...
                initial_displacement = data0.u, initial_velocity = data0.u_t, ...
                source_term = data0.f, datum_neumann = neumann);
            states = {initial_state(pde, mesh, k), random_state(pde, mesh, k, dt)};
            pde.Data.source_term = data1.f;
            pde.Data.datum_dirichlet = data1.u;
            pde.Data.datum_neumann = {data1.sigma{1, 1}; data1.sigma{2, 1}};
            stiffness_only = pde;
            stiffness_only.Data.source_term = {@(x, y) 0*x; @(x, y) 0*x};
            [e_b1, e_b2, e_b1r, e_b2r] = deal(0);
            for s = 1:numel(states)
                pde.Data.state = states{s};
                next = get_kinematic_state(pde, solve(pde, mesh, k), states{s});
                a1 = next.solid.acceleration;
                % the scale of the terms of the momentum balance: |a¹| and
                % |M⁻¹ K u¹| (the acceleration of u¹ without source)
                scale = max(volume_norm(mesh, k, a1), volume_norm(mesh, k, ...
                    stiffness_only.balanced_acceleration(mesh, k, next.solid.displacement, ...
                    boundary_trace = next.solid.displacement_skeleton)));
                a = pde.balanced_acceleration(mesh, k, next.solid.displacement);
                e_b1r = max(e_b1r, volume_norm(mesh, k, a - a1)/volume_norm(mesh, k, a1));
                e_b1 = max(e_b1, volume_norm(mesh, k, a - a1)/scale);
                a = pde.balanced_acceleration(mesh, k, next.solid.displacement, ...
                    boundary_trace = next.solid.displacement_skeleton);
                e_b2r = max(e_b2r, volume_norm(mesh, k, a - a1)/volume_norm(mesh, k, a1));
                e_b2 = max(e_b2, volume_norm(mesh, k, a - a1)/scale);
            end

            %% (c) energy identity, zero data, b0 = 0
            zero_data = LinearElasticityProblem(k, young_modulus = E, poisson_ratio = nu, ...
                density = rho, time_scheme = "newmark", time_step = dt);
            mesh_c = mesh;   % zero Neumann data on the square's Neumann side
            u1 = randn(size(u_s));
            u2 = randn(size(u_s));
            [a1, sigma1, ~, uh1] = zero_data.balanced_acceleration(mesh_c, k, u1);
            a2 = zero_data.balanced_acceleration(mesh_c, k, u2);
            work = -rho*volume_inner(mesh_c, k, u1, a1);
            stored = complementary_energy(zero_data, mesh_c, k, sigma1) ...
                + stabilization_energy(mesh_c, k, u1, uh1);
            e_c1 = abs(work - stored)/stored;
            e_c2 = abs(volume_inner(mesh_c, k, u1, a2) - volume_inner(mesh_c, k, u2, a1)) ...
                / sqrt(volume_inner(mesh_c, k, u1, a1)*volume_inner(mesh_c, k, u2, a2));
            positive = work > 0 && stored > 0;

            %% (d) pseudo-step route on the all-Dirichlet mesh
            mesh_d = all_dirichlet(mesh);
            a_dir = pde.balanced_acceleration(mesh_d, k, u1);
            e_d2 = volume_norm(mesh_d, k, a_dir - pseudo_step(pde, mesh_d, k, u1)) ...
                / volume_norm(mesh_d, k, a_dir);
            hom = pde;
            hom.Data.source_term = {@(x, y) 0*x; @(x, y) 0*x};
            hom.Data.datum_dirichlet = {@(x, y) 0*x; @(x, y) 0*x};
            a_hom = hom.balanced_acceleration(mesh_d, k, u1);
            e_d1 = volume_norm(mesh_d, k, a_hom - pseudo_step(hom, mesh_d, k, u1)) ...
                / volume_norm(mesh_d, k, a_hom);

            fprintf("%-6s  %d  %-5g  %.3e    %.3e    %.3e [%.3e]  %.3e [%.3e]  %.3e    %.3e    %.3e    %.3e\n", ...
                mesh_names(im), k, E, e_a1, e_a2, e_b1, e_b1r, e_b2, e_b2r, e_c1, e_c2, e_d1, e_d2);
            n_failed = n_failed + check(e_a1 <= 1e-10, tag + "(a) static lift, data from handles");
            n_failed = n_failed + check(e_a2 <= 1e-10, tag + "(a) static lift, boundary_trace");
            n_failed = n_failed + check(e_b1 <= 1e-9, tag + "(b) Newmark acceleration, handles");
            n_failed = n_failed + check(e_b2 <= 1e-9, tag + "(b) Newmark acceleration, boundary_trace");
            n_failed = n_failed + check(e_c1 <= 1e-10 && positive, tag + "(c) energy identity");
            n_failed = n_failed + check(e_c2 <= 1e-10, tag + "(c) symmetry");
            n_failed = n_failed + check(e_d1 <= 1e-4, tag + "(d) pseudo-step route, homogeneous");
            n_failed = n_failed + check(e_d2 <= 1e-4, tag + "(d) pseudo-step route, data");
        end
    end
end

fprintf("\n%d check(s) failed.\n", n_failed);

%=========================================================================%
function mesh = square_mesh(N)
% unit square, Dirichlet except the side x = 1 (Neumann)
mesh = all_dirichlet(StructuredSimplicialGrid(N, N).add_stabilization(1));
is_boundary = accumarray(mesh.boundaries(:), 1, [size(mesh.faces, 2), 1]) == 1;
x_mid = (mesh.coordinates(1, mesh.faces(1, :)) + mesh.coordinates(1, mesh.faces(2, :)))/2;
mesh.faces(3, is_boundary(:).' & abs(x_mid - 1) < 1e-12) = mesh.boundary_tags("Neumann");
end

function mesh = disc_mesh(N, curve)
% the solid of eulerian_solid_march's carve_both (snapped), all Dirichlet
[~, mesh] = StructuredSimplicialGrid(N, N, "crossed").split(curve, max_normal_deviation = 90);
mesh = all_dirichlet(snap_to_curve(mesh, curve, normal_sign = -1));
end

function mesh = all_dirichlet(mesh)
is_boundary = accumarray(mesh.boundaries(:), 1, [size(mesh.faces, 2), 1]) == 1;
mesh.faces(3, is_boundary) = mesh.boundary_tags("Dirichlet");
end

function [sigma, u, gamma, u_hat] = unknowns(solution)
U = solution.Unknowns;
sigma = U.field(1).dof;
u = U.field(2).dof;
gamma = U.field(3).dof;
if iscell(gamma), gamma = cat(1, gamma{:}); end
u_hat = U.hybrid(1).dof;
end

function state = random_state(pde, mesh, k, dt)
% initial_state's layout with random volume and skeleton histories
state = initial_state(pde, mesh, k);
s = state.solid;
s.displacement = randn(size(s.displacement));
s.velocity = randn(size(s.displacement))/dt;
s.acceleration = randn(size(s.displacement))/dt^2;
state.solid = s;
end

function a = pseudo_step(pde, mesh, k, u)
% a Newmark step of length 1e-6 from (u, 0, 0), then consistent_acceleration
% of (u, σ, γ, û) of that step
q = pde;
q.Data.time_step = 1e-6;
q.Data.newmark_beta = 0.25;
q.Data.newmark_gamma = 0.5;
state = struct("solid", struct("displacement", u, "velocity", zeros(size(u)), ...
    "acceleration", zeros(size(u))));
q.Data.state = state;
[sigma, ~, gamma, u_hat] = unknowns(solve(q, mesh, k));
state.solid.stress = sigma;
state.solid.rotation = gamma;
state.solid.displacement_skeleton = u_hat;
a = pde.consistent_acceleration(mesh, k, state);
end

function e = rel(a, b)
e = max(abs(a - b), [], 'all')/max(abs(b), [], 'all');
end

function v = volume_inner(mesh, k, a, b)
% L² inner product of two 2*d2 x N volume dof arrays (x block over y block)
M = mass_matrix(mesh, k);
d2 = size(M, 1);
v = 0;
for c = 1:2
    a_c = reshape(a((c-1)*d2 + (1:d2), :), d2, 1, []);
    b_c = reshape(b((c-1)*d2 + (1:d2), :), d2, 1, []);
    v = v + sum(pagemtimes(pagemtimes(a_c, 'transpose', M, 'none'), b_c), 'all');
end
end

function n = volume_norm(mesh, k, e)
n = sqrt(volume_inner(mesh, k, e, e));
end

function w = complementary_energy(pde, mesh, k, sigma)
% Σ_K c1 (σ, σ)_K - c2 (tr σ, tr σ)_K: the stress-stress block of get_matrices' A1
static = pde;
static.Data.density = 0;
local_operators = get_matrices(static, mesh, k);
d2peers = size(sigma, 1);
t_s = local_operators{1}(1:d2peers, 1:d2peers, :);
s = reshape(sigma, d2peers, 1, []);
w = sum(pagemtimes(pagemtimes(s, 'transpose', t_s, 'none'), s), 'all');
end

function w = stabilization_energy(mesh, k, u, u_hat)
% Σ_K <S (u - û), u - û>_∂K
[tPP, ~, ~, tLP, tLL] = boundary_matrix(mesh, k);
d1 = k + 1;
d2 = size(tPP, 1);
N = size(mesh.elements, 2);
w = 0;
for c = 1:2
    u_c = reshape(u((c-1)*d2 + (1:d2), :), d2, 1, N);
    uh_c = u_hat((c-1)*d1 + (1:d1), :);
    uh_K = reshape(uh_c(:, mesh.boundaries), 3*d1, 1, N);
    S = [tPP, -permute(tLP, [2 1 3]); -tLP, tLL];
    x = [u_c; uh_K];
    w = w + sum(pagemtimes(pagemtimes(x, 'transpose', S, 'none'), x), 'all');
end
end
