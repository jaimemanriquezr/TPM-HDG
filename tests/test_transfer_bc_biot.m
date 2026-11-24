% test_transfer_bc_biot.m
% TPM transfer boundary conditions for the HDG Biot element (BiotMatrices,
% trace (eta_hat_x, eta_hat_y, p_hat)) with the biot_path_gradient
% reconstruction, one Euler step from rest (Delta t = 1).
%
% Stage 1 (consistency): fitted 8x8 square, identity mapping -> transfer
%   reduces to weak Dirichlet; all fields must match dof for dof.
% Stage 2a (reconstruction): carved disk, n = 8, alpha = s0 = 0,
%   anisotropic kappa: equal dof for dof to the elasticity and Darcy
%   transfer factories solved separately (nonzero transfer paths).
% Stage 2 (unfitted): disk R = 0.75 carved from [-1,1]^2 by immerse_mesh,
%   data transferred radially from the circle, k = 2, n = 8..64,
%   alpha = 1, s0 = 0.1, kappa = I: rate over all levels >= k + 1 - 0.25
%   for eta and p. When BiotProblem is on the path, its transfer branch
%   (transfer_mapping) is compared dof for dof with the hand-built problem
%   on the coarsest level.
% Stage 3 (reported, not asserted): the same with kappa = 1e-4 I.
% Stage 4: the disk meshed to the linear interpolant of the circle
%   (interface_meshes "fitted", inner mesh), vertex gap delta = 0 and
%   delta = h^2, h = 1/4..1/32, alpha = 1, s0 = 0.1, kappa = I: the same
%   rate threshold for eta and p.
clc; clear;
n_fail = 0;

k = 2;
E = 1.0;
nu = 0.3;
alpha = 1.0;
s0 = 0.1;
dt = 1.0;
rate_min = k + 1 - 0.25;

R = 0.75;
phi = @(X) sqrt(X(1, :, :).^2 + X(2, :, :).^2) - R;
radial_map = @(X) R * X ./ vecnorm(X);

%% Stage 1: fitted square, identity mapping == weak Dirichlet
data = biot_data(E, nu, alpha, s0, eye(2), dt);
mesh = StructuredSimplicialGrid(8, 8);
problem = build_problem(mesh, k, E, nu, alpha, s0, eye(2), dt, data);
problem.BoundaryConditions = WeakEssentialBoundaryConditions(data.dirichlet, "Dirichlet");
sol_weak = solve(problem, mesh, k);
problem.BoundaryConditions = TransferEssentialBoundaryConditions(data.dirichlet, "Dirichlet", ...
    @(X) X, path_gradient = biot_path_gradient(k, E, nu, alpha, eye(2)));
sol_transfer = solve(problem, mesh, k);
d = 0;
for f = 1:5
    d = max(d, max(abs(sol_weak.Unknowns.field(f).dof - sol_transfer.Unknowns.field(f).dof), [], 'all'));
end
d = max(d, max(abs(sol_weak.Unknowns.hybrid(1).dof - sol_transfer.Unknowns.hybrid(1).dof), [], 'all'));
fprintf("Stage 1 (identity mapping vs weak BC): max dof diff = %.3e\n", d);
if ~(d < 1e-9)
    fprintf("Stage 1 FAILED\n"); n_fail = n_fail + 1;
end

%% Stage 2a: decoupled oracle at nonzero transfer paths
% alpha = 0, s0 = 0: BiotMatrices is the direct sum of the elasticity and
% Darcy elements, and biot_path_gradient must reproduce the .elasticity
% and .darcy transfer factories (anisotropic kappa exercises kappa^{-T} m)
kappa = [2 0.5; 0.5 1];
data = biot_data(E, nu, 0, 0, kappa, dt);
mesh = immerse_mesh(StructuredSimplicialGrid(8, 8, -1, 1, -1, 1), phi);
problem = build_problem(mesh, k, E, nu, 0, 0, kappa, dt, data);
problem.BoundaryConditions = TransferEssentialBoundaryConditions(data.dirichlet, ...
    "Dirichlet", radial_map, path_gradient = biot_path_gradient(k, E, nu, 0, kappa));
sol_biot = solve(problem, mesh, k);
Pk = DubinerBasis(k);
w = TestFunction(Pk);
el = FEProblem(Spaces = [FESpace([PEERSBasis(k), Pk.to_vector(), Pk], Condense = true), ...
    FESpace(LegendreBasis(k).to_vector(2))]);
[lhs, flux] = LinearElasticityMatrices(k, mesh, E, nu);
F1 = double(data.f{1} * w * Measure(mesh, 'dx'));
F2 = double(data.f{2} * w * Measure(mesh, 'dx'));
el.Operator = AffineOperator(lhs, flux) + IndexedMatrix(IndexedArray(cat(1, F1, F2), [-1, 2]));
el.BoundaryConditions = TransferEssentialBoundaryConditions.elasticity( ...
    data.u, "Dirichlet", radial_map, E, nu);
sol_el = solve(el, mesh, k);
da = FEProblem(Spaces = [FESpace([Pk.to_vector(), Pk], Condense = true), FESpace(LegendreBasis(k))]);
% the Biot element's Darcy stabilization: tau_p = ||kappa|| tau
mesh_p = mesh;  mesh_p.tau = norm(kappa)*mesh.tau;
[lhs, flux] = DarcyMatrices(k, mesh_p, kappa);
G = double(data.g * w * Measure(mesh, 'dx'));
da.Operator = AffineOperator(lhs, flux) + IndexedMatrix(IndexedArray(G, [-1, 2]));
da.BoundaryConditions = TransferEssentialBoundaryConditions.darcy( ...
    data.p, "Dirichlet", radial_map, kappa);
sol_da = solve(da, mesh_p, k);
d = 0;
for f = 1:3
    d = max(d, max(abs(sol_biot.Unknowns.field(f).dof - sol_el.Unknowns.field(f).dof), [], 'all'));
end
for f = 1:2
    d = max(d, max(abs(sol_biot.Unknowns.field(3 + f).dof - sol_da.Unknowns.field(f).dof), [], 'all'));
end
fprintf("Stage 2a (alpha = s0 = 0 vs elasticity + Darcy transfer factories): max dof diff = %.3e\n", d);
if ~(d < 1e-9)
    fprintf("Stage 2a FAILED\n"); n_fail = n_fail + 1;
end

%% Stages 2 and 3: carved disk, radial transfer
for kappa0 = [1, 1e-4]
    kappa = kappa0*eye(2);
    data = biot_data(E, nu, alpha, s0, kappa, dt);
    n_ref = 4;
    [err_u, err_p, err_q, h_list] = deal(zeros(n_ref, 1));
    fprintf("\n=== carved disk, kappa = %.0e I ===\n", kappa0);
    for i = 1:n_ref
        n = 8 * 2^(i-1);
        h_list(i) = 2/n;
        mesh = immerse_mesh(StructuredSimplicialGrid(n, n, -1, 1, -1, 1), phi);
        problem = build_problem(mesh, k, E, nu, alpha, s0, kappa, dt, data);
        problem.BoundaryConditions = TransferEssentialBoundaryConditions(data.dirichlet, ...
            "Dirichlet", radial_map, path_gradient = biot_path_gradient(k, E, nu, alpha, kappa));
        solution = solve(problem, mesh, k);
        [err_u(i), err_p(i), err_q(i)] = biot_errors(solution, mesh, data);
        fprintf("n = %3d, h = %.4f, L2 err eta = %.4e, p = %.4e, q = %.4e\n", ...
            n, h_list(i), err_u(i), err_p(i), err_q(i));
        if i == 1 && kappa0 == 1
            n_fail = n_fail + class_route_check(solution, mesh, k, E, nu, alpha, s0, kappa, dt, data, radial_map);
        end
    end
    [rate_u, rate_p] = report_rates(err_u, err_p, err_q, h_list, k);
    if kappa0 == 1 && min(rate_u, rate_p) < rate_min
        fprintf("Stage 2 rate FAILED\n"); n_fail = n_fail + 1;
    end
end

%% Stage 4: disk fitted to the linear interpolant, gap 0 and h^2
data = biot_data(E, nu, alpha, s0, eye(2), dt);
curve = CircularCurve([0; 0], R);
for gap_kind = ["zero", "h2"]
    n_ref = 4;
    [err_u, err_p, err_q, h_list] = deal(zeros(n_ref, 1));
    fprintf("\n=== fitted disk, gap %s ===\n", gap_kind);
    for i = 1:n_ref
        h = 0.25 * 2^(1-i);
        gap = (gap_kind == "h2") * h^2;
        [~, mesh] = interface_meshes(curve, "fitted", h, box = [-1 1 -1 1], gap = gap);
        h_list(i) = h;
        problem = build_problem(mesh, k, E, nu, alpha, s0, eye(2), dt, data);
        problem.BoundaryConditions = TransferEssentialBoundaryConditions(data.dirichlet, ...
            "Interface", radial_map, path_gradient = biot_path_gradient(k, E, nu, alpha, eye(2)));
        solution = solve(problem, mesh, k);
        [err_u(i), err_p(i), err_q(i)] = biot_errors(solution, mesh, data);
        fprintf("h = %.5f, n_el = %5d, L2 err eta = %.4e, p = %.4e, q = %.4e\n", ...
            h, size(mesh.elements, 2), err_u(i), err_p(i), err_q(i));
    end
    [rate_u, rate_p] = report_rates(err_u, err_p, err_q, h_list, k);
    if min(rate_u, rate_p) < rate_min
        fprintf("Stage 4 rate FAILED (gap %s)\n", gap_kind); n_fail = n_fail + 1;
    end
end

fprintf("%d check(s) failed.\n", n_fail);

%% ------------------------------------------------------------------ %%
function data = biot_data(E, nu, alpha, s0, kappa, dt)
% one Euler step from rest: d_t w = w/dt; total stress, Darcy flux and the
% sources f = -div sigma, g = (s0 p + alpha div eta)/dt + div q
[lambda, mu] = LinearElasticityProblem.get_lame_parameters(E, nu);
syms x y
u = [sin(pi*x) .* sin(pi*y); cos(pi*x) .* cos(pi*y)];
p = cos(pi*x) .* sin(pi*y) + x;
gu = jacobian(u, [x y]);
div_u = trace(gu);
sig = mu*(gu + gu.') + lambda*div_u*eye(2) - alpha*p*eye(2);
f = -[diff(sig(1,1), x) + diff(sig(1,2), y); diff(sig(2,1), x) + diff(sig(2,2), y)];
q = -kappa * gradient(p, [x y]);
g = (s0*p + alpha*div_u)/dt + diff(q(1), x) + diff(q(2), y);
data.u = Datum(u);
data.p = Datum(p);
data.q = Datum(q);
data.f = {Datum(f(1)), Datum(f(2))};
data.g = Datum(g);
data.dirichlet = Datum([u; p]);
end

function problem = build_problem(mesh, k, E, nu, alpha, s0, kappa, dt, data)
Pk = DubinerBasis(k);
problem = FEProblem(Spaces = [FESpace([PEERSBasis(k), Pk.to_vector(), Pk, Pk.to_vector(), Pk], ...
    Condense = true), FESpace(LegendreBasis(k).to_vector(3))]);
[lhs, flux] = BiotMatrices(k, mesh, E, nu, kappa, biot_alpha = alpha, storage = s0, ...
    effective_time_step = dt);
w = TestFunction(Pk);
F1 = double(data.f{1} * w * Measure(mesh, 'dx'));
F2 = double(data.f{2} * w * Measure(mesh, 'dx'));
G = double(data.g * w * Measure(mesh, 'dx'));
problem.Operator = AffineOperator(lhs, flux) ...
    + IndexedMatrix(IndexedArray(cat(1, F1, F2), [-1, 2])) ...
    + IndexedMatrix(IndexedArray(G, [-1, 5]));
end

function [err_u, err_p, err_q] = biot_errors(solution, mesh, data)
err_u = compute_l2_error(solution.Unknowns.field(2), mesh, data.u);
err_q = compute_l2_error(solution.Unknowns.field(4), mesh, data.q);
err_p = compute_l2_error(solution.Unknowns.field(5), mesh, data.p);
end

function [rate_u, rate_p] = report_rates(err_u, err_p, err_q, h_list, k)
% rates between the coarsest and finest levels (the mean EOC): the
% per-level EOCs of a carved disk oscillate around k+1
fprintf("EOC(eta): %s\n", num2str(compute_eoc(err_u, h_list).', "%.2f "));
fprintf("EOC(p):   %s\n", num2str(compute_eoc(err_p, h_list).', "%.2f "));
fprintf("EOC(q):   %s\n", num2str(compute_eoc(err_q, h_list).', "%.2f "));
r = @(e) log(e(1)/e(end)) / log(h_list(1)/h_list(end));
rate_u = r(err_u);  rate_p = r(err_p);
fprintf("rate over all levels: eta %.2f, p %.2f, q %.2f (pass >= %.2f for eta, p)\n", ...
    rate_u, rate_p, r(err_q), k + 1 - 0.25);
end

function n_fail = class_route_check(solution, mesh, k, E, nu, alpha, s0, kappa, dt, data, mapping)
% BiotProblem's transfer branch (Data.transfer_mapping) against the
% hand-built problem; skipped while the class is not on the path
n_fail = 0;
if exist("BiotProblem", "class") ~= 8
    fprintf("BiotProblem route: SKIPPED (class not on the path)\n");
    return
end
pb = BiotProblem(k, young_modulus = E, poisson_ratio = nu, biot_alpha = alpha, ...
    storage = s0, diffusion = kappa, time_step = dt, ...
    source_term = {data.f{1}.value; data.f{2}.value}, source_mass = data.g.value, ...
    datum_dirichlet = data.dirichlet.value, transfer_mapping = mapping);
sol = solve(pb, mesh, k);
d = 0;
for f = 1:5
    d = max(d, max(abs(sol.Unknowns.field(f).dof - solution.Unknowns.field(f).dof), [], 'all'));
end
fprintf("BiotProblem route vs hand-built problem: max dof diff = %.3e\n", d);
if ~(d < 1e-9)
    fprintf("BiotProblem route FAILED\n"); n_fail = 1;
end
end
