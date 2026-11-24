% test_biot_convergence.m
% The single-domain quasi-static Biot problem (BiotProblem, BiotMatrices):
%   (a) with biot_alpha = 0 and storage = 0 the solution is the direct sum
%       of the LinearElasticityMatrices and DarcyMatrices solves, dof for
%       dof;
%   (b) polynomial data of degree <= k (alpha = 1, s0 = 1, one step from
%       rest) are reproduced exactly;
%   (c) h-convergence of one step from rest: rate between the coarsest
%       and the finest level >= k + 1 - 0.25 for sigma, eta, q, p (gamma
%       printed).
clc; clear;
n_fail = 0;

%% (a) direct sum
k = 2; E = 1; nu = 0.3; kappa = [2, 0.5; 0.5, 1];
mesh = StructuredSimplicialGrid(8, 8);
mms = BiotProblem.get_manufactured_solution(["sin(pi*x)*y"; "x*exp(y)"], ...
    "cos(pi*x)*sin(pi*y)", young_modulus = E, poisson_ratio = nu, ...
    biot_alpha = 0, storage = 0, diffusion = kappa);
data = mms.at(1);
biot = BiotProblem(k, young_modulus = E, poisson_ratio = nu, biot_alpha = 0, ...
    storage = 0, diffusion = kappa, source_term = data.f_P, source_mass = data.g_P, ...
    datum_dirichlet = [data.u_P; {data.p_P}]);
sol = solve(biot, mesh, k);

Pk = DubinerBasis(k);
w = TestFunction(Pk);
dx = Measure(mesh, 'dx');
% elasticity alone
pe = FEProblem(Order = k, Spaces = [ ...
    FESpace([PEERSBasis(k), Pk.to_vector(), Pk], Condense = true), ...
    FESpace(LegendreBasis(k).to_vector(2))]);
pe.Mesh = mesh;
[lhs, flux] = LinearElasticityMatrices(k, mesh, E, nu);
pe.Operator = AffineOperator(lhs) + flux + IndexedMatrix(IndexedArray(cat(1, ...
    double(Datum(data.f_P{1}) * w * dx), double(Datum(data.f_P{2}) * w * dx)), [-1, 2]));
pe.BoundaryConditions = StrongEssentialBoundaryConditions(Datum(data.u_P), ...
    "Dirichlet", projection_order = k);
sol_e = solve(pe, mesh, k);
% Darcy alone
pd = FEProblem(Order = k, Spaces = [FESpace([Pk.to_vector(), Pk], Condense = true), ...
    FESpace(LegendreBasis(k))]);
% the Biot element's Darcy stabilization: tau_p = ||kappa|| tau
mesh_p = mesh;  mesh_p.tau = norm(kappa)*mesh.tau;
pd.Mesh = mesh_p;
[lhs, flux] = DarcyMatrices(k, mesh_p, kappa);
pd.Operator = AffineOperator(lhs) + flux + IndexedMatrix(IndexedArray( ...
    double(Datum(data.g_P) * w * dx), [-1, 2]));
pd.BoundaryConditions = StrongEssentialBoundaryConditions(Datum(data.p_P), ...
    "Dirichlet", projection_order = k);
sol_d = solve(pd, mesh_p, k);

d1 = k + 1;
diffs = [ ...
    max(abs(sol.Unknowns.field(1).dof - sol_e.Unknowns.field(1).dof), [], 'all'), ...
    max(abs(sol.Unknowns.field(2).dof - sol_e.Unknowns.field(2).dof), [], 'all'), ...
    max(abs(sol.Unknowns.field(3).dof - sol_e.Unknowns.field(3).dof), [], 'all'), ...
    max(abs(sol.Unknowns.field(4).dof - sol_d.Unknowns.field(1).dof), [], 'all'), ...
    max(abs(sol.Unknowns.field(5).dof - sol_d.Unknowns.field(2).dof), [], 'all'), ...
    max(abs(sol.Unknowns.hybrid(1).dof(1:2*d1, :) - sol_e.Unknowns.hybrid(1).dof), [], 'all'), ...
    max(abs(sol.Unknowns.hybrid(1).dof(2*d1 + (1:d1), :) - sol_d.Unknowns.hybrid(1).dof), [], 'all')];
fprintf("(a) direct sum, k = %d: max dof diff sigma %.2e eta %.2e gamma %.2e q %.2e p %.2e eta_hat %.2e p_hat %.2e\n", ...
    k, diffs);
n_fail = n_fail + check(all(diffs < 1e-9), "(a) direct sum");

% the same with a non-unit effective time step: still the direct sum
biot.Data.time_step = 0.37;
sol_dt = solve(biot, mesh, k);
d_dt = max(abs(sol_dt.Unknowns.field(5).dof - sol.Unknowns.field(5).dof), [], 'all');
fprintf("(a) direct sum independent of the time step: max p dof diff %.2e\n", d_dt);
n_fail = n_fail + check(d_dt < 1e-9, "(a) time-step independence");

%% (b) polynomial exactness
names = ["sigma", "eta", "gamma", "q", "p"];
for k = 1:3
    eta = ["(1 + t)*(x^" + k + " + 2*x*y^" + (k-1) + " + y + 1)"; ...
           "y^" + k + " - x^" + (k-1) + "*y + 3*x - t"];
    p = "(2 - t)*(1 + x*y^" + (k-1) + " + 2*y^" + k + " - x)";
    dt = 0.5;
    mms = BiotProblem.get_manufactured_solution(eta, p, young_modulus = 3, ...
        poisson_ratio = 0.3, biot_alpha = 1, storage = 1, time_step = dt, ...
        rate = "from_rest");
    data = mms.at(dt);
    biot = BiotProblem(k, young_modulus = 3, poisson_ratio = 0.3, biot_alpha = 1, ...
        storage = 1, time_step = dt, source_term = data.f_P, source_mass = data.g_P, ...
        datum_dirichlet = [data.u_P; {data.p_P}]);
    mesh = StructuredSimplicialGrid(4, 4);
    sol = solve(biot, mesh, k);
    err = biot_errors(sol, mesh, data);
    fprintf("(b) k = %d exactness: L2 errors %s = %s\n", k, strjoin(names, "/"), ...
        strjoin(compose("%.2e", err), " "));
    n_fail = n_fail + check(all(err < 1e-9), "(b) polynomial exactness k = " + k);
end

%% (b') the history block: one step from an exact state, data linear in t
% (the backward difference is then exact), assembled from BiotMatrices
% with history = the discrete (sigma, p) at t0
k = 2; dt = 0.25; t0 = 0.5; E = 3; nu = 0.3; s0 = 0.7;
eta = ["(1 + t)*(x^2 + 2*x*y + y + 1)"; "y^2 - x*y + 3*x*t"];
p = "(2 - t)*(1 + x*y + 2*y^2 - x)";
opts = {"young_modulus", E, "poisson_ratio", nu, "biot_alpha", 1, "storage", s0};
mesh = StructuredSimplicialGrid(4, 4);
data0 = BiotProblem.get_manufactured_solution(eta, p, opts{:}, time_step = dt, ...
    rate = "from_rest").at(t0);
sol0 = solve(BiotProblem(k, opts{:}, time_step = dt, source_term = data0.f_P, ...
    source_mass = data0.g_P, datum_dirichlet = [data0.u_P; {data0.p_P}]), mesh, k);
data1 = BiotProblem.get_manufactured_solution(eta, p, opts{:}).at(t0 + dt);
fe = get_declarative_problem(BiotProblem(k, opts{:}, time_step = dt, ...
    source_term = data1.f_P, source_mass = data1.g_P, ...
    datum_dirichlet = [data1.u_P; {data1.p_P}]), mesh, k);
[lhs, flux, rhs] = BiotMatrices(k, mesh, E, nu, eye(2), biot_alpha = 1, storage = s0, ...
    effective_time_step = dt, history = struct("stress", sol0.Unknowns.field(1).dof, ...
    "pressure", sol0.Unknowns.field(5).dof));
fe.Operator.Terms(1:2) = {lhs, flux};      % the zero-history operator had no rhs
fe.Operator = fe.Operator + rhs;
err = biot_errors(solve(fe, mesh, k), mesh, data1);
fprintf("(b') k = %d step from an exact state: L2 errors %s = %s\n", k, strjoin(names, "/"), ...
    strjoin(compose("%.2e", err), " "));
n_fail = n_fail + check(all(err < 1e-9), "(b') history block");

%% (c) convergence
E = 1; nu = 0.3; alpha = 1; s0 = 0.1; dt = 1;
mms = BiotProblem.get_manufactured_solution(["sin(pi*x)*cos(pi*y)"; "x*exp(y)*(1 + t)"], ...
    "cos(pi*x)*sin(pi*y) + t*x*y", young_modulus = E, poisson_ratio = nu, ...
    biot_alpha = alpha, storage = s0, time_step = dt, rate = "from_rest");
data = mms.at(dt);
biot = BiotProblem([], young_modulus = E, poisson_ratio = nu, biot_alpha = alpha, ...
    storage = s0, time_step = dt, source_term = data.f_P, source_mass = data.g_P, ...
    datum_dirichlet = [data.u_P; {data.p_P}]);
for k = 1:3
    if k < 3, levels = [4 8 16 32]; else, levels = [4 8 16]; end
    err = zeros(numel(levels), 5);
    for i = 1:numel(levels)
        mesh = StructuredSimplicialGrid(levels(i), levels(i));
        sol = solve(biot, mesh, k);
        err(i, :) = biot_errors(sol, mesh, data);
    end
    eoc = [nan(1, 5); log2(err(1:end-1, :)./err(2:end, :))];
    fprintf("(c) k = %d\n  %4s  %s\n", k, "n", ...
        strjoin(compose("%10s  %5s", names.', repmat("eoc", 5, 1)).', "  "));
    for i = 1:numel(levels)
        fprintf("  %4d  %s\n", levels(i), ...
            strjoin(compose("%10.3e  %5.2f", err(i, :).', eoc(i, :).').', "  "));
    end
    rate = log(err(1, :)./err(end, :))/log(levels(end)/levels(1));
    fprintf("  rate coarsest-finest: %s\n", strjoin(names + " " + compose("%.2f", rate), ", "));
    asserted = [1 2 4 5];    % sigma, eta, q, p; gamma printed only
    n_fail = n_fail + check(all(rate(asserted) >= k + 1 - 0.25), "(c) rates k = " + k);
end

fprintf("%d check(s) failed.\n", n_fail);
if n_fail > 0, error("test_biot_convergence:failed", "%d checks failed", n_fail); end

function err = biot_errors(sol, mesh, data)
% L2 errors of (sigma, eta, gamma, q, p)
err = [compute_l2_error(sol.Unknowns.field(1), mesh, Datum(data.sigma_P)), ...
       compute_l2_error(sol.Unknowns.field(2), mesh, Datum(data.u_P)), ...
       compute_l2_error(sol.Unknowns.field(3), mesh, Datum(data.gamma_P)), ...
       compute_l2_error(sol.Unknowns.field(4), mesh, Datum(data.q_P)), ...
       compute_l2_error(sol.Unknowns.field(5), mesh, Datum(data.p_P))];
end
