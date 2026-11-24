% test_interface_stokes_biot.m
% Stokes-Biot transmission rows (InterfaceProblem.stokes_biot) on exact
% states that need no manufactured data. Flat Sigma: y = 1/2, fluid above,
% fitted_split meshes of sizes (1/10, 1/8) (fluid, Biot) with gap h^2, k = 2.
%
% (a) Hydrostatic state, alpha = 1, s0 = 0, zero sources and state,
%     friction_bjs = 1, fluid u = 0 and Biot eta = 0, p = P on the outer
%     boundaries: u_F = 0, sigma_F = -P I, eta = 0, q = 0, p = P,
%     sigma_P = -P I must be reproduced (R_p fixes the fluid pressure
%     level); the free block is invertible (condest finite). The same
%     over 3 steps of "bdf2" (Euler start, then BDF2).
% (b) Steady state with a transported history, alpha = 0: sigma_P = -P I,
%     p = P and the linear displacement eta = -c x, c = P/(2(lambda+mu)),
%     with the state itself as the history. The interface rows see
%     eta(x_bar) through Th(eta_hat) on one side and through
%     h* = eta_hat* + Ext(eta*) on the other, so the state is reproduced
%     only if the history transport onto Sigma is exact.
% (c) The fluid "constrained" closure is rejected.
% (d) The hydrostatic state of (a) on a closed Sigma (circle, R = 1/4, Biot
%     inside with no boundary of its own, fitted meshes (1/10, 1/8),
%     gap 1/64): with s0 = 0 and alpha = 1 the pore-pressure level is
%     free, the coupling keeps its level scalar and pore_pressure_mean = P
%     must reproduce the state (< 1e-8), with a finite condest and a
%     level scalar ~ 0. In (a) the Biot boundary fixes the level and the
%     scalar must be eliminated.
% (e) Flat meshes (1/8, 1/10) (Biot finer at Sigma): interface_problem
%     warns StokesBiotCoupling:biotFiner, and errors with the same id when
%     friction_bjs = 0.
clc; clear;
k = 2; P = 3; h = [1/10, 1/8]; gap = h(2)^2;
E_p = 1.0; nu_p = 0.25; mu_f = 1.0;
interface = FlatCurveHorizontal(0.5);
[fluid_mesh, biot_mesh] = fitted_split(interface, h, gap);
d1 = k + 1;  d2 = nchoosek(k + 2, 2);
phi0 = dubiner_basis(0, 0, 0);          % the constant Dubiner function
l0 = legendre_basis(0, 0.5);            % the constant Legendre function
zero = @(x,y) 0*x.*y;
n_fail = 0;

fluid = IncompressibleElasticityProblem(k, viscosity = mu_f, density = 0, ...
    pressure_mean = "pinned", time_step = 1.0);

%% (a) hydrostatic, Euler from rest
biot = BiotProblem(k, young_modulus = E_p, poisson_ratio = nu_p, biot_alpha = 1, ...
    storage = 0, time_step = 1.0, time_scheme = "euler", ...
    datum_dirichlet = {zero; zero; @(x,y) P + 0*x.*y});
[fs, bs, ~, gs] = StokesBiotCoupling.solve(fluid, biot, fluid_mesh, biot_mesh, ...
    interface, struct.empty, friction_bjs = 1.0);
err = hydrostatic_errors(fs, bs, P, k, phi0, l0, [0, 0]);
fprintf("(a) euler: max dof errors sigma_F %.1e u_F %.1e gamma_F %.1e u_hat_F %.1e | " + ...
    "sigma_P %.1e eta %.1e gamma_P %.1e q %.1e p %.1e eta_hat %.1e p_hat %.1e\n", err);
n_fail = n_fail + check(all(err < 1e-8), "hydrostatic state (euler)");
kappa_free = condest(gs.matrix(gs.dof, gs.dof));
fprintf("(a) condest of the free block %.3e, residue %.1e\n", kappa_free, gs.residue);
n_fail = n_fail + check(isfinite(kappa_free) && kappa_free < 1e16, ...
    "free block invertible (condest finite)");
n_fail = n_fail + check(~ismember(numel(gs.solution), gs.dof), ...
    "the Biot boundary fixes the pore-pressure level: its scalar is eliminated");

%% (a') hydrostatic, 3 steps of bdf2 (the first an Euler step)
biot.Data.time_scheme = "bdf2";
state = struct.empty;
for n = 1:3
    [fs, bs, state] = StokesBiotCoupling.solve(fluid, biot, fluid_mesh, biot_mesh, ...
        interface, struct.empty, friction_bjs = 1.0, state = state);
    err = hydrostatic_errors(fs, bs, P, k, phi0, l0, [0, 0]);
    fprintf("(a) bdf2 step %d: worst dof error %.1e\n", n, max(err));
    n_fail = n_fail + check(all(err < 1e-8), sprintf("hydrostatic state (bdf2 step %d)", n));
end

%% (b) steady state with a transported history, alpha = 0
[lambda_p, mu_p] = LinearElasticityProblem.get_lame_parameters(E_p, nu_p);
c = P/(2*(lambda_p + mu_p));
eta_x = @(x,y) -c*x + 0*y;  eta_y = @(x,y) -c*y + 0*x;
biot_b = BiotProblem(k, young_modulus = E_p, poisson_ratio = nu_p, biot_alpha = 0, ...
    storage = 0, time_step = 1.0, time_scheme = "euler", ...
    datum_dirichlet = {eta_x; eta_y; @(x,y) P + 0*x.*y});
n_el = size(biot_mesh.elements, 2);
dS = PEERSBasis(k).Dimension;
sigma_ex = zeros(dS, n_el);
sigma_ex([1, 3*d2 + 1], :) = -P/phi0;
p_ex = zeros(d2, n_el);  p_ex(1, :) = P/phi0;
skeleton = [project(Datum(eta_x), biot_mesh, LegendreBasis(k), order = k).dof; ...
            project(Datum(eta_y), biot_mesh, LegendreBasis(k), order = k).dof];
state_b = struct("time", 0, "step", 1, "fluid", struct.empty, "diagnostics", struct.empty, ...
    "biot", struct("displacement_skeleton", skeleton, "stress", sigma_ex, ...
        "rotation", zeros(d2, n_el), "pressure", p_ex, ...
        "pressure_skeleton", [], "displacement_skeleton_previous", [], ...
        "stress_previous", [], "rotation_previous", [], "pressure_previous", []));
[fs, bs] = StokesBiotCoupling.solve(fluid, biot_b, fluid_mesh, biot_mesh, ...
    interface, struct.empty, friction_bjs = 1.0, state = state_b);
err = hydrostatic_errors(fs, bs, P, k, phi0, l0, [c, c]);
fprintf("(b) steady history: worst dof error %.1e (eta %.1e, eta_hat %.1e, u_F %.1e)\n", ...
    max(err), err(6), err(10), err(2));
n_fail = n_fail + check(all(err < 1e-8), "steady state with transported history");

%% (d) closed Sigma, hydrostatic, alpha = 1, s0 = 0: no boundary fixes the
% pore-pressure level, so the level scalar is kept and pore_pressure_mean
% = P selects the state
circle = CircularCurve([0.5; 0.5], 0.25);
[fluid_d, biot_d] = interface_meshes(circle, "fitted", 1/10, gap = 1/64, h_solid = 1/8);
biot = BiotProblem(k, young_modulus = E_p, poisson_ratio = nu_p, biot_alpha = 1, ...
    storage = 0, time_step = 1.0, time_scheme = "euler");
[fs, bs, ~, gs] = StokesBiotCoupling.solve(fluid, biot, fluid_d, biot_d, ...
    circle, struct.empty, friction_bjs = 1.0, pore_pressure_mean = P);
err = hydrostatic_errors(fs, bs, P, k, phi0, l0, [0, 0]);
kappa_free = condest(gs.matrix(gs.dof, gs.dof));
beta = gs.solution(end);
fprintf("(d) closed Sigma, s0 = 0: worst dof error %.1e, condest %.3e, level scalar %.1e\n", ...
    max(err), kappa_free, beta);
n_fail = n_fail + check(ismember(numel(gs.solution), gs.dof), ...
    "closed Sigma, s0 = 0: the pore-pressure level scalar is kept");
n_fail = n_fail + check(all(err < 1e-8), "hydrostatic state on a closed Sigma, s0 = 0");
n_fail = n_fail + check(isfinite(kappa_free) && kappa_free < 1e16, ...
    "closed Sigma, s0 = 0: free block invertible (condest finite)");
n_fail = n_fail + check(abs(beta) < 1e-8, "closed Sigma, s0 = 0: level scalar ~ 0 for compatible data");

%% (e) a Biot mesh finer at Sigma: a warning, and an error with friction_bjs = 0
[fluid_e, biot_e] = fitted_split(interface, [1/8, 1/10], gap);
lastwarn("");
StokesBiotCoupling.interface_problem(fluid, biot, fluid_e, biot_e, interface, ...
    struct.empty, friction_bjs = 1.0);
[~, warn_id] = lastwarn();
n_fail = n_fail + check(warn_id == "StokesBiotCoupling:biotFiner", ...
    "finer Biot mesh: StokesBiotCoupling:biotFiner warning");
n_fail = n_fail + check(throws(@() StokesBiotCoupling.interface_problem(fluid, biot, ...
    fluid_e, biot_e, interface, struct.empty, friction_bjs = 0.0), ...
    "StokesBiotCoupling:biotFiner"), "finer Biot mesh with friction_bjs = 0 rejected");
fprintf("(e) finer Biot mesh: warning %s\n", warn_id);

%% (c) the constrained fluid closure is rejected
Pk = DubinerBasis(k);
fluid_c = FEProblem(Order = k, Spaces = [ ...
    FESpace([PEERSBasis(k), Pk.to_vector(), Pk, DubinerBasis(0)], Condense = true), ...
    FESpace([LegendreBasis(k).to_vector(2), DubinerBasis(0), ConstantBasis()])]);
biot_pb = get_declarative_problem(biot, biot_mesh, k);
n_fail = n_fail + check(throws(@() InterfaceProblem.stokes_biot(fluid_c, biot_pb, interface), ...
    "InterfaceProblem:stokesBiotClosure"), "constrained fluid closure rejected");

fprintf("%d check(s) failed.\n", n_fail);
if n_fail > 0
    error("test_interface_stokes_biot:failed", "%d checks failed", n_fail);
end

function err = hydrostatic_errors(fs, bs, P, k, phi0, l0, grad)
% max dof errors against sigma_F = sigma_P = -P I, u_F = 0, q = 0, p = P,
% eta = -(grad(1) x, grad(2) y) (zero rotation, linear displacement)
d2 = nchoosek(k + 2, 2);  d1 = k + 1;
iso = @(F) F.dof - full(sparse([1, 3*d2 + 1], 1, -P/phi0, size(F.dof, 1), 1));
f = fs.Unknowns.field;  b = bs.Unknowns.field;
uh = fs.Unknowns.hybrid(1).dof;
bh = reshape(bs.Unknowns.hybrid(1).dof, 3*d1, []);
p_ex = zeros(d2, 1);  p_ex(1) = P/phi0;
% the linear displacement: eta_x = -grad(1) x, eta_y = -grad(2) y
eta = b(2);
eta_ex = [interpolate_linear(eta.mesh, k, -grad(1), 1); ...
          interpolate_linear(eta.mesh, k, -grad(2), 2)];
eta_hat_ex = [face_linear(eta.mesh, k, -grad(1), 1); face_linear(eta.mesh, k, -grad(2), 2)];
p_hat_ex = zeros(d1, 1);  p_hat_ex(1) = P/l0;
err = [max(abs(iso(f(1))), [], 'all'), max(abs(f(2).dof), [], 'all'), ...
       max(abs(f(3).dof), [], 'all'), max(abs(uh), [], 'all'), ...
       max(abs(iso(b(1))), [], 'all'), max(abs(b(2).dof - eta_ex), [], 'all'), ...
       max(abs(b(3).dof), [], 'all'), max(abs(b(4).dof), [], 'all'), ...
       max(abs(b(5).dof - p_ex), [], 'all'), ...
       max(abs(bh(1:2*d1, :) - eta_hat_ex), [], 'all'), ...
       max(abs(bh(2*d1 + (1:d1), :) - p_hat_ex), [], 'all')];
end

function dof = interpolate_linear(mesh, k, slope, comp)
% Dubiner dofs (d2 x n_el) of slope * x_comp, by L2 projection
if slope == 0
    dof = zeros(nchoosek(k + 2, 2), size(mesh.elements, 2));
    return
end
proj = project(Datum(@(x,y) slope*(comp == 1)*x + slope*(comp == 2)*y), ...
    mesh, DubinerBasis(k), order = k + 1);
dof = proj.dof;
end

function dof = face_linear(mesh, k, slope, comp)
% Legendre dofs (d1 x n_faces) of slope * x_comp, by L2 projection
if slope == 0
    dof = zeros(k + 1, size(mesh.faces, 2));
    return
end
proj = project(Datum(@(x,y) slope*(comp == 1)*x + slope*(comp == 2)*y), ...
    mesh, LegendreBasis(k), order = k + 1);
dof = proj.dof;
end
