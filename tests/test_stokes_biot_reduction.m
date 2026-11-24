% test_stokes_biot_reduction.m
% Stokes-Biot reduces to the FSI rows as kappa -> 0 (InterfaceProblem path).
%
% With alpha = 0 the Biot momentum is linear elasticity, and with s0 = 0
% and kappa = kappa0 I the numerical Darcy flux E(q).nu + tau_p (p - p_hat)
% in R_n is O(kappa0) because tau_p = stabilization_darcy ||kappa|| tau
% scales with kappa0 (with tau_p = O(1) its jump part would not vanish):
% the rows R_f, R_n, R_t then become the fsi traction row and the normal /
% tangential rows of the "navier" law, with slip_coefficient =
% friction_bjs / cdt. Flat interface, k = 2, h = 1/8, meshes 0.8 h (fluid) / h (Biot),
% gap h^2, Delta t = 1, friction_bjs = 10; the data are those of
% test_interface_fsi_slip.m (FluidStructureCoupling.load_problems,
% "pinned"): T_1 = T1, T_2 = T2, T_4 = 10*T2 (the fsi Navier row carries
% a*T2 on its right-hand side), T_3 = 0; Biot eta = u_S and p = 0 on the
% outer boundary. Reference: the same fluid and solid on
% FluidStructureCoupling.solve, "steady", interface_law = "navier",
% slip_coefficient = 10 (= InterfaceProblem.fsi with time_step = 1,
% pressure_mean = "pinned").
%
% For kappa0 = 1e-2, 1e-4, 1e-6, 1e-8 the max dof difference of (u_F, eta)
% must decrease with kappa0 (or both values sit at round-off, < 1e-10) and
% be < 1e-6 at kappa0 = 1e-8. The observed slope against kappa0 (expected
% ~1) and max |q| are printed; a slope well below 1 is to be reported.
%
% Skips (prints SKIPPED, 0 checks) while the Stokes-Biot code is not on
% the path.
clc; clear;

k = 2; mu_f = 1.0; E_s = 1.0; nu_s = 0.25; dt = 1.0; h = 1/8;
friction = 10;
kappas = [1e-2, 1e-4, 1e-6, 1e-8];
floor_roundoff = 1e-10;

n_fail = 0;
need = ["BiotProblem", "StokesBiotCoupling.solve"];
missing = need(arrayfun(@(f) isempty(which(f)), need));
if ~isempty(missing)
    fprintf("SKIPPED test_stokes_biot_reduction: missing %s\n", strjoin(missing, ", "));
    fprintf("%d check(s) failed.\n", n_fail);
    return
end

interface = FlatCurveHorizontal(0.5);
[fluid_problem, solid_problem, interface_sources, pd] = FluidStructureCoupling.load_problems( ...
    k, interface, 0.0, dt, "x^2*y^2", "x^2 + y^2", ["x^2*y"; "x*y^2"], ...
    0.0, mu_f, 0.0, E_s, nu_s, 0.0, pressure_mean = "pinned");
[fluid_mesh, solid_mesh] = fitted_split(interface, [0.8*h, h], h^2);

%% reference: FSI, Navier law with a = friction
[fs_ref, ss_ref] = FluidStructureCoupling.solve(fluid_problem, solid_problem, ...
    fluid_mesh, solid_mesh, interface, interface_sources{1}, interface_sources{2}, ...
    "steady", Order = k, interface_law = "navier", slip_coefficient = friction);
uF_ref = fs_ref.Unknowns.field(2).dof;
eta_ref = ss_ref.Unknowns.field(2).dof;

%% Stokes-Biot with alpha = s0 = 0
zero = @(x, y) 0*x.*y;
T2 = interface_sources{2};
sources = struct("T_1", {interface_sources{1}}, "T_2", {T2}, "T_3", zero, ...
    "T_4", {cellfun(@(g) @(x, y) friction*g(x, y), T2, UniformOutput = false)});
d = zeros(size(kappas));
for i = 1:numel(kappas)
    biot_problem = BiotProblem(k, young_modulus = E_s, poisson_ratio = nu_s, ...
        biot_alpha = 0.0, storage = 0.0, diffusion = kappas(i)*eye(2), ...
        time = 0.0, time_step = dt, time_scheme = "euler", ...
        source_term = pd.f_S, ...
        datum_dirichlet = {pd.u_S{1}; pd.u_S{2}; zero});
    [fs, bs] = StokesBiotCoupling.solve(fluid_problem, biot_problem, ...
        fluid_mesh, solid_mesh, interface, sources, Order = k, ...
        friction_bjs = friction);
    d(i) = max(max(abs(fs.Unknowns.field(2).dof - uF_ref), [], "all"), ...
               max(abs(bs.Unknowns.field(2).dof - eta_ref), [], "all"));
    q_max = max(abs(bs.Unknowns.field(4).dof), [], "all");
    fprintf("kappa0 %.0e: max|(u_F, eta) - fsi| %.3e   max|q| %.3e\n", kappas(i), d(i), q_max);
    if i > 1
        slope = log(d(i - 1)/d(i)) / log(kappas(i - 1)/kappas(i));
        fprintf("              slope against kappa0 %.2f (expected ~1)\n", slope);
        n_fail = n_fail + check(d(i) < d(i - 1) || max(d(i - 1), d(i)) < floor_roundoff, ...
            sprintf("kappa0 = %g does not decrease the difference (%.3e >= %.3e)", ...
                kappas(i), d(i), d(i - 1)));
    end
end
n_fail = n_fail + check(d(end) < 1e-6, ...
    sprintf("kappa0 = %g: difference to the fsi rows %.3e >= 1e-6", kappas(end), d(end)));

fprintf("%d check(s) failed.\n", n_fail);
if n_fail > 0
    error("test_stokes_biot_reduction:failed", "%d checks failed", n_fail);
end
