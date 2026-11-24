% test_interface_fsi_postprocess.m
% divergence_free_velocity on the fluid solution of the coupled
% InterfaceProblem path (FluidStructureCoupling.interface_problem): the
% declarative fluid problem must carry the coefficients it needs
% (Data.viscosity for the PEERS stress reconstruction). On the flat fitted
% FSI problem, a Stokes coupled solve gives u* in P_{k+1} with
% max_K ||div u*||_{L2(K)} <= 1e-10; a second, Oseen coupled solve
% convected by that u* (history = the HDG velocity) must give a finite
% solution whose u* is again divergence-free to 1e-10.
clc; clear;
n_fail = 0;
k = 2;
mu_f = 1.0; E_s = 1.0; nu_s = 0.25; dt = 1e-2;
interface = FlatCurveHorizontal(0.5);
h = 1/8;
tol = 1e-10;

[fp, sp, srcs] = FluidStructureCoupling.load_problems(k, interface, 0.0, dt, ...
    "x^2*y^2", "x^2 + y^2", ["x^2*y"; "x*y^2"], ...
    0.0, mu_f, 0.0, E_s, nu_s, 0.0, pressure_mean = "pinned");
fp.Data.density = 1.0;
% fitted_split: mesh 1 (upper) is the fluid, mesh 2 (lower) the solid
[fluid_mesh, solid_mesh] = fitted_split(interface, [h, h], h^2);

%% Stokes
fp.Data.convective_velocity = {};
coupled = FluidStructureCoupling.interface_problem(fp, sp, fluid_mesh, solid_mesh, interface, ...
    srcs{:}, "steady", Order = k);
fs = solve(coupled, k);
[u_star, div_star] = divergence_free_velocity(fs);
fprintf("Stokes: u* %d x %d per component (P_%d: %d), max div u* = %.3e\n", ...
    size(u_star{1}), k + 1, nchoosek(k + 3, 2), div_star);
if size(u_star{1}, 1) ~= nchoosek(k + 3, 2) || size(u_star{1}, 2) ~= size(fluid_mesh.elements, 2)
    fprintf(2, "FAILED: u* is not P_%d per element\n", k + 1); n_fail = n_fail + 1;
end
if ~(div_star <= tol)
    fprintf(2, "FAILED: Stokes u* not divergence-free (%.3e)\n", div_star); n_fail = n_fail + 1;
end

%% Oseen, convected by u*
fp.Data.convective_velocity = u_star;
coupled = FluidStructureCoupling.interface_problem(fp, sp, fluid_mesh, solid_mesh, interface, ...
    srcs{:}, "steady", Order = k, state = struct("fluid", struct("velocity", fs.Unknowns.field(2).dof)));
[fs2, ss2] = solve(coupled, k);
dofs = {fs2.Unknowns.field.dof, fs2.Unknowns.hybrid.dof, ss2.Unknowns.field.dof};
finite = all(cellfun(@(a) all(isfinite(a), "all"), dofs));
[~, div_star2] = divergence_free_velocity(fs2);
u2 = fs2.Unknowns.field(2).dof;
fprintf("Oseen: finite %d, ||u_h|| = %.6e, ||u_h - u_h(Stokes)||/||u_h|| = %.3e, max div u* = %.3e\n", ...
    finite, norm(u2, "fro"), norm(u2 - fs.Unknowns.field(2).dof, "fro")/norm(u2, "fro"), div_star2);
if ~finite
    fprintf(2, "FAILED: the convected solve is not finite\n"); n_fail = n_fail + 1;
end
if ~(div_star2 <= tol)
    fprintf(2, "FAILED: Oseen u* not divergence-free (%.3e)\n", div_star2); n_fail = n_fail + 1;
end

fprintf("%d check(s) failed.\n", n_fail);
if n_fail > 0, error("test_interface_fsi_postprocess:failed", "%d checks failed", n_fail); end
