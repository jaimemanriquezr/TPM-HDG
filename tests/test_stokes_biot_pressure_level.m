% test_stokes_biot_pressure_level.m
% The pore-pressure level scalar of the Stokes-Biot coupling
% (InterfaceProblem.stokes_biot, InterfaceProblem.solve) against the
% boundary conditions off Sigma.
%
% Hydrostatic state, P = 3: u_F = 0, sigma_F = sigma_P = -P I, eta = 0,
% q = 0, p = P; zero sources and state, alpha = 1, s0 = 0 unless stated,
% friction_bjs = 1, k = 2. Flat Sigma: y = 1/2, fitted_split meshes
% (0.8 h, h) (fluid above, Biot below), h = 1/8, gap h^2. Fluid u = 0 and
% Biot eta = 0 on the faces tagged "Dirichlet"; the Biot p_hat is natural
% (zero flux) on the whole outer boundary (componentwise specs).
%
% Level free (the scalar must be kept; pore_pressure_mean = P must give
% the state to 1e-8 with a finite condest and a level scalar ~ 0):
%   roller     Biot y = 0: eta_x natural, (sigma n)_x = 0 (n_x = 0)
%   slip       fluid y = 1: u_x natural, (sigma n)_x = 0 (n_x = 0)
%   symmetry   fluid x = 0: u_y natural, (sigma n)_y = 0 (n_y = 0)
%   all three at once
%   s0 = 1e-14 and alpha = 1 - 1e-13 (below the activation tolerance)
% Level fixed (the scalar must be eliminated; pore_pressure_mean = 0 is
% then ignored and the state is reproduced from the traction datum):
%   Biot x = 1: eta_x natural, (sigma n)_x = -P (n_x = 1)
%   fluid x = 1: u_x and u_y natural, (sigma n) = (-P, 0)
%   s0 = 1e-6 (above the tolerance; zero datum: the state p = 0)
clc; clear;
k = 2; h = 1/8; P = 3;
d1 = k + 1;  d2 = nchoosek(k + 2, 2);
phi0 = dubiner_basis(0, 0, 0);
zero = @(x, y) 0*x;
minus_P = @(x, y) -P + 0*x;
curve = FlatCurveHorizontal(0.5);
[fluid_mesh, biot_mesh] = fitted_split(curve, [0.8*h, h], h^2);
side = @(m, tag, f) apply_tag(m, tag, f);
on_bottom = @(X) X(2, :) < 1e-12;
on_top    = @(X) X(2, :) > 1 - 1e-12;
on_left   = @(X) X(1, :) < 1e-12;
on_right  = @(X) X(1, :) > 1 - 1e-12;

% name, kept, biot specs, fluid specs, s0, alpha, pore_pressure_mean, level P
cases = {
  "roller",            true,  "roller",  "",         0,     1,         P, P
  "slip",              true,  "",        "slip",     0,     1,         P, P
  "symmetry",          true,  "",        "symmetry", 0,     1,         P, P
  "roller+slip+symm.", true,  "roller",  "slip+symmetry", 0, 1,      P, P
  "s0 = 1e-14",        true,  "",        "",         1e-14, 1,         P, P
  "alpha = 1 - 1e-13", true,  "",        "",         0,     1 - 1e-13, P, P
  "Biot traction x=1", false, "traction", "",        0,     1,         0, P
  "fluid outflow x=1", false, "",        "outflow",  0,     1,         0, P
  "s0 = 1e-6",         false, "",        "",         1e-6,  1,         P, 0
};

n_fail = 0;
for i = 1:size(cases, 1)
    [name, kept, biot_kind, fluid_kind, s0, alpha, level_mean, p_level] = cases{i, :};
    fm = fluid_mesh;  bm = biot_mesh;
    fluid_specs = boundary_condition("Dirichlet");
    if contains(fluid_kind, "slip")
        fm = side(fm, "Top", on_top);
        fluid_specs = [fluid_specs; boundary_condition("Top", x = "neumann")];
    end
    if contains(fluid_kind, "symmetry")
        fm = side(fm, "Left", on_left);
        fluid_specs = [fluid_specs; boundary_condition("Left", y = "neumann")];
    end
    if fluid_kind == "outflow"
        fm = side(fm, "Outflow", on_right);
        fluid_specs = [fluid_specs; boundary_condition("Outflow", x = "neumann", ...
            y = "neumann", x_data = minus_P)];
    end
    biot_specs = boundary_condition("Dirichlet");
    if biot_kind == "roller"
        bm = side(bm, "Bottom", on_bottom);
        biot_specs = [biot_specs; boundary_condition("Bottom", x = "neumann")];
    elseif biot_kind == "traction"
        bm = side(bm, "Right", on_right);
        biot_specs = [biot_specs; boundary_condition("Right", x = "neumann", x_data = minus_P)];
    end
    fluid = IncompressibleElasticityProblem(k, viscosity = 1, density = 0, time_step = 1, ...
        source_term = {zero; zero}, datum_dirichlet = {zero; zero}, pressure_mean = "pinned");
    fluid.Data.boundary_conditions = fluid_specs;
    biot = BiotProblem(k, biot_alpha = alpha, storage = s0, time_step = 1);
    biot.Data.boundary_conditions = biot_specs;

    [fs, bs, ~, gs] = StokesBiotCoupling.solve(fluid, biot, fm, bm, curve, struct.empty, ...
        Order = k, friction_bjs = 1.0, pore_pressure_mean = level_mean);
    is_kept = ismember(numel(gs.solution), gs.dof);
    err = state_errors(fs, bs, p_level, d2, phi0);
    c = condest(gs.matrix(gs.dof, gs.dof));
    fprintf("%-20s scalar kept %d (expected %d)  worst dof error %.1e  condest %.2e  level scalar %.1e\n", ...
        name, is_kept, kept, max(err), c, gs.level_scalar);
    n_fail = n_fail + check(is_kept == kept, sprintf("%s: level scalar kept = %d", name, is_kept));
    n_fail = n_fail + check(all(err < 1e-8), sprintf("%s: hydrostatic state (worst %.1e)", name, max(err)));
    if kept
        n_fail = n_fail + check(c < 1e12, sprintf("%s: condest %.1e", name, c));
        n_fail = n_fail + check(abs(gs.level_scalar) < 1e-8, ...
            sprintf("%s: level scalar %.1e", name, gs.level_scalar));
    else
        n_fail = n_fail + check(gs.level_scalar == 0, sprintf("%s: eliminated scalar not 0", name));
    end
end

fprintf("%d check(s) failed.\n", n_fail);
if n_fail > 0
    error("test_stokes_biot_pressure_level:failed", "%d checks failed", n_fail);
end

function err = state_errors(fs, bs, P, d2, phi0)
% max dof errors against u_F = 0, sigma_F = sigma_P = -P I, eta = 0, q = 0,
% p = P
iso = @(F) F.dof - full(sparse([1, 3*d2 + 1], 1, -P/phi0, size(F.dof, 1), 1));
f = fs.Unknowns.field;  b = bs.Unknowns.field;
p_ex = zeros(d2, 1);  p_ex(1) = P/phi0;
err = [max(abs(iso(f(1))), [], 'all'), max(abs(f(2).dof), [], 'all'), ...
       max(abs(iso(b(1))), [], 'all'), max(abs(b(2).dof), [], 'all'), ...
       max(abs(b(4).dof), [], 'all'), max(abs(b(5).dof - p_ex), [], 'all')];
end
