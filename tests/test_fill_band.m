% FILL_BAND: a carved mesh with its cut band meshed back in.
%
% 1. Geometry, on the stenosis wall at h = 0.05:
%    every interior edge is shared by exactly two elements, the boundary
%    lies only on Sigma and the three walls, Sigma carries the Interface
%    tag, the area equals the integral under Sigma, and the strip holds no
%    sliver (min angle > 15 deg with the default drop 0.3h).
% 2. Coupling, on the unit-square polynomial of test_fsi_curved_carve: the
%    fluid carved then filled (fitted to Sigma), the solid carved, the march
%    from FluidStructureCoupling.initial_state with the state returned by
%    each step (the coupling transports the T2 history). Both integrators
%    reproduce this solution exactly, so the error must stay at round-off
%    over the steps.
%
%See also: test_fsi_curved_carve, SimplicialGrid/fill_band
clc; clear;
n_fail = 0;

%% 1. Geometry -----------------------------------------------------------
wall = SmoothCurveHorizontal("0.5 - 0.05*exp(-((t - 3)/0.5)^2)", [0 6]);
h = 0.05;
bg = StructuredSimplicialGrid(round(6/h), 13, 0, 6, 0, 0.6, style = "left");
[~, carved] = bg.split(wall, "cut_elements", "discard");
mesh = carved.fill_band(wall, h = h);

X = mesh.coordinates;  T = mesh.elements;
edges = sort([T([1 2], :), T([2 3], :), T([3 1], :)], 1).';
[unique_edges, ~, ic] = unique(edges, "rows");
count = accumarray(ic, 1);
boundary = unique_edges(count == 1, :);
bx = reshape(X(1, boundary.'), 2, []);  by = reshape(X(2, boundary.'), 2, []);
on_sigma = all(abs(by - wall.y_coord(bx)) < 1e-9, 1);
on_wall = all(abs(bx) < 1e-9, 1) | all(abs(bx - 6) < 1e-9, 1) | all(abs(by) < 1e-9, 1);
xq = linspace(0, 6, 60001);
area_error = abs(sum(abs(mesh.volume)) - trapz(xq, wall.y_coord(xq)));
P1 = X(:, T(1, :));  P2 = X(:, T(2, :));  P3 = X(:, T(3, :));
angle = @(A, B, C) acosd(sum((B - A).*(C - A), 1)./(vecnorm(B - A).*vecnorm(C - A)));
min_angle = min(min([angle(P1, P2, P3); angle(P2, P3, P1); angle(P3, P1, P2)], [], 1));
n_interface = numel(get_face_by_tag(mesh, "Interface"));

checks = { ...
    "no edge shared by more than two elements", all(count <= 2)
    "boundary only on Sigma and the walls",     all(on_sigma | on_wall)
    "Interface tag covers Sigma",               n_interface == nnz(on_sigma) && n_interface > 0
    "area equals the integral under Sigma",     area_error < 1e-12
    "min angle above 15 deg",                   min_angle > 15};
fprintf("geometry: %d elements, %d Interface faces, area error %.1e, min angle %.1f deg\n", ...
    size(T, 2), n_interface, area_error, min_angle);
for i = 1:size(checks, 1)
    if ~checks{i, 2}
        fprintf("  FAILED: %s\n", checks{i, 1});
        n_fail = n_fail + 1;
    end
end

%% 2. Coupling -------------------------------------------------------------
k = 2; E = 1.0; nu = 0.25; visc = 0.4; rho_f = 1.0; rho_s = 1.0;
N = 9; dt = 1e-3; n_steps = 12;
tol_exact = 1e-8;

uu_fluid = "x*y*(1 + t)";
p_fluid  = "(x + y)*(1 + t)";
u_solid  = ["x*(1 + t + t^2/2)"; "y*(1 + t + t^2/2)"];

interface = SmoothCurveHorizontal();
[fluid_mesh, solid_mesh] = StructuredSimplicialGrid(N, N, "left") ...
    .split(interface, cut_elements = "discard");
fluid_mesh = fluid_mesh.fill_band(interface, h = 1/N);

[pde_f, pde_s, src, d0] = FluidStructureCoupling.load_problems(k, interface, 0.0, dt, ...
    uu_fluid, p_fluid, u_solid, rho_f, visc, rho_s, E, nu, 0);
pde_s.Data.time_scheme = "newmark";
pde_f.Data.time_scheme = "cranknicolson";
state = FluidStructureCoupling.initial_state(pde_f, pde_s, fluid_mesh, solid_mesh, ...
    interface, src{:}, source_dirichlet_rate = d0.T_2_t);
worst = NaN; blew = 0;
for n = 1:n_steps
    t = n*dt;
    [pde_f, pde_s, src, pd] = FluidStructureCoupling.load_problems(k, interface, t, dt, ...
        uu_fluid, p_fluid, u_solid, rho_f, visc, rho_s, E, nu, 0);
    pde_s.Data.time_scheme = "newmark";
    pde_f.Data.time_scheme = "cranknicolson";
    [sol_f, sol_s, state] = FluidStructureCoupling.solve(pde_f, pde_s, fluid_mesh, solid_mesh, ...
        interface, src{:}, "newmark", quiet = true, state = state);

    es = sol_s.compute_l2_error({pd.sigma_S, pd.u_S, pd.gamma_S}, {pd.u_S});
    ef = sol_f.compute_l2_error({pd.sigma_F, pd.u_F, pd.gamma_F}, {pd.u_F});
    worst = max([es(:); ef(:)]);
    if ~isfinite(worst) || worst > 1e6, blew = n; break; end
end
if blew
    fprintf("coupling: DIVERGED at step %d\n", blew);
    n_fail = n_fail + 1;
else
    fprintf("coupling: worst L2 error after %d steps %.3e\n", n_steps, worst);
    if worst > tol_exact
        fprintf("  FAILED: filled fluid + carved solid is not exact\n");
        n_fail = n_fail + 1;
    end
end

fprintf("%d check(s) failed.\n", n_fail);
if n_fail > 0
    error("test_fill_band:failed", "%d check(s) failed", n_fail);
end
