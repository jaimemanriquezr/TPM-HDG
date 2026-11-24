function out = eulerian_carved_march(k, background, solid_mesh, amplitude, dt, n_steps, ...
    uu_fluid, p_fluid, u_solid, options)
%EULERIAN_CARVED_MARCH  FSI march with an Eulerian fluid re-carved every step.
%
%   out = eulerian_carved_march(k, background, solid_mesh, A, dt, n_steps, ...
%       uu_fluid, p_fluid, u_solid, ...)
%
%   Sigma(t) = Sigma_0 + d(t) e_x, d(t) = A sin(2 pi t), Sigma_0 the circle
%   of radius 0.4 at (0.5, 0.5); data from translating_circle_problems.
%   At t_n the fluid mesh is the part of the fixed background grid outside
%   Sigma(t_n) (SimplicialGrid.split, cut elements discarded, so the fluid
%   mesh lies inside Omega_F(t_n) and its transfer paths are at most about
%   one element long); the solid mesh translates rigidly, X + d(t_n) e_x
%   (as in ale_translating_march). The fluid is Eulerian: its dofs live on
%   the background elements, so no convective mesh-velocity term enters.
%   Its history (velocity, acceleration) is copied on the elements active
%   at both steps (identified by their background element; re-expressed
%   from their previous shape if snapping moved them); an element that
%   switches on gets the L2 projection of
%   the polynomial of the nearest element active at the previous step,
%   evaluated on it (the same extrapolation the transfer paths use).
%   The fluid is backward Euler on the first step and Crank-Nicolson
%   after; the solid is Newmark.
%
%   out: errors at t = n_steps*dt (eF, eS: compute_l2_error of
%   (sigma, u, gamma)), the transfer path lengths [min, max] over the run,
%   the final solutions, and per step (history) the largest |u_h| dof of
%   the fluid and of the solid and the number of fluid elements switched
%   on and off.
%
%See also: ale_translating_march, translating_circle_problems, SimplicialGrid.split
arguments
    k, background SimplicialGrid, solid_mesh, amplitude, dt, n_steps, uu_fluid, p_fluid, u_solid
    options.rho_f = 1; options.mu = 1; options.rho_s = 1; options.E = 10; options.nu = 0.3;
    options.patch_weight (1,1) string = "jacobian"
    options.n_quadpoints double {mustBeScalarOrEmpty} = []
    options.verbose (1,1) logical = false
    % "flat": Sigma(t) = {y = 0.5 + d(t)}, fluid above, the solid mesh
    % (below) translating by d(t) e_y; data from translating_flat_problems
    options.geometry (1,1) string {mustBeMember(options.geometry, ["circle", "flat"])} = "circle"
    % circle only: after carving, move the fluid interface vertices just
    % outside Sigma(t_n) (snap_outside); the reshaped first-layer elements
    % then get their history re-expressed from their previous shape
    options.snap (1,1) logical = false
    % every snapshot_every steps, append the solution to snapshot_file
    % (save_march_snapshots) for plots/movies
    options.snapshot_file (1,1) string = ""
    options.snapshot_every (1,1) double {mustBeInteger, mustBePositive} = 1
end
if options.geometry == "flat"
    data = @(t) translating_flat_problems(k, 0.5, amplitude, t, dt, uu_fluid, p_fluid, ...
        u_solid, options.rho_f, options.mu, options.rho_s, options.E, options.nu);
    direction = [0; 1];
else
    curve_0 = CircularCurve([0.5; 0.5], 0.4);
    data = @(t) translating_circle_problems(k, curve_0, amplitude, t, dt, uu_fluid, p_fluid, ...
        u_solid, options.rho_f, options.mu, options.rho_s, options.E, options.nu);
    direction = [1; 0];
end
d = @(t) amplitude*sin(2*pi*t);
% [fluid mesh, its background element keys]: the carved fluid (cut
% elements discarded; snapped when asked) and the centroids of the
% unsnapped carve, which identify each element's background element
if options.snap
    assert(options.geometry == "circle", "eulerian_carved_march:snap", "snap needs the circle.");
end
carve = @(curve) carve_fluid(background, curve, options.snap);

[fp, sp, src, ~, curve] = data(0);
sp.Data.time_scheme = "newmark";
fp.Data.time_scheme = "euler";
quad = {"n_quadpoints", options.n_quadpoints, "n_quadpoints_path", options.n_quadpoints + 1};
if isempty(options.n_quadpoints), quad = {"n_quadpoints", k + 2, "n_quadpoints_path", k + 3}; end
[fm, keys] = carve(curve);
state = FluidStructureCoupling.initial_state(fp, sp, fm, solid_mesh, curve, src{:}, "Order", k, quad{:});
paths = [Inf, 0];
history = zeros(n_steps, 4);
for n = 1:n_steps
    t = n*dt;
    [fp, sp, src, pd, curve] = data(t);
    [fm_old, keys_old] = deal(fm, keys);
    [fm, keys] = carve(curve);
    [state.fluid, switched] = remap_fluid_state(state.fluid, fm_old, fm, k, keys_old, keys);
    sm = solid_mesh + d(t)*direction;
    sp.Data.time_scheme = "newmark";
    fp.Data.time_scheme = "cranknicolson";
    if n == 1, fp.Data.time_scheme = "euler"; end
    [fs, ss, state, gs] = FluidStructureCoupling.solve(fp, sp, fm, sm, curve, src{:}, ...
        "newmark", "quiet", true, "state", state, "patch_weight", options.patch_weight, quad{:});
    for part = {gs.partition_1, gs.partition_2}
        l = vecnorm(part{1}.quadrature.interface - part{1}.quadrature.physical, 2, 1);
        paths = [min(paths(1), min(l(:))), max(paths(2), max(l(:)))];
    end
    history(n, :) = [max(abs(fs.Unknowns.field(2).dof), [], "all"), ...
        max(abs(ss.Unknowns.field(2).dof), [], "all"), switched];
    if options.snapshot_file ~= "" && (mod(n, options.snapshot_every) == 0 || n == n_steps)
        save_march_snapshots(options.snapshot_file, k, t, fs, ss, state, ...
            curve.get_2d_points(linspace(curve.t_range(1), curve.t_range(2), 401)));
    end
    if options.verbose
        eF = fs.compute_l2_error({pd.sigma_F, pd.u_F, pd.gamma_F});
        eS = ss.compute_l2_error({pd.sigma_S, pd.u_S, pd.gamma_S});
        fprintf("  step %d t %.4f  uF %.3e  uS %.3e  on/off %d/%d\n", n, t, eF(2), eS(2), switched);
    end
end
eF = fs.compute_l2_error({pd.sigma_F, pd.u_F, pd.gamma_F});
eS = ss.compute_l2_error({pd.sigma_S, pd.u_S, pd.gamma_S});
out = struct("eF", eF(:).', "eS", eS(:).', "paths", paths, "fs", fs, "ss", ss, ...
    "h", max([fm.area(:); solid_mesh.area(:)]), "history", history);
end

function [fluid, keys] = carve_fluid(background, curve, snap)
fluid = background.split(curve);
keys = squeeze(mean(reshape(fluid.coordinates(:, fluid.elements), 2, 3, []), 2)).';
if snap
    fluid = snap_outside(fluid, curve);
end
end

function [s, switched] = remap_fluid_state(s, old_mesh, new_mesh, k, old_keys, new_keys)
% the volume histories (velocity, acceleration; Dubiner dofs, component
% blocks stacked, one column per element) from old_mesh to new_mesh, both
% carved from the same background grid, elements identified by their
% background element (the keys): an element present at both steps keeps
% its polynomial, re-expressed on its new shape if snapping moved it
% (exact); an element that switches on gets the L2 projection of the
% nearest old element's polynomial. switched counts the background
% elements switched on and off.
centroid = @(m) squeeze(mean(reshape(m.coordinates(:, m.elements), 2, 3, []), 2)).';
[present, from] = ismembertol(new_keys, old_keys, 1e-10, "ByRows", true);
switched = [nnz(~present), size(old_mesh.elements, 2) - nnz(present)];
c_old = centroid(old_mesh);
c_new = centroid(new_mesh);
same_shape = present;
same_shape(present) = vecnorm(c_new(present, :) - c_old(from(present), :), 2, 2) < 1e-14;
kept = find(same_shape);
added = find(~same_shape);
source = zeros(size(added));
for j = 1:numel(added)
    if present(added(j))
        source(j) = from(added(j));
    else
        [~, source(j)] = min(vecnorm(c_old - c_new(added(j), :), 2, 2));
    end
end
d2 = nchoosek(k + 2, 2);
[rs, w] = get_2d_reference_quadrature(2*k + 2);
P = dubiner_basis(k, rs);
fit = (P.'*(w(:).*P)) \ (P.'.*w(:).');          % d2 x n_q: values -> dofs
for name = ["velocity", "acceleration"]
    if ~isfield(s, name) || ~isnumeric(s.(name)) || isempty(s.(name)) || isequal(s.(name), 0)
        continue
    end
    D = s.(name);
    E = zeros(size(D, 1), size(new_mesh.elements, 2));
    E(:, kept) = D(:, from(kept));
    for j = 1:numel(added)
        xy = new_mesh.get_2d_physical(rs, added(j));
        rs_old = old_mesh.get_2d_local(xy, source(j));
        Q = dubiner_basis(k, rs_old);
        for c = 0:size(D, 1)/d2 - 1
            rows = c*d2 + (1:d2);
            E(rows, added(j)) = fit*(Q*D(rows, source(j)));
        end
    end
    s.(name) = E;
end
end
