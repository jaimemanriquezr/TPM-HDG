% test_fsi_lagrangian_circle_march.m
% The closed-curve path of lagrangian_carved_march: Sigma_0 the circle of
% radius 0.4 about (0.5, 0.5) (lagrangian_march_problems, geometry
% "circle"), Sigma_h a periodic DeformedCurve, the fluid carved from a
% "crossed" background grid (N = 24, the coarsest level of
% lagrangian_march_convergence) and snapped (snap_to_curve on a periodic
% curve), the disc meshed by mesh2d with its boundary vertices on the
% circle (as cavity_smoke meshes its disc; mesh2d must be on the path),
% k = 1, A = 0.2, 4 steps of dt = 0.1, linear and hyperelastic (Saint
% Venant-Kirchhoff) solid:
%   1. the march completes, the solution stays bounded and the velocity
%      and displacement errors are at the level of this coarse mesh
%      (bounds per solid model);
%   2. Sigma_h follows the motion: the geometry error stays below half the
%      exact interface displacement;
%   3. the transfer paths stay O(h^2): max path below 1.5 h^2;
%   4. only the interface vertices of the carve move, onto the final
%      Sigma_h or just outside it (offset below h^2), and all of them are
%      snapped;
%   5. hyperelastic: the coupled Newton loop converges at every step.
clc; clear;
warning("off", "MATLAB:nearlySingularMatrix");   % trace extension on sliver patches

assert(exist("refine2", "file") == 2, "test_fsi_lagrangian_circle_march:mesh2d", ...
    "mesh2d (refine2) is not on the path: add src/mesh-generation/utilities/mesh2d " + ...
    "(populated in the main checkout) before running this test.");

n_fail = 0;
k = 1;  N = 24;  A = 0.2;
center = [0.5; 0.5];  R = 0.4;
background = StructuredSimplicialGrid(N, N, "crossed");
solid_mesh = disc_mesh(center, R, 0.05);
fprintf("circle march: N = %d, k = %d, A = %g, disc %d elements\n", N, k, A, size(solid_mesh.elements, 2));

% bounds: [|u_F - u|, |u_S - u|]
runs = struct("name", {"linear", "hyperelastic"}, "dt", {0.1, 0.1}, "n_steps", {4, 4}, ...
    "options", {{}, {"solid_model", "hyperelastic", "material", "svk"}}, ...
    "bounds", {[9e-5, 1.1e-3], [1e-4, 1.1e-3]});
for run = runs
    timer = tic;
    out = lagrangian_carved_march(k, "circle", background, solid_mesh, A, run.dt, run.n_steps, ...
        "snap", true, run.options{:});
    T = run.dt*run.n_steps;
    tag = run.name + ": ";

    %% 1. completes, bounded, errors at the mesh level
    max_u = max(out.history(:, 4:5), [], "all");
    fprintf("%sT = %g (%.1f s): |sF - s| %.3e, |uF - u| %.3e, |sS - s| %.3e, |uS - u| %.3e, max |u_h| %.3g\n", ...
        tag, T, toc(timer), out.eF(1:2), out.eS(1:2), max_u);
    fprintf("%sswitched on/off per step: %s\n", tag, mat2str(out.history(:, 6:7).'));
    n_fail = n_fail + check(abs(out.history(end, 1) - T) < 1e-12, tag + "the march did not reach T");
    n_fail = n_fail + check(all(isfinite([out.eF, out.eS])) && max_u < 10, tag + "the march is unbounded");
    n_fail = n_fail + check(out.eF(2) < run.bounds(1) && out.eS(2) < run.bounds(2), ...
        tag + "errors above the level of N = 24");

    %% 2. Sigma_h follows the exact interface
    X = center + R*[cos(2*pi*(0:399)/400); sin(2*pi*(0:399)/400)];
    motion = max(vecnorm(A*T^4*[cos(pi*X(2, :)) + X(1, :); sin(pi*X(1, :)) + X(2, :)]));
    geometry = max(out.history(:, 2));
    fprintf("%sgeometry error %.3e (predictor %.3e), lag %.3e, exact interface displacement %.3e\n", ...
        tag, geometry, max(out.history(:, 3)), max(out.history(:, 8)), motion);
    n_fail = n_fail + check(geometry < 0.5*motion, tag + "Sigma_h does not follow the interface");

    %% 3. transfer paths O(h^2)
    fprintf("%stransfer paths [%.2e, %.2e], max/h^2 = %.2f\n", tag, out.paths, out.paths(2)*N^2);
    n_fail = n_fail + check(out.paths(2)*N^2 < 1.5, tag + "transfer paths above 1.5 h^2");

    %% 4. only the interface vertices move, onto Sigma_h
    carved = background.split(out.curve);
    snapped = out.fluid_mesh;
    same_mesh = isequal(size(carved.coordinates), size(snapped.coordinates)) ...
        && isequal(carved.elements, snapped.elements);
    n_fail = n_fail + check(same_mesh, tag + "the final fluid mesh is not the snapped carve of out.curve");
    if same_mesh
        v = unique(carved.faces(1:2, carved.intfaces));
        moved = find(any(snapped.coordinates ~= carved.coordinates, 1));
        P = snapped.coordinates(:, moved);
        offset = [vecnorm(P - out.curve.map_points_to_interface(P, coordinates = "physical")), 0];
        fprintf("%s%d of %d interface vertices snapped, %d other vertices moved; offset from Sigma_h up to %.2e (h^2 = %.2e)\n", ...
            tag, numel(intersect(moved, v)), numel(v), numel(setdiff(moved, v)), max(offset), 1/N^2);
        n_fail = n_fail + check(isempty(setdiff(moved, v)), tag + "a vertex off the interface moved");
        n_fail = n_fail + check(max(offset) < 1/N^2, tag + "a snapped vertex is not on Sigma_h");
        n_fail = n_fail + check(numel(intersect(moved, v)) == numel(v), ...
            tag + "an interface vertex is not snapped");
    end

    %% 5. Newton
    if isfield(out, "newton")
        fprintf("%sNewton: iterations %s, stop %s\n", tag, mat2str([out.newton.iterations]), ...
            strjoin([out.newton.stop_reason], ", "));
        n_fail = n_fail + check(all([out.newton.converged]), tag + "the coupled Newton loop did not converge");
    end
end

fprintf("%d check(s) failed.\n", n_fail);
if n_fail > 0
    error("test_fsi_lagrangian_circle_march:failed", "%d check(s) failed", n_fail);
end

function mesh = disc_mesh(center, R, h)
% the disc meshed by mesh2d, its boundary polygon inscribed in the circle
% (vertices on it, no boundary refinement), boundary faces "Interface"
% (cavity_smoke's disc)
n = ceil(2*pi*R/h);
theta = 2*pi*(0:n - 1)/n;
mesh = SimplicialGrid(center + R*[cos(theta); sin(theta)], [], element_diameter = h, ...
    refine_boundaries = false);
mesh = apply_tag(mesh, "Interface", @(P) abs(vecnorm(P - center) - R) < 1e-10);
end
