% test_interface_fsi_cache.m
% InterfaceProblem.solve with a Cache (side 2, the solid, frozen) must
% reproduce the uncached solve exactly (isequal of the dofs) over a
% three-step Picard loop (Stokes, then two convected steps) on the flat
% fitted FSI problem, with transferred essential conditions on both sides
% (fluid "Wall" faces on y = 1, solid "Clamp" faces on y = 0), for
% swap_transmission false and true; the cached loop also passes the cache
% to FluidStructureCoupling.interface_problem (fluid geometric blocks,
% frozen solid problem). The fluid factory with a geometry cache must
% reproduce its uncached blocks. A cache built on another mesh must be
% rejected (InterfaceProblem:cacheKey, IncompressibleElasticityMatrices:
% cacheKey, FluidStructureCoupling:cacheKey), and so must a cache built
% for another cdt, viscosity or patch weight (InterfaceProblem:cacheKey) or
% one that does not declare frozen = 2 (InterfaceProblem:cacheFrozen).
clc; clear;
n_fail = 0;
k = 2;
mu_f = 1.0; E_s = 1.0; nu_s = 0.25; dt = 1e-2;
interface = FlatCurveHorizontal(0.5);
h = 1/8;
n_steps = 3;

[fp, sp, srcs, pd] = FluidStructureCoupling.load_problems(k, interface, 0.0, dt, ...
    "x^2*y^2", "x^2 + y^2", ["x^2*y"; "x*y^2"], ...
    0.0, mu_f, 0.0, E_s, nu_s, 0.0, pressure_mean = "pinned");
fp.Data.density = 1.0;
wall_bc = TransferEssentialBoundaryConditions.incompressible_stress(Datum(pd.u_F), ...
    "Wall", @(X) [X(1, :); 1.05 + 0*X(2, :)], mu_f);
clamp_bc = TransferEssentialBoundaryConditions.elasticity(Datum(pd.u_S), ...
    "Clamp", @(X) [X(1, :); -0.05 + 0*X(2, :)], E_s, nu_s);

% fitted_split: mesh 1 (upper) is the fluid, mesh 2 (lower) the solid
[fluid_mesh, solid_mesh] = tagged_meshes(interface, h);
n_wall = numel(get_face_by_tag(fluid_mesh, "Wall"));
n_clamp = numel(get_face_by_tag(solid_mesh, "Clamp"));
fprintf("transfer faces: Wall %d (fluid), Clamp %d (solid)\n", n_wall, n_clamp);
if n_wall == 0 || n_clamp == 0
    fprintf(2, "FAILED: no faces carry a transfer condition\n"); n_fail = n_fail + 1;
end

%% cached vs uncached Picard loops
for swap = [false, true]
    cache = struct("frozen", 2);
    u = 0;
    for step = 1:n_steps
        state = struct("fluid", struct("velocity", u));   % the pseudo-time history
        fp.Data.convective_velocity = {};
        if ~isequal(u, 0)
            d2 = size(u, 1)/2;
            fp.Data.convective_velocity = {u(1:d2, :); u(d2+1:end, :)};
        end
        coupled = coupled_problem(fp, sp, fluid_mesh, solid_mesh, interface, srcs, k, swap, wall_bc, clamp_bc, struct.empty, state);
        [fs_ref, ss_ref] = solve(coupled, k);
        [coupled, cache] = coupled_problem(fp, sp, fluid_mesh, solid_mesh, interface, srcs, k, swap, wall_bc, clamp_bc, cache, state);
        [fs, ss, ~, cache] = solve(coupled, k, Cache = cache);
        ref = {fs_ref.Unknowns.field(1).dof, fs_ref.Unknowns.field(2).dof, ...
               fs_ref.Unknowns.hybrid(1).dof, ss_ref.Unknowns.field(2).dof};
        got = {fs.Unknowns.field(1).dof, fs.Unknowns.field(2).dof, ...
               fs.Unknowns.hybrid(1).dof, ss.Unknowns.field(2).dof};
        d = max(cellfun(@(a, b) max(abs(a - b), [], "all"), ref, got));
        fprintf("swap %d step %d: max|cached - uncached| = %.3e\n", swap, step, d);
        if ~isequal(ref, got)
            fprintf(2, "FAILED: swap %d step %d: cached solve differs\n", swap, step);
            n_fail = n_fail + 1;
        end
        u = fs_ref.Unknowns.field(2).dof;
    end
end

%% fluid factory with and without a geometry cache (Oseen, transient)
d2 = size(u, 1)/2;
fluid_opts = {"pressure_mean", "pinned", "density", 1.0, "time_step", dt, ...
    "state", struct("velocity", u), "convective_velocity", {u(1:d2, :); u(d2+1:end, :)}};
[lhs, flux, rhs, geometry] = IncompressibleElasticityMatrices(k, fluid_mesh, fp.Data.viscosity, ...
    fp.Data.source_term, fluid_opts{:}, Cache = struct.empty);
[lhs_c, flux_c, rhs_c] = IncompressibleElasticityMatrices(k, fluid_mesh, fp.Data.viscosity, ...
    fp.Data.source_term, fluid_opts{:}, Cache = geometry);
ref = [lhs.Arrays; flux.Arrays; rhs.Arrays];
got = [lhs_c.Arrays; flux_c.Arrays; rhs_c.Arrays];
d = max(cellfun(@(a, b) max(abs(double(a) - double(b)), [], "all"), ref, got));
fprintf("factory: %d blocks, max|cached - uncached| = %.3e\n", numel(ref), d);
for i = 1:numel(ref)
    if ~isequal(ref{i}.id, got{i}.id) || ~isequal(size(ref{i}), size(got{i})) ...
            || max(abs(double(ref{i}) - double(got{i})), [], "all") > 1e-13*max(abs(double(ref{i})), [], "all")
        fprintf(2, "FAILED: factory block %d [%d %d] differs with a cache\n", i, ref{i}.id);
        n_fail = n_fail + 1;
    end
end

%% a cache built on another mesh is rejected
[fluid_fine, solid_fine] = tagged_meshes(interface, h/2);
fp.Data.convective_velocity = {};
coupled_fine = coupled_problem(fp, sp, fluid_fine, solid_fine, interface, srcs, k, true, wall_bc, clamp_bc, struct.empty, struct.empty);
n_fail = n_fail + expect_error(@() solve(coupled_fine, k, Cache = cache), "InterfaceProblem:cacheKey");
n_fail = n_fail + expect_error(@() IncompressibleElasticityMatrices(k, fluid_fine, fp.Data.viscosity, ...
    fp.Data.source_term, fluid_opts{:}, Cache = geometry), "IncompressibleElasticityMatrices:cacheKey");
[~, cache_fine] = coupled_problem(fp, sp, fluid_fine, solid_fine, interface, srcs, k, true, wall_bc, clamp_bc, struct("frozen", 2), struct.empty);
n_fail = n_fail + expect_error(@() coupled_problem(fp, sp, fluid_mesh, solid_mesh, interface, srcs, k, true, ...
    wall_bc, clamp_bc, rmfield(cache_fine, "geometry"), struct.empty), "FluidStructureCoupling:cacheKey");

%% a cache built for another cdt is rejected
coupled_dt = coupled;  coupled_dt.coupling.cdt = 2*coupled.coupling.cdt;
n_fail = n_fail + expect_error(@() solve(coupled_dt, k, Cache = cache), "InterfaceProblem:cacheKey");

%% a cache built for another viscosity or patch weight is rejected
coupled_mu = coupled;  coupled_mu.coupling.viscosity = 2*coupled.coupling.viscosity;
n_fail = n_fail + expect_error(@() solve(coupled_mu, k, Cache = cache), "InterfaceProblem:cacheKey");
coupled_pw = coupled;  coupled_pw.coupling.patch_weight = "jacobian";
if coupled.coupling.patch_weight == "jacobian", coupled_pw.coupling.patch_weight = "rectangular"; end
n_fail = n_fail + expect_error(@() solve(coupled_pw, k, Cache = cache), "InterfaceProblem:cacheKey");

%% the cache must declare frozen = 2
n_fail = n_fail + expect_error(@() solve(coupled, k, Cache = struct()), "InterfaceProblem:cacheFrozen");
n_fail = n_fail + expect_error(@() solve(coupled, k, Cache = struct("frozen", 1)), "InterfaceProblem:cacheFrozen");

fprintf("%d check(s) failed.\n", n_fail);
if n_fail > 0, error("test_interface_fsi_cache:failed", "%d checks failed", n_fail); end

%% ----------------------------------------------------------------------
function [fluid_mesh, solid_mesh] = tagged_meshes(interface, h)
[fluid_mesh, solid_mesh] = fitted_split(interface, [h, h], h^2);
fluid_mesh = apply_tag(fluid_mesh, "Wall", @(X) isapprox(X(2, :), 1));
solid_mesh = apply_tag(solid_mesh, "Clamp", @(X) isapprox(X(2, :), 0));
end

function [coupled, cache] = coupled_problem(fp, sp, fluid_mesh, solid_mesh, interface, srcs, k, swap, wall_bc, clamp_bc, cache, state)
[coupled, ~, ~, cache] = FluidStructureCoupling.interface_problem(fp, sp, fluid_mesh, solid_mesh, interface, ...
    srcs{:}, "steady", Order = k, swap_transmission = swap, Cache = cache, state = state);
coupled.problems(1).BoundaryConditions = [coupled.problems(1).BoundaryConditions, wall_bc];
coupled.problems(2).BoundaryConditions = [coupled.problems(2).BoundaryConditions, clamp_bc];
end

function failed = expect_error(f, id)
failed = 0;
try
    f();
    fprintf(2, "FAILED: no error, expected %s\n", id);
    failed = 1;
catch err
    fprintf("%s raised: %s\n", id, err.identifier);
    if err.identifier ~= id
        fprintf(2, "FAILED: expected %s, got %s (%s)\n", id, err.identifier, err.message);
        failed = 1;
    end
end
end
