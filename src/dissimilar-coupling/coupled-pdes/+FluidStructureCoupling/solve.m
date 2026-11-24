function [fluid_solution, solid_solution, state, global_system, cache] = solve(...
    fluid_problem, solid_problem, fluid_mesh, solid_mesh, ...
    interface, source_neumann, source_dirichlet, scheme, ...
    options)
%SOLVE  TPM-coupled FSI solve: the declarative InterfaceProblem coupling.
%
%   [fluid_solution, solid_solution, state] = FluidStructureCoupling.solve( ...
%       fluid_problem, solid_problem, fluid_mesh, solid_mesh, interface, ...
%       T1, T2, scheme, state = state, ...)
%
%   Takes the two low-level problem objects load_problems builds, the two
%   meshes and the curve Sigma (or the per-side pair {curve_fluid,
%   curve_solid} of a Lagrangian solid, see InterfaceProblem.fsi), builds
%   the declarative coupling (interface_problem) and solves it.
%   scheme = "steady" (default) |
%   "onestep" | "newmark"; the three differ only in the T2 kinematic row.
%   The interface law is "noslip" by default; "navier" (with
%   slip_coefficient a) and "full" (a = 0) replace the tangential half of
%   the kinematic row by the fluid Navier law (see InterfaceProblem.fsi).
%   swap_transmission = true tests the traction row T1 with the solid's
%   interface traces and the kinematic row T2 with the fluid's ("noslip"
%   only; see InterfaceProblem.fsi).
%
%   Time stepping: state is the time-stepping state of both sides
%   (state.time, state.step, state.solid, state.fluid, state.diagnostics;
%   see LinearElasticityProblem and IncompressibleElasticityProblem). Build
%   it at t = 0 with FluidStructureCoupling.initial_state and pass each
%   step's third output to the next step:
%     state = FluidStructureCoupling.initial_state(fp, sp, fm, sm, curve, T1, T2);
%     for n = 1:n_steps
%         % set the sources of fp, sp at t_n
%         [fs, ss, state] = FluidStructureCoupling.solve(fp, sp, fm, sm, curve, ...
%             T1, T2, "newmark", state = state);
%     end
%   The coupling forms the T2 history and its transport onto Sigma from
%   state.solid (see interface_problem). Under "onestep" / "newmark" the
%   returned state is advanced by each side's get_kinematic_state and
%   merged (state.solid from the solid, state.fluid from the fluid,
%   state.time and state.step from the solid, state.diagnostics carried).
%   Under "steady" the input state is returned unchanged. An
%   empty state is the zero history on every scheme (no initial data are
%   read); the steady problem with time = Inf and an empty state is the
%   stationary-with-reaction convention.
%
%   Cache: start a loop in which only the fluid data change (e.g. Picard)
%   from struct("frozen", 2) and pass each call's fifth output back; see
%   interface_problem and InterfaceProblem.solve. Without a cache the
%   fifth output is struct.empty.
%
%See also: FluidStructureCoupling.initial_state,
%          FluidStructureCoupling.interface_problem, InterfaceProblem.fsi,
%          FluidStructureCoupling.load_problems
arguments
    fluid_problem   IncompressibleElasticityProblem
    solid_problem   LinearElasticityProblem
    fluid_mesh      SimplicialGrid
    solid_mesh      SimplicialGrid
    interface       {mustBeInterfaceCurves}   % one curve, or {curve_fluid, curve_solid} (InterfaceProblem.fsi)
    source_neumann  cell = {@(x,y) 0*x.*y, @(x,y) 0*x.*y}
    source_dirichlet cell = {@(x,y) 0*x.*y, @(x,y) 0*x.*y}
    scheme          (1,1) string {mustBeMember(scheme, ["steady", "onestep", "newmark"])} = "steady"
    options.Order double {mustBeScalarOrEmpty} = fluid_problem.Order
    options.n_quadpoints      double {mustBeScalarOrEmpty} = []
    options.n_quadpoints_path double {mustBeScalarOrEmpty} = []
    options.quiet             (1,1) logical = false   % accepted for the callers; the solve is silent
    options.patch_weight      (1,1) string {mustBePatchWeight} = "rectangular"
    options.sliver_correction (1,1) logical = true
    options.vertex_tol        (1,1) double = 1e-8
    options.vertex_tol_relative double = []
    options.mapping_1 (1,1) string = ""   % fluid / solid point-to-interface maps
    options.mapping_2 (1,1) string = ""   % ("" keeps the curve's own)
    options.state struct = struct.empty
    options.interface_law (1,1) string ...
        {mustBeMember(options.interface_law, ["noslip", "navier", "full"])} = "noslip"
    options.slip_coefficient (1,1) double {mustBeNonnegative} = 0.0
    options.source_slip cell = {@(x,y) 0*x.*y, @(x,y) 0*x.*y}
    % see InterfaceProblem.fsi and get_partitions
    % (tip_mapping = "collapse" is a fallback where "spread" fails with
    % orphanFace; where "spread" works it is more accurate, and collapse
    % is not meant for convergence studies)
    options.swap_transmission (1,1) logical = false
    options.tip_mapping (1,1) string ...
        {mustBeMember(options.tip_mapping, ["spread", "collapse"])} = "spread"
    options.collapsed_normal (1,1) string ...
        {mustBeMember(options.collapsed_normal, ["path", "face", "bisector"])} = "path"
    options.Cache struct = struct.empty
end
[coupled, ~, ~, cache] = FluidStructureCoupling.interface_problem(fluid_problem, solid_problem, ...
    fluid_mesh, solid_mesh, interface, source_neumann, source_dirichlet, scheme, ...
    Order = options.Order, patch_weight = options.patch_weight, ...
    sliver_correction = options.sliver_correction, ...
    vertex_tol = options.vertex_tol, vertex_tol_relative = options.vertex_tol_relative, ...
    mapping_1 = options.mapping_1, mapping_2 = options.mapping_2, ...
    n_quadpoints = options.n_quadpoints, n_quadpoints_path = options.n_quadpoints_path, ...
    state = options.state, ...
    interface_law = options.interface_law, slip_coefficient = options.slip_coefficient, ...
    source_slip = options.source_slip, swap_transmission = options.swap_transmission, ...
    tip_mapping = options.tip_mapping, collapsed_normal = options.collapsed_normal, ...
    Cache = options.Cache);
[fluid_solution, solid_solution, global_system, cache] = ...
    solve(coupled, options.Order, Cache = cache);
state = options.state;
if scheme ~= "steady"
    if isempty(state)
        % the zero history at the solid's time
        state = struct("time", solid_problem.Data.time, "step", 0, "solid", struct(), ...
            "fluid", struct.empty, "diagnostics", struct.empty);
    end
    solid_state = get_kinematic_state(solid_problem, solid_solution, state);
    fluid_state = get_kinematic_state(fluid_problem, fluid_solution, state);
    state.solid = solid_state.solid;
    state.fluid = fluid_state.fluid;
    state.time  = solid_state.time;
    state.step  = solid_state.step;
end
end
