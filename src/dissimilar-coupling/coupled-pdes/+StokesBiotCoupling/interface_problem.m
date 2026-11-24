function [coupled, fluid_pb, biot_pb, cache] = interface_problem(fluid_problem, biot_problem, ...
    fluid_mesh, biot_mesh, interface, sources, options)
%INTERFACE_PROBLEM  The TPM-coupled Stokes-Biot problem of one time step.
%
%   [coupled, fluid_pb, biot_pb] = StokesBiotCoupling.interface_problem( ...
%       fluid_problem, biot_problem, fluid_mesh, biot_mesh, interface, ...
%       sources, state = state, friction_bjs = g, ...)
%
%   Builds the declarative fluid FEProblem as FluidStructureCoupling
%   .interface_problem does (IncompressibleElasticityMatrices, "pinned"
%   closure; a convective_velocity in the fluid Data passes through), the
%   Biot FEProblem by BiotProblem.get_declarative_problem with Data.state
%   = state, and returns InterfaceProblem.stokes_biot; solve it with
%   [fs, bs] = solve(coupled, k). The time scheme and step are read from
%   biot_problem.Data only: cdt and the history combinations come from
%   biot_step_history(time_scheme, time_step, state.biot), and the history
%   is transported onto Sigma by InterfaceProblem.solve.
%
%   Fluid inertia: with fluid Data.density = rho_F > 0 the fluid step is
%   the same scheme (Euler or BDF2) as the Biot step: the fluid's
%   time_step must equal the Biot one (its own time_scheme / theta is not
%   read), (cdt, u*) come from fluid_step_history(time_scheme, time_step,
%   state.fluid), and the fluid matrices get rho_F M (u - u*)/cdt (theta =
%   1, time_step = cdt). Both media take the same step: under BDF2 a step
%   from a state without the level n-1 (the initial state) is an Euler
%   step on both, and a state where only one medium has its level n-1 is
%   the error StokesBiotCoupling:historyLevels (StokesBiotCoupling
%   .initial_state and advance_state keep the two levels together). With
%   rho_F = 0 the fluid is quasi-static and state.fluid is not read.
%   sources: struct with fields T_1, T_2, T_3, T_4 (the interface data of
%   StokesBiotCoupling.get_manufactured_solution, whose mms.at(t) may be
%   passed as is); a missing or empty field is zero. An empty state is the
%   zero history (one step from rest).
%
%   pore_pressure_mean: the mean of p over the Biot mesh, imposed only
%   when the pore-pressure level is otherwise free (s0 = 0, alpha = 1 and
%   no boundary off Sigma that fixes it; see InterfaceProblem.stokes_biot).
%   The fourth output of StokesBiotCoupling.solve carries the level
%   scalar as global_system.level_scalar.
%
%   Convection: a convective_velocity in the fluid Data is passed to
%   IncompressibleElasticityMatrices. The transmission rows read the
%   stabilization from fluid_mesh.tau (S = mu tau), so the element must
%   use the same S: the fluid Data.convective_stabilization must then be
%   false (error StokesBiotCoupling:stabilization), and the caller puts
%   tau(beta) on the fluid mesh (S - rho_c beta . n/2 > 0 on every element
%   face; see StokesBiotCoupling.solve_navier_stokes).
%
%   Cache: the fourth output of IncompressibleElasticityMatrices (the
%   iterate-independent fluid element blocks) from a previous call on the
%   same fluid mesh, degree and fluid_mesh.tau; empty builds them. The
%   fourth output cache returns them for the next call (e.g. the next
%   Picard iterate or time step).
%
%   Mesh sizes at Sigma: the Biot mesh must not be finer than the fluid
%   mesh. The data the Biot side receives come from the fluid faces, and
%   on finer Biot faces their error is not superclose to the face
%   projection: the Biot stress then converges at about O(h) for kappa =
%   O(1) (the Darcy flux in R_n amplifies the loss by a factor ~ kappa),
%   and with friction_bjs = 0 the fine-scale part of eta_hat . t has no
%   equation (a singular system). A Biot side with more interface faces
%   than the fluid side gives the warning StokesBiotCoupling:biotFiner,
%   and with friction_bjs = 0 the error StokesBiotCoupling:biotFiner. The
%   check compares the total face counts: a Biot mesh finer on part of
%   Sigma only, with fewer faces overall, is not detected.
%
%See also: StokesBiotCoupling.solve, InterfaceProblem.stokes_biot,
%          BiotProblem, FluidStructureCoupling.interface_problem
arguments
    fluid_problem IncompressibleElasticityProblem
    biot_problem  BiotProblem
    fluid_mesh    SimplicialGrid
    biot_mesh     SimplicialGrid
    interface     AbstractCurve
    sources       struct = struct.empty
    options.Order double {mustBeScalarOrEmpty} = fluid_problem.Order
    options.friction_bjs (1,1) double {mustBeNonnegative} = 0.0
    options.pore_pressure_mean (1,1) double = 0.0
    options.state struct = struct.empty
    options.n_quadpoints      double {mustBeScalarOrEmpty} = []
    options.n_quadpoints_path double {mustBeScalarOrEmpty} = []
    options.vertex_tol        (1,1) double = 1e-8
    options.vertex_tol_relative double = []
    options.Cache struct = struct.empty
end
k = options.Order;
assert(fluid_problem.Order == biot_problem.Order, ...
    "StokesBiotCoupling:order", "PDE problems must have equal polynomial degree!");
fd = fluid_problem.Data;
bd = biot_problem.Data;
% assert(fd.density == 0, "StokesBiotCoupling:fluidInertia", ...
%     "The Stokes-Biot coupling is quasi-static: the fluid density must be 0 (got %g).", ...
%     fd.density);
assert(fd.density >= 0, "StokesBiotCoupling:fluidInertia", ...
    "The fluid density must be nonnegative (got %g).", fd.density);
if fd.density > 0
    assert(fd.time_step == bd.time_step, "StokesBiotCoupling:timeStep", ...
        "The fluid time_step (%g) must equal the Biot time_step (%g): one step " + ...
        "for both media.", fd.time_step, bd.time_step);
end
assert(fd.pressure_mean == "pinned", "StokesBiotCoupling:closure", ...
    "The fluid pressure level is fixed by the R_p row: pressure_mean must be ""pinned"".");
Pk = DubinerBasis(k);
n_fluid = numel(fluid_mesh.intfaces);
n_biot = numel(biot_mesh.intfaces);
if n_biot > n_fluid
    assert(options.friction_bjs > 0, "StokesBiotCoupling:biotFiner", ...
        "friction_bjs = 0 with a Biot mesh finer than the fluid mesh at Sigma " + ...
        "(%d > %d interface faces) is singular: mesh the Biot side no finer " + ...
        "than the fluid side.", n_biot, n_fluid);
    warning("StokesBiotCoupling:biotFiner", ...
        "The Biot mesh is finer than the fluid mesh at Sigma (%d > %d interface " + ...
        "faces): the Biot stress converges suboptimally. Mesh the Biot side no " + ...
        "finer than the fluid side.", n_biot, n_fluid);
end

%% ---- fluid: incompressible elasticity (PEERS stress / velocity), rho_F >= 0
biot_state = struct.empty;
if ~isempty(options.state) && isfield(options.state, "biot")
    biot_state = options.state.biot;
end
[cdt, history] = biot_step_history(bd.time_scheme, bd.time_step, biot_state);

fluid_bc = boundary_conditions_of(fd);
if fd.density > 0
    % the Biot scheme and step, with the fluid's own history (u* = w*)
    fluid_state = struct.empty;
    if ~isempty(options.state) && isfield(options.state, "fluid")
        fluid_state = options.state.fluid;
    end
    [cdt_fluid, history_fluid] = fluid_step_history(bd.time_scheme, bd.time_step, fluid_state);
    assert(cdt_fluid == cdt, "StokesBiotCoupling:historyLevels", ...
        "The fluid and Biot states must both have, or both lack, the level n-1 " + ...
        "(state.fluid.velocity_previous, state.biot.displacement_skeleton_previous): " + ...
        "advance them together with StokesBiotCoupling.advance_state.");
    fluid_opts = {"density", fd.density, "theta", 1.0, "time_step", cdt_fluid, ...
                  "state", struct("velocity", {history_fluid.velocity}), ...
                  "natural_faces", fluid_bc.natural_faces(fluid_mesh)};
else
    fluid_opts = {"density", fd.density, "time_step", fd.time_step, ...
                  "natural_faces", fluid_bc.natural_faces(fluid_mesh)};
end
if isfield(fd, "convective_velocity") && ~isempty(fd.convective_velocity)
    assert(~(isfield(fd, "convective_stabilization") && fd.convective_stabilization), ...
        "StokesBiotCoupling:stabilization", ...
        "With a convective_velocity the fluid Data.convective_stabilization must be " + ...
        "false: the interface rows use S = mu tau from the fluid mesh, so put " + ...
        "tau(beta) on fluid_mesh.tau instead.");
    fluid_opts = [fluid_opts, {"convective_velocity", fd.convective_velocity, ...
                               "convective_stabilization", fd.convective_stabilization}];
end
fluid_pb = FEProblem(Order = k, Spaces = [ ...
    FESpace([PEERSBasis(k), Pk.to_vector(), Pk, DubinerBasis(0)], Condense = true), ...
    FESpace([LegendreBasis(k).to_vector(2), DubinerBasis(0)])]);
fluid_pb.Mesh = fluid_mesh;
fluid_pb.Data = fd;
[lhs, flux, rhs, cache] = IncompressibleElasticityMatrices(k, fluid_mesh, fd.viscosity, ...
    fd.source_term, "pressure_mean", "pinned", fluid_opts{:}, "Cache", options.Cache);
fluid_pb.Operator = AffineOperator(lhs + rhs) + flux;
fluid_pb.BoundaryConditions = fluid_bc;

%% ---- Biot -----------------------------------------------------------------
biot_problem.Data.state = options.state;
biot_pb = get_declarative_problem(biot_problem, biot_mesh, k);

%% ---- the coupling -----------------------------------------------------------
coupled = InterfaceProblem.stokes_biot(fluid_pb, biot_pb, interface, ...
    viscosity = fd.viscosity, young_modulus = bd.young_modulus, ...
    poisson_ratio = bd.poisson_ratio, biot_alpha = bd.biot_alpha, ...
    diffusion = bd.diffusion, stabilization_darcy = bd.stabilization_darcy, ...
    stabilization_solid = bd.stabilization_solid, ...
    storage = bd.storage, pore_pressure_mean = options.pore_pressure_mean, ...
    friction_bjs = options.friction_bjs, effective_time_step = cdt, ...
    history = history, ...
    source_neumann = vector_source(sources, "T_1"), ...
    source_dirichlet = vector_source(sources, "T_2"), ...
    source_pressure = scalar_source(sources, "T_3"), ...
    source_slip = vector_source(sources, "T_4"), ...
    quadrature_patch = options.n_quadpoints, quadrature_path = options.n_quadpoints_path, ...
    vertex_tol = options.vertex_tol, vertex_tol_relative = options.vertex_tol_relative);
end

function f = vector_source(sources, name)
% a {2,1} cell of handles as a 1 x 2 cell; missing or empty = zero
f = {@(x,y) 0*x.*y, @(x,y) 0*x.*y};
if ~isempty(sources) && isfield(sources, name) && ~isempty(sources.(name))
    f = reshape(sources.(name), 1, []);
end
end

function f = scalar_source(sources, name)
f = @(x,y) 0*x.*y;
if ~isempty(sources) && isfield(sources, name) && ~isempty(sources.(name))
    f = sources.(name);
end
end

function bc = boundary_conditions_of(data)
% the low-level get_boundary_conditions convention (as in
% FluidStructureCoupling.interface_problem)
if isfield(data, "boundary_conditions") && ~isempty(data.boundary_conditions)
    specs = data.boundary_conditions;
else
    gD = data.datum_dirichlet;  gN = data.datum_neumann;
    specs = [ ...
        boundary_condition("Dirichlet", x_data = gD{1}, y_data = gD{2});
        boundary_condition("Neumann", x = "neumann", y = "neumann", ...
            x_data = gN{1}, y_data = gN{2})];
end
bc = ComponentwiseBoundaryConditions(specs);
end
