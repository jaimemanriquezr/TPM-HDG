function [coupled, fluid_pb, solid_pb, cache] = interface_problem(fluid_problem, solid_problem, ...
    fluid_mesh, solid_mesh, interface, source_neumann, source_dirichlet, scheme, options)
%INTERFACE_PROBLEM  The TPM-coupled FSI problem, from low-level problem objects.
%
%   Same positional signature as FluidStructureCoupling.solve, which is a
%   thin wrapper over this: builds the two declarative FEProblems that
%   InterfaceProblem.fsi needs from the low-level objects' Data structs
%   and returns the coupled problem; solve it with
%   [fs, ss] = solve(coupled, k).
%
%   Repeated builds in which only the fluid data change (e.g. a Picard
%   loop) pass a cache struct and get it back,
%     [coupled, ~, ~, cache] = interface_problem(..., Cache = cache);
%     [fs, ss, ~, cache] = solve(coupled, k, Cache = cache);
%   starting from cache = struct("frozen", 2). The fluid's geometric
%   blocks are kept in cache.geometry (see IncompressibleElasticityMatrices)
%   and, under frozen = 2, the solid problem is built once and kept in
%   cache.solid_problem (its data are assumed unchanged). Without a cache
%   the fourth output is struct.empty. The solid state (the state option)
%   is part of those data.
%
%   The state option carries the time-stepping state (see
%   FluidStructureCoupling.initial_state and FluidStructureCoupling.solve):
%   state.fluid and state.solid are the factories' volume histories, and
%   under scheme = "onestep" / "newmark" state.solid is also passed as
%   the T2 history of InterfaceProblem.fsi: its skeleton traces give h0,
%   and its stress and rotation fields the transport Ext(h0) onto Sigma
%   that InterfaceProblem.solve builds on the solid side of its
%   partition (see the history option of InterfaceProblem.fsi). On a
%   fitted solid Ext(h0) vanishes to round-off. An empty state (or an
%   empty side, or [] fields) is the zero history; a steady scheme reads
%   no T2 history.
%
%   Two optional solid Data flags reformulate the solid for large lambda_S
%   without changing its discrete solution; they are alternatives:
%     Data.factored_condensation  keeps each element's constant isotropic
%       PEERS stress coefficient global (FESpace.UncondensedModes,
%       LinearElasticityProblem.isotropic_stress_dofs); a
%       LinearElasticityProblem option.
%     Data.pressure_split  keeps the element mean of tr(sigma_S) global
%       (LinearElasticityMatrices pressure_split: local xi, global rho per
%       element); set it on the Data struct,
%       solid_problem.Data.pressure_split = true.
%
%   An optional solid Data.constitutive_source (PEERS dofs per element,
%   PEERSBasis(k).Dimension x n_elements, as in LinearElasticityProblem)
%   holds the moments (r, tau)_K of a right-hand side r of the solid's
%   constitutive row, A sigma - eps(u) + gamma = r; it is added to that
%   row's right-hand side. Absent or empty, the solid operator has no such
%   term.
%
%   An optional solid Data.eulerian_blocks replaces the solid operator by
%   given declarative blocks (a struct with fields lhs, flux and rhs,
%   sums of IndexedArray blocks), on the spaces local
%   [sigma (PEERS), w (P_k^2), gamma (P_k), u (P_k^2)] (ids 1-4) and
%   global [w_hat, u_hat] (ids 5, 6; Legendre, x over y per face): a fully
%   Eulerian solid with the velocity unknown w (scaled as the T2 row
%   needs, w = time_step * velocity under "onestep") and the displacement
%   u. The transmission rows act on (sigma, w, w_hat); the path gradient
%   of the T2 row is the lifting L_h of (w, w_hat) (as for a
%   HyperelasticProblem iterate below); the rows of u_hat on Sigma stay
%   the solid's own. Sigma must be the solid's whole boundary (no solid
%   boundary conditions). Absent or empty, the solid is built as
%   described here.
%
%   A HyperelasticProblem solid with Data.iterate set is assembled as its
%   Newton step at that iterate (HyperelasticProblem.newton_blocks, the
%   declarative form of its get_matrices); the caller runs the Newton
%   loop. The traction of the T1 row is then sigma_S n_S with sigma_S the
%   first Piola stress unknown and n_S the normal of the reference curve.
%
%   Boundary conditions follow the low-level convention exactly:
%   Data.boundary_conditions (componentwise specs) when set, otherwise the
%   "Dirichlet"/"Neumann" tags with datum_dirichlet / datum_neumann.
%
%See also: FluidStructureCoupling.solve, InterfaceProblem.fsi,
%          ComponentwiseBoundaryConditions
arguments
    fluid_problem    IncompressibleElasticityProblem
    solid_problem    LinearElasticityProblem
    fluid_mesh       SimplicialGrid
    solid_mesh       SimplicialGrid
    interface        {mustBeInterfaceCurves}   % one curve, or {curve_fluid, curve_solid} (InterfaceProblem.fsi)
    source_neumann   cell = {@(x,y) 0*x.*y, @(x,y) 0*x.*y}
    source_dirichlet cell = {@(x,y) 0*x.*y, @(x,y) 0*x.*y}
    scheme           (1,1) string {mustBeMember(scheme, ["steady", "onestep", "newmark"])} = "steady"
    options.Order double {mustBeScalarOrEmpty} = fluid_problem.Order
    options.patch_weight      (1,1) string {mustBePatchWeight} = "rectangular"
    options.sliver_correction (1,1) logical = true
    options.vertex_tol        (1,1) double = 1e-8
    options.vertex_tol_relative double = []
    % per-side point-to-interface maps (get_partitions); "" keeps the curve's own
    options.mapping_1 (1,1) string = ""
    options.mapping_2 (1,1) string = ""
    options.n_quadpoints      double {mustBeScalarOrEmpty} = []
    options.n_quadpoints_path double {mustBeScalarOrEmpty} = []
    % the time-stepping state, see above
    options.state struct = struct.empty
    % interface law; see InterfaceProblem.fsi
    options.interface_law (1,1) string ...
        {mustBeMember(options.interface_law, ["noslip", "navier", "full"])} = "noslip"
    options.slip_coefficient (1,1) double {mustBeNonnegative} = 0.0
    options.source_slip cell = {@(x,y) 0*x.*y, @(x,y) 0*x.*y}
    % see InterfaceProblem.fsi
    options.swap_transmission (1,1) logical = false
    % see get_partitions and InterfaceProblem.fsi
    % (tip_mapping = "collapse" is a fallback where "spread" fails with
    % orphanFace; where "spread" works it is more accurate, and collapse
    % is not meant for convergence studies)
    options.tip_mapping (1,1) string ...
        {mustBeMember(options.tip_mapping, ["spread", "collapse"])} = "spread"
    options.collapsed_normal (1,1) string ...
        {mustBeMember(options.collapsed_normal, ["path", "face", "bisector"])} = "path"
    % see above
    options.Cache struct = struct.empty
end
k = options.Order;
assert(fluid_problem.Order == solid_problem.Order, ...
    "FluidStructureCoupling:order", "PDE problems must have equal polynomial degree!");
fd = fluid_problem.Data;
sd = solid_problem.Data;
no_dof_histories(fd);
no_dof_histories(sd);
% the time step is read from the fluid problem
dt = fd.time_step;
if fd.density > 0 && sd.density > 0
    assert(sd.time_step == dt, "FluidStructureCoupling:timeStep", ...
        "The fluid and solid problems must share the time step (%g vs %g).", dt, sd.time_step);
end
closure = fd.pressure_mean;
Pk = DubinerBasis(k);

%% ---- fluid: incompressible elasticity (PEERS stress / velocity) --------
fluid_bc = boundary_conditions_of(fd);
fluid_opts = {"density", fd.density, "time_step", dt, ...
              "state", side_of(options.state, "fluid"), ...
              "natural_faces", fluid_bc.natural_faces(fluid_mesh)};
if fd.density > 0
    fluid_opts = [fluid_opts, {"theta", get_theta(fluid_problem)}];
end
if isfield(fd, "convective_velocity") && ~isempty(fd.convective_velocity)
    fluid_opts = [fluid_opts, {"convective_velocity", fd.convective_velocity, ...
                               "convective_stabilization", fd.convective_stabilization}];
    if isfield(fd, "convective_form")
        fluid_opts = [fluid_opts, {"convective_form", fd.convective_form}];
    end
end
fluid_global = [LegendreBasis(k).to_vector(2), DubinerBasis(0)];
if closure == "constrained"
    fluid_global = [fluid_global, ConstantBasis()];
end
fluid_pb = FEProblem(Order = k, Spaces = [ ...
    FESpace([PEERSBasis(k), Pk.to_vector(), Pk, DubinerBasis(0)], Condense = true), ...
    FESpace(fluid_global)]);
fluid_pb.Mesh = fluid_mesh;
% the coefficients, for postprocessing the solution (divergence_free_velocity)
fluid_pb.Data = fd;
cache = options.Cache;
if isempty(cache)
    [lhs, flux, rhs] = IncompressibleElasticityMatrices(k, fluid_mesh, fd.viscosity, ...
        fd.source_term, "pressure_mean", closure, fluid_opts{:});
else
    geometry = struct.empty;
    if isfield(cache, "geometry"), geometry = cache.geometry; end
    [lhs, flux, rhs, cache.geometry] = IncompressibleElasticityMatrices(k, fluid_mesh, ...
        fd.viscosity, fd.source_term, "pressure_mean", closure, fluid_opts{:}, Cache = geometry);
end
fluid_pb.Operator = AffineOperator(lhs + rhs) + flux;
fluid_pb.BoundaryConditions = fluid_bc;

%% ---- solid: linear elasticity (+ Winkler) ------------------------------
frozen = ~isempty(cache) && isfield(cache, "frozen") && isequal(cache.frozen, 2);
solid_key = struct("k", k, "n_el", size(solid_mesh.elements, 2), ...
    "n_faces", size(solid_mesh.faces, 2));
eulerian = isfield(sd, "eulerian_blocks") && ~isempty(sd.eulerian_blocks);
if eulerian
    % an Eulerian velocity solid: its operator is given (see above)
    assert(~frozen, "FluidStructureCoupling:eulerian", ...
        "The operator of an Eulerian solid changes with its iterate and is not cached.");
    blocks = sd.eulerian_blocks;
    solid_pb = FEProblem(Order = k, Spaces = [ ...
        FESpace([PEERSBasis(k), Pk.to_vector(), Pk, Pk.to_vector()], Condense = true), ...
        FESpace([LegendreBasis(k).to_vector(2), LegendreBasis(k).to_vector(2)])]);
    solid_pb.Mesh = solid_mesh;
    solid_pb.Data = sd;
    solid_pb.Operator = AffineOperator(blocks.lhs) + blocks.flux + blocks.rhs;
    % Sigma is its whole boundary: no boundary conditions of its own
    solid_pb.BoundaryConditions = BoundaryConditions.empty;
elseif frozen && isfield(cache, "solid_problem")
    if ~isequal(cache.solid_key, solid_key)
        error("FluidStructureCoupling:cacheKey", ...
            "The cached solid problem was built for another mesh or degree.");
    end
    solid_pb = cache.solid_problem;
else
    solid_opts = {"density", sd.density, "time_step", sd.time_step, ...
                  "time_scheme", sd.time_scheme, "newmark_beta", sd.newmark_beta, ...
                  "state", side_of(options.state, "solid"), ...
                  "hooke_modulus", sd.hooke_modulus};
    % optional Data.pressure_split: the element mean of tr(sigma_S) kept global
    % (LinearElasticityMatrices pressure_split; robust as lambda_S -> infinity)
    split = isfield(sd, "pressure_split") && sd.pressure_split;
    solid_opts = [solid_opts, {"pressure_split", split}];
    if split
        solid_pb = FEProblem(Order = k, Spaces = [ ...
            FESpace([PEERSBasis(k), Pk.to_vector(), Pk, DubinerBasis(0)], Condense = true), ...
            FESpace([LegendreBasis(k).to_vector(2), DubinerBasis(0)])]);
    else
        solid_pb = FEProblem(Order = k, Spaces = [ ...
            FESpace([PEERSBasis(k), Pk.to_vector(), Pk], Condense = true), ...
            FESpace(LegendreBasis(k).to_vector(2))]);
    end
    % optional Data.factored_condensation: keep the constant isotropic
    % stress coefficient global (FESpace.UncondensedModes on the stress
    % space); an exact re-elimination, no new spaces or blocks
    factored = isfield(sd, "factored_condensation") && sd.factored_condensation;
    assert(~(factored && split), "FluidStructureCoupling:lambdaForms", ...
        "pressure_split and factored_condensation are alternatives; set one.");
    if factored
        solid_pb.Spaces(1).UncondensedModes = {LinearElasticityProblem.isotropic_stress_dofs(k)};
    end
    solid_pb.Mesh = solid_mesh;
    solid_pb.Data = sd;
    [slhs, sflux, srhs] = LinearElasticityMatrices(k, solid_mesh, ...
        sd.young_modulus, sd.poisson_ratio, solid_opts{:});
    solid_pb.Operator = AffineOperator(slhs) + sflux;
    if ~isempty(srhs), solid_pb.Operator = solid_pb.Operator + srhs; end
    % optional Data.constitutive_source: the moments (r, tau)_K of a
    % right-hand side r of the constitutive row (see above)
    if isfield(sd, "constitutive_source") && ~isempty(sd.constitutive_source)
        solid_pb.Operator = solid_pb.Operator + IndexedMatrix(IndexedArray(reshape( ...
            sd.constitutive_source, [], 1, 1, 1, size(solid_mesh.elements, 2)), [-1, 1]));
    end
    w = TestFunction(Pk);
    solid_pb.Operator = solid_pb.Operator + IndexedMatrix(IndexedArray(cat(1, ...
        double(Datum(sd.source_term{1}) * w * Measure(solid_mesh, 'dx')), ...
        double(Datum(sd.source_term{2}) * w * Measure(solid_mesh, 'dx'))), [-1, 2]));
    solid_pb.BoundaryConditions = boundary_conditions_of(sd);
    % a HyperelasticProblem with Data.iterate: its Newton step at the
    % iterate (HyperelasticProblem.newton_blocks); without one, the
    % linear solid
    if isa(solid_problem, "HyperelasticProblem") && ~isempty(sd.iterate)
        assert(~frozen && ~split, "FluidStructureCoupling:hyperelastic", ...
            "The Newton step of a HyperelasticProblem is not cached and has no pressure_split form.");
        % the step's state carries the history of a viscous solid's strain rate
        newton_problem = solid_problem;
        newton_problem.Data.state = options.state;
        [newton_lhs, newton_rhs] = newton_blocks(newton_problem, solid_mesh, k);
        solid_pb.Operator = solid_pb.Operator + newton_lhs + newton_rhs;
    end
    if frozen
        cache.solid_problem = solid_pb;
        cache.solid_key = solid_key;
    end
end

%% ---- the coupling ------------------------------------------------------
% the solid state is the T2 history; InterfaceProblem.solve transports it
history = struct.empty;
if scheme ~= "steady"
    history = side_of(options.state, "solid");
end
coupled = InterfaceProblem.fsi(fluid_pb, solid_pb, interface, ...
    viscosity = fd.viscosity, young_modulus = sd.young_modulus, ...
    poisson_ratio = sd.poisson_ratio, domain_measure = fd.domain_measure, ...
    scheme = scheme, time_step = dt, ...
    newmark_beta = sd.newmark_beta, newmark_gamma = sd.newmark_gamma, ...
    history = history, ...
    source_neumann = source_neumann, source_dirichlet = source_dirichlet, ...
    pressure_mean = closure, ...
    patch_weight = options.patch_weight, sliver_correction = options.sliver_correction, ...
    quadrature_patch = options.n_quadpoints, quadrature_path = options.n_quadpoints_path, ...
    vertex_tol = options.vertex_tol, vertex_tol_relative = options.vertex_tol_relative, ...
    mapping_1 = options.mapping_1, mapping_2 = options.mapping_2, ...
    interface_law = options.interface_law, slip_coefficient = options.slip_coefficient, ...
    source_slip = options.source_slip, swap_transmission = options.swap_transmission, ...
    tip_mapping = options.tip_mapping, collapsed_normal = options.collapsed_normal);
end

function no_dof_histories(data)
% the histories are read from the state option only: a dof array in a
% Data.initial_* field would be ignored, so it is an error
names = "initial_" + ["displacement", "velocity", "acceleration", "displacement_previous", ...
    "displacement_skeleton", "velocity_skeleton", "acceleration_skeleton"];
for name = names
    if isfield(data, name) && isnumeric(data.(name))
        error("FluidStructureCoupling:dofHistory", ...
            "Data.%s holds dof values, which the coupling does not read: " + ...
            "pass the histories in the state option (Data.initial_displacement " + ...
            "and Data.initial_velocity are cells of handles, read by initial_state).", name);
    end
end
end

function s = side_of(state, side)
% state.fluid / state.solid; struct.empty when the state or the side is empty
s = struct.empty;
if ~isempty(state) && isfield(state, side), s = state.(side); end
end

function bc = boundary_conditions_of(data)
% the low-level get_boundary_conditions convention, verbatim
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
