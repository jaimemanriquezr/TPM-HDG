function state = initial_state(fluid_problem, solid_problem, fluid_mesh, solid_mesh, ...
    interface, source_neumann, source_dirichlet, options)
%INITIAL_STATE  Time-stepping state of the TPM-coupled FSI problem at t = 0.
%
%   state = FluidStructureCoupling.initial_state(fluid_problem, solid_problem, ...
%       fluid_mesh, solid_mesh, interface, T1, T2, ...)
%
%   Positional arguments as FluidStructureCoupling.solve, including the
%   per-side pair {curve_fluid, curve_solid} (InterfaceProblem.fsi). Builds the state
%   that solve(..., "onestep" | "newmark", state = state) marches from,
%   with the physical initial data of the two problems (solid
%   Data.initial_displacement / initial_velocity, fluid
%   Data.initial_velocity, cells of handles) and their sources set to
%   their t = 0 values, as a sequence:
%     1. the solid state, LinearElasticityProblem.initial_state: the
%        projections (the HDG elliptic projections under newmark), their
%        (sigma, gamma) and, under newmark, the consistent acceleration
%        a_S and its skeleton trace;
%     2. the fluid velocity u_F = Pi u_F(0) and its skeleton trace,
%        IncompressibleElasticityProblem.initial_state;
%     3. when the fluid step reads its acceleration (theta < 1) or
%        with_pressure = true: the fluid acceleration problem with the
%        kinematic datum on Sigma, a_F = a_S + d/dt T2 (the kinematic
%        condition u_F = d/dt u_S + T2 differentiated in time). a_S on
%        Sigma is the solid skeleton acceleration transported along the
%        transfer paths (trace_extension_coefficients of the lifted
%        acceleration stress), as the T2 history is in each step; it is
%        imposed on the fluid faces tagged "Interface" through
%        TransferEssentialBoundaryConditions.field_gradient with the TPM
%        map of those faces onto Sigma (interface_datum). This needs a
%        newmark solid (FluidStructureCoupling:onestepAcceleration
%        otherwise).
%   Pressure level of step 3: with a natural (traction) face on the fluid
%   boundary the fluid's own closure fixes it; without one the problem is
%   solved under pressure_mean = "constrained" and the level is then
%   shifted so that the normal component of the traction residual below
%   has zero mean on Sigma, which is how the traction balance fixes the
%   level in the coupled step.
%
%   state.diagnostics (step 3 only; [] fields otherwise):
%     traction_residual           ||sigma_F n_F + sigma_S n_S - T1||_L2(Sigma)
%                                 at t = 0, the stresses of the state
%                                 extrapolated to the interface quadrature
%                                 points as the T1 row does
%     traction_scale              max of ||sigma_F n_F||, ||sigma_S n_S||,
%                                 ||T1|| in the same norm
%     traction_residual_relative  their ratio, the scale floored at 1e-10
%                                 (zero-traction data compare round-off
%                                 against that absolute level)
%   With a curve pair the fluid traction of the residual carries the
%   Nanson factor, as in the T1 row, and the step-3 datum on Sigma(t) is
%   a_S(X(s)) + d/dt T2: the transport term -(v_S . grad) u_F of the
%   moving interface is not included, so it is exact for an interface
%   at rest at t = 0 (v_S = 0 on Sigma_0).
%   The kinematic route leaves the dynamic condition as a compatibility
%   check of the initial data. When the relative residual exceeds 0.13
%   the warning FluidStructureCoupling:incompatibleInitialData is raised.
%   For compatible data the residual is a discretization error that
%   decreases with h; on smooth manufactured data it stays at or below
%   0.11 down to the coarsest meshes tried (k = 1, h = 1/4 fitted and
%   h = 1/9 with both meshes retreated by h/2). A solid with twice the
%   stiffness of those data gives 0.16 on the coarsest mesh and a rest
%   start under a solid load gives 1. On coarse unfitted meshes whose
%   gap is a large fraction of h the transported solid acceleration is
%   less accurate and compatible data can exceed the threshold.
%
%   Options:
%     Order                  polynomial degree (fluid_problem.Order)
%     with_pressure          solve step 3 also when theta = 1 (false)
%     source_dirichlet_rate  {d/dt T2_x; d/dt T2_y} at t = 0 (zero)
%     n_quadpoints, n_quadpoints_path, vertex_tol, vertex_tol_relative,
%     tip_mapping, collapsed_normal, mapping_1, mapping_2
%                            the partition, as in FluidStructureCoupling.solve
%                            (see get_partitions for the accuracy limitation
%                            of tip_mapping = "collapse"); collapsed patches
%                            carry no weight in the pressure-level shift and
%                            the traction diagnostics
%
%See also: FluidStructureCoupling.solve, LinearElasticityProblem.initial_state,
%          IncompressibleElasticityProblem.initial_state, interface_datum
arguments
    fluid_problem    IncompressibleElasticityProblem
    solid_problem    LinearElasticityProblem
    fluid_mesh       SimplicialGrid
    solid_mesh       SimplicialGrid
    interface        {mustBeInterfaceCurves}
    source_neumann   cell = {@(x,y) 0*x.*y, @(x,y) 0*x.*y}
    source_dirichlet cell = {@(x,y) 0*x.*y, @(x,y) 0*x.*y} %#ok<INUSA> positional, as solve
    options.Order double {mustBeScalarOrEmpty} = fluid_problem.Order
    options.with_pressure (1,1) logical = false
    options.source_dirichlet_rate cell = {@(x,y) 0*x.*y, @(x,y) 0*x.*y}
    options.n_quadpoints      double {mustBeScalarOrEmpty} = []
    options.n_quadpoints_path double {mustBeScalarOrEmpty} = []
    options.vertex_tol        (1,1) double = 1e-8
    options.vertex_tol_relative double = []
    options.tip_mapping (1,1) string ...
        {mustBeMember(options.tip_mapping, ["spread", "collapse"])} = "spread"
    options.collapsed_normal (1,1) string ...
        {mustBeMember(options.collapsed_normal, ["path", "face", "bisector"])} = "path"
    options.mapping_1 (1,1) string = ""
    options.mapping_2 (1,1) string = ""
end
% relative traction residual above which the initial data are flagged
% as incompatible (see the help)
INCOMPATIBLE_TOLERANCE = 0.13;
% absolute floor of the traction scale
TRACTION_SCALE_FLOOR = 1e-10;
k = options.Order;

state = initial_state(solid_problem, solid_mesh, k);
state.diagnostics = struct("traction_residual", [], "traction_scale", [], ...
    "traction_residual_relative", []);
if get_theta(fluid_problem) == 1 && ~options.with_pressure
    fluid = initial_state(fluid_problem, fluid_mesh, k);
    state.fluid = fluid.fluid;
    return
end
if solid_problem.Data.time_scheme ~= "newmark"
    error("FluidStructureCoupling:onestepAcceleration", ...
        "The fluid acceleration at t = 0 needs the solid acceleration, which " + ...
        "only a newmark solid carries (time_scheme = ""%s"").", solid_problem.Data.time_scheme);
end
s = state.solid;

%% T: the solid acceleration on Sigma, plus the rate of T2
sq = options.n_quadpoints;       if isempty(sq), sq = k + 1; end
tq = options.n_quadpoints_path;  if isempty(tq), tq = k + 2; end
[curve_fluid, curve_solid] = interface_curves(interface);
[part_1, part_2] = curve_fluid.get_partitions(fluid_mesh, solid_mesh, ...
    quadrature_patch = sq, quadrature_path = tq, ...
    vertex_tol = options.vertex_tol, vertex_tol_relative = options.vertex_tol_relative, ...
    tip_mapping = options.tip_mapping, collapsed_normal = options.collapsed_normal, ...
    mapping_1 = options.mapping_1, mapping_2 = options.mapping_2, ...
    interface_2 = curve_solid);
ext_a = trace_extension_coefficients(solid_problem, solid_mesh, part_2, ...
    s.stress_acceleration, s.rotation_acceleration, Order = k);
a_sigma = s.acceleration_skeleton(:, part_2.faces(:)) + ext_a;
[a_solid, mapping] = interface_datum(curve_fluid, fluid_mesh, part_2, a_sigma, ...
    tip_mapping = options.tip_mapping, mapping = options.mapping_1);
rate = options.source_dirichlet_rate;
datum = {@(x, y) a_solid{1}(x, y) + rate{1}(x, y) + 0*x; ...
         @(x, y) a_solid{2}(x, y) + rate{2}(x, y) + 0*x};
sigma_bc = TransferEssentialBoundaryConditions.field_gradient(Datum(datum), ...
    "Interface", mapping);

%% F: the fluid acceleration problem
fluid_pb = fluid_problem;
shift_level = isempty(natural_faces(fluid_problem.Data, fluid_mesh));
if shift_level
    fluid_pb.Data.pressure_mean = "constrained";
end
fluid = initial_state(fluid_pb, fluid_mesh, k, with_pressure = true, ...
    interface_conditions = sigma_bc);
state.fluid = fluid.fluid;

%% the traction residual at t = 0
% the two partitions must pair their quadrature nodes: the same
% interface parameters, and with one curve the same points
same_points = isempty(curve_solid) || ...
    max(abs(part_1.quadrature.interface - part_2.quadrature.interface), [], "all") ...
    <= 1e-12*max(1, max(abs(part_1.quadrature.interface), [], "all"));
assert(isequal(part_1.vertices_interface, part_2.vertices_interface) && same_points, ...
    "FluidStructureCoupling:interfacePoints", ...
    "The two partitions must share their interface quadrature parameters " + ...
    "(and, with one curve, their points).");
[~, w] = lgwt(sq, 0, 1);
% A collapsed patch (tip_mapping = "collapse") has the fluid face's
% length on side 1 but a single point of Sigma, where the normal is
% undefined (quadrature.normals there is the arbitrary collapsed_normal);
% it carries no traction balance and gets zero weight here.
patch_weight = reshape(part_1.area, 1, 1, []);
if isfield(part_1, "collapsed")
    patch_weight = patch_weight .* reshape(~part_1.collapsed, 1, 1, []);
end
weights = reshape(w, 1, sq) .* patch_weight;   % 1 x sq x n_patches
norm_sigma = @(v) sqrt(sum(weights .* sum(v.^2, 1), "all"));
n_F = -part_1.quadrature.normals;
t_S = interface_traction(s.stress, solid_mesh, part_2, k, part_2.quadrature.normals);
T1 = cellfun(@(f) f(part_1.quadrature.interface(1, :, :), part_1.quadrature.interface(2, :, :)) ...
    + 0*part_1.quadrature.interface(1, :, :), source_neumann(:), "UniformOutput", false);
T1 = cat(1, T1{:});
% the fluid traction per unit reference length: Nanson factor of a pair
nanson = 1;
if ~isempty(curve_solid)
    nanson = part_1.quadrature.speed ./ part_2.quadrature.speed;
end
t_F = nanson .* interface_traction(state.fluid.stress, fluid_mesh, part_1, k, n_F);
if shift_level
    % the level that zeroes the mean of the normal residual: p + c shifts
    % sigma_F n_F by -c n_F
    c = sum(weights .* sum((t_F + t_S - T1) .* n_F, 1), "all") / sum(weights .* nanson, "all");
    phi_1 = dubiner_basis(0, 0, 0);   % the constant Dubiner function
    d2 = nchoosek(k + 2, 2);
    state.fluid.pressure(1, :) = state.fluid.pressure(1, :) + c/phi_1;
    state.fluid.stress([1, 3*d2 + 1], :) = state.fluid.stress([1, 3*d2 + 1], :) - c/phi_1;
    t_F = nanson .* interface_traction(state.fluid.stress, fluid_mesh, part_1, k, n_F);
end
residual = norm_sigma(t_F + t_S - T1);
scale = max([norm_sigma(t_F), norm_sigma(t_S), norm_sigma(T1)]);
relative = residual/max(scale, TRACTION_SCALE_FLOOR);
state.diagnostics = struct("traction_residual", residual, "traction_scale", scale, ...
    "traction_residual_relative", relative);
if relative > INCOMPATIBLE_TOLERANCE
    warning("FluidStructureCoupling:incompatibleInitialData", ...
        "The initial data violate the traction balance on Sigma at t = 0: " + ...
        "residual %.3e, %.3e relative to the traction scale %.3e.", ...
        residual, relative, scale);
end
end

function t = interface_traction(sigma_dof, mesh, partition, k, normals)
% sigma n at the interface quadrature points (2 x sq x n_patches), the
% PEERS stress of each patch's element extrapolated there, as
% transfer_traction_to_interface evaluates it
s_quad = size(partition.quadrature.local, 2);
n_of_patches = size(partition.faces, 3);
G = get_affine_map(mesh, partition.elements, inverse = true);
t = zeros(2, s_quad, n_of_patches);
for e = 1:n_of_patches
    K = partition.elements(e);
    s_iface = partition.quadrature.interface_2d(1, :, e);
    t_iface = partition.quadrature.interface_2d(2, :, e);
    P = PEERS_basis(k, s_iface(:), t_iface(:), G(:, :, e));   % s x nb x 2 x 2
    sigma = tensorprod(P, sigma_dof(:, K), 2, 1);             % s x 2 x 2
    n = normals(:, :, e);
    t(1, :, e) = squeeze(sigma(:, 1, 1)).' .* n(1, :) + squeeze(sigma(:, 1, 2)).' .* n(2, :);
    t(2, :, e) = squeeze(sigma(:, 2, 1)).' .* n(1, :) + squeeze(sigma(:, 2, 2)).' .* n(2, :);
end
end

function faces = natural_faces(data, mesh)
% the faces with a traction component, by the low-level boundary-condition
% convention of interface_problem
if isfield(data, "boundary_conditions") && ~isempty(data.boundary_conditions)
    specs = data.boundary_conditions;
else
    gD = data.datum_dirichlet;  gN = data.datum_neumann;
    specs = [ ...
        boundary_condition("Dirichlet", x_data = gD{1}, y_data = gD{2});
        boundary_condition("Neumann", x = "neumann", y = "neumann", ...
            x_data = gN{1}, y_data = gN{2})];
end
faces = ComponentwiseBoundaryConditions(specs).natural_faces(mesh);
end
