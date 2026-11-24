function [fluid_solution, solid_solution, global_system] = solve(...
    fluid_problem, solid_problem, fluid_mesh, solid_mesh, ...
    interface, source_neumann, source_dirichlet, scheme, ...
    options)
% TPM-coupled FSI solve, all schemes (scheme = "steady" (default) |
% "onestep" | "newmark"). The three schemes differ only in the T2 kinematic
% interface row; kinematic state updates are postprocesses of the
% solutions (the problems' get_kinematic_state methods):
%   steady:  <cdt*û₁ - û₂, μ₂> = <cdt*T₂, μ₂>,          cdt = dt, no history
%   onestep: same row with history H⁰ = u⁰ (backward Euler of ∂t u_S = u_F)
%   newmark: cdt = (β/γ)Δt and H⁰ = u⁰ + Δt(1-β/γ) v⁰
%            + Δt²[(1/2-β) - β(1-γ)/γ] a⁰ (elimination of a^{n+1}).
% Volume time terms live in the problems' get_matrices, not here.
% Frozen low-level reference coupler, used by the tests only (the toolbox
% entry point is FluidStructureCoupling.solve). Also carries the
% circulation, solid_rigid_modes and traction_resultants_{fluid,solid}
% anchors, which the toolbox does not expose. No slip law.
% test_interface_fsi and test_interface_oracle_parity compare
% FluidStructureCoupling.solve (built on InterfaceProblem) against this
% coupler dof for dof. It is not part of the toolbox API; do not extend it.
arguments
    fluid_problem   IncompressibleElasticityProblem
    solid_problem   LinearElasticityProblem
    fluid_mesh      SimplicialGrid
    solid_mesh      SimplicialGrid
    interface       AbstractCurve
    source_neumann  cell = {@(x,y) 0*x.*y, @(x,y) 0*x.*y}
    source_dirichlet cell = {@(x,y) 0*x.*y, @(x,y) 0*x.*y}
    scheme          (1,1) string {mustBeMember(scheme, ["steady", "onestep", "newmark"])} = "steady"
    options.Order double {mustBeScalarOrEmpty} = fluid_problem.Order
    options.n_quadpoints      double {mustBeScalarOrEmpty} = []
    options.n_quadpoints_path double {mustBeScalarOrEmpty} = []
    options.quiet             (1,1) logical = false
    options.patch_weight      (1,1) string {mustBePatchWeight} = "rectangular"
    % Ablation switch: false zeroes the sliver weight only (alpha row,
    % p0 column and domain measure stay in place)
    options.sliver_correction (1,1) logical = true
    % Tolerance for merging the two meshes' interface breakpoints into
    % the patch list
    options.vertex_tol        (1,1) double = 1e-8
    % Merge tolerance relative to the typical interface face; overrides
    % vertex_tol when set. See get_partitions.
    options.vertex_tol_relative double = []
    % EXPERIMENTAL history transport for the T2 row. The kinematic row
    % compares a PATH-EXTENDED trace at n+1 against the RAW face trace of
    % the history h0, which on a retreated mesh live at different places;
    % the mismatch is an O(delta/dt) datum error absent only when the mesh
    % sits on Sigma. Supply Ext(u^0) as (2*d1) x n_patches coefficients
    % (see trace_extension_coefficients) and it is added to h0, putting
    % both sides on Sigma. Empty = off, which is the historical behaviour.
    options.history_extension double = []
    % EXPERIMENTAL circulation anchor for closed interface loops: one
    % multiplier row enforcing the fluid trace circulation. Empty = off.
    options.circulation double {mustBeScalarOrEmpty} = []
    % EXPERIMENTAL rigid-mode anchor for an all-interface solid:
    % [t_x; t_y; rot] multiplier rows over the solid interface faces.
    % Empty = off.
    options.solid_rigid_modes (:, 1) double = []
    % EXPERIMENTAL traction-resultant anchor: [F_x; F_y; M] multiplier
    % rows enforcing the loop resultants of the extrapolated traction.
    % Empty = off.
    options.traction_resultants_fluid (:, 1) double = []
    options.traction_resultants_solid (:, 1) double = []
end
assert(fluid_problem.Order == solid_problem.Order, ...
    "PDE problems must have equal polynomial degree!")

% Pressure-level closure (see the pressure_mean option of
% IncompressibleElasticityProblem)
pressure_mean = "constrained";
if isfield(fluid_problem.Data, "pressure_mean")
    pressure_mean = fluid_problem.Data.pressure_mean;
end
mustBeMember(pressure_mean, ["constrained", "pinned"])
mesh_fluid = fluid_mesh;
mesh_solid = solid_mesh;

if options.quiet
    print_message = @(s) disp('');
else
    print_message = @(s) disp(s);
end

%=========================== SOME ASSIGNMENTS  ===========================%
k = options.Order;
fluid_problem.Order = k;
solid_problem.Order = k;
dt = fluid_problem.Data.time_step;

if scheme == "newmark"
    assert(solid_problem.Data.time_scheme == "newmark", ...
        "scheme=""newmark"" requires the solid problem to be constructed with time_scheme=""newmark"".");
    beta_nm  = solid_problem.Data.newmark_beta;
    gamma_nm = solid_problem.Data.newmark_gamma;
    % Interface kinematics: eliminating a^{n+1} from the Newmark updates and
    % imposing v^{n+1} = û₁ gives
    %   û₂ = (β/γ)Δt û₁ + u⁰ + Δt(1-β/γ) v⁰ + Δt²[(1/2-β) - β(1-γ)/γ] a⁰
    cdt = (beta_nm/gamma_nm)*dt;                                  % û₁ coupling factor
    c_v0 = dt*(1 - beta_nm/gamma_nm);                             % v⁰ history factor
    c_a0 = dt^2*((1/2 - beta_nm) - beta_nm*(1 - gamma_nm)/gamma_nm); % a⁰ history factor
else
    cdt = dt;
end

% Degree-adaptive quadrature; n_quadpoints=8, n_quadpoints_path=9
% recovers the fixed rule
s_quad = k + 1;  tp_quad = k + 2;
if ~isempty(options.n_quadpoints),      s_quad  = options.n_quadpoints;      end
if ~isempty(options.n_quadpoints_path), tp_quad = options.n_quadpoints_path; end
[~, wS]  = lgwt(s_quad, 0, 1);
[~, wTP] = lgwt(tp_quad, 0, 1);
%% I. LOADING PARAMETERS
print_message("Loading parameters...")
%======================= LOAD EQUATION PARAMETERS ========================%
viscosity_fluid = fluid_problem.Data.viscosity;
mu_f = viscosity_fluid;
c1_1 = 1/(2*mu_f);
c2_1 = 1/(4*mu_f);
measure_stokes = fluid_problem.Data.domain_measure;
young_modulus = solid_problem.Data.young_modulus;
poisson_ratio = solid_problem.Data.poisson_ratio;
[lambda_s, mu_s] = LinearElasticityProblem.get_lame_parameters(young_modulus, poisson_ratio);
c1_2 = 1/(2*mu_s);
c2_2 = lambda_s/(2*mu_s*(2*lambda_s + 2*mu_s));
%=========================== USEFUL QUANTITIES ===========================%
d1 = k + 1;                             % dim P_k(dK)
% d2peers = 4*nchoosek(k + 2, 2) + min(1, k)*(k + 1);                 % dim P_k(K)

N_F_1 = size(mesh_fluid.faces, 2);        % # of faces of first mesh
N_DOF_1 = 2*d1*N_F_1 + size(mesh_fluid.elements, 2); % # of DOF belonging to 1st mesh

N_F_2 = size(mesh_solid.faces, 2);         % # of faces of second mesh
N_DOF_2 = 2*d1*N_F_2;                          % # of DOF belonging to 2nd mesh

A_DOF = N_DOF_1 + N_DOF_2 + 1;

print_message("Parameters loaded.")
%% II. INDIVIDUAL SYSTEMS
print_message('Generating individual systems...')
[localsolvers_inner, val_inner, rhs_inner] = fluid_problem.get_matrices(mesh_fluid);
[localsolvers_outer, val_outer, rhs_outer] = solid_problem.get_matrices(mesh_solid);
print_message('Individual global systems generated.')
%% III. TRANSMISSION CONDITIONS (OFFLINE)
print_message('Computing induced partition...')
[partition_fluid, partition_solid] = ...
    get_partitions(interface, ...
                   mesh_fluid, mesh_solid, ...
                   quadrature_patch = s_quad, ...
                   quadrature_path = tp_quad, ...
                   vertex_tol = options.vertex_tol, ...
                   vertex_tol_relative = options.vertex_tol_relative);
n_of_patches = size(partition_fluid.faces, 3); % or partition_solid
    
% the + 0*x term broadcasts spatially constant source handles (scalar
% return) to the quadrature-point array shape
evaluate_on_points = @(f, x) cellfun(@(f) reshape(f(x(1,:,:), x(2,:,:)) + 0*x(1,:,:), s_quad, 1, 1, []), f, 'UniformOutput', 0);

T1 = evaluate_on_points(source_neumann, partition_fluid.quadrature.interface);
T1 = cat(3, T1{:});

T2 = evaluate_on_points(source_dirichlet, partition_solid.quadrature.interface);
T2 = cat(3, T2{:});
%=========================== LOCAL SOLVERS ===============================%
A_fluid = cellfun(...
    @(A) A(:, :, partition_fluid.elements(:)), ...
    localsolvers_inner, 'UniformOutput', false);
[A1_fluid, A2_fluid, Af_fluid] = deal(A_fluid{:});
% The traction balance is the Cauchy traction plus the HDG stabilization,
% also for a Navier-Stokes fluid (no momentum is advected through the
% material interface).


A_solid = cellfun(...
    @(A) A(:, :, partition_solid.elements(:)), ...
    localsolvers_outer, 'UniformOutput', false);
[A1_solid, A2_solid, Af_solid] = deal(A_solid{:});

%======================= PRE-ALLOCATED ARRAYS ============================%
[val_T1_1, row_T1_1, col_T1_1] = deal(zeros(2*d1, 3*2*d1 + 1, n_of_patches));
[val_T1_2, row_T1_2, col_T1_2] = deal(zeros(2*d1, 3*2*d1, n_of_patches));

[val_T1_3, row_T1_3, col_T1_3] = deal(zeros(2*d1, 2*d1, n_of_patches));
[val_T1_4, row_T1_4, col_T1_4] = deal(zeros(2*d1, 2*d1, n_of_patches));
[val_T1_5, row_T1_5, col_T1_5] = deal(zeros(2*d1, 1, n_of_patches));

[rhs_T1, rowRHS_T1] = deal(zeros(2*d1, 1, n_of_patches));

% Traction-resultant functionals (3 per side), appended to the system
% only when the options are set
val_RES1 = zeros(3, 3*2*d1 + 1, n_of_patches);
val_RES2 = zeros(3, 3*2*d1, n_of_patches);
[rhs_RES1, rhs_RES2] = deal(zeros(3, 1, n_of_patches));

[val_T2_1, row_T2_1, col_T2_1] = deal(zeros(2*d1, 3*2*d1 + 1, n_of_patches));
[val_T2_2, row_T2_2, col_T2_2] = deal(zeros(2*d1, 3*2*d1, n_of_patches));

[val_T2_3, row_T2_3, col_T2_3] = deal(zeros(2*d1, 2*d1, n_of_patches));
[val_T2_4, row_T2_4, col_T2_4] = deal(zeros(2*d1, 2*d1, n_of_patches));
[rhs_T2, rowRHS_T2] = deal(zeros(2*d1, 1, n_of_patches));
    
[val_A_1, row_A_1, col_A_1] = deal(zeros(1, 3*2*d1 + 1, n_of_patches));
[rhs_A, rowRHS_A] = deal(zeros(1, 1, n_of_patches));

%% III-2. TRANSMISSION CONDITIONS (PARALLEL)
print_message("Computing transmission conditions (parallel)...")

[DOF_1, SC_DOF_1] = get_partition_dof(mesh_fluid, partition_fluid, d1);
SC_DOF_1 = [SC_DOF_1; partition_fluid.elements(:).' + 2*d1*N_F_1];
[DOF_2, SC_DOF_2] = get_partition_dof(mesh_solid, partition_solid, d1, N_DOF_1);
    
partitions = [partition_fluid; partition_solid];
% patch_normals = cat(4, partitions.normal);
patch_area = cat(3, partitions.area);
patch_tau = cat(4, partitions.tau);

quadratures = [partition_fluid.quadrature; partition_solid.quadrature];
patch_locals = cat(4, quadratures.local); 
patch_locals_2d = cat(4, quadratures.local_2d);
patch_interface_2d = cat(4, quadratures.interface_2d);
patch_path_2d = cat(5, quadratures.path_2d);
patch_normals = cat(4, quadratures.normals);
patch_iface_phys = cat(4, quadratures.interface);

m_fluid = partition_fluid.quadrature.interface - partition_fluid.quadrature.physical;
m_solid = partition_solid.quadrature.interface - partition_solid.quadrature.physical;
patch_tpaths = cat(4, m_fluid, m_solid);
    
G_fluid = get_affine_map(mesh_fluid, partition_fluid.elements, inverse=true);
G_solid = get_affine_map(mesh_solid, partition_solid.elements, inverse=true);

D_A = dubiner_basis(0, 0.0, 0.0);

% Signed sliver area element for the gap-sliver integral (alpha row),
% per options.patch_weight; see get_sliver_jacobian. Under
% pressure_mean = "pinned" the alpha row is eliminated, so the weight
% is not built.
if pressure_mean == "pinned" || ~options.sliver_correction
    if pressure_mean == "pinned" && ...
            (options.patch_weight ~= "rectangular" || ~options.sliver_correction)
        warning("FluidStructureCoupling:inertOption", ...
            "patch_weight/sliver_correction have no effect under " + ...
            "pressure_mean = ""pinned"": the alpha row is eliminated.")
    end
    sliver_jacobian = zeros(s_quad, tp_quad, n_of_patches);
else
    sliver_jacobian = interface.get_sliver_jacobian(mesh_fluid, partition_fluid, ...
        patch_weight = options.patch_weight, ...
        quadrature_patch = s_quad, quadrature_path = tp_quad);
end

% Interface history H⁰ on the skeleton trace DOFs of the interface faces:
% steady none, onestep u⁰, newmark u⁰ + Δt(1-β/γ) v⁰ + Δt²[(1/2-β)-β(1-γ)/γ] a⁰.
if scheme == "steady"
    h0 = zeros(2*d1, n_of_patches);
else
    h0 = solid_problem.Data.initial_displacement_skeleton;
    if scheme == "newmark"
        if isfield(solid_problem.Data, "initial_velocity_skeleton") && ...
                ~isequal(solid_problem.Data.initial_velocity_skeleton, 0.0)
            h0 = h0 + c_v0 * solid_problem.Data.initial_velocity_skeleton;
        end
        if isfield(solid_problem.Data, "initial_acceleration_skeleton") && ...
                ~isequal(solid_problem.Data.initial_acceleration_skeleton, 0.0)
            h0 = h0 + c_a0 * solid_problem.Data.initial_acceleration_skeleton;
        end
    end
    if isscalar(h0)
        % scalar-0 placeholder history (at-rest start) expands to zeros
        h0 = h0 + zeros(2*d1, n_of_patches);
    else
        h0 = h0(:, partition_solid.faces(:));
    end
end
% Transport the history onto Sigma, so that both sides of the T2 row are
% measured at the same place (see the history_extension option).
if ~isempty(options.history_extension)
    assert(isequal(size(options.history_extension), [2*d1, n_of_patches]), ...
        "FluidStructureCoupling:solve:badHistoryExtension", ...
        "history_extension must be %d x %d, got %s.", ...
        2*d1, n_of_patches, mat2str(size(options.history_extension)));
    h0 = h0 + options.history_extension;
end
parfor e = 1:n_of_patches
    DOF1 = DOF_1(:, e);
    DOF2 = DOF_2(:, e);

    rowRHS_T1(:,:,e) = DOF1;
    rowRHS_T2(:,:,e) = DOF2;
    rowRHS_A(:, :, e) = A_DOF;

    faceDOF1 = SC_DOF_1(:, e);
    faceDOF2 = SC_DOF_2(:, e);

    LS_inner = A1_fluid(:,:,e) \ A2_fluid(:,:,e);
    LS_RHS_inner = A1_fluid(:,:,e) \ Af_fluid(:,:,e);
    G_inner = G_fluid(:, :, e);

    LS_outer = A1_solid(:,:,e) \ A2_solid(:,:,e);
    LS_RHS_outer = A1_solid(:,:,e) \ Af_solid(:,:,e);
    G_outer = G_solid(:, :, e);

    s_patch = patch_locals(:, :, e, :);
    st_boundary = patch_locals_2d(:, :, e, :);
    st_path = patch_path_2d(:, :, e, :, :);
    st_interface = patch_interface_2d(:, :, e, :);
    m_tpath = patch_tpaths(:, :, e, :);
    x_iface = patch_iface_phys(:, :, e, :);
    n_interface = reshape(patch_normals(:, :, e, :), 2, 1, [], 2);

    nF = -n_interface(:, :, :, 1);
    nS = n_interface(:, :, :, 2);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% POLYNOMIALS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %============================ COLLOCATION ============================%
    %---------------------------- LEGENDRE (1D) --------------------------% 
    % On triangle side facing interface
    s_fluid = s_patch(:, :, :, 1);
    l1 = legendre_basis(k, s_fluid); % p_gauss x d1
    L1 = vector_legendre_basis(k, s_fluid); % p_gauss x d1
    s_solid = s_patch(:, :, :, 2);
    l2 = legendre_basis(k, s_solid); % p_gauss x d1
    L2 = vector_legendre_basis(k, s_solid); % p_gauss x d1

    %----------------------------- DUBINER (2D) --------------------------%
    % On triangle side facing interface
    s_bdy_inner = st_boundary(1, :, :, 1);
    t_bdy_inner = st_boundary(2, :, :, 1);
    D1 = vector_dubiner_basis(k, s_bdy_inner, t_bdy_inner);

    s_bdy_outer = st_boundary(1, :, :, 2);
    t_bdy_outer = st_boundary(2, :, :, 2);
    D2 = vector_dubiner_basis(k, s_bdy_outer, t_bdy_outer);

    % On connecting segment (EXTRAPOLATED)
    s_conn_inner = reshape(st_path(1, :, :, :, 1), s_quad, tp_quad).';
    t_conn_inner = reshape(st_path(2, :, :, :, 1), s_quad, tp_quad).';
    TP_P1 = PEERS_basis(k, s_conn_inner(:), t_conn_inner(:), G_inner);
    TP_P1 = reshape(TP_P1, [tp_quad s_quad size(TP_P1, 2:4)]);
    TP_A1 = asym_tensor_dubiner_basis(k, s_conn_inner(:), t_conn_inner(:));
    TP_A1 = reshape(TP_A1, [tp_quad s_quad size(TP_A1, 2:4)]);

    s_conn_outer = reshape(st_path(1, :, :, :, 2), s_quad, tp_quad).';
    t_conn_outer = reshape(st_path(2, :, :, :, 2), s_quad, tp_quad).';
    TP_P2 = PEERS_basis(k, s_conn_outer(:), t_conn_outer(:), G_outer);
    TP_P2 = reshape(TP_P2, [tp_quad s_quad size(TP_P2, 2:4)]);
    TP_A2 = asym_tensor_dubiner_basis(k, s_conn_outer(:), t_conn_outer(:));
    TP_A2 = reshape(TP_A2, [tp_quad s_quad size(TP_A2, 2:4)]);

    % On interface (EXTRAPOLATED)
    s_iface_inner = st_interface(1, :, :, 1);
    t_iface_inner = st_interface(2, :, :, 1);
    E_P1 = PEERS_basis(k, s_iface_inner, t_iface_inner, G_inner);
    E_P1 = permute(E_P1, [3 4 1 2]);
    E_P1_N = permute(pagemtimes(E_P1, nF), [3 4 1 2]);
    D_A_N = permute(D_A * nF, [3 2 1]);
    E_A1 = asym_tensor_dubiner_basis(k, s_iface_inner, t_iface_inner);
    E_A1 = permute(E_A1, [3 4 1 2]);
    E_A1_N = permute(pagemtimes(E_A1, nF), [3 4 1 2]);

    s_iface_outer = st_interface(1, :, :, 2);
    t_iface_outer = st_interface(2, :, :, 2);
    E_P2 = PEERS_basis(k, s_iface_outer, t_iface_outer, G_outer);
    E_P2 = permute(E_P2, [3 4 1 2]);
    E_P2_N = permute(pagemtimes(E_P2, nS), [3 4 1 2]);
    E_A2 = asym_tensor_dubiner_basis(k, s_iface_outer, t_iface_outer);
    E_A2 = permute(E_A2, [3 4 1 2]);
    E_A2_N = permute(pagemtimes(E_A2, nS), [3 4 1 2]);
        

    % Integration along connecting segment (tangent vector is NOT normalized = includes length of segment)
    mF = reshape(m_tpath(:, :, :, 1), 2, 1, []);
    TP_P1 = permute(TP_P1, [4 5 2 1 3]);
    TP_P1_tr = TP_P1(1, 1, :, :, :) + TP_P1(2, 2, :, :, :);
    TP_P1_tr_Id = pagemtimes(TP_P1_tr, eye(2));
    TP_P1_grad_M = permute(pagemtimes(c1_1*TP_P1 - c2_1*TP_P1_tr_Id, mF), [3 4 5 1 2]);
    TP_INT_P1 = tensorprod(TP_P1_grad_M, wTP, 2, 1);
    TP_A1 = permute(TP_A1, [4 5 2 1 3]);
    TP_A1_M = permute(pagemtimes(TP_A1, mF), [3 4 5 1 2]);
    TP_INT_A1 = tensorprod(TP_A1_M, wTP, 2, 1);

    mS = reshape(m_tpath(:, :, :, 2), 2, 1, []);
    TP_P2 = permute(TP_P2, [4 5 2 1 3]);
    TP_P2_tr_Id = pagemtimes(TP_P2(1, 1, :, :, :) + TP_P2(2, 2, :, :, :), eye(2));
    TP_P2_grad_M = permute(pagemtimes(c1_2*TP_P2 - c2_2*TP_P2_tr_Id, mS), [3 4 5 1 2]);
    TP_INT_P2 = tensorprod(TP_P2_grad_M, wTP, 2, 1);
    TP_A2 = permute(TP_A2, [4 5 2 1 3]);
    TP_A2_M = permute(pagemtimes(TP_A2, mS), [3 4 5 1 2]);
    TP_INT_A2 = tensorprod(TP_A2_M, wTP, 2, 1);
     
    % Explicit reshape, not squeeze: at k = 0 the adaptive rule gives
    % s_quad = 1 and squeeze would drop the patch-quadrature axis, so the
    % product with W_sliver (s_quad x tp_quad) below silently mismatches.
    TP_P1_tr = reshape(TP_P1_tr, s_quad, tp_quad, []);
    W_sliver = (wS(:) .* sliver_jacobian(:, :, e)) .* wTP(:)';
    EP_INT_INT_P1 = reshape(sum(TP_P1_tr .* W_sliver, [1 2]), 1, []);
    
    tau = patch_tau(:, :, e, :);
    area = patch_area(:, e, :);
    %% <E(σ₁)n₁−μτu₁+μτû₁+E(σ₂)n₂−τu₂+τû₂,μ₁>_I₁ = <T₁ o φ,μ₁>_I₁
    tau_1 = tau(:, :, :, 1);
    tau_2 = tau(:, :, :, 2);
    patch_len = area(:, :, 1);
    integrate = @(D, L) patch_len * reshape(permute(pagemtimes(wS.*L, 'transpose', D, 'none'), [1 3 2]), 2*d1, []);

    sigma_term = integrate(E_P1_N, l1);
    u_term = -tau_1*mu_f*integrate(D1, l1);
    gamma_term = integrate(0*E_A1_N, l1);
    xi_term = integrate(0*D1(:, 1, :), l1);
    S1_inner = [sigma_term, u_term, gamma_term, xi_term];
    val_T1_1(:,:,e) = S1_inner*LS_inner;
    row_T1_1(:,:,e) = kron(ones(1,3*2*d1 + 1), DOF1);
    col_T1_1(:,:,e) = kron(ones(2*d1,1), faceDOF1');
    rhs_T1(:,:,e) = rhs_T1(:,:,e) - S1_inner*LS_RHS_inner;

    sigma_term = integrate(E_P2_N, l1);
    u_term = -tau_2*integrate(D2, l1);
    gamma_term = integrate(0*E_A2_N, l1);
    S1_outer = [sigma_term, u_term, gamma_term];
    val_T1_2(:,:,e) = S1_outer*LS_outer;
    row_T1_2(:,:,e) = kron(ones(1,3*2*d1), DOF1);
    col_T1_2(:,:,e) = kron(ones(2*d1,1), faceDOF2'); 
    rhs_T1(:,:,e) = rhs_T1(:,:,e) - S1_outer*LS_RHS_outer;
    
    val_T1_3(:,:,e) = tau_1*mu_f*integrate(L1, l1);
    row_T1_3(:,:,e) = kron(ones(1,2*d1), DOF1);
    col_T1_3(:,:,e) = kron(ones(2*d1,1), DOF1');

    val_T1_4(:,:,e) = tau_2*integrate(L2, l1);
    row_T1_4(:,:,e) = kron(ones(1,2*d1), DOF1);
    col_T1_4(:,:,e) = kron(ones(2*d1,1), DOF2');

    val_T1_5(:,:,e) = integrate(D_A_N, l1);
    row_T1_5(:,:,e) = DOF1;
    col_T1_5(:,:,e) = A_DOF;

    rhs_T1(:,:,e) = rhs_T1(:,:,e) + integrate(T1(:, :, :, e), l1);

    % Resultant functionals of the extrapolated traction, pushed through
    % the local solvers: net force x/y and torque about (.5, .5)
    rF = x_iface(:, :, 1, 1) - 0.5;
    RES_sig = [area(:, :, 1) * (wS' * E_P1_N(:, :, 1));
               area(:, :, 1) * (wS' * E_P1_N(:, :, 2));
               area(:, :, 1) * ((wS .* rF(1, :).')' * E_P1_N(:, :, 2) ...
                              - (wS .* rF(2, :).')' * E_P1_N(:, :, 1))];
    R_inner = [RES_sig, zeros(3, size(D1, 2) + size(E_A1_N, 2) + 1)];
    val_RES1(:,:,e) = R_inner*LS_inner;
    rhs_RES1(:,:,e) = R_inner*LS_RHS_inner;

    rS = x_iface(:, :, 1, 2) - 0.5;
    RES_sig = [area(:, :, 2) * (wS' * E_P2_N(:, :, 1));
               area(:, :, 2) * (wS' * E_P2_N(:, :, 2));
               area(:, :, 2) * ((wS .* rS(1, :).')' * E_P2_N(:, :, 2) ...
                              - (wS .* rS(2, :).')' * E_P2_N(:, :, 1))];
    R_outer = [RES_sig, zeros(3, size(D2, 2) + size(E_A2_N, 2))];
    val_RES2(:,:,e) = R_outer*LS_outer;
    rhs_RES2(:,:,e) = R_outer*LS_RHS_outer;

    %%
    sigma_term = D_A'*EP_INT_INT_P1;
    u_term = zeros(1, size(D1, 2));
    gamma_term = zeros(1, size(E_A1_N, 2));
    xi_term = 0;
    SA_inner = [sigma_term, u_term, gamma_term, xi_term];
    val_A_1(:,:,e) = SA_inner*LS_inner;
    row_A_1(:,:,e) = A_DOF;
    col_A_1(:,:,e) = faceDOF1';
    rhs_A(:,:,e) = rhs_A(:,:,e) - SA_inner*LS_RHS_inner;
    

    %% <dt*û₁ - û₂, μ₂>_I₂ = <T₂ o φ,μ₂>_I₂
    patch_len = area(:, :, 2);
    integrate = @(D, L) patch_len * reshape(permute(pagemtimes(wS.*L, 'transpose', D, 'none'), [1 3 2]), 2*d1, []);

    sigma_term = integrate(TP_INT_P1, l2);
    u_term = integrate(0*D1, l2);
    gamma_term = integrate(TP_INT_A1, l2);
    xi_term = integrate(0*D1(:, 1, :), l2);
    S2_inner = cdt*[sigma_term, u_term, gamma_term, xi_term];

    val_T2_1(:,:,e) = S2_inner*LS_inner;
    row_T2_1(:,:,e) = kron(ones(1,3*2*d1 + 1),DOF2);
    col_T2_1(:,:,e) = kron(ones(2*d1,1),faceDOF1');
    rhs_T2(:,:,e) = rhs_T2(:,:,e) - S2_inner*LS_RHS_inner;

    sigma_term = -integrate(TP_INT_P2, l2);
    u_term = -integrate(0*D2, l2);
    gamma_term = -integrate(TP_INT_A2, l2);
    S2_outer = [sigma_term, u_term, gamma_term];
    val_T2_2(:,:,e) = S2_outer*LS_outer;
    row_T2_2(:,:,e) = kron(ones(1,3*2*d1), DOF2);
    col_T2_2(:,:,e) = kron(ones(2*d1,1), faceDOF2');
    rhs_T2(:,:,e) = rhs_T2(:,:,e) - S2_outer*LS_RHS_outer;
    
    val_T2_3(:,:,e) = cdt*integrate(L1, l2);
    row_T2_3(:,:,e) = kron(ones(1,2*d1),DOF2);
    col_T2_3(:,:,e) = kron(ones(2*d1,1),DOF1');

    val_T2_4(:,:,e) = -integrate(L2, l2);
    row_T2_4(:,:,e) = kron(ones(1,2*d1),DOF2);
    col_T2_4(:,:,e) = kron(ones(2*d1,1),DOF2');
   
    rhs_T2(:,:,e) = rhs_T2(:,:,e) - integrate(L2, l2)*h0(:, e) + cdt*integrate(T2(:, :, :, e), l2);
end
val_A_2 = 2 * measure_stokes * (D_A')*D_A;
row_A_2 = A_DOF;
col_A_2 = A_DOF;
print_message('Transmission conditions ready.')
%% IV. ASSEMBLING GLOBAL SYSTEM
print_message('Assembling global system...')
[row_inner, col_inner, rowRHS_inner] = fluid_problem.get_assemble_indices(mesh_fluid, k);
[row_outer, col_outer, rowRHS_outer] = solid_problem.get_assemble_indices(mesh_solid, k);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% DOF RE-NUMBERING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Delete entries that belong to interface DOF from original systems
val_inner(ismember(row_inner, DOF_1)) = [];
col_inner(ismember(row_inner, DOF_1)) = [];
row_inner(ismember(row_inner, DOF_1)) = [];
rhs_inner(ismember(rowRHS_inner, DOF_1)) = [];
rowRHS_inner(ismember(rowRHS_inner, DOF_1)) = [];

val_outer(ismember(row_outer, DOF_2 - N_DOF_1)) = [];
col_outer(ismember(row_outer, DOF_2 - N_DOF_1)) = [];
row_outer(ismember(row_outer, DOF_2 - N_DOF_1)) = [];
rhs_outer(ismember(rowRHS_outer, DOF_2 - N_DOF_1)) = [];
rowRHS_outer(ismember(rowRHS_outer, DOF_2 - N_DOF_1)) = [];

%%% Shift entries belonging to second mesh
row_outer = row_outer + N_DOF_1;
col_outer = col_outer + N_DOF_1;
rowRHS_outer = rowRHS_outer + N_DOF_1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% SYSTEM ASSEMBLY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
row = [row_inner(:); row_T1_1(:); row_T1_2(:); row_T1_3(:); row_T1_4(:); 
       row_outer(:); row_T2_1(:); row_T2_2(:); row_T2_3(:); row_T2_4(:);
                row_T1_5(:); row_A_1(:);  row_A_2(:)];
col = [col_inner(:); col_T1_1(:); col_T1_2(:); col_T1_3(:); col_T1_4(:);
       col_outer(:); col_T2_1(:); col_T2_2(:); col_T2_3(:); col_T2_4(:);
                 col_T1_5(:); col_A_1(:);  col_A_2(:)];
val = [val_inner(:); val_T1_1(:); val_T1_2(:); val_T1_3(:); val_T1_4(:);
       val_outer(:); val_T2_1(:); val_T2_2(:); val_T2_3(:); val_T2_4(:); 
                 val_T1_5(:); val_A_1(:);  val_A_2(:)];
H = sparse(row, col,val);

F = [rhs_inner(:); rhs_T1(:);
     rhs_outer(:); rhs_T2(:);
              rhs_A(:)];

RowRHS = [rowRHS_inner(:); rowRHS_T1(:);
          rowRHS_outer(:); rowRHS_T2(:);
                           rowRHS_A(:)];
RHS = accumarray(RowRHS,F);

% \sum_K \rho_K = 0
RHS(end+1) = 0;
H(end+1, :) = [zeros(1, 2*d1*N_F_1), mesh_fluid.volume, zeros(1, N_DOF_2), 0];
H(:, end+1) = [zeros(1, 2*d1*N_F_1), mesh_fluid.volume, zeros(1, N_DOF_2), 0, 0]';

% EXPERIMENTAL circulation anchor (see the option's doc)
if ~isempty(options.circulation)
    ifaces_F = mesh_fluid.intfaces;
    [tA, wA] = lgwt(s_quad, 0, 1);
    int_basis = (wA' * legendre_basis(k, tA)).';    % \int_0^1 l_i, per mode
    a = zeros(1, size(H, 2));
    for f_if = reshape(ifaces_F, 1, [])
        va = mesh_fluid.coordinates(1:2, mesh_fluid.faces(1, f_if));
        vb = mesh_fluid.coordinates(1:2, mesh_fluid.faces(2, f_if));
        tau = (vb - va) / norm(vb - va);
        mid = (va + vb) / 2 - [.5; .5];             % about the loop center
        if tau(1)*mid(2) - tau(2)*mid(1) > 0        % orient counterclockwise
            tau = -tau;
        end
        len = mesh_fluid.area(f_if);
        a(d1*(f_if - 1) + (1:d1))            = len * tau(1) * int_basis;
        a(d1*N_F_1 + d1*(f_if - 1) + (1:d1)) = len * tau(2) * int_basis;
    end
    H(end+1, :) = a;
    H(:, end+1) = [a, 0]';
    RHS(end+1) = options.circulation;
end

% EXPERIMENTAL solid rigid-mode anchor (see the option's doc)
if ~isempty(options.solid_rigid_modes)
    assert(numel(options.solid_rigid_modes) == 3, ...
        "solid_rigid_modes must be [t_x; t_y; rot].");
    [tR, wR] = lgwt(s_quad, 0, 1);
    lR = legendre_basis(k, tR);                     % s x d1
    int_basis_R = (wR' * lR).';
    A3 = zeros(3, size(H, 2));
    for f_if = reshape(mesh_solid.intfaces, 1, [])
        va = mesh_solid.coordinates(1:2, mesh_solid.faces(1, f_if));
        vb = mesh_solid.coordinates(1:2, mesh_solid.faces(2, f_if));
        len = mesh_solid.area(f_if);
        dx = N_DOF_1 + d1*(f_if - 1) + (1:d1);      % solid u_hat, x comp
        dy = dx + d1*N_F_2;                         % solid u_hat, y comp
        % translations: mean coefficient per component
        A3(1, dx) = len * int_basis_R;
        A3(2, dy) = len * int_basis_R;
        % rotation weight along the face: w = [-(y - .5); x - .5]
        X = va(1) + tR * (vb(1) - va(1));
        Y = va(2) + tR * (vb(2) - va(2));
        A3(3, dx) = len * ((wR .* (-(Y - .5)))' * lR);
        A3(3, dy) = len * ((wR .* (X - .5))' * lR);
    end
    for j = 1:3
        H(end+1, :) = [A3(j, :), zeros(1, j - 1)]; %#ok<AGROW>
        H(:, end+1) = [A3(j, :), zeros(1, j - 1), 0]'; %#ok<AGROW>
        RHS(end+1) = options.solid_rigid_modes(j); %#ok<AGROW>
    end
end

% EXPERIMENTAL traction-resultant anchor (see the option's doc)
resultant_sides = { ...
    options.traction_resultants_fluid, val_RES1, rhs_RES1, SC_DOF_1;
    options.traction_resultants_solid, val_RES2, rhs_RES2, SC_DOF_2};
for side = 1:2
    [targets, vals, rhss, sc_dof] = resultant_sides{side, :};
    if isempty(targets), continue; end
    assert(numel(targets) == 3, ...
        "traction_resultants must be [F_x; F_y; M].");
    B = zeros(3, size(H, 2));
    c0 = zeros(3, 1);
    for e = 1:n_of_patches
        cols = sc_dof(:, e);
        B(:, cols) = B(:, cols) + vals(:, :, e);
        c0 = c0 + rhss(:, :, e);
    end
    for j = 1:3
        H(end+1, :) = [B(j, :), zeros(1, j - 1)]; %#ok<AGROW>
        H(:, end+1) = [B(j, :), zeros(1, j - 1), 0]'; %#ok<AGROW>
        RHS(end+1) = targets(j) - c0(j); %#ok<AGROW>
    end
end

%% V. BOUNDARY CONDITIONS
print_message('Implementing boundary conditions...')
[x1_hat, dirfaces_1, g1_hat] = fluid_problem.get_boundary_conditions(mesh_fluid);
[x2_hat, dirfaces_2, g2_hat] = solid_problem.get_boundary_conditions(mesh_solid);
dirfaces = [dirfaces_1(:); dirfaces_2(:) + N_DOF_1];
if pressure_mean == "pinned"
    % Eliminate p0 and the mean-pressure constraint (p0_h = 0)
    dirfaces = [dirfaces; A_DOF; A_DOF + 1];
end

x_hat = [x1_hat; x2_hat; 0; 0];
if ~isempty(options.circulation)
    x_hat = [x_hat; 0];   % the circulation multiplier is always free
end
if ~isempty(options.solid_rigid_modes)
    x_hat = [x_hat; 0; 0; 0];   % rigid-mode multipliers are always free
end
if ~isempty(options.traction_resultants_fluid)
    x_hat = [x_hat; 0; 0; 0];   % resultant multipliers are always free
end
if ~isempty(options.traction_resultants_solid)
    x_hat = [x_hat; 0; 0; 0];
end
RHS(1:N_DOF_1) = RHS(1:N_DOF_1) + g1_hat;
RHS(N_DOF_1 + (1:N_DOF_2)) = RHS(N_DOF_1 + (1:N_DOF_2)) + g2_hat;
RHS = RHS - H*x_hat;

%% VI-1. SOLVING GLOBAL SYSTEM
print_message('Solving linear systems...')
freefaces = setdiff((1:size(H,1)).',dirfaces);
H_   = H(freefaces, freefaces);
RHS_ = RHS(freefaces);
x_hat(freefaces) = mldivide(H_, RHS_);

if nargout > 2
    global_system = struct( ...
        "matrix",   H, ...
        "rhs",      RHS, ...
        "dof",      freefaces, ...
        "solution", x_hat, ...
        "residue",  norm(H_*x_hat(freefaces) - RHS_));
end

x1_hat = x_hat(1:N_DOF_1);
x2_hat = x_hat(N_DOF_1 + (1:N_DOF_2));
p0_h   = x_hat(A_DOF);
print_message('Global system solved.')
%% VI-2. SOLVING LOCAL SYSTEMS
print_message('Recovering local solutions...')
[fluid_sols, fluid_hybrid_sols] = fluid_problem.get_local_solutions(x1_hat, localsolvers_inner, mesh_fluid);
[solid_sols, solid_hybrid_sols] = solid_problem.get_local_solutions(x2_hat, localsolvers_outer, mesh_solid);
print_message('Local solutions computed.')

%% VII. OUTPUT
% pressure post-process: add mean-pressure correction to the constant
% coefficient of the diagonal components L11 (cell 1) and L22 (cell 4)
fluid_sols{1}{1}(1, :) = fluid_sols{1}{1}(1, :) + p0_h;
fluid_sols{1}{4}(1, :) = fluid_sols{1}{4}(1, :) + p0_h;

fluid_solution = FESolution(fluid_problem, mesh_fluid);
fluid_solution = fluid_solution.assign_unknowns(fluid_sols, fluid_hybrid_sols);

solid_solution = FESolution(solid_problem, mesh_solid);
solid_solution = solid_solution.assign_unknowns(solid_sols, solid_hybrid_sols);

end

function [DOF, SC_DOF] = get_partition_dof(mesh, partition, d1, shift)
    arguments
        mesh
        partition
        d1
        shift = 0
    end
    N_PATCHES = length(partition.faces);
    N_FACES = size(mesh.faces, 2);

    get_1d_dof = @(e) d1*(e(:).' - 1) + (1:d1)';
    get_DOF = @(e) reshape(get_1d_dof(e) + permute([0, d1*N_FACES] + shift, ...
                            [1 3 2]), [], 1);

    DOF = cell2mat(cellfun(get_DOF, num2cell(partition.faces), ...
                           'UniformOutput', 0));

    partition_FBE = mesh.boundaries(:, partition.elements);
    SC_DOF = cell2mat(cellfun(get_DOF, ...
                              mat2cell(partition_FBE, 3, ones(N_PATCHES, 1)), ...
                              'UniformOutput', 0));
    % SC_DOF = [SC_DOF; partition.elements(:).' + 2*d1*N_FACES];
end