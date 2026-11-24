function [stokes_solution, darcy_solution, global_system] = solve(...
    stokes_problem, darcy_problem, upper_mesh, lower_mesh, ...
    interface, wbjs, source_neumann, source_dirichlet, ...
    options)
%SOLVE  Steady Stokes-Darcy coupling (Beavers-Joseph-Saffman) on dissimilar meshes by TPM.
%
%   [stokes_solution, darcy_solution, global_system] = solve( ...
%       stokes_problem, darcy_problem, upper_mesh, lower_mesh, interface, ...
%       wbjs, source_neumann, source_dirichlet, Order = ..., ...
%       n_quadpoints = ..., n_quadpoints_path = ..., quiet = ...)
%
%   solves the Stokes problem on upper_mesh and the Darcy problem on
%   lower_mesh, coupled across the AbstractCurve interface by the
%   transfer path method. The two meshes need not match at the interface.
%
%   Positional arguments:
%     stokes_problem    StokesProblem
%     darcy_problem     DarcyProblem
%     upper_mesh        SimplicialGrid of the Stokes domain
%     lower_mesh        SimplicialGrid of the Darcy domain
%     interface         AbstractCurve
%     wbjs              BJS friction coefficient (0.0)
%     source_neumann    cell of two function handles of (x, y), the
%                       components of the datum of the traction
%                       transmission condition, evaluated on the
%                       interface ({@(x,y) 0*x.*y, @(x,y) 0*x.*y})
%     source_dirichlet  function handle of (x, y), the datum of the
%                       normal-velocity transmission condition, evaluated
%                       on the interface (@(x,y) 0*x.*y)
%
%   Name-value options:
%     Order              polynomial degree k, set on both problems
%                        (default: stokes_problem.Order)
%     n_quadpoints       quadrature points per patch on the interface
%                        (default: k + 1)
%     n_quadpoints_path  quadrature points per transfer path
%                        (default: k + 1)
%     quiet              suppress progress messages (false)
%
%   Outputs:
%     stokes_solution, darcy_solution  FESolution of each side
%     global_system                    struct with fields matrix, rhs, dof,
%                                      solution and residue (computed only
%                                      when requested)
%
%See also: InterfaceProblem
arguments
    stokes_problem  StokesProblem
    darcy_problem   DarcyProblem
    upper_mesh      SimplicialGrid
    lower_mesh      SimplicialGrid
    interface       AbstractCurve
    wbjs            double                  = 0.0
    source_neumann  cell                    = {@(x,y) 0*x.*y, @(x,y) 0*x.*y}
    source_dirichlet function_handle        = @(x,y) 0*x.*y
    options.Order double {mustBeScalarOrEmpty} = stokes_problem.Order
    options.n_quadpoints        = []   % defaults to k+1 below
    options.n_quadpoints_path   = []   % defaults to k+1 below
    options.quiet               = false
end

if options.quiet
    print_message = @(s) disp('');
else
    print_message = @(s) disp(s);
end

k = options.Order;
stokes_problem.Order = k;
darcy_problem.Order  = k;

assert(stokes_problem.Order == darcy_problem.Order, ...
    "PDE problems must have equal polynomial degree!")

%=========================== SOME ASSIGNMENTS  ===========================%
mesh_stokes = upper_mesh;
mesh_darcy  = lower_mesh;

s_quad  = k + 1; if ~isempty(options.n_quadpoints),      s_quad  = options.n_quadpoints;      end
tp_quad = k + 1; if ~isempty(options.n_quadpoints_path), tp_quad = options.n_quadpoints_path; end
[~, wS]  = lgwt(s_quad,  0, 1);
[~, wTP] = lgwt(tp_quad, 0, 1);

%% I. LOADING PARAMETERS
print_message("Loading parameters...")
%======================= LOAD EQUATION PARAMETERS ========================%
viscosity_coefficient = stokes_problem.Data.viscosity;
diffusion_tensor = darcy_problem.Data.diffusion;
measure_stokes = stokes_problem.Data.domain_measure;

%=========================== USEFUL QUANTITIES ===========================%
d1 = k + 1;                             % dim P_k(dK)
d2 = nchoosek(k + 2, 2);                 % dim P_k(K)

N_F_1 = size(mesh_stokes.faces, 2);        % # of faces of first mesh
N_DOF_1 = 2*d1*N_F_1 + size(mesh_stokes.elements, 2); % # of DOF belonging to 1st mesh

N_F_2 = size(mesh_darcy.faces, 2);         % # of faces of second mesh
N_DOF_2 = d1*N_F_2;                          % # of DOF belonging to 2nd mesh

A_DOF = N_DOF_1 + N_DOF_2 + 1;             % DOF of extrapolated pressure

print_message("Parameters loaded.")
%% II. INDIVIDUAL SYSTEMS
print_message('Generating individual systems...')
[localsolvers_inner, val_inner, rhs_inner] = stokes_problem.get_matrices(mesh_stokes);
[localsolvers_outer, val_outer, rhs_outer] = darcy_problem.get_matrices(mesh_darcy);
print_message('Individual global systems generated.')
%% III. TRANSMISSION CONDITIONS (OFFLINE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% INTERFACE PARTITION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
print_message('Computing induced partition...')
[partition_stokes, partition_darcy] = ...
    get_partitions(interface, mesh_stokes, mesh_darcy, ...
                        quadrature_patch = s_quad, quadrature_path = tp_quad);

n_of_patches = size(partition_stokes.faces, 3); % or partition_darcy

%========= EVALUATE SOURCES AND PHYSICAL FUNCTIONS ON INTERFACE ==========%
evaluate_on_points = @(f, x) cellfun(@(f) reshape(f(x(1,:,:), x(2,:,:)), s_quad, 1, 1, []), f, 'UniformOutput', 0);

T1 = evaluate_on_points(source_neumann, partition_stokes.quadrature.interface);
T1 = cat(3, T1{:});

T2 = reshape(source_dirichlet(partition_darcy.quadrature.interface(1,:,:), ...
             partition_darcy.quadrature.interface(2,:,:)), s_quad, 1, []);

nu = viscosity_coefficient;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% ARRAY RESHAPING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%=========================== LOCAL SOLVERS ===============================%
A_stokes = cellfun(...
    @(A) A(:, :, partition_stokes.elements(:)), ...
    localsolvers_inner, 'UniformOutput', false);
[A1_stokes, A2_stokes, Af_stokes] = deal(A_stokes{:});

A_darcy = cellfun(...
    @(A) A(:, :, partition_darcy.elements(:)), ...
    localsolvers_outer, 'UniformOutput', false);
[A1_darcy, A2_darcy, Af_darcy] = deal(A_darcy{:});

%======================= PRE-ALLOCATED ARRAYS ============================%
% Matrices for "Neumann" transfer
% E(νL₁ - p₁I)n₁ - ντu₁
[val_T1_1, row_T1_1, col_T1_1] = deal(zeros(2*d1, 3*2*d1 + 1, n_of_patches));
% (∫κ⁻¹Eu₂·m)n₂
[val_T1_2, row_T1_2, col_T1_2] = deal(zeros(2*d1, 3*d1, n_of_patches));
% ντû₁
[val_T1_3, row_T1_3, col_T1_3] = deal(zeros(2*d1, 2*d1, n_of_patches));
% -p̂₂n₂
[val_T1_4, row_T1_4, col_T1_4] = deal(zeros(2*d1, d1, n_of_patches));
% -αn₁
[val_T1_5, row_T1_5, col_T1_5] = deal(zeros(2*d1, 1, n_of_patches));

% = RHS
[rhs_T1, rowRHS_T1] = deal(zeros(2*d1, 1, n_of_patches));

% Matrices for "Dirichlet" transfer
% ∫EL₁m·n₁
[val_T2_1, row_T2_1, col_T2_1] = deal(zeros(d1, 3*2*d1 + 1, n_of_patches));
% Eu₂·n₂ + τp₂
[val_T2_2, row_T2_2, col_T2_2] = deal(zeros(d1, 3*d1, n_of_patches));
% û₁·n₁
[val_T2_3, row_T2_3, col_T2_3] = deal(zeros(d1, 2*d1, n_of_patches));
% -τp̂₂
[val_T2_4, row_T2_4, col_T2_4] = deal(zeros(d1, d1, n_of_patches));
% = RHS
[rhs_T2, rowRHS_T2] = deal(zeros(d1, 1, n_of_patches));

% Matrices for extrapolated pressure
% ∫p₁
[val_A_1, row_A_1, col_A_1] = deal(zeros(1, 3*2*d1 + 1,n_of_patches));
% α
[val_A_2, row_A_2, col_A_2] = deal(zeros(1, 1, n_of_patches));
% = RHS
[rhs_A, rowRHS_A] = deal(zeros(1, 1, n_of_patches));

%% III-2. TRANSMISSION CONDITIONS (PARALLEL)
print_message("Computing transmission conditions (parallel)...")
[DOF_1, SC_DOF_1] = get_partition_dof(mesh_stokes, partition_stokes, d1, 2);
SC_DOF_1 = [SC_DOF_1; partition_stokes.elements(:).' + 2*d1*N_F_1];
[DOF_2, SC_DOF_2] = get_partition_dof(mesh_darcy, partition_darcy, d1, 1, N_DOF_1);
    
partitions = [partition_stokes; partition_darcy];
patch_area = cat(3, partitions.area);
patch_tau = cat(4, partitions.tau);

quadratures = [partition_stokes.quadrature; partition_darcy.quadrature];
patch_locals = cat(4, quadratures.local); 
patch_locals_2d = cat(4, quadratures.local_2d);
patch_interface_2d = cat(4, quadratures.interface_2d);
patch_path_2d = cat(5, quadratures.path_2d);   
patch_normals = cat(4, quadratures.normals);

m_stokes = partition_stokes.quadrature.interface - partition_stokes.quadrature.physical;
m_darcy = partition_darcy.quadrature.interface - partition_darcy.quadrature.physical;
patch_tpaths = cat(4, m_stokes, m_darcy);
    
% Integration helpers (defined once, used in each patch)
matfun = @(f, C) cell2mat(cellfun(f, C, 'UniformOutput', false));
wP   = wS;                             % 1D patch quadrature weights
wCS  = repmat(wTP(:), s_quad, 1);      % (tp_quad*s_quad,1) path weights
wI   = kron(1:s_quad, ones(1, tp_quad)); % index expansion: patch→path ordering

for e = 1:n_of_patches
    %% DOF numbering
    DOF1 = DOF_1(:, e);
    DOF2 = DOF_2(:, e);
    rowRHS_T1(:,:,e) = DOF1;
    rowRHS_T2(:,:,e) = DOF2;
    rowRHS_A(:, :, e) = A_DOF;
    faceDOF1 = SC_DOF_1(:, e);
    faceDOF2 = SC_DOF_2(:, e);

    LS_inner     = A1_stokes(:,:,e) \ A2_stokes(:,:,e);
    LS_RHS_inner = A1_stokes(:,:,e) \ Af_stokes(:,:,e);
    LS_outer     = A1_darcy(:,:,e)  \ A2_darcy(:,:,e);
    LS_RHS_outer = A1_darcy(:,:,e)  \ Af_darcy(:,:,e);

    m_tpath    = patch_tpaths(:, :, e, :);
    n_interface = reshape(patch_normals(:, :, e, :), 2, 1, [], 2);

    % Interface normals (varying with quadrature point for curved interfaces)
    % nF = Stokes outward normal (pointing into Darcy), nS = Darcy outward normal
    nF = -n_interface(:, :, :, 1);           % (2,1,s_quad)
    tF =  pagemtimes([0,-1;1,0], nF);
    nS =  n_interface(:, :, :, 2);           % (2,1,s_quad)

    % Scalar normal/tangent components as (s_quad,1) column vectors
    nS1 = reshape(nF(1,:,:), s_quad, 1);
    nS2 = reshape(nF(2,:,:), s_quad, 1);
    nD1 = reshape(nS(1,:,:), s_quad, 1);
    nD2 = reshape(nS(2,:,:), s_quad, 1);
    tS1 = reshape(tF(1,:,:), s_quad, 1);
    tS2 = reshape(tF(2,:,:), s_quad, 1);

    % Transfer path vectors and Stokes path lengths
    m_i = reshape(m_tpath(:,:,:,1), 2, s_quad);  % (2,s_quad) Stokes side
    m_o = reshape(m_tpath(:,:,:,2), 2, s_quad);  % (2,s_quad) Darcy side
    len_cs_inner = vecnorm(m_i).';               % (s_quad,1)

    %% Basis functions
    % Legendre 1D on each face
    l1 = legendre_basis(k, reshape(patch_locals(:,:,e,1), s_quad, 1));  % (s_quad,d1)
    l2 = legendre_basis(k, reshape(patch_locals(:,:,e,2), s_quad, 1));  % (s_quad,d1)

    % Dubiner 2D at FACE coordinates (for u_h: mixed-point evaluation)
    s_bdy_i = patch_locals_2d(1,:,e,1); t_bdy_i = patch_locals_2d(2,:,e,1);
    P1 = dubiner_basis(k, s_bdy_i(:), t_bdy_i(:));   % (s_quad,d2)
    s_bdy_o = patch_locals_2d(1,:,e,2); t_bdy_o = patch_locals_2d(2,:,e,2);
    P2 = dubiner_basis(k, s_bdy_o(:), t_bdy_o(:));   % (s_quad,d2)

    % Dubiner 2D at INTERFACE coordinates (for L_h, u_D: extrapolated)
    s_if_i = patch_interface_2d(1,:,e,1); t_if_i = patch_interface_2d(2,:,e,1);
    Ep1 = dubiner_basis(k, s_if_i(:), t_if_i(:));    % (s_quad,d2) — L_h at interface
    s_if_o = patch_interface_2d(1,:,e,2); t_if_o = patch_interface_2d(2,:,e,2);
    Ep2 = dubiner_basis(k, s_if_o(:), t_if_o(:));    % (s_quad,d2) — u_D at interface

    % Dubiner 2D along PATH (tp_quad*s_quad points, tp_quad varies fastest)
    s_path_i = reshape(patch_path_2d(1,:,e,:,1), s_quad, tp_quad).';  % (tp_quad,s_quad)
    t_path_i = reshape(patch_path_2d(2,:,e,:,1), s_quad, tp_quad).';
    CSp1 = dubiner_basis(k, s_path_i(:), t_path_i(:));  % (tp_quad*s_quad,d2)
    s_path_o = reshape(patch_path_2d(1,:,e,:,2), s_quad, tp_quad).';
    t_path_o = reshape(patch_path_2d(2,:,e,:,2), s_quad, tp_quad).';
    CSp2 = dubiner_basis(k, s_path_o(:), t_path_o(:));  % (tp_quad*s_quad,d2)

    % Path tangent vectors expanded to (tp_quad*s_quad,1): each value repeated tp_quad times
    mx_i = m_i(1, wI)';  my_i = m_i(2, wI)';
    mx_o = m_o(1, wI)';  my_o = m_o(2, wI)';

    % Stokes transfer path integrals: ∫ψ_β * m_α dl  →  (s_quad,d2)
    p1m_1 = reshape(sum(reshape(wCS .* CSp1 .* mx_i, [tp_quad d2*s_quad]), 1), [s_quad d2]);
    p1m_2 = reshape(sum(reshape(wCS .* CSp1 .* my_i, [tp_quad d2*s_quad]), 1), [s_quad d2]);

    % κ⁻¹ at Darcy path quadrature points
    x_path = reshape(partition_darcy.quadrature.path(1,:,e,:), s_quad, tp_quad).';  % (tp_quad,s_quad)
    y_path = reshape(partition_darcy.quadrature.path(2,:,e,:), s_quad, tp_quad).';
    kI_path = DarcyProblem.get_inverse_tensor(diffusion_tensor, x_path, y_path);
    [kI11, kI21, kI12, kI22] = deal(kI_path{:});  % each (tp_quad,s_quad) → (:) is (tp_quad*s_quad,1)

    % Darcy transfer path integrals weighted by κ⁻¹: ∫(κ⁻¹)_{αβ} ψ_γ * m_β dl
    Kp2m_1 = reshape(sum(reshape(wCS .* CSp2 .* (kI11(:) .* mx_o + kI21(:) .* my_o), [tp_quad d2*s_quad]), 1), [s_quad d2]);
    Kp2m_2 = reshape(sum(reshape(wCS .* CSp2 .* (kI12(:) .* mx_o + kI22(:) .* my_o), [tp_quad d2*s_quad]), 1), [s_quad d2]);

    % Stokes path integral for pressure mean: ∫ψ_β dl weighted by path length
    I_CS_1 = len_cs_inner .* reshape(sum(reshape(wCS .* CSp1, [tp_quad d2*s_quad]), 1), [s_quad d2]);

    O = zeros(s_quad, d2);

    %% T1: <E(νL₁-p₁I)n₁ - αn₁ - ντu₁ + ντû₁ + BJS(û₁+∫L₁m)·t - p̂₂n₂ + ∫κ⁻¹u₂m·n₂, μ₁>
    tau = patch_tau(:, :, e, 1);
    patch_len = abs(patch_area(:, e, 1));
    wphi = wP .* l1;
    test_inner = @(C) matfun(@(M) patch_len * wphi'*M, C);

    % Static condensation: volume DOFs → (2*d1, 7*d2+1) → (2*d1, 6*d1+1)
    L_term = nu*test_inner({Ep1.*nS1, O, Ep1.*nS2, O;
                             O, Ep1.*nS1, O, Ep1.*nS2}) ...
           + wbjs*test_inner({p1m_1.*(tS1.*tS1), p1m_1.*(tS2.*tS1), p1m_2.*(tS1.*tS1), p1m_2.*(tS2.*tS1);
                               p1m_1.*(tS1.*tS2), p1m_1.*(tS2.*tS2), p1m_2.*(tS1.*tS2), p1m_2.*(tS2.*tS2)});
    u_term = -nu*tau*test_inner({P1, O; O, P1});
    p_term = -test_inner({Ep1.*nS1; Ep1.*nS2});
    a_term = test_inner({O(:,1); O(:,1)});
    S1_inner = [L_term, u_term, p_term, a_term];

    val_T1_1(:,:,e) = S1_inner*LS_inner;
    row_T1_1(:,:,e) = kron(ones(1,3*2*d1+1), DOF1);
    col_T1_1(:,:,e) = kron(ones(2*d1,1), faceDOF1');

    u_term = test_inner({Kp2m_1.*nD1, Kp2m_2.*nD1;
                         Kp2m_1.*nD2, Kp2m_2.*nD2});
    p_term = test_inner({O; O});
    S1_outer = [u_term, p_term];

    val_T1_2(:,:,e) = S1_outer*LS_outer;
    row_T1_2(:,:,e) = kron(ones(1,3*d1), DOF1);
    col_T1_2(:,:,e) = kron(ones(2*d1,1), faceDOF2');
    rhs_T1(:,:,e)   = -S1_inner*LS_RHS_inner - S1_outer*LS_RHS_outer;

    % Trace terms
    val_T1_3(:,:,e) = test_inner({(nu*tau + wbjs*(tS1.*tS1)).*l1, (wbjs*(tS1.*tS2)).*l1;
                                   (wbjs*(tS2.*tS1)).*l1, (nu*tau + wbjs*(tS2.*tS2)).*l1});
    row_T1_3(:,:,e) = kron(ones(1,2*d1), DOF1);
    col_T1_3(:,:,e) = kron(ones(2*d1,1), DOF1');

    val_T1_4(:,:,e) = -test_inner({l2.*nD1; l2.*nD2});
    row_T1_4(:,:,e) = kron(ones(1,d1), DOF1);
    col_T1_4(:,:,e) = kron(ones(2*d1,1), DOF2');

    val_T1_5(:,:,e) = -test_inner({P1(:,1).*nS1; P1(:,1).*nS2});
    row_T1_5(:,:,e) = DOF1;
    col_T1_5(:,:,e) = A_DOF;

    % Source
    T1_1_e = reshape(T1(:,1,1,e), s_quad, 1);
    T1_2_e = reshape(T1(:,1,2,e), s_quad, 1);
    rhs_T1(:,:,e) = rhs_T1(:,:,e) + test_inner({T1_1_e; T1_2_e});

    %% T2: <û₁·n₁ + ∫L₁m·n₁ + u₂·n₂ + τ(p₂ - p̂₂), μ₂>
    tau = patch_tau(:, :, e, 2);
    patch_len = abs(patch_area(:, e, 2));
    wphi = wP .* l2;
    test_outer = @(C) matfun(@(M) patch_len * wphi'*M, C);

    % Static condensation: (d1, 7*d2+1) → (d1, 6*d1+1)
    L_term = test_outer({p1m_1.*nS1, p1m_1.*nS2, p1m_2.*nS1, p1m_2.*nS2});
    u_term = test_outer({O, O});
    p_term = test_outer({O});
    a_term = test_outer({O(:,1)});
    S2_inner = [L_term, u_term, p_term, a_term];

    val_T2_1(:,:,e) = S2_inner*LS_inner;
    row_T2_1(:,:,e) = kron(ones(1,3*2*d1+1), DOF2);
    col_T2_1(:,:,e) = kron(ones(d1,1), faceDOF1');

    u_term = test_outer({Ep2.*nD1, Ep2.*nD2});
    p_term = tau*test_outer({P2});
    S2_outer = [u_term, p_term];

    val_T2_2(:,:,e) = S2_outer*LS_outer;
    row_T2_2(:,:,e) = kron(ones(1,3*d1), DOF2);
    col_T2_2(:,:,e) = kron(ones(d1,1), faceDOF2');
    rhs_T2(:,:,e)   = -S2_inner*LS_RHS_inner - S2_outer*LS_RHS_outer;

    % Trace terms
    val_T2_3(:,:,e) = test_outer({l1.*nS1, l1.*nS2});
    row_T2_3(:,:,e) = kron(ones(1,2*d1), DOF2);
    col_T2_3(:,:,e) = kron(ones(d1,1), DOF1');

    val_T2_4(:,:,e) = -tau*test_outer({l2});
    row_T2_4(:,:,e) = kron(ones(1,d1), DOF2);
    col_T2_4(:,:,e) = kron(ones(d1,1), DOF2');

    % Source
    T2_e = T2(:, :, e);
    rhs_T2(:,:,e) = rhs_T2(:,:,e) + test_outer({T2_e});

    %% Alpha: <∫_Kext p₁ + |Ωₛ|α, ξ>  (pressure mean constraint)
    patch_len = patch_area(:, e, 1);
    wphi = wP .* P1(:,1);
    test_alpha = @(C) matfun(@(M) patch_len * wphi'*M, C);

    L_term  = test_alpha({O, O, O, O});
    u_term  = test_alpha({O, O});
    p_term  = test_alpha({I_CS_1});
    a_term  = test_alpha({O(:,1)});
    S_alpha = [L_term, u_term, p_term, a_term];

    val_A_1(:,:,e) = S_alpha*LS_inner;
    row_A_1(:,:,e) = A_DOF;
    col_A_1(:,:,e) = faceDOF1;
    rhs_A(:,:,e)   = -S_alpha*LS_RHS_inner;

    val_A_2(:,:,e) = measure_stokes * test_alpha({P1(:,1)});
    row_A_2(:,:,e) = A_DOF;
    col_A_2(:,:,e) = A_DOF;
end
print_message('Transmission conditions ready.')
%% IV. ASSEMBLING GLOBAL SYSTEM
print_message('Assembling global system...')
[row_inner, col_inner, rowRHS_inner] = stokes_problem.get_assemble_indices(mesh_stokes, k);
[row_outer, col_outer, rowRHS_outer] = darcy_problem.get_assemble_indices(mesh_darcy, k);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% DOF RE-NUMBERING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Delete entries that belong to interface DOF from original systems
innerDOF = [d1*((partition_stokes.faces(:))-1) + (1:d1);
          d1*((partition_stokes.faces(:))-1) + (1:d1) + d1*N_F_1];
val_inner(ismember(row_inner,innerDOF)) = [];
col_inner(ismember(row_inner,innerDOF)) = [];
row_inner(ismember(row_inner,innerDOF)) = [];
rhs_inner(ismember(rowRHS_inner,innerDOF)) = [];
rowRHS_inner(ismember(rowRHS_inner,innerDOF)) = [];

outerDOF = d1*((partition_darcy.faces(:))-1) + (1:d1);
val_outer(ismember(row_outer,outerDOF)) = [];
col_outer(ismember(row_outer,outerDOF)) = [];
row_outer(ismember(row_outer,outerDOF)) = [];
rhs_outer(ismember(rowRHS_outer,outerDOF)) = [];
rowRHS_outer(ismember(rowRHS_outer,outerDOF)) = [];

%%% Shift entries belonging to second meshRow1
row_outer = row_outer + N_DOF_1;
col_outer = col_outer + N_DOF_1;
rowRHS_outer = rowRHS_outer + N_DOF_1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% SYSTEM ASSEMBLY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
row = [row_inner(:); row_T1_1(:); row_T1_2(:); row_T1_3(:); row_T1_4(:); row_T1_5(:);
       row_outer(:); row_T2_1(:); row_T2_2(:); row_T2_3(:); row_T2_4(:);
                row_A_1(:);  row_A_2(:)];
col = [col_inner(:); col_T1_1(:); col_T1_2(:); col_T1_3(:); col_T1_4(:); col_T1_5(:);
       col_outer(:); col_T2_1(:); col_T2_2(:); col_T2_3(:); col_T2_4(:);
                col_A_1(:);  col_A_2(:)];
val = [val_inner(:); val_T1_1(:); val_T1_2(:); val_T1_3(:); val_T1_4(:); val_T1_5(:);
       val_outer(:); val_T2_1(:); val_T2_2(:); val_T2_3(:); val_T2_4(:); 
                val_A_1(:);  val_A_2(:)];
H = sparse(row,col,val);

    F = [rhs_inner(:); rhs_T1(:);
         rhs_outer(:); rhs_T2(:);
                  rhs_A(:)];
    
    RowRHS = [rowRHS_inner(:); rowRHS_T1(:);
              rowRHS_outer(:); rowRHS_T2(:);
                               rowRHS_A(:)];
RHS = accumarray(RowRHS,F);

% \sum_K \rho_K = 0
RHS(end+1) = 0;
H(end+1,:) = [zeros(1,2*d1*N_F_1), mesh_stokes.volume, zeros(1,d1*N_F_2), 0];
H(:,end+1) = [zeros(1,2*d1*N_F_1), mesh_stokes.volume, zeros(1,d1*N_F_2), 0, 0]';

%% V. BOUNDARY CONDITIONS
print_message('Implementing boundary conditions...')
[x1_hat, dirfaces_1, g1_hat] = stokes_problem.get_boundary_conditions(mesh_stokes);
[x2_hat, dirfaces_2, g2_hat] = darcy_problem.get_boundary_conditions(mesh_darcy);
dirfaces = [dirfaces_1(:); dirfaces_2(:) + N_DOF_1];

x_hat = [x1_hat; x2_hat; 0; 0];
RHS(1:N_DOF_1) = RHS(1:N_DOF_1) + g1_hat;
RHS(N_DOF_1 + (1:N_DOF_2)) = RHS(N_DOF_1 + (1:N_DOF_2)) + g2_hat;
RHS = RHS - H*x_hat;

%% VI-1. SOLVING GLOBAL SYSTEM
print_message('Solving linear systems...')
freefaces = setdiff((1:size(H,1)).',dirfaces);
x_hat(freefaces) = H(freefaces,freefaces)\RHS(freefaces);

if nargout > 2
    global_system = struct( ...
        "matrix",   H, ...
        "rhs",      RHS, ...
        "dof",      freefaces, ...
        "solution", x_hat, ...
        "residue",  H*x_hat - RHS);
end

x1_hat = x_hat(1:N_DOF_1);
x2_hat = x_hat(N_DOF_1 + (1:N_DOF_2));
alpha_h = x_hat(A_DOF);
print_message('Global system solved.')
%% VI-2. SOLVING LOCAL SYSTEMS
print_message('Recovering local solutions...')
%%% RECOVERING LOCAL SOLUTIONS
[stokes_sols, stokes_hybrid_sols] = stokes_problem.get_local_solutions(x1_hat, localsolvers_inner, mesh_stokes);
[darcy_sols, darcy_hybrid_sols] = darcy_problem.get_local_solutions(x2_hat, localsolvers_outer, mesh_darcy);
print_message('Local solutions computed.')

%% VII. OUTPUT
% pressure post-process: add mean-pressure correction to constant coefficient
stokes_sols{end}(1, :) = stokes_sols{end}(1, :) + alpha_h;

stokes_solution = FESolution(stokes_problem, mesh_stokes);
stokes_solution = stokes_solution.assign_unknowns(stokes_sols, stokes_hybrid_sols);

darcy_solution = FESolution(darcy_problem, mesh_darcy);
darcy_solution = darcy_solution.assign_unknowns(darcy_sols, darcy_hybrid_sols);

end

function [DOF, SC_DOF] = get_partition_dof(mesh, partition, d1, nd, shift)
    arguments
        mesh
        partition
        d1
        nd
        shift = 0
    end
    N_PATCHES = length(partition.faces);
    N_FACES = size(mesh.faces, 2);

    get_1d_dof = @(e) d1*(e(:).' - 1) + (1:d1)';
    if nd == 2
        vec = permute([0, d1*N_FACES] + shift, [1 3 2]);
    elseif nd == 1
        vec = shift;
    end
    get_DOF = @(e) reshape(get_1d_dof(e) + vec, [], 1);

    DOF = cell2mat(cellfun(get_DOF, num2cell(partition.faces), ...
                           'UniformOutput', 0));

    partition_FBE = mesh.boundaries(:, partition.elements);
    SC_DOF = cell2mat(cellfun(get_DOF, ...
                              mat2cell(partition_FBE, 3, ones(N_PATCHES, 1)), ...
                              'UniformOutput', 0));
end