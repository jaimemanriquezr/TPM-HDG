function [solution_1, solution_2, global_system, cache] = solve(obj, k, options)
% Coupled TPM solve: condense both subproblems (condense_global_system),
% concatenate their global-support layouts, delete each side's interface-face
% flux rows and substitute the coupling's transmission rows over the
% induced patch partition, apply each side's boundary conditions, solve,
% and recover both local solutions.
%
% Under the fsi coupling, the T2 history h0 is transported onto Sigma on
% the partition this solve builds: the coupling's explicit Ext(h0) when
% given, otherwise trace_extension_coefficients of the solid stress
% history (see the history option of InterfaceProblem.fsi); without a
% stress history there is no transport.
%
% Under the stokes_biot coupling, the pore-pressure level scalar is kept
% only when the applied boundary conditions leave that level free (see
% InterfaceProblem.stokes_biot); otherwise it is eliminated.
% global_system.level_scalar reports its value.
%
% Cache (fsi coupling without symbolic transmission rows only): for
% repeated solves in which only side 1 (the fluid) changes, e.g. a Picard
% loop. Pass struct("frozen", 2) on the first call to declare side 2
% iterate-independent; the returned cache then holds the partitions, the
% layout, side 2's condensed system and boundary-condition outputs, the
% transmission functionals (fsi_functionals) and side 1's transfer blocks,
% and is passed back on later calls. Its key (k, mesh sizes, quadrature
% orders, coupling options, cdt, slip coefficient and interface
% histories) must match the call, or solve errors. The sources and side 2
% are assumed unchanged. Results equal the uncached solve's exactly.
% Without a cache the fourth output is struct.empty.
arguments
    obj (1,1) InterfaceProblem
    k   (1,1) double = obj.problems(1).Order
    options.Cache struct = struct.empty   % see above
end
cache = options.Cache;
use_cache = ~isempty(cache);
if use_cache
    if ~isfield(cache, "frozen") || ~isequal(cache.frozen, 2)
        error("InterfaceProblem:cacheFrozen", ...
            "A cache must declare frozen = 2 (side 2 iterate-independent).");
    end
    if obj.coupling.kind ~= "fsi" || ~isempty(obj.transmission)
        error("InterfaceProblem:cacheUnsupported", ...
            "A cache is supported for the built-in fsi coupling only.");
    end
end
% per-side curves: curve_2 = [] for one shared curve (interface_curves)
[curve_1, curve_2] = interface_curves(obj.interface);
if ~isempty(curve_2) && (obj.coupling.kind ~= "fsi" || ~isempty(obj.transmission))
    error("InterfaceProblem:curvePair", ...
        "A pair of interface curves is supported by the built-in fsi coupling only.");
end
p1 = obj.problems(1);
p2 = obj.problems(2);

mesh_1 = p1.Mesh;      mesh_2 = p2.Mesh;
assert(~isempty(mesh_1) && ~isempty(mesh_2), ...
    "InterfaceProblem:mesh", "both subproblems must carry their meshes");
% bind each side's spaces to its mesh (enumeration needs bound spaces);
% assemble_transmission reads obj.problems, so rebind there too
[p1.Spaces.Mesh] = deal(mesh_1);
[p2.Spaces.Mesh] = deal(mesh_2);
obj.problems = [p1, p2];

sq = obj.options.quadrature_patch;  if isempty(sq), sq = k + 1; end
tq = obj.options.quadrature_path;   if isempty(tq), tq = k + 2; end

filled = false;
if use_cache
    key = struct("k", k, ...
        "elements", [size(mesh_1.elements, 2), size(mesh_2.elements, 2)], ...
        "faces", [size(mesh_1.faces, 2), size(mesh_2.faces, 2)], ...
        "interface_faces", [numel(mesh_1.intfaces), numel(mesh_2.intfaces)], ...
        "sq", sq, "tq", tq, ...
        "mappings", [obj.options.mapping_1, obj.options.mapping_2], ...
        "swap_transmission", obj.coupling.swap_transmission, ...
        "interface_law", obj.coupling.interface_law, ...
        "pin_extra", obj.coupling.pin_extra, ...
        "kind", obj.coupling.kind, ...
        "cdt", obj.coupling.cdt, ...
        "viscosity", obj.coupling.viscosity, ...
        "young_modulus", obj.coupling.young_modulus, ...
        "poisson_ratio", obj.coupling.poisson_ratio, ...
        "domain_measure", obj.coupling.domain_measure, ...
        "solid_weight", obj.coupling.solid_weight, ...
        "patch_weight", obj.coupling.patch_weight, ...
        "sliver_correction", obj.coupling.sliver_correction, ...
        "slip_coefficient", obj.coupling.slip_coefficient, ...
        "history_skeleton", obj.coupling.history_skeleton, ...
        "history_extension", obj.coupling.history_extension, ...
        "history_stress", obj.coupling.history_stress, ...
        "history_rotation", obj.coupling.history_rotation, ...
        "tip_mapping", obj.options.tip_mapping, ...
        "collapsed_normal", obj.options.collapsed_normal, ...
        "frozen", cache.frozen);
    filled = isfield(cache, "key");
    if filled && ~isequal(cache.key, key)
        error("InterfaceProblem:cacheKey", ...
            "The cache was built for another problem (k, meshes, quadrature or coupling options differ).");
    end
end

%% condensed subsystems and coupled layout
S1 = condense_global_system(p1, mesh_1);
if filled
    S2 = cache.S2;
    [off, n_dof, n_extra, extra_base, del] = deal(cache.layout.off, ...
        cache.layout.n_dof, cache.layout.n_extra, cache.layout.extra_base, ...
        cache.layout.del);
    [part_1, part_2] = cache.partitions{:};
else
    S2 = condense_global_system(p2, mesh_2);
    off   = S1.L.n_dof;
    n_dof = off + S2.L.n_dof;
    % coupling-owned extra scalars (e.g. the Stokes-Darcy extrapolated
    % pressure alpha) live after both subproblem layouts
    n_extra = 0;
    if isfield(obj.coupling, "n_extra"), n_extra = obj.coupling.n_extra; end
    extra_base = n_dof;
    n_dof = n_dof + n_extra;

    %% induced patch partition
    [part_1, part_2] = get_partitions(curve_1, mesh_1, mesh_2, ...
        quadrature_patch = sq, quadrature_path = tq, ...
        vertex_tol = obj.options.vertex_tol, ...
        vertex_tol_relative = obj.options.vertex_tol_relative, ...
        tip_mapping = obj.options.tip_mapping, ...
        collapsed_normal = obj.options.collapsed_normal, ...
        mapping_1 = obj.options.mapping_1, mapping_2 = obj.options.mapping_2, ...
        interface_2 = curve_2);

    check_collapsed_rows(obj, part_1);

    %% base system: block-diagonal, interface flux rows deleted
    tr1 = p1.GlobalSpaces([p1.GlobalSpaces.Support] == "face");
    tr2 = p2.GlobalSpaces([p2.GlobalSpaces.Support] == "face");
    if obj.coupling.kind == "fsi", tr2 = transmitted_trace(p2); end
    del_1 = get_global_dof(tr1, subset=mesh_1.intfaces);
    del_2 = get_global_dof(tr2, subset=mesh_2.intfaces);
    del   = [del_1(:); del_2(:) + off];
end

H   = blkdiag(S1.H, S2.H, sparse(n_extra, n_extra));
RHS = [S1.RHS; S2.RHS; zeros(n_extra, 1)];
H(del, :) = 0;
RHS(del)  = 0;

%% transmission rows
if ~isempty(obj.transmission)
    [Hc, Fc] = assemble_transmission(obj, k, S1, S2, mesh_1, mesh_2, ...
        part_1, part_2, off, n_dof, sq, tq);
    H   = H + Hc;
    RHS = RHS + Fc;
else
    switch obj.coupling.kind
        case "darcy_darcy"
            [Hc, Fc] = darcy_darcy_rows(obj, k, S1, S2, ...
                part_1, part_2, off, n_dof, sq, tq);
        case "stokes_darcy"
            [Hc, Fc] = stokes_darcy_rows(obj, k, S1, S2, ...
                part_1, part_2, off, n_dof, extra_base, sq, tq);
        case "fsi"
            if filled
                functionals = cache.functionals;
            else
                functionals = fsi_functionals(obj, k, mesh_1, mesh_2, ...
                    part_1, part_2, off, sq, tq, curve_1);
            end
            [Hc, Fc] = fsi_rows(k, S1, S2, off, n_dof, extra_base, functionals);
        case "stokes_biot"
            [Hc, Fc] = stokes_biot_rows(obj, k, S1, S2, mesh_1, mesh_2, ...
                part_1, part_2, off, n_dof, extra_base, sq, tq);
        otherwise
            error("InterfaceProblem:coupling", ...
                "Unknown coupling kind ""%s"".", obj.coupling.kind);
    end
    H   = H + Hc;
    RHS = RHS + Fc;
end

%% boundary conditions per side (skeleton prefixes only, as in the BVPs)
if use_cache
    % side 1 reuses its condensed system and transfer blocks; side 2's
    % outputs are condensed and applied once
    bc_1 = struct.empty;
    if filled, bc_1 = cache.bc_1; end
    [x1, g1, Hh1, dof1, bc_1] = p1.BoundaryConditions.apply(p1, S = S1, cache = bc_1);
    if filled
        [x2, g2, Hh2, dof2] = deal(cache.bc_2.x_hat, cache.bc_2.g_hat, ...
            cache.bc_2.H_hat, cache.bc_2.dof_list);
    else
        [x2, g2, Hh2, dof2] = p2.BoundaryConditions.apply(p2, S = S2);
        cache.key = key;
        cache.partitions = {part_1, part_2};
        cache.layout = struct("off", off, "n_dof", n_dof, "n_extra", n_extra, ...
            "extra_base", extra_base, "del", del);
        cache.S2 = S2;
        cache.bc_2 = struct("x_hat", x2, "g_hat", g2, "H_hat", Hh2, "dof_list", dof2);
        cache.functionals = functionals;
    end
    cache.bc_1 = bc_1;
else
    [x1, g1, Hh1, dof1] = p1.BoundaryConditions.apply(p1);
    [x2, g2, Hh2, dof2] = p2.BoundaryConditions.apply(p2);
end
% uncondensed local modes (FESpace.UncondensedModes) append
% element-supported unknowns after a side's global spaces; they carry
% no boundary condition
[x1, g1, Hh1] = pad_side(S1.L.n_dof, x1, g1, Hh1);
[x2, g2, Hh2] = pad_side(S2.L.n_dof, x2, g2, Hh2);
x_hat = [x1; x2; zeros(n_extra, 1)];
g     = [g1; g2; zeros(n_extra, 1)];
weak   = [reshape(dof1.essential.weak, [], 1); reshape(dof2.essential.weak, [], 1) + off];
strong = [reshape(dof1.essential.strong, [], 1); reshape(dof2.essential.strong, [], 1) + off];
if isfield(obj.coupling, "pin_extra") && obj.coupling.pin_extra
    % e.g. FSI under pressure_mean = "pinned": the alpha scalar is
    % eliminated (p0_h = 0)
    strong = [strong; extra_base + (1:n_extra).'];
end
if obj.coupling.kind == "stokes_biot" && ~stokes_biot_pressure_mode(obj, k, ...
        mesh_1, mesh_2, tr1, tr2, dof1, dof2)
    % the boundary conditions fix the pore-pressure level: the level
    % scalar beta is eliminated (beta = 0)
    strong = [strong; extra_base + (1:n_extra).'];
end
Hh = blkdiag(Hh1, Hh2, sparse(n_extra, n_extra));
H(weak, :)  = Hh(weak, :);
RHS(weak)   = 0;
RHS = RHS - H*x_hat + g;
free = setdiff((1:n_dof).', strong);
x_hat(free) = H(free, free) \ RHS(free);

%% recovery per side
solution_1 = recover(p1, mesh_1, S1, x_hat(1:off), k);
solution_2 = recover(p2, mesh_2, S2, x_hat(off + (1:S2.L.n_dof)), k);
if obj.coupling.kind == "stokes_darcy"
    % extrapolated-pressure post-process: add alpha to the constant
    % coefficient of the Stokes pressure
    alpha_h = x_hat(extra_base + 1);
    solution_1.Unknowns.field(3).dof(1, :) = ...
        solution_1.Unknowns.field(3).dof(1, :) + alpha_h;
elseif obj.coupling.kind == "fsi"
    % pressure-level post-process: sigma <- sigma + p0 * phi_0 * I on the
    % fluid stress (constant coefficients of the sigma_11/sigma_22
    % Dubiner blocks)
    p0_h = x_hat(extra_base + 1);
    d2 = nchoosek(k + 2, 2);
    solution_1.Unknowns.field(1).dof(1, :) = ...
        solution_1.Unknowns.field(1).dof(1, :) + p0_h;
    solution_1.Unknowns.field(1).dof(3*d2 + 1, :) = ...
        solution_1.Unknowns.field(1).dof(3*d2 + 1, :) + p0_h;
end
if nargout > 2
    global_system = struct("matrix", H, "rhs", RHS, "dof", free, ...
        "solution", x_hat, ...
        "residue", norm(H(free, free)*x_hat(free) - RHS(free)), ...
        "partition_1", part_1, "partition_2", part_2);   % as used, for history transport
    if obj.coupling.kind == "stokes_biot"
        % the pore-pressure level scalar beta (0 when eliminated; for
        % consistent data it tends to 0 with h, a large value flags
        % incompatible data)
        global_system.level_scalar = x_hat(extra_base + 1);
    end
end
end

%% ======================================================================
function tf = stokes_biot_pressure_mode(obj, k, mesh_1, mesh_2, tr1, tr2, dof1, dof2)
% True when the Stokes-Biot step has the pore-pressure null mode
% (p, p_hat, sigma_P, sigma_F) -> (p + c, p_hat + c, sigma_P - c I,
% sigma_F - c I). The Biot equations are blind to it when s0 = 0 and
% alpha = 1: the shift leaves residuals (alpha - 1)/(2(lambda + mu)) in
% the constitutive and s0 + alpha (alpha - 1)/(lambda + mu) in the mass
% equation, and below the relative tolerance 1e-10 (s0 (lambda + mu) and
% |1 - alpha|) they cannot fix the level in double precision. Off Sigma, a
% prescribed p_hat fixes it, and so does a natural displacement or velocity
% component c on a face with n_c ~= 0 (the shift changes the traction
% datum (sigma n)_c by -c n_c); a natural component with n_c = 0 (slip and
% symmetry walls, rollers on axis-aligned faces) does not.
c = obj.coupling;
[lambda, mu] = LinearElasticityProblem.get_lame_parameters(c.young_modulus, c.poisson_ratio);
level_tol = 1e-10;
tf = c.storage*(lambda + mu) <= level_tol && abs(1 - c.biot_alpha) <= level_tol;
if ~tf, return; end
d1 = k + 1;
fixed_1 = [reshape(dof1.essential.strong, [], 1); reshape(dof1.essential.weak, [], 1)];
fixed_2 = [reshape(dof2.essential.strong, [], 1); reshape(dof2.essential.weak, [], 1)];
outer_1 = outer_faces(mesh_1);
if ~isempty(outer_1)
    u_1 = reshape(get_global_dof(tr1, subset=outer_1), 2*d1, []);
    tf = ~loads_level(u_1, fixed_1, mesh_1, outer_1, d1);
end
outer_2 = outer_faces(mesh_2);
if tf && ~isempty(outer_2)
    d_2 = reshape(get_global_dof(tr2, subset=outer_2), 3*d1, []);
    tf = ~loads_level(d_2(1:2*d1, :), fixed_2, mesh_2, outer_2, d1) ...
        && ~any(ismember(reshape(d_2(2*d1 + (1:d1), :), [], 1), fixed_2));
end
end

function tf = loads_level(dofs, fixed, mesh, faces, d1)
% True when a natural component c of a two-component trace ([c_1; c_2]
% blocks of d1 dofs per face, one column per face) sits on a face whose
% unit normal has |n_c| > 1e-8
X = mesh.coordinates(1:2, mesh.faces(1:2, faces));
t = X(:, 2:2:end) - X(:, 1:2:end);
n_abs = abs([t(2, :); t(1, :)]) ./ vecnorm(t);          % |n_x|; |n_y|
natural = [any(~ismember(dofs(1:d1, :), fixed), 1); ...
           any(~ismember(dofs(d1 + (1:d1), :), fixed), 1)];
tf = any(natural & n_abs > 1e-8, "all");
end

function f = outer_faces(mesh)
% the boundary faces of the mesh that are not on the interface
n_faces = size(mesh.faces, 2);
on_boundary = accumarray(mesh.boundaries(:), 1, [n_faces, 1]) == 1;
f = setdiff(find(on_boundary).', mesh.intfaces);
end

%% ======================================================================
function check_collapsed_rows(obj, part_1)
% A face collapsed to one point of Sigma (get_partitions, tip_mapping =
% "collapse") is normal-free only in the kinematic row, and only when
% that row is tested on the collapsed face's own side: the row tested on
% the other side has zero measure there, and the traction row on the
% collapsed face would need the normal of Sigma at a point where it is
% undefined. So the kinematic row must be tested on the side that
% collapsed: fluid (side 1) faces need swap_transmission = true, solid
% (side 2) faces the default rows; both sides at once admit no
% orientation.
if ~isfield(part_1, "collapsed") || ~any(part_1.collapsed)
    return
end
if obj.coupling.kind ~= "fsi" || ~isempty(obj.transmission)
    error("InterfaceProblem:collapsedCoupling", ...
        "tip_mapping = ""collapse"" is supported for the built-in fsi " + ...
        "coupling only (coupling ""%s"").", obj.coupling.kind);
end
sides = unique(part_1.collapsed_side(part_1.collapsed));
if numel(sides) > 1
    error("InterfaceProblem:collapsedBothSides", ...
        "Faces of both meshes collapse to points of Sigma; the kinematic " + ...
        "row can be tested on one side only. Use tip_mapping = ""spread"".");
end
if obj.coupling.interface_law ~= "noslip"
    error("InterfaceProblem:collapsedSlip", ...
        "tip_mapping = ""collapse"" needs interface_law = ""noslip"": the " + ...
        "%s rows project onto the normal of Sigma, undefined at a collapse point.", ...
        obj.coupling.interface_law);
end
swap = obj.coupling.swap_transmission;
if sides == 1 && ~swap
    error("InterfaceProblem:collapsedRows", ...
        "Fluid (side 1) faces collapse to a point of Sigma: the kinematic " + ...
        "row must be tested on them. Set swap_transmission = true.");
elseif sides == 2 && swap
    error("InterfaceProblem:collapsedRows", ...
        "Solid (side 2) faces collapse to a point of Sigma: the kinematic " + ...
        "row must be tested on them. Set swap_transmission = false.");
end
end

%% ======================================================================
function sol = recover(pde, mesh, S, x_side, k)
% local recovery: x_loc = A11 \ A1f - sum_t (A11 \ A12_t) x_t
n_el = size(mesh.elements, 2);
x_loc = reshape(S.S1f, [], 1, n_el);
% every condensed block: the global spaces, plus the uncondensed local modes
for t = 1:numel(S.S12)
    S12flat = reshape(permute(S.S12{t}, [1 2 7 5 3 4 6]), ...
        size(S.S12{t}, 1), S.d_g(t)*S.f_g(t), n_el);
    u_t = reshape(x_side(S.L.footprint{t}), S.d_g(t)*S.f_g(t), 1, n_el);
    x_loc = x_loc - pagemtimes(S12flat, u_t);
end

sol = FESolution(pde, mesh, k);
U = get_field_functions(pde.Space, pde.Id.field);
[U.mesh] = deal(mesh);
dims_loc = arrayfun(@(sp) sp.Dimension, pde.LocalSpaces);
x_loc = reshape(x_loc, [], n_el);
field_dofs = mat2cell(x_loc, dims_loc(:), n_el);
[U.dof] = field_dofs{:};

U_hat = get_field_functions(pde.Space, pde.Id.hybrid);
[U_hat.mesh] = deal(mesh);
for t = 1:numel(pde.GlobalSpaces)
    U_hat(t).dof = x_side(S.L.offset(t) + (1:S.L.count(t)));
end
sol.Unknowns = struct('field', U, 'hybrid', U_hat);
end

%% ======================================================================
function [Hc, Fc] = darcy_darcy_rows(obj, k, S1, S2, ...
    part_1, part_2, off, n_dof, sq, tq)
% Scalar TPM transmission over the induced patches, with the T1/T2
% structure of fsi_rows (local vector per side: [q_x; q_y; p]):
%   T1 (tested on side-1 interface faces): the extrapolated numerical
%       fluxes balance,  sum_i < E(q_i).n_i + tau_i (E(p_i) - p_hat_i), mu_1 > = 0
%   T2 (tested on side-2 faces): the transferred traces agree,
%       < (p_hat_1 - PI_1) - (p_hat_2 - PI_2), mu_2 > = 0,
%       PI_i = int_path (kappa_i^{-1} E(q_i)) . m_i  (since grad p = -kappa^{-1} q)
d1 = k + 1;
d2 = nchoosek(k + 2, 2);
n_loc = 3*d2;
kappa = obj.coupling.kappa;
[~, wS]  = lgwt(sq, 0, 1);
[~, wTP] = lgwt(tq, 0, 1);
n_p = size(part_1.faces, 3);

S12flat = cell(1, 2);
S1f_pages = {S1.S1f, S2.S1f};
Ss = {S1, S2};
parts  = {part_1, part_2};
offs   = [0, off];
assert(numel(S1.S12) == numel(obj.problems(1).GlobalSpaces) && numel(S2.S12) == numel(obj.problems(2).GlobalSpaces), ...
    "InterfaceProblem:uncondensedModes", "uncondensed local modes are not supported by this coupling");
for i = 1:2
    S12flat{i} = reshape(permute(Ss{i}.S12{1}, [1 2 7 5 3 4 6]), ...
        n_loc, Ss{i}.d_g(1)*Ss{i}.f_g(1), []);
end

[rows_acc, cols_acc, vals_acc] = deal({});
rhs_acc = zeros(n_dof, 1);

for e = 1:n_p
    % per-side patch data
    [K, F, tau_e, s_loc, l, D_iface, D_face, n_iface, m, a_kappa, PI, foot, fdof] = ...
        deal(cell(1, 2));
    for i = 1:2
        P = parts{i};
        K{i} = P.elements(e);
        F{i} = P.faces(e);
        tau_e{i} = P.tau(e);
        s_loc{i} = reshape(P.quadrature.local(1, :, e), [], 1);
        l{i} = legendre_basis(k, s_loc{i});                       % sq x d1
        rs_if = P.quadrature.interface_2d(:, :, e);
        D_iface{i} = dubiner_basis(k, rs_if(1, :).', rs_if(2, :).');  % sq x d2 (extrapolated)
        rs_face = P.quadrature.local_2d(:, :, e);
        D_face{i} = dubiner_basis(k, rs_face(1, :).', rs_face(2, :).');  % sq x d2 (on the face)
        n_iface{i} = P.quadrature.normals(:, :, e);               % 2 x sq
        m{i} = P.quadrature.interface(:, :, e) - P.quadrature.physical(:, :, e);
        a_kappa{i} = kappa{i}.' \ m{i};                           % 2 x sq
        % path evaluation for the transferred trace
        rs_path = reshape(P.quadrature.path_2d(:, :, e, :), 2, []);   % 2 x (sq*tq)
        Dp = dubiner_basis(k, rs_path(1, :).', rs_path(2, :).');      % (sq*tq) x d2
        Dp = reshape(Dp, sq, tq, d2);
        % PI matrix vs q columns: PI(j, :) = sum_t wTP(t) * a.m-contraction
        PIx = reshape(sum(Dp .* reshape(wTP, 1, tq), 2), sq, d2) .* a_kappa{i}(1, :).';
        PIy = reshape(sum(Dp .* reshape(wTP, 1, tq), 2), sq, d2) .* a_kappa{i}(2, :).';
        PI{i} = [PIx, PIy];                                       % sq x 2d2
        foot{i} = Ss{i}.L.footprint{1}(:, K{i}) + offs(i);        % 3*d1 x 1
        fdof{i} = (F{i} - 1)*d1 + (1:d1).' + offs(i);             % d1 x 1
    end
    % sign convention as in fsi_rows: side 1's normal negated
    n_iface{1} = -n_iface{1};

    area_1 = part_1.area(e);
    area_2 = part_2.area(e);
    int_1 = @(B) area_1 * (l{1}.' * (wS .* B));   % test with mu_1
    int_2 = @(B) area_2 * (l{2}.' * (wS .* B));   % test with mu_2

    for i = 1:2
        % ---- T1: numerical-flux balance, tested on side 1 ----
        % flux extrapolated to Sigma; the stabilization pair stays on the
        % face (extrapolating it would be inconsistent by O(gap))
        Rq = [n_iface{i}(1, :).' .* D_iface{i}, n_iface{i}(2, :).' .* D_iface{i}];
        R1 = int_1([Rq, tau_e{i} * D_face{i}]);           % d1 x n_loc
        % x_loc = S1f - S12 x_hat, so the condensed part enters with a minus
        V  = -R1 * S12flat{i}(:, :, K{i});                % d1 x 3d1
        rows_acc{end+1} = repmat(fdof{1}, 1, 3*d1);       %#ok<AGROW>
        cols_acc{end+1} = repmat(foot{i}.', d1, 1);       %#ok<AGROW>
        vals_acc{end+1} = V;                              %#ok<AGROW>
        rhs_acc(fdof{1}) = rhs_acc(fdof{1}) - R1 * S1f_pages{i}(:, :, 1, 1, K{i});
        % trace term -tau_i <p_hat_i, mu_1>
        M1 = -tau_e{i} * int_1(l{i});                     % d1 x d1
        rows_acc{end+1} = repmat(fdof{1}, 1, d1);         %#ok<AGROW>
        cols_acc{end+1} = repmat(fdof{i}.', d1, 1);       %#ok<AGROW>
        vals_acc{end+1} = M1;                             %#ok<AGROW>

        % ---- T2: transferred-trace agreement, tested on side 2 ----
        sgn = 3 - 2*i;                                    % +1 side 1, -1 side 2
        R2 = -sgn * int_2([PI{i}, zeros(sq, d2)]);        % d1 x n_loc
        V  = -R2 * S12flat{i}(:, :, K{i});
        rows_acc{end+1} = repmat(fdof{2}, 1, 3*d1);       %#ok<AGROW>
        cols_acc{end+1} = repmat(foot{i}.', d1, 1);       %#ok<AGROW>
        vals_acc{end+1} = V;                              %#ok<AGROW>
        rhs_acc(fdof{2}) = rhs_acc(fdof{2}) - R2 * S1f_pages{i}(:, :, 1, 1, K{i});
        % trace term +/- <p_hat_i, mu_2>
        M2 = sgn * int_2(l{i});
        rows_acc{end+1} = repmat(fdof{2}, 1, d1);         %#ok<AGROW>
        cols_acc{end+1} = repmat(fdof{i}.', d1, 1);       %#ok<AGROW>
        vals_acc{end+1} = M2;                             %#ok<AGROW>
    end
end

rows = cellfun(@(x) x(:), rows_acc.', 'UniformOutput', false);
cols = cellfun(@(x) x(:), cols_acc.', 'UniformOutput', false);
vals = cellfun(@(x) x(:), vals_acc.', 'UniformOutput', false);
Hc = sparse(vertcat(rows{:}), vertcat(cols{:}), vertcat(vals{:}), n_dof, n_dof);
Fc = rhs_acc;
end

%% ======================================================================
function [Hc, Fc] = stokes_darcy_rows(obj, k, S1, S2, ...
    part_1, part_2, off, n_dof, extra_base, sq, tq)
% Stokes-Darcy transmission with BJS friction: side 1 = velocity-gradient
% Stokes (local [L(4d2); u(2d2); p(d2); xi]), side 2 = Darcy (local
% [q(2d2); p(d2)], scalar trace p_hat), plus the coupling-owned
% extrapolated-pressure scalar alpha at extra_base+1. Local solves are
% x = S1f - S12 x_hat, so a functional R on the local unknowns enters as
% -R*S12 in the matrix and -R*S1f in the rhs.
d1 = k + 1;
d2 = nchoosek(k + 2, 2);
nu    = obj.coupling.viscosity;
kappa = obj.coupling.diffusion;
wbjs  = obj.coupling.friction_bjs;
measure_stokes = obj.coupling.domain_measure;
T1_data = obj.coupling.source_neumann;
T2_data = obj.coupling.source_dirichlet;
[~, wS]  = lgwt(sq, 0, 1);
[~, wTP] = lgwt(tq, 0, 1);
n_p = size(part_1.faces, 3);
n_loc1 = 7*d2 + 1;
n_loc2 = 3*d2;
a_dof = extra_base + 1;

% the condensed part includes every global space of each side (faces and
% element-supported unknowns), not just the trace
S12f1 = cell(1, numel(S1.S12));
for t = 1:numel(S1.S12)
    S12f1{t} = reshape(permute(S1.S12{t}, [1 2 7 5 3 4 6]), n_loc1, S1.d_g(t)*S1.f_g(t), []);
end
S12f2 = cell(1, numel(S2.S12));
for t = 1:numel(S2.S12)
    S12f2{t} = reshape(permute(S2.S12{t}, [1 2 7 5 3 4 6]), n_loc2, S2.d_g(t)*S2.f_g(t), []);
end
tr1 = obj.problems(1).GlobalSpaces([obj.problems(1).GlobalSpaces.Support] == "face");
tr2 = obj.problems(2).GlobalSpaces([obj.problems(2).GlobalSpaces.Support] == "face");

matfun = @(f, C) cell2mat(cellfun(f, C, 'UniformOutput', false));
wP  = wS;
wCS = repmat(wTP(:), sq, 1);            % (tq*sq, 1), tq fastest
wI  = kron(1:sq, ones(1, tq));          % patch -> path index expansion

[rows_acc, cols_acc, vals_acc] = deal({});
rhs_acc = zeros(n_dof, 1);

for e = 1:n_p
    K1 = part_1.elements(e);  F1 = part_1.faces(e);
    K2 = part_2.elements(e);  F2 = part_2.faces(e);
    fd1 = get_global_dof(tr1, subset=F1);  DOF1 = fd1(:);            % 2d1
    fd2 = get_global_dof(tr2, subset=F2);  DOF2 = fd2(:) + off;      % d1
    S12K1 = cell2mat(cellfun(@(A) A(:, :, K1), S12f1, 'UniformOutput', false));
    S12K2 = cell2mat(cellfun(@(A) A(:, :, K2), S12f2, 'UniformOutput', false));
    foot1 = cell2mat(cellfun(@(fp) fp(:, K1), S1.L.footprint(:).', 'UniformOutput', false).');
    foot2 = cell2mat(cellfun(@(fp) fp(:, K2), S2.L.footprint(:).', 'UniformOutput', false).') + off;

    % normals / tangents (Stokes outward = -side1 normals)
    nF = -part_1.quadrature.normals(:, :, e);          % 2 x sq
    tF = [0, -1; 1, 0] * nF;
    nS = part_2.quadrature.normals(:, :, e);
    nS1v = nF(1, :).';  nS2v = nF(2, :).';
    nD1v = nS(1, :).';  nD2v = nS(2, :).';
    tS1v = tF(1, :).';  tS2v = tF(2, :).';

    m_i = part_1.quadrature.interface(:, :, e) - part_1.quadrature.physical(:, :, e);
    m_o = part_2.quadrature.interface(:, :, e) - part_2.quadrature.physical(:, :, e);
    len_cs = vecnorm(m_i).';

    l1 = legendre_basis(k, reshape(part_1.quadrature.local(1, :, e), [], 1));
    l2 = legendre_basis(k, reshape(part_2.quadrature.local(1, :, e), [], 1));
    P1f = dubiner_basis(k, part_1.quadrature.local_2d(1, :, e).', part_1.quadrature.local_2d(2, :, e).');
    P2f = dubiner_basis(k, part_2.quadrature.local_2d(1, :, e).', part_2.quadrature.local_2d(2, :, e).');
    Ep1 = dubiner_basis(k, part_1.quadrature.interface_2d(1, :, e).', part_1.quadrature.interface_2d(2, :, e).');
    Ep2 = dubiner_basis(k, part_2.quadrature.interface_2d(1, :, e).', part_2.quadrature.interface_2d(2, :, e).');

    % path evaluations, tq fastest within each patch node
    sp = reshape(part_1.quadrature.path_2d(1, :, e, :), sq, tq).';
    tp = reshape(part_1.quadrature.path_2d(2, :, e, :), sq, tq).';
    CSp1 = dubiner_basis(k, sp(:), tp(:));
    sp = reshape(part_2.quadrature.path_2d(1, :, e, :), sq, tq).';
    tp = reshape(part_2.quadrature.path_2d(2, :, e, :), sq, tq).';
    CSp2 = dubiner_basis(k, sp(:), tp(:));

    mx_i = m_i(1, wI).';  my_i = m_i(2, wI).';
    mx_o = m_o(1, wI).';  my_o = m_o(2, wI).';

    p1m_1 = reshape(sum(reshape(wCS .* CSp1 .* mx_i, [tq, d2*sq]), 1), [sq, d2]);
    p1m_2 = reshape(sum(reshape(wCS .* CSp1 .* my_i, [tq, d2*sq]), 1), [sq, d2]);

    x_path = reshape(part_2.quadrature.path(1, :, e, :), sq, tq).';
    y_path = reshape(part_2.quadrature.path(2, :, e, :), sq, tq).';
    kI = DarcyProblem.get_inverse_tensor(kappa, x_path, y_path);
    [kI11, kI21, kI12, kI22] = deal(kI{:});
    Kp2m_1 = reshape(sum(reshape(wCS .* CSp2 .* (kI11(:).*mx_o + kI21(:).*my_o), [tq, d2*sq]), 1), [sq, d2]);
    Kp2m_2 = reshape(sum(reshape(wCS .* CSp2 .* (kI12(:).*mx_o + kI22(:).*my_o), [tq, d2*sq]), 1), [sq, d2]);

    I_CS_1 = len_cs .* reshape(sum(reshape(wCS .* CSp1, [tq, d2*sq]), 1), [sq, d2]);

    O = zeros(sq, d2);

    %% T1 rows (2d1, tested with mu_1 on side 1)
    tau = part_1.tau(e);
    patch_len = abs(part_1.area(e));
    wphi = wP .* l1;
    test_inner = @(C) matfun(@(M) patch_len * (wphi.' * M), C);

    L_term = nu*test_inner({Ep1.*nS1v, O, Ep1.*nS2v, O;
        O, Ep1.*nS1v, O, Ep1.*nS2v}) ...
        + wbjs*test_inner({p1m_1.*(tS1v.*tS1v), p1m_1.*(tS2v.*tS1v), p1m_2.*(tS1v.*tS1v), p1m_2.*(tS2v.*tS1v);
        p1m_1.*(tS1v.*tS2v), p1m_1.*(tS2v.*tS2v), p1m_2.*(tS1v.*tS2v), p1m_2.*(tS2v.*tS2v)});
    u_term = -nu*tau*test_inner({P1f, O; O, P1f});
    p_term = -test_inner({Ep1.*nS1v; Ep1.*nS2v});
    a_term = test_inner({O(:, 1); O(:, 1)});
    R1_inner = [L_term, u_term, p_term, a_term];
    V = -R1_inner * S12K1;
    rows_acc{end+1} = repmat(DOF1, 1, size(V, 2));  cols_acc{end+1} = repmat(foot1.', 2*d1, 1);  vals_acc{end+1} = V; %#ok<AGROW>
    rhs_acc(DOF1) = rhs_acc(DOF1) - R1_inner * S1.S1f(:, :, 1, 1, K1);

    R1_outer = [test_inner({Kp2m_1.*nD1v, Kp2m_2.*nD1v; Kp2m_1.*nD2v, Kp2m_2.*nD2v}), ...
        test_inner({O; O})];
    V = -R1_outer * S12K2;
    rows_acc{end+1} = repmat(DOF1, 1, size(V, 2));  cols_acc{end+1} = repmat(foot2.', 2*d1, 1);  vals_acc{end+1} = V; %#ok<AGROW>
    rhs_acc(DOF1) = rhs_acc(DOF1) - R1_outer * S2.S1f(:, :, 1, 1, K2);

    M = test_inner({(nu*tau + wbjs*(tS1v.*tS1v)).*l1, (wbjs*(tS1v.*tS2v)).*l1;
        (wbjs*(tS2v.*tS1v)).*l1, (nu*tau + wbjs*(tS2v.*tS2v)).*l1});
    rows_acc{end+1} = repmat(DOF1, 1, 2*d1);  cols_acc{end+1} = repmat(DOF1.', 2*d1, 1);  vals_acc{end+1} = M; %#ok<AGROW>

    M = -test_inner({l2.*nD1v; l2.*nD2v});
    rows_acc{end+1} = repmat(DOF1, 1, d1);  cols_acc{end+1} = repmat(DOF2.', 2*d1, 1);  vals_acc{end+1} = M; %#ok<AGROW>

    M = -test_inner({P1f(:, 1).*nS1v; P1f(:, 1).*nS2v});
    rows_acc{end+1} = DOF1;  cols_acc{end+1} = a_dof*ones(2*d1, 1);  vals_acc{end+1} = M; %#ok<AGROW>

    T1v = {reshape(T1_data{1}(part_1.quadrature.interface(1, :, e), part_1.quadrature.interface(2, :, e)) ...
        + 0*part_1.quadrature.interface(1, :, e), [], 1); ...
        reshape(T1_data{2}(part_1.quadrature.interface(1, :, e), part_1.quadrature.interface(2, :, e)) ...
        + 0*part_1.quadrature.interface(1, :, e), [], 1)};
    rhs_acc(DOF1) = rhs_acc(DOF1) + test_inner(T1v);

    %% T2 rows (d1, tested with mu_2 on side 2)
    tau = part_2.tau(e);
    patch_len = abs(part_2.area(e));
    wphi = wP .* l2;
    test_outer = @(C) matfun(@(M) patch_len * (wphi.' * M), C);

    R2_inner = [test_outer({p1m_1.*nS1v, p1m_1.*nS2v, p1m_2.*nS1v, p1m_2.*nS2v}), ...
        test_outer({O, O}), test_outer({O}), test_outer({O(:, 1)})];
    V = -R2_inner * S12K1;
    rows_acc{end+1} = repmat(DOF2, 1, size(V, 2));  cols_acc{end+1} = repmat(foot1.', d1, 1);  vals_acc{end+1} = V; %#ok<AGROW>
    rhs_acc(DOF2) = rhs_acc(DOF2) - R2_inner * S1.S1f(:, :, 1, 1, K1);

    R2_outer = [test_outer({Ep2.*nD1v, Ep2.*nD2v}), tau*test_outer({P2f})];
    V = -R2_outer * S12K2;
    rows_acc{end+1} = repmat(DOF2, 1, size(V, 2));  cols_acc{end+1} = repmat(foot2.', d1, 1);  vals_acc{end+1} = V; %#ok<AGROW>
    rhs_acc(DOF2) = rhs_acc(DOF2) - R2_outer * S2.S1f(:, :, 1, 1, K2);

    M = test_outer({l1.*nS1v, l1.*nS2v});
    rows_acc{end+1} = repmat(DOF2, 1, 2*d1);  cols_acc{end+1} = repmat(DOF1.', d1, 1);  vals_acc{end+1} = M; %#ok<AGROW>

    M = -tau*test_outer({l2});
    rows_acc{end+1} = repmat(DOF2, 1, d1);  cols_acc{end+1} = repmat(DOF2.', d1, 1);  vals_acc{end+1} = M; %#ok<AGROW>

    T2v = reshape(T2_data(part_2.quadrature.interface(1, :, e), part_2.quadrature.interface(2, :, e)) ...
        + 0*part_2.quadrature.interface(1, :, e), [], 1);
    rhs_acc(DOF2) = rhs_acc(DOF2) + test_outer({T2v});

    %% alpha row (extrapolated pressure mean over the sliver)
    patch_len = part_1.area(e);
    wphi = wP .* P1f(:, 1);
    test_alpha = @(C) matfun(@(M) patch_len * (wphi.' * M), C);
    R_a = [test_alpha({O, O, O, O}), test_alpha({O, O}), test_alpha({I_CS_1}), test_alpha({O(:, 1)})];
    V = -R_a * S12K1;
    rows_acc{end+1} = a_dof*ones(1, size(V, 2));  cols_acc{end+1} = foot1.';  vals_acc{end+1} = V; %#ok<AGROW>
    rhs_acc(a_dof) = rhs_acc(a_dof) - R_a * S1.S1f(:, :, 1, 1, K1);

    rows_acc{end+1} = a_dof;  cols_acc{end+1} = a_dof; %#ok<AGROW>
    vals_acc{end+1} = measure_stokes * test_alpha({P1f(:, 1)}); %#ok<AGROW>
end

rows = cellfun(@(v) v(:), rows_acc.', 'UniformOutput', false);
cols = cellfun(@(v) v(:), cols_acc.', 'UniformOutput', false);
vals = cellfun(@(v) v(:), vals_acc.', 'UniformOutput', false);
Hc = sparse(vertcat(rows{:}), vertcat(cols{:}), vertcat(vals{:}), n_dof, n_dof);
Fc = rhs_acc;
end

%% ======================================================================
function F = fsi_functionals(obj, k, mesh_1, mesh_2, part_1, part_2, off, sq, tq, curve_fluid)
% The iterate-independent part of the fluid-structure transmission rows
% (see fsi_rows): per patch, the row functionals on the local unknowns,
% the trace (M) blocks, the interface source terms and the index vectors.
% None of it involves the condensed systems, so it can be kept across
% solves on the same meshes and coupling. curve_fluid is the fluid side's
% curve (the interface itself, or the first of a per-side pair); with a
% pair the fluid terms of T1 carry the Nanson factor (InterfaceProblem.fsi).
d1 = k + 1;
% local solid fields after (sigma, u, gamma) (the pressure split's xi, an
% Eulerian solid's displacement) are absent from the transmission
% functionals: pad them with zero columns
dims_2 = arrayfun(@(sp) sp.Dimension, obj.problems(2).LocalSpaces);
n_xi2 = sum(dims_2(4:end));
mu_f = obj.coupling.viscosity;
c1_1 = 1/(2*mu_f);
c2_1 = 1/(4*mu_f);
[lambda_s, mu_s] = LinearElasticityProblem.get_lame_parameters( ...
    obj.coupling.young_modulus, obj.coupling.poisson_ratio);
c1_2 = 1/(2*mu_s);
c2_2 = lambda_s/(2*mu_s*(2*lambda_s + 2*mu_s));
measure_stokes = obj.coupling.domain_measure;
cdt = obj.coupling.cdt;
% the kinematic row <cdt*u_F - c_S*u_S, mu>; c_S = 0 only for the
% stationary limit time_step = Inf (InterfaceProblem.fsi), <u_F, mu> = <T2, mu>
c_S = obj.coupling.solid_weight;
T1_data = obj.coupling.source_neumann;
T2_data = obj.coupling.source_dirichlet;
[~, wS]  = lgwt(sq, 0, 1);
[~, wTP] = lgwt(tq, 0, 1);
n_p = size(part_1.faces, 3);

tr1 = obj.problems(1).GlobalSpaces([obj.problems(1).GlobalSpaces.Support] == "face");
tr2 = transmitted_trace(obj.problems(2));

G_fluid = get_affine_map(mesh_1, part_1.elements, inverse = true);
D_A = dubiner_basis(0, 0.0, 0.0);

% signed sliver area element for the alpha row (zeroed when the row is
% eliminated or the correction is disabled)
if obj.coupling.pin_extra || ~obj.coupling.sliver_correction
    sliver_jacobian = zeros(sq, tq, n_p);
else
    sliver_jacobian = curve_fluid.get_sliver_jacobian(mesh_1, part_1, ...
        patch_weight = obj.coupling.patch_weight, ...
        quadrature_patch = sq, quadrature_path = tq);
end

% interface history H0 on the solid interface-face trace dofs
h0 = obj.coupling.history_skeleton;
if isscalar(h0)
    h0 = h0 + zeros(2*d1, n_p);
else
    h0 = h0(:, part_2.faces(:));
end
% A HyperelasticProblem solid at a Newton iterate (Data.iterate) under
% the no-slip law: its path gradient in the kinematic row is the discrete
% displacement gradient L_h of (u_S, u_hat_S) (discrete_gradient_operators),
% linear in the unknowns, since grad u = zeta^{-1}(sigma - N(L_h)) + gamma
% differs from the linear compliance; the history h0 is transported with
% L_h of its Newmark combination (solid_gradient_history), on the patch
% points, in place of Ext(h0)
% The same path gradient serves an Eulerian velocity solid
% (Data.eulerian_blocks, see FluidStructureCoupling.interface_problem): L_h
% of its velocity unknown (field 2) and the velocity skeleton.
solid_data = obj.problems(2).Data;
solid_gradient = obj.coupling.interface_law == "noslip" && isstruct(solid_data) ...
    && ((isfield(solid_data, "iterate") && ~isempty(solid_data.iterate)) ...
    || (isfield(solid_data, "eulerian_blocks") && ~isempty(solid_data.eulerian_blocks)));
if solid_gradient
    [R_u, R_hat] = HyperelasticProblem.discrete_gradient_operators(mesh_2, k);
    G_solid = get_affine_map(mesh_2, part_2.elements, inverse = true);
    history_gradient = [];
    if isempty(obj.coupling.history_extension)
        history_gradient = solid_gradient_history(obj, k, mesh_2, R_u, R_hat);
    end
end
% its transport onto Sigma, Ext(h0): given, or built from the solid
% stress history on this partition
history_extension = obj.coupling.history_extension;
if isempty(history_extension) && ~isempty(obj.coupling.history_stress) && ~solid_gradient
    history_extension = history_transport(obj, k, mesh_2, part_2);
end
if ~isempty(history_extension)
    % per patch, like h0 after the gather above: a face may carry
    % several patches
    assert(isequal(size(history_extension), size(h0)), ...
        "InterfaceProblem:badHistoryExtension", ...
        "history_extension must be %d x %d, got %s.", ...
        size(h0, 1), size(h0, 2), mat2str(size(history_extension)));
    h0 = h0 + history_extension;
end

% interface sources at each side's own interface quadrature points
evaluate_on_points = @(f, x) cellfun(@(f) ...
    reshape(f(x(1, :, :), x(2, :, :)) + 0*x(1, :, :), sq, 1, 1, []), ...
    f, 'UniformOutput', 0);
T1 = evaluate_on_points(T1_data, part_1.quadrature.interface);
T1 = cat(3, T1{:});
T2 = evaluate_on_points(T2_data, part_2.quadrature.interface);
T2 = cat(3, T2{:});
% Interface law: "noslip" is the classical coupling; "navier"/"full"
% replace the tangential half of the kinematic row by the fluid Navier
% law of Cesik, Kampschulte & Schwarzacher (arXiv:2603.12030), a = 0
% being full slip.
is_slip = obj.coupling.interface_law ~= "noslip";
a_slip = obj.coupling.slip_coefficient;
% The traction balance is the Cauchy traction plus the HDG stabilization,
% also for a Navier-Stokes fluid: Sigma is a material surface, so no
% momentum is advected through it, and the convective part of the fluid's
% interior-face numerical flux does not belong in this row.
Gs = 0;
if is_slip
    Gs = evaluate_on_points(obj.coupling.source_slip, part_2.quadrature.interface);
    Gs = cat(3, Gs{:});
end
% Nanson factor J(s) = |x'(s)|/|X'(s)| of a per-side curve pair (fluid
% over solid parametric speed at the shared patch nodes); [] for one curve
nanson = [];
if isfield(part_1.quadrature, "speed")
    if is_slip
        error("InterfaceProblem:curvePairSlip", ...
            "A pair of interface curves supports interface_law = ""noslip"" only.");
    end
    nanson = part_1.quadrature.speed ./ part_2.quadrature.speed;   % 1 x sq x n_p
end

patches = cell(1, n_p);
for e = 1:n_p
    K1 = part_1.elements(e);  F1 = part_1.faces(e);
    K2 = part_2.elements(e);  F2 = part_2.faces(e);
    fd1 = get_global_dof(tr1, subset=F1);  DOF1 = fd1(:);
    fd2 = get_global_dof(tr2, subset=F2);  DOF2 = fd2(:) + off;
    G_inner = G_fluid(:, :, e);

    n_interface = cat(4, ...
        reshape(part_1.quadrature.normals(:, :, e), 2, 1, []), ...
        reshape(part_2.quadrature.normals(:, :, e), 2, 1, []));
    nF = -n_interface(:, :, :, 1);
    nS = n_interface(:, :, :, 2);

    %% polynomials
    s_fluid = part_1.quadrature.local(:, :, e);
    l1 = legendre_basis(k, s_fluid);
    L1 = vector_legendre_basis(k, s_fluid);
    s_solid = part_2.quadrature.local(:, :, e);
    l2 = legendre_basis(k, s_solid);
    L2 = vector_legendre_basis(k, s_solid);

    D1 = vector_dubiner_basis(k, part_1.quadrature.local_2d(1, :, e), ...
        part_1.quadrature.local_2d(2, :, e));
    D2 = vector_dubiner_basis(k, part_2.quadrature.local_2d(1, :, e), ...
        part_2.quadrature.local_2d(2, :, e));

    s_conn_inner = reshape(part_1.quadrature.path_2d(1, :, e, :), sq, tq).';
    t_conn_inner = reshape(part_1.quadrature.path_2d(2, :, e, :), sq, tq).';
    TP_P1 = PEERS_basis(k, s_conn_inner(:), t_conn_inner(:), G_inner);
    TP_P1 = reshape(TP_P1, [tq sq size(TP_P1, 2:4)]);
    TP_A1 = asym_tensor_dubiner_basis(k, s_conn_inner(:), t_conn_inner(:));
    TP_A1 = reshape(TP_A1, [tq sq size(TP_A1, 2:4)]);

    % the solid stress is evaluated with the problem's own declared basis
    % (it coincides with PEERS_basis, used on the fluid side)
    solid_sigma_basis = obj.problems(2).LocalSpaces(1).Basis;
    s_conn_outer = reshape(part_2.quadrature.path_2d(1, :, e, :), sq, tq).';
    t_conn_outer = reshape(part_2.quadrature.path_2d(2, :, e, :), sq, tq).';
    TP_P2 = evaluate(solid_sigma_basis, [s_conn_outer(:).'; t_conn_outer(:).'], ...
        mesh_2, element_list = K2);
    TP_P2 = reshape(TP_P2, [tq sq size(TP_P2, 2:4)]);
    TP_A2 = asym_tensor_dubiner_basis(k, s_conn_outer(:), t_conn_outer(:));
    TP_A2 = reshape(TP_A2, [tq sq size(TP_A2, 2:4)]);

    s_iface_inner = part_1.quadrature.interface_2d(1, :, e);
    t_iface_inner = part_1.quadrature.interface_2d(2, :, e);
    E_P1 = PEERS_basis(k, s_iface_inner, t_iface_inner, G_inner);
    E_P1 = permute(E_P1, [3 4 1 2]);
    E_P1_N = permute(pagemtimes(E_P1, nF), [3 4 1 2]);
    D_A_N = permute(D_A * nF, [3 2 1]);
    E_A1 = asym_tensor_dubiner_basis(k, s_iface_inner, t_iface_inner);
    E_A1 = permute(E_A1, [3 4 1 2]);
    E_A1_N = permute(pagemtimes(E_A1, nF), [3 4 1 2]);

    s_iface_outer = part_2.quadrature.interface_2d(1, :, e);
    t_iface_outer = part_2.quadrature.interface_2d(2, :, e);
    E_P2 = evaluate(solid_sigma_basis, [s_iface_outer; t_iface_outer], ...
        mesh_2, element_list = K2);
    E_P2 = permute(E_P2, [3 4 1 2]);
    E_P2_N = permute(pagemtimes(E_P2, nS), [3 4 1 2]);
    E_A2 = asym_tensor_dubiner_basis(k, s_iface_outer, t_iface_outer);
    E_A2 = permute(E_A2, [3 4 1 2]);
    E_A2_N = permute(pagemtimes(E_A2, nS), [3 4 1 2]);

    % path integrals of the velocity/displacement reconstructions
    mF = reshape(part_1.quadrature.interface(:, :, e) ...
        - part_1.quadrature.physical(:, :, e), 2, 1, []);
    TP_P1 = permute(TP_P1, [4 5 2 1 3]);
    TP_P1_tr = TP_P1(1, 1, :, :, :) + TP_P1(2, 2, :, :, :);
    TP_P1_tr_Id = pagemtimes(TP_P1_tr, eye(2));
    TP_P1_grad_M = permute(pagemtimes(c1_1*TP_P1 - c2_1*TP_P1_tr_Id, mF), [3 4 5 1 2]);
    TP_INT_P1 = tensorprod(TP_P1_grad_M, wTP, 2, 1);
    TP_A1 = permute(TP_A1, [4 5 2 1 3]);
    TP_A1_M = permute(pagemtimes(TP_A1, mF), [3 4 5 1 2]);
    TP_INT_A1 = tensorprod(TP_A1_M, wTP, 2, 1);

    mS = reshape(part_2.quadrature.interface(:, :, e) ...
        - part_2.quadrature.physical(:, :, e), 2, 1, []);
    TP_P2 = permute(TP_P2, [4 5 2 1 3]);
    TP_P2_tr_Id = pagemtimes(TP_P2(1, 1, :, :, :) + TP_P2(2, 2, :, :, :), eye(2));
    TP_P2_grad_M = permute(pagemtimes(c1_2*TP_P2 - c2_2*TP_P2_tr_Id, mS), [3 4 5 1 2]);
    TP_INT_P2 = tensorprod(TP_P2_grad_M, wTP, 2, 1);
    TP_A2 = permute(TP_A2, [4 5 2 1 3]);
    TP_A2_M = permute(pagemtimes(TP_A2, mS), [3 4 5 1 2]);
    TP_INT_A2 = tensorprod(TP_A2_M, wTP, 2, 1);

    TP_P1_tr = reshape(TP_P1_tr, sq, tq, []);
    W_sliver = (wS(:) .* sliver_jacobian(:, :, e)) .* wTP(:)';
    EP_INT_INT_P1 = reshape(sum(TP_P1_tr .* W_sliver, [1 2]), 1, []);

    tau_1 = part_1.tau(e);
    tau_2 = part_2.tau(e);

    % swap_transmission: the traction balance is tested on side 2's faces
    % and the kinematic row on side 1's (default: T1 on side 1, T2 on 2)
    if obj.coupling.swap_transmission
        R1 = DOF2;  t1 = l2;  T1_len = part_2.area(e);
        R2 = DOF1;  t2 = l1;  T2_len = part_1.area(e);
    else
        R1 = DOF1;  t1 = l1;  T1_len = part_1.area(e);
        R2 = DOF2;  t2 = l2;  T2_len = part_2.area(e);
    end
    P = struct("K1", K1, "K2", K2, "DOF1", DOF1, "DOF2", DOF2, "R1", R1, "R2", R2);

    %% T1: <E(sigma_1)n_1 - mu tau u_1 + mu tau u_hat_1 + E(sigma_2)n_2 - tau u_2 + tau u_hat_2, mu_1>
    patch_len = T1_len;
    integrate = @(D, L) patch_len * reshape(permute(pagemtimes(wS.*L, 'transpose', D, 'none'), [1 3 2]), 2*d1, []);
    % test functions of the fluid terms: weighted by the Nanson factor
    % under a curve pair
    t1F = t1;
    if ~isempty(nanson)
        t1F = reshape(nanson(1, :, e), [], 1) .* t1;
    end

    sigma_term = integrate(E_P1_N, t1F);
    u_term = -tau_1*mu_f*integrate(D1, t1F);
    gamma_term = integrate(0*E_A1_N, t1);
    xi_term = integrate(0*D1(:, 1, :), t1);
    P.S1_inner = [sigma_term, u_term, gamma_term, xi_term];

    sigma_term = integrate(E_P2_N, t1);
    u_term = -tau_2*integrate(D2, t1);
    gamma_term = integrate(0*E_A2_N, t1);
    P.S1_outer = [sigma_term, u_term, gamma_term, zeros(size(sigma_term, 1), n_xi2)];

    % the fluid u_hat weight
    w_hat1 = tau_1*mu_f;
    P.M_hat_1 = integrate(L1, w_hat1 .* t1F);
    P.M_hat_2 = tau_2*integrate(L2, t1);
    P.M_alpha = integrate(D_A_N, t1F);

    % T1 is the full traction balance: the interface force is internal, so
    % the traction jump is zero and Navier slip prescribes only the fluid
    % traction.

    P.T1_term = integrate(T1(:, :, :, e), t1);

    %% alpha row: sliver-weighted mean of the extrapolated tr(sigma_1)
    sigma_term = D_A'*EP_INT_INT_P1;
    u_term = zeros(1, size(D1, 2));
    gamma_term = zeros(1, size(E_A1_N, 2));
    xi_term = 0;
    P.SA_inner = [sigma_term, u_term, gamma_term, xi_term];

    %% T2: <cdt*u_F - c_S*u_S, mu_2> on the transferred reconstructions
    %  interface_law "navier"/"full" (Cesik-Kampschulte-Schwarzacher (1.4)):
    %  only the normal component is coupled, and the freed d1 rows carry the
    %  fluid's own tangential traction law (their (1.5) with the balance).
    patch_len = part_2.area(e);
    % scalar projections of a vector-valued integrand onto a direction field
    project = @(D, L, v) patch_len * ((wS(:).*L).' * ...
        (reshape(v(1, 1, :), [], 1).*reshape(D(:, :, 1), size(D, 1), []) ...
       + reshape(v(2, 1, :), [], 1).*reshape(D(:, :, 2), size(D, 1), [])));
    % nS is 2 x 1 x nq; normalize over dim 1 to keep that shape
    nS_hat = nS ./ vecnorm(nS, 2, 1);
    tS_hat = [-nS_hat(2, 1, :); nS_hat(1, 1, :)];

    if is_slip
        % --- (1.4) impermeability: the normal component only, d1 rows ---
        P.DOF2n = DOF2(1:d1);
        proj_n = @(D, L) project(D, L, nS_hat);
        P.S2_inner = cdt*[proj_n(TP_INT_P1, l2), proj_n(0*D1, l2), ...
                          proj_n(TP_INT_A1, l2), proj_n(0*D1(:, 1, :), l2)];
        P.S2_outer = -c_S*[proj_n(TP_INT_P2, l2), proj_n(0*D2, l2), proj_n(TP_INT_A2, l2), zeros(d1, n_xi2)];
        P.M_kin_1 = cdt*proj_n(L1, l2);
        P.M_kin_2 = -c_S*proj_n(L2, l2);
        P.h0_term = c_S*proj_n(L2, l2)*h0(:, e);
        P.T2_term = cdt*proj_n(T2(:, :, :, e), l2);

        % --- the freed d1 rows: the fluid's tangential Navier law (3.11),
        %     t.(sigma_F n_F) + a t.(v - u_S) = 0, with the velocity
        %     difference taken on the transferred reconstructions of the
        %     kinematic row. a = 0 is free slip; a -> infinity recovers the
        %     no-slip coupling. The interface force is internal, so T1 is
        %     unchanged and the slip law prescribes only the fluid traction
        %     (reading (1.5) as a traction-jump law instead double-counts
        %     the friction on the solid and breaks the no-slip limit).
        P.DOF2t = DOF2(d1 + (1:d1));
        proj_t = @(D, L) project(D, L, tS_hat);
        S1t_inner = [proj_t(E_P1_N, l2), -tau_1*mu_f*proj_t(D1, l2), ...
                     proj_t(0*E_A1_N, l2), proj_t(0*D1(:, 1, :), l2)];
        T2t_inner = a_slip*cdt*[proj_t(TP_INT_P1, l2), proj_t(0*D1, l2), ...
                                proj_t(TP_INT_A1, l2), proj_t(0*D1(:, 1, :), l2)];
        P.St_inner = S1t_inner + T2t_inner;
        P.T2t_outer = -a_slip*c_S*[proj_t(TP_INT_P2, l2), proj_t(0*D2, l2), proj_t(TP_INT_A2, l2), zeros(d1, n_xi2)];

        % the same fluid-flux weight as the T1 row
        P.M_slip_1 = proj_t(L1, w_hat1 .* l2) + a_slip*cdt*proj_t(L1, l2);
        P.M_slip_2 = -a_slip*c_S*proj_t(L2, l2);
        P.slip_term = a_slip*(-c_S*proj_t(L2, l2)*h0(:, e) ...
            + cdt*proj_t(T2(:, :, :, e), l2));
        P.Gs_term = proj_t(Gs(:, :, :, e), l2);
    else
        patch_len = T2_len;
        integrate = @(D, L) patch_len * reshape(permute(pagemtimes(wS.*L, 'transpose', D, 'none'), [1 3 2]), 2*d1, []);
        sigma_term = integrate(TP_INT_P1, t2);
        u_term = integrate(0*D1, t2);
        gamma_term = integrate(TP_INT_A1, t2);
        xi_term = integrate(0*D1(:, 1, :), t2);
        P.S2_inner = cdt*[sigma_term, u_term, gamma_term, xi_term];

        if solid_gradient
            % int_path L_h m: (sq x columns x component) on u_K and u_hat_K
            TP_L = reshape(PEERS_basis(k, s_conn_outer(:), t_conn_outer(:), G_solid(:, :, e)), ...
                tq, sq, [], 2, 2);
            TP_L = reshape(tensorprod(wTP(:), TP_L, 1, 1), sq, [], 2, 2);
            TP_L_M = TP_L(:, :, :, 1) .* reshape(mS(1, 1, :), [], 1) ...
                + TP_L(:, :, :, 2) .* reshape(mS(2, 1, :), [], 1);
            sigma_term = zeros(2*d1, size(TP_INT_P2, 2));
            u_term = -integrate(pagemtimes(TP_L_M, R_u(:, :, K2)), t2);
            gamma_term = zeros(2*d1, size(TP_INT_A2, 2));
            % the u_hat_K columns: the faces of K2 (element_traces order)
            face_dofs = arrayfun(@(f) reshape(get_global_dof(tr2, subset = f), [], 2), ...
                mesh_2.boundaries(:, K2), "UniformOutput", false);
            face_dofs = cat(1, face_dofs{:});
            P.DOF2_hat = face_dofs(:) + off;
            P.S2_hat = -c_S*integrate(pagemtimes(TP_L_M, R_hat(:, :, K2)), t2);
        else
            sigma_term = -integrate(TP_INT_P2, t2);
            u_term = -integrate(0*D2, t2);
            gamma_term = -integrate(TP_INT_A2, t2);
        end
        P.S2_outer = c_S*[sigma_term, u_term, gamma_term, zeros(size(sigma_term, 1), n_xi2)];

        P.M_kin_1 = cdt*integrate(L1, t2);
        P.M_kin_2 = -c_S*integrate(L2, t2);
        P.h0_term = c_S*integrate(L2, t2)*h0(:, e);
        if solid_gradient && ~isempty(history_gradient)
            P.h0_term = P.h0_term ...
                + c_S*integrate(pagemtimes(TP_L_M, history_gradient(:, K2)), t2);
        end
        P.T2_term = cdt*integrate(T2(:, :, :, e), t2);
    end
    patches{e} = P;
end

% mean-value diagonal of the alpha unknown (the zero-mean rho row is
% provided by the fluid factory's constrained closure blocks)
F = struct("is_slip", is_slip, "patches", {patches}, ...
    "alpha_diagonal", 2 * measure_stokes * (D_A') * D_A);
end

%% ======================================================================
function ext = history_transport(obj, k, mesh_2, part_2)
% Ext(h0) of the fsi T2 history, (2*(k+1)) x n_patches, on the solid
% side of the partition: trace_extension_coefficients of the history's
% (sigma, gamma) (InterfaceProblem.fsi; a scalar 0 is a zero field),
% with the coupling's solid coefficients.
n_el = size(mesh_2.elements, 2);
sigma_h0 = obj.coupling.history_stress + zeros(PEERSBasis(k).Dimension, n_el);
gamma_h0 = obj.coupling.history_rotation + zeros(nchoosek(k + 2, 2), n_el);
solid = LinearElasticityProblem(k, young_modulus = obj.coupling.young_modulus, ...
    poisson_ratio = obj.coupling.poisson_ratio);
ext = trace_extension_coefficients(solid, mesh_2, part_2, sigma_h0, gamma_h0, Order = k);
end

%% ======================================================================
function l = solid_gradient_history(obj, k, mesh_2, R_u, R_hat)
% L_h dofs (d2peers x n_elements) of the T2 history of a HyperelasticProblem
% solid, the Newmark combination u0 + c_v v0 + c_a a0 of the volume and
% skeleton displacements of its Data.state.solid (beta, gamma of its
% Data; dt = cdt*gamma/beta, InterfaceProblem.fsi); [] without a history
% (steady row, or a zero history and no state)
has_history = ~isempty(obj.coupling.history_stress) || ~isequal(obj.coupling.history_skeleton, 0);
data = obj.problems(2).Data;
has_state = isfield(data, "state") && ~isempty(data.state) && isfield(data.state, "solid") ...
    && isfield(data.state.solid, "displacement");
l = [];
if ~has_state
    assert(~has_history, "InterfaceProblem:hyperelasticHistory", ...
        "A HyperelasticProblem solid with a T2 history needs its state in Data.state.");
    return
end
if ~has_history, return; end
s = data.state.solid;
beta_nm = data.newmark_beta;
gamma_nm = data.newmark_gamma;
dt = obj.coupling.cdt*gamma_nm/beta_nm;
c_v = dt*(1 - beta_nm/gamma_nm);
c_a = dt^2*((1/2 - beta_nm) - beta_nm*(1 - gamma_nm)/gamma_nm);
combine = @(name) newmark_field(s, name) + c_v*newmark_field(s, replace(name, "displacement", "velocity")) ...
    + c_a*newmark_field(s, replace(name, "displacement", "acceleration"));
u = combine("displacement");
u_hat = combine("displacement_skeleton");
n_el = size(mesh_2.elements, 2);
d1 = k + 1;
u = reshape(u + zeros(2*nchoosek(k + 2, 2), n_el), [], 1, n_el);
u_hat = u_hat + zeros(2*d1, size(mesh_2.faces, 2));
% the trace dofs of each element's faces (x block over y block)
u_K = reshape([reshape(u_hat(1:d1, mesh_2.boundaries), 3*d1, n_el);
               reshape(u_hat(d1 + 1:end, mesh_2.boundaries), 3*d1, n_el)], [], 1, n_el);
l = reshape(pagemtimes(R_u, u) + pagemtimes(R_hat, u_K), [], n_el);
end

function A = newmark_field(s, name)
% a state.solid dof array (cells stacked), 0 when absent or empty
A = 0;
if isfield(s, name) && ~isempty(s.(name))
    A = s.(name);
    if iscell(A), A = cat(1, A{:}); end
end
end

%% ======================================================================
function [Hc, Fc] = fsi_rows(k, S1, S2, off, n_dof, extra_base, F)
% Fluid-structure transmission: side 1 = PEERS stress-velocity fluid
% (local [sigma(dS); u(2d2); gamma(d2); xi]), side 2 = PEERS linear
% elasticity (local [sigma(dS); u(2d2); gamma(d2)]), plus the fluid
% pressure-level scalar alpha at extra_base+1 (sliver-weighted mean row).
% Local solves are x = S1f - S12 x_hat, so a functional R on the local
% unknowns enters as -R*S12 in the matrix and -R*S1f in the rhs.
% The functionals, trace blocks and sources come from fsi_functionals (F);
% this function forms only their products with S1, S2 and assembles.
% Closed-loop anchors (circulation, rigid modes, traction resultants) are
% not supported.
d1 = k + 1;
n_p = numel(F.patches);
a_dof = extra_base + 1;

nt1 = numel(S1.S12);  nt2 = numel(S2.S12);
S12f1 = cell(1, nt1);
for t = 1:nt1
    S12f1{t} = reshape(permute(S1.S12{t}, [1 2 7 5 3 4 6]), ...
        size(S1.S12{t}, 1), S1.d_g(t)*S1.f_g(t), []);
end
S12f2 = cell(1, nt2);
for t = 1:nt2
    S12f2{t} = reshape(permute(S2.S12{t}, [1 2 7 5 3 4 6]), ...
        size(S2.S12{t}, 1), S2.d_g(t)*S2.f_g(t), []);
end

[rows_acc, cols_acc, vals_acc] = deal({});
rhs_acc = zeros(n_dof, 1);

for e = 1:n_p
    P = F.patches{e};
    K1 = P.K1;  K2 = P.K2;
    DOF1 = P.DOF1;  DOF2 = P.DOF2;
    R1 = P.R1;  R2 = P.R2;
    S12K1 = cell2mat(cellfun(@(A) A(:, :, K1), S12f1, 'UniformOutput', false));
    S12K2 = cell2mat(cellfun(@(A) A(:, :, K2), S12f2, 'UniformOutput', false));
    foot1 = cell2mat(cellfun(@(fp) fp(:, K1), S1.L.footprint(:).', 'UniformOutput', false).');
    foot2 = cell2mat(cellfun(@(fp) fp(:, K2), S2.L.footprint(:).', 'UniformOutput', false).') + off;

    %% T1: <E(sigma_1)n_1 - mu tau u_1 + mu tau u_hat_1 + E(sigma_2)n_2 - tau u_2 + tau u_hat_2, mu_1>
    S1_inner = P.S1_inner;
    V = -S1_inner * S12K1;
    rows_acc{end+1} = repmat(R1, 1, size(V, 2));  cols_acc{end+1} = repmat(foot1.', 2*d1, 1);  vals_acc{end+1} = V; %#ok<AGROW>
    rhs_acc(R1) = rhs_acc(R1) - S1_inner * S1.S1f(:, :, 1, 1, K1);

    S1_outer = P.S1_outer;
    V = -S1_outer * S12K2;
    rows_acc{end+1} = repmat(R1, 1, size(V, 2));  cols_acc{end+1} = repmat(foot2.', 2*d1, 1);  vals_acc{end+1} = V; %#ok<AGROW>
    rhs_acc(R1) = rhs_acc(R1) - S1_outer * S2.S1f(:, :, 1, 1, K2);

    % the fluid u_hat weight
    M = P.M_hat_1;
    rows_acc{end+1} = repmat(R1, 1, 2*d1);  cols_acc{end+1} = repmat(DOF1.', 2*d1, 1);  vals_acc{end+1} = M; %#ok<AGROW>

    M = P.M_hat_2;
    rows_acc{end+1} = repmat(R1, 1, 2*d1);  cols_acc{end+1} = repmat(DOF2.', 2*d1, 1);  vals_acc{end+1} = M; %#ok<AGROW>

    M = P.M_alpha;
    rows_acc{end+1} = R1;  cols_acc{end+1} = a_dof*ones(2*d1, 1);  vals_acc{end+1} = M; %#ok<AGROW>

    rhs_acc(R1) = rhs_acc(R1) + P.T1_term;

    %% alpha row: sliver-weighted mean of the extrapolated tr(sigma_1)
    SA_inner = P.SA_inner;
    V = -SA_inner * S12K1;
    rows_acc{end+1} = a_dof*ones(1, size(V, 2));  cols_acc{end+1} = foot1.';  vals_acc{end+1} = V; %#ok<AGROW>
    rhs_acc(a_dof) = rhs_acc(a_dof) - SA_inner * S1.S1f(:, :, 1, 1, K1);

    %% T2: <cdt*u_F - c_S*u_S, mu_2> on the transferred reconstructions
    if F.is_slip
        % --- (1.4) impermeability: the normal component only, d1 rows ---
        DOF2n = P.DOF2n;
        S2_inner = P.S2_inner;
        V = -S2_inner * S12K1;
        rows_acc{end+1} = repmat(DOF2n, 1, size(V, 2));  cols_acc{end+1} = repmat(foot1.', d1, 1);  vals_acc{end+1} = V; %#ok<AGROW>
        rhs_acc(DOF2n) = rhs_acc(DOF2n) - S2_inner * S1.S1f(:, :, 1, 1, K1);

        S2_outer = P.S2_outer;
        V = -S2_outer * S12K2;
        rows_acc{end+1} = repmat(DOF2n, 1, size(V, 2));  cols_acc{end+1} = repmat(foot2.', d1, 1);  vals_acc{end+1} = V; %#ok<AGROW>
        rhs_acc(DOF2n) = rhs_acc(DOF2n) - S2_outer * S2.S1f(:, :, 1, 1, K2);

        M = P.M_kin_1;
        rows_acc{end+1} = repmat(DOF2n, 1, 2*d1);  cols_acc{end+1} = repmat(DOF1.', d1, 1);  vals_acc{end+1} = M; %#ok<AGROW>
        M = P.M_kin_2;
        rows_acc{end+1} = repmat(DOF2n, 1, 2*d1);  cols_acc{end+1} = repmat(DOF2.', d1, 1);  vals_acc{end+1} = M; %#ok<AGROW>
        rhs_acc(DOF2n) = rhs_acc(DOF2n) - P.h0_term + P.T2_term;

        % --- the freed d1 rows: the fluid's tangential Navier law ---
        DOF2t = P.DOF2t;
        V = -P.St_inner * S12K1;
        rows_acc{end+1} = repmat(DOF2t, 1, size(V, 2));  cols_acc{end+1} = repmat(foot1.', d1, 1);  vals_acc{end+1} = V; %#ok<AGROW>
        rhs_acc(DOF2t) = rhs_acc(DOF2t) - P.St_inner * S1.S1f(:, :, 1, 1, K1);

        T2t_outer = P.T2t_outer;
        V = -T2t_outer * S12K2;
        rows_acc{end+1} = repmat(DOF2t, 1, size(V, 2));  cols_acc{end+1} = repmat(foot2.', d1, 1);  vals_acc{end+1} = V; %#ok<AGROW>
        rhs_acc(DOF2t) = rhs_acc(DOF2t) - T2t_outer * S2.S1f(:, :, 1, 1, K2);

        % the same fluid-flux weight as the T1 row
        M = P.M_slip_1;
        rows_acc{end+1} = repmat(DOF2t, 1, 2*d1);  cols_acc{end+1} = repmat(DOF1.', d1, 1);  vals_acc{end+1} = M; %#ok<AGROW>
        M = P.M_slip_2;
        rows_acc{end+1} = repmat(DOF2t, 1, 2*d1);  cols_acc{end+1} = repmat(DOF2.', d1, 1);  vals_acc{end+1} = M; %#ok<AGROW>
        rhs_acc(DOF2t) = rhs_acc(DOF2t) + P.slip_term + P.Gs_term;
    else
        S2_inner = P.S2_inner;
        V = -S2_inner * S12K1;
        rows_acc{end+1} = repmat(R2, 1, size(V, 2));  cols_acc{end+1} = repmat(foot1.', 2*d1, 1);  vals_acc{end+1} = V; %#ok<AGROW>
        rhs_acc(R2) = rhs_acc(R2) - S2_inner * S1.S1f(:, :, 1, 1, K1);

        S2_outer = P.S2_outer;
        V = -S2_outer * S12K2;
        rows_acc{end+1} = repmat(R2, 1, size(V, 2));  cols_acc{end+1} = repmat(foot2.', 2*d1, 1);  vals_acc{end+1} = V; %#ok<AGROW>
        rhs_acc(R2) = rhs_acc(R2) - S2_outer * S2.S1f(:, :, 1, 1, K2);

        M = P.M_kin_1;
        rows_acc{end+1} = repmat(R2, 1, 2*d1);  cols_acc{end+1} = repmat(DOF1.', 2*d1, 1);  vals_acc{end+1} = M; %#ok<AGROW>

        M = P.M_kin_2;
        rows_acc{end+1} = repmat(R2, 1, 2*d1);  cols_acc{end+1} = repmat(DOF2.', 2*d1, 1);  vals_acc{end+1} = M; %#ok<AGROW>

        if isfield(P, "S2_hat")
            % the solid path gradient's trace columns (fsi_functionals)
            M = P.S2_hat;
            rows_acc{end+1} = repmat(R2, 1, size(M, 2));  cols_acc{end+1} = repmat(P.DOF2_hat.', 2*d1, 1);  vals_acc{end+1} = M; %#ok<AGROW>
        end

        rhs_acc(R2) = rhs_acc(R2) - P.h0_term + P.T2_term;
    end
end

% mean-value diagonal of the alpha unknown (the zero-mean rho row is
% provided by the fluid factory's constrained closure blocks)
rows_acc{end+1} = a_dof;  cols_acc{end+1} = a_dof;
vals_acc{end+1} = F.alpha_diagonal;

rows = cellfun(@(v) v(:), rows_acc.', 'UniformOutput', false);
cols = cellfun(@(v) v(:), cols_acc.', 'UniformOutput', false);
vals = cellfun(@(v) v(:), vals_acc.', 'UniformOutput', false);
Hc = sparse(vertcat(rows{:}), vertcat(cols{:}), vertcat(vals{:}), n_dof, n_dof);
Fc = rhs_acc;
end

%% ======================================================================
function tr = transmitted_trace(problem)
% the face space of side 2 that the fsi transmission rows act on: its first
% one (a solid may carry further face spaces, e.g. an Eulerian solid's
% displacement skeleton beside its velocity skeleton; their rows stay)
tr = problem.GlobalSpaces(find([problem.GlobalSpaces.Support] == "face", 1));
end

%% ======================================================================
function [x, g, Hh] = pad_side(n, x, g, Hh)
n_pad = n - numel(x);
if n_pad > 0
    x = [x; zeros(n_pad, 1)];  g = [g; zeros(n_pad, 1)];
    Hh = blkdiag(Hh, sparse(n_pad, n_pad));
end
end
