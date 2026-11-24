function [induced_partition_1, induced_partition_2, interface_partition] = ...
    get_partitions(interface, mesh_1, mesh_2, options)
%-------------------------------------------------------------------------
%Description: Builds the TPM induced partitions of two dissimilar meshes
%             against a common interface Sigma.
%
%             Each mesh carries a map phi_i from its discrete interface
%             (the faces tagged "Interface") onto Sigma. The construction
%             maintains three properties, each asserted below:
%
%               total       every interface face owns at least one patch;
%                           otherwise its trace dofs are unconstrained
%                           (a singular system, or a spurious Neumann
%                           condition imposed by omission)
%               complete    each face's patches tile its whole [0, 1], so
%                           the integral over that face is well defined
%               equal range both maps cover the same part of Sigma
%
%             Injectivity is not required in principle (a folded discrete
%             interface may cover a stretch of Sigma twice; see the
%             temporary restriction below), nor is surjectivity: Sigma may
%             have stretches covered by neither mesh.
%
%             Patches are the elementary intervals of the merged endpoint
%             list, each paired with its covering face on either side.
%
%             Collapsed faces (tip_mapping = "collapse", see
%             map_mesh_to_interface): a face whose vertices share one
%             closest point t* (a sharp tip or cusp of Sigma) is not
%             spread along Sigma but sent whole to xbar = gamma(t*). It
%             owns one extra patch, appended after the regular ones:
%               own side   the whole face, s in [0, 1] (its area is the
%                          face length), paths fanning from each point
%                          of the face to xbar
%               other side the single point paired with t*, s = [s*, s*]
%                          on the face covering t* (area 0)
%               Sigma      vertices_interface = [t*, t*]
%             A row integrated over the other side's patch measure (the
%             kinematic row of the default fsi coupling, tested on the
%             solid) gets nothing from it; a row integrated over the own
%             side's measure sees the whole face. At xbar the normal of
%             Sigma is undefined, so quadrature.normals on these patches
%             follows collapsed_normal:
%               "path"     (x - xbar)/|x - xbar|, x the own-side point of
%                          each path, oriented like the curve's normal
%               "face"     the own face's unit normal, oriented likewise
%               "bisector" the normalized mean of the curve's one-sided
%                          normals at t* +- 1e-6*period; an error where
%                          they cancel (a cusp)
%             The patch's .collapsed flag (1 x n_patches, on both
%             partitions) marks these patches and .collapsed_side
%             (1 x n_patches, 0 on regular patches) names the mesh whose
%             face collapsed. A collapsed face is normal-free only in the
%             row tested on its own side: InterfaceProblem.solve requires
%             the kinematic row there and refuses the other orientation.
%             "collapse" needs a periodic interface (see
%             map_mesh_to_interface).
%
%             Accuracy: "collapse" is a fallback that makes meshes solvable where
%             "spread" fails (orphanFace when the tied run straddles the parameter
%             seam, as on an exact cusp). Where "spread" works it is the more accurate
%             choice: at k = 2 on meshes with six tied faces at a cusp, collapse gave
%             errors about 4-5 times larger, and with collapse the error near a cusp
%             depends on the tie pattern and is not observed to converge for k >= 2.
%             Do not use it for convergence studies.
%
%             Per-side curves: with interface_2 set, mesh_1 is mapped onto
%             `interface` and mesh_2 onto interface_2, two curves on one
%             parameter t (equal t_range and periodicity). Patches are
%             still the merged intervals in t, so a patch quadrature node t
%             pairs the point interface(t) on side 1 with interface_2(t)
%             on side 2; each side's quadrature.interface, .normals and
%             transfer paths refer to its own curve. Each partition then
%             also carries quadrature.speed (1 x nq x n), the parametric
%             speed |gamma_j'(t)| of its curve at the patch nodes, so that
%             a coupling can form the ratio of the two area elements.
%-------------------------------------------------------------------------
arguments
    interface AbstractCurve
    mesh_1 SimplicialGrid
    mesh_2 SimplicialGrid
    % Absolute merge tolerance in the interface parameter. Two interface
    % vertices closer than this are treated as one breakpoint.
    options.vertex_tol double = 1e-8

    % Merge tolerance relative to the typical face, as a fraction in
    % (0, 1). When set it overrides vertex_tol with
    %       tol = vertex_tol_relative * median(segment length in t),
    % so the tolerance scales with h. Keep it well below 1: the tolerance
    % must not swallow a whole face, or the orphanFace assertion fires.
    options.vertex_tol_relative double = []

    options.quadrature_patch = 1;
    options.quadrature_path = 2;

    % Faces with a zero-length image (see the description and
    % map_mesh_to_interface): "spread" along Sigma, or "collapse" to the
    % common closest point with the normal rule collapsed_normal.
    % "collapse" is a fallback for meshes where "spread" fails; it is less
    % accurate near a cusp (see the description).
    options.tip_mapping (1,1) string ...
        {mustBeMember(options.tip_mapping, ["spread", "collapse"])} = "spread"
    options.collapsed_normal (1,1) string ...
        {mustBeMember(options.collapsed_normal, ["path", "face", "bisector"])} = "path"

    % Point-to-interface map of each mesh's vertices (map_mesh_to_interface
    % option `mapping`); "" keeps the curve's own. Inside a face the map is
    % affine in the interface parameter whatever the vertex rule, so a
    % patch's quadrature nodes, Gauss in t, are Gauss in either face-local
    % coordinate and pair the same points from either side.
    options.mapping_1 (1,1) string = ""
    options.mapping_2 (1,1) string = ""

    % Side 2's own curve (see the description); [] maps both meshes onto
    % `interface`.
    options.interface_2 = []
end
mesh_list = {mesh_1; mesh_2};
is_pair = ~isempty(options.interface_2);
if is_pair
    mustBeInterfaceCurves({interface, options.interface_2});
    if options.tip_mapping == "collapse"
        % collapsed patches take their normal from `interface` alone
        error("get_partitions:collapsePerSideCurves", ...
            "tip_mapping = ""collapse"" is not supported with interface_2.");
    end
    curves = {interface; options.interface_2};
else
    curves = {interface; interface};
end
segments = {map_mesh_to_interface(curves{1}, mesh_1, ...
    mapping = options.mapping_1, tip_mapping = options.tip_mapping); ...
            map_mesh_to_interface(curves{2}, mesh_2, ...
    mapping = options.mapping_2, tip_mapping = options.tip_mapping)};

%==================== MERGED BREAKPOINTS ALONG SIGMA =====================%
n_segments = cellfun(@(s) size(s.t, 2), segments).';

% Resolve the merge tolerance (see the argument block)
vertex_tol = options.vertex_tol;
if ~isempty(options.vertex_tol_relative)
    segment_length = abs([diff(segments{1}.t), diff(segments{2}.t)]);
    vertex_tol = options.vertex_tol_relative * median(segment_length);
end

% DataScale 1: uniquetol otherwise scales the tolerance by max(|t|), which
% loosens it on interfaces whose t_range is not [0, 1].
[breaks, ~, break_of] = uniquetol( ...
    [segments{1}.t(:); segments{2}.t(:)], vertex_tol, DataScale = 1);
breaks = reshape(breaks, 1, []);
break_of = mat2cell(reshape(break_of, 2, []), 2, n_segments);
n_intervals = numel(breaks) - 1;

% Snap each segment's endpoints onto the merged representatives. uniquetol
% reports a representative per cluster, not the cluster mean, so a segment
% endpoint can sit up to the merge tolerance away from its own breakpoint;
% evaluating the face-local coordinate there would land just outside
% [0, 1] and trip the assertions in get_1d_physical.
for j = 1:2
    segments{j}.t = reshape(breaks(break_of{j}), 2, []);
end

%=========================== COVERAGE PER MESH ===========================%
% owner{j}(i) is the segment of mesh j covering elementary interval i, or
% 0 where mesh j does not reach it.
[owner, n_covering] = deal(cell(2, 1));
for j = 1:2
    [owner{j}, n_covering{j}] = deal(zeros(1, n_intervals));
    for k = 1:n_segments(j)
        covered = break_of{j}(1, k) : break_of{j}(2, k) - 1;
        n_covering{j}(covered) = n_covering{j}(covered) + 1;
        owner{j}(covered) = k;
    end
end

% Temporary restriction: folded discrete interfaces (two faces of one
% mesh covering the same stretch of Sigma) are rejected until their
% patch weighting is settled.
for j = 1:2
    folded = find(n_covering{j} > 1, 1);
    if ~isempty(folded)
        error("get_partitions:foldedProjection", ...
            "Mesh %d folds over the interface: more than one interface " + ...
            "face covers t in [%g, %g]. Folded discrete interfaces are " + ...
            "not supported yet.", j, breaks(folded), breaks(folded + 1));
    end
end

reached = cellfun(@(c) c > 0, n_covering, 'UniformOutput', false);
lopsided = find(xor(reached{1}, reached{2}), 1);
if ~isempty(lopsided)
    error("get_partitions:unequalRange", ...
        "Only mesh %d reaches t in [%g, %g]. The two interface maps must " + ...
        "have equal range, otherwise the face on the other side loses " + ...
        "part of its integral.", ...
        1 + reached{2}(lopsided), breaks(lopsided), breaks(lopsided + 1));
end

patch_interval = find(reached{1} & reached{2});
n_of_patches = numel(patch_interval);
t_lo = breaks(patch_interval);
t_hi = breaks(patch_interval + 1);

%============================== PATCH DATA ===============================%
% Face-local coordinates are affine in the interface parameter across a
% segment, also when a face spans several merged vertices.
induced = cell(2, 1);
for j = 1:2
    k = owner{j}(patch_interval);
    seg_t = segments{j}.t(:, k);
    seg_s = segments{j}.s(:, k);
    local = @(t) seg_s(1, :) + (t - seg_t(1, :)) ./ ...
        (seg_t(2, :) - seg_t(1, :)) .* (seg_s(2, :) - seg_s(1, :));
    induced{j} = [local(t_lo); local(t_hi); segments{j}.face(k)];
end

% Collapsed faces: one appended patch each (see the description)
[t_lo, t_hi, collapsed] = append_collapsed_patches(interface, segments, ...
    induced, t_lo, t_hi);
[induced, collapsed_patch] = deal(collapsed.induced, collapsed.flag);
n_of_patches = numel(t_lo);

for j = 1:2
    orphans = setdiff(mesh_list{j}.intfaces(:).', induced{j}(3, :));
    if ~isempty(orphans)
        error("get_partitions:orphanFace", ...
            "Mesh %d: interface face(s) %s own no patch, so their trace " + ...
            "degrees of freedom would be left unconstrained.", ...
            j, mat2str(orphans));
    end
end

induced_partition_1 = enrich_partition(mesh_1, induced{1});
induced_partition_2 = enrich_partition(mesh_2, induced{2});

% On a closed interface the first and last breakpoints name the same
% physical point twice; deduplicating them would shift every patch index.
interface_partition.coordinates = interface.get_2d_points(breaks);
interface_partition.connectivity = [patch_interval; patch_interval + 1];
if any(collapsed_patch)
    % a collapsed patch is the single breakpoint nearest to t*
    t_star = interface.t_range(1) + mod(t_lo(collapsed_patch) - interface.t_range(1), ...
        diff(interface.t_range));
    [~, nearest] = min(abs(breaks(:) - t_star), [], 1);
    interface_partition.connectivity = [interface_partition.connectivity, ...
        [nearest; nearest]];
end

% Global interface-parameter endpoints of each patch, shared by both
% meshes so their quadrature rules evaluate identical interface points.
patch_endpoints_interface = reshape([t_lo; t_hi], 1, 2, n_of_patches);
induced_partition_1.vertices_interface = patch_endpoints_interface;
induced_partition_2.vertices_interface = patch_endpoints_interface;

get_quadrature = @(curve, mesh, partition) get_induced_partition_quadrature(curve, mesh, partition, "quadrature_patch", options.quadrature_patch, "quadrature_path", options.quadrature_path);

induced_partition_1.quadrature = get_quadrature(curves{1}, mesh_1, induced_partition_1);
induced_partition_2.quadrature = get_quadrature(curves{2}, mesh_2, induced_partition_2);

if is_pair
    % parametric speed at the patch nodes, which are Gauss in t (the
    % rule of get_induced_partition_quadrature)
    t_nodes = patch_endpoints_interface(:, 1, :) + lgwt(options.quadrature_patch, 0, 1).' ...
        .* (patch_endpoints_interface(:, 2, :) - patch_endpoints_interface(:, 1, :));
    speed = @(curve) reshape(vecnorm(curve.tangent(reshape(t_nodes, 1, [])), 2, 1), size(t_nodes));
    induced_partition_1.quadrature.speed = speed(curves{1});
    induced_partition_2.quadrature.speed = speed(curves{2});
end

if any(collapsed_patch)
    [induced_partition_1.quadrature.normals, induced_partition_2.quadrature.normals] = ...
        deal(collapsed_normals(interface, mesh_list, ...
        {induced_partition_1, induced_partition_2}, collapsed, options.collapsed_normal));
end
induced_partition_1.collapsed = collapsed_patch;
induced_partition_2.collapsed = collapsed_patch;
collapsed_side = zeros(1, n_of_patches);
collapsed_side(collapsed_patch) = collapsed.side;
induced_partition_1.collapsed_side = collapsed_side;
induced_partition_2.collapsed_side = collapsed_side;

to_boundary = @(mesh, partition) get_2d_local_boundary(mesh, reshape(partition.quadrature.local, [], n_of_patches), partition.elements, partition.side_number);
induced_partition_1.quadrature.local_2d = to_boundary(mesh_1, induced_partition_1);
induced_partition_2.quadrature.local_2d = to_boundary(mesh_2, induced_partition_2);

to_interface = @(mesh, partition) mesh.get_2d_local(partition.quadrature.interface, partition.elements);
induced_partition_1.quadrature.interface_2d = to_interface(mesh_1, induced_partition_1);
induced_partition_2.quadrature.interface_2d = to_interface(mesh_2, induced_partition_2);

to_path = @(mesh, partition) permute(reshape(mesh.get_2d_local(reshape(permute(partition.quadrature.path, [1 2 4 3]), 2, [], n_of_patches), partition.elements), 2, options.quadrature_patch, options.quadrature_path, n_of_patches), [1 2 4 3]);
induced_partition_1.quadrature.path_2d = to_path(mesh_1, induced_partition_1);
induced_partition_2.quadrature.path_2d = to_path(mesh_2, induced_partition_2);
end

function [t_lo, t_hi, collapsed] = append_collapsed_patches(interface, segments, induced, t_lo, t_hi)
% Appends one patch per collapsed face (either mesh): the whole face on
% its own side, the point paired with t* on the other. collapsed.flag
% marks the patches, collapsed.side the mesh whose face collapsed.
n_regular = numel(t_lo);
period = diff(interface.t_range);
side = zeros(1, 0);
for j = 1:2
    other = 3 - j;
    for c = 1:numel(segments{j}.collapsed_face)
        t_star = segments{j}.collapsed_t(c);
        % the other mesh's segment covering t* (or its periodic copy)
        seg_t = segments{other}.t;
        k = [];
        for shift = [0, period, -period]
            k = find(seg_t(1, :) <= t_star + shift & t_star + shift <= seg_t(2, :), 1);
            if ~isempty(k), t_star = t_star + shift; break; end
        end
        if isempty(k)
            error("get_partitions:collapsedUnpaired", ...
                "Mesh %d: face %d collapses to t = %g, which mesh %d does " + ...
                "not reach.", j, segments{j}.collapsed_face(c), t_star, other);
        end
        seg_s = segments{other}.s(:, k);
        s_star = seg_s(1) + (t_star - seg_t(1, k)) / (seg_t(2, k) - seg_t(1, k)) ...
            * (seg_s(2) - seg_s(1));
        induced{j}(:, end + 1) = [0; 1; segments{j}.collapsed_face(c)];
        induced{other}(:, end + 1) = [s_star; s_star; segments{other}.face(k)];
        t_lo(end + 1) = t_star; %#ok<AGROW>
        t_hi(end + 1) = t_star; %#ok<AGROW>
        side(end + 1) = j; %#ok<AGROW>
    end
end
collapsed.induced = induced;
collapsed.flag = [false(1, n_regular), true(1, numel(side))];
collapsed.side = side;
end

function N = collapsed_normals(interface, mesh_list, partitions, collapsed, rule)
% quadrature.normals with the collapsed patches' normals replaced by the
% collapsed_normal rule (see the description). The curve's orientation
% relative to each mesh is read off that mesh's regular patches: s(j) =
% +1 when the curve normal points from Sigma into mesh j.
N = partitions{1}.quadrature.normals;
regular = ~collapsed.flag;
orientation = zeros(1, 2);
for j = 1:2
    q = partitions{j}.quadrature;
    into_mesh = q.physical(:, :, regular) - q.interface(:, :, regular);
    orientation(j) = sign(sum(q.normals(:, :, regular) .* into_mesh, "all", "omitnan"));
end
patches = find(collapsed.flag);
period = diff(interface.t_range);
for c = 1:numel(patches)
    e = patches(c);
    j = collapsed.side(c);
    q = partitions{j}.quadrature;
    switch rule
        case "path"
            v = q.physical(:, :, e) - q.interface(:, :, e);
            v = v ./ vecnorm(v, 2, 1);
        case "face"
            f = partitions{j}.faces(1, 1, e);
            P = mesh_list{j}.coordinates(1:2, mesh_list{j}.faces(1:2, f));
            v = [P(2, 2) - P(2, 1); P(1, 1) - P(1, 2)];
            v = v / norm(v);
            % from Sigma into mesh j
            v = v * sign(sum(v .* (mean(q.physical(:, :, e), 2) - q.interface(:, 1, e))));
            v = repmat(v, 1, size(q.physical, 2));
        case "bisector"
            t_star = partitions{j}.vertices_interface(1, 1, e);
            n_pm = interface.unit_normal(reshape(t_star + [-1, 1] * 1e-6 * period, 1, []));
            v = sum(reshape(n_pm, 2, 2), 2);
            if norm(v) < 1e-3
                error("get_partitions:bisectorUndefined", ...
                    "The one-sided normals of Sigma at t = %g cancel (a cusp); " + ...
                    "use collapsed_normal = ""path"" or ""face"".", t_star);
            end
            N(:, :, e) = repmat(v / norm(v), 1, size(N, 2));
            continue
    end
    N(:, :, e) = orientation(j) * v;
end
end

function enriched_partition = enrich_partition(mesh, induced_partition)
n_patches = size(induced_partition, 2);

S_f = induced_partition(end, :);

[elements, side_numbers] = mesh.get_element_by_face(S_f);

S_K = reshape(elements, 1, 1, n_patches);
S_sn = reshape(side_numbers, 1, 1, n_patches);
S_t = zeros(1, 1, n_patches);

for e = 1:n_patches
    S_t(e) = mesh.tau(S_sn(e), S_K(e));
end

enriched_partition.vertices = reshape(induced_partition(1:end-1, :), 1, [], n_patches);
enriched_partition.faces = reshape(induced_partition(end, :), 1, [], n_patches);
enriched_partition.elements = S_K;
enriched_partition.side_number = S_sn;
enriched_partition.tau = S_t;
enriched_partition.area = reshape(abs(diff(enriched_partition.vertices)), 1, []) .* mesh.area(enriched_partition.faces);
end
