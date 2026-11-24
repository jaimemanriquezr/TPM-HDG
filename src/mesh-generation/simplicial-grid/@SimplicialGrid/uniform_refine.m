function [fine, info] = uniform_refine(obj, n_times, options)
%-------------------------------------------------------------------------
%Description: Uniform red refinement: every element is split into four by
%             joining its edge midpoints. Returns the refined mesh together
%             with the bookkeeping that makes a coarse field evaluable on
%             it without any point location.
%
%             Each child lies inside exactly one parent, so a coarse field
%             restricted to a child is the same polynomial and coarse/fine
%             comparisons need no interpolation. `info` carries the maps:
%
%               parent_element  which coarse element each child came from,
%               child_type      which of the four slots it occupies,
%               A, b            r_coarse = A*r_fine + b, the reference-to-
%                               reference map, identical for every element
%                               of a given type, so a coarse field is
%                               evaluated for all of them in one broadcast
%                               call.
%
%             Without `interface`, the geometry does not move: Omega_h, its
%             boundary and (on an unfitted mesh) the staircase and the gap
%             are unchanged, as an unfitted refinement ladder requires. Do
%             not re-run `split`, `immerse_mesh` or the prune passes on the
%             result: they select elements of a background grid and are
%             calibrated in units of the face length.
%
%             With `interface`, new interface vertices are projected onto
%             Sigma, so a fitted ladder's interpolation refines with the
%             mesh. This is required there: otherwise the TPM gap stays at
%             the base mesh's sagitta and the convergence rate collapses
%             once h approaches it.
%
%Usage:       [fine, info] = uniform_refine(mesh)
%             [fine, info] = uniform_refine(mesh, 3)     % 4^3 children each
%
%See also: refinement_l2_error, SimplicialGrid.split
%-------------------------------------------------------------------------
arguments
    obj     SimplicialGrid
    n_times (1,1) double {mustBeInteger, mustBeNonnegative} = 1

    % Curve to project new interface vertices onto. Without it a new
    % vertex sits at the chord midpoint, the discrete interface never
    % approaches Sigma and the TPM gap stays O(h0^2); with it the gap
    % shrinks like h^2. A projected child is no longer exactly its
    % parent's quarter (area changes of O(h^3)), so the maps in `info` are
    % built from the actual child and parent geometry; a fine quadrature
    % point just outside its parent evaluates the coarse polynomial by
    % continuation.
    options.interface AbstractCurve {mustBeScalarOrEmpty} = AbstractCurve.empty
end

fine = obj;
info = identity_info(obj);
for level = 1:n_times
    [fine, step] = refine_once(fine, options.interface);
    info = compose(info, step);
end
end

%=========================================================================%
function [fine, info] = refine_once(obj, interface)
n_vertices = size(obj.coordinates, 2);
n_faces    = size(obj.faces, 2);
n_elements = size(obj.elements, 2);

% Midpoints are numbered by face, so neighbours share them by construction
% and the result is conforming without any merge step. get_boundaries
% asserts that every element edge appears exactly once in `faces`, which is
% what makes this enumeration global and unique.
midpoints = 0.5 * (obj.coordinates(:, obj.faces(1, :)) ...
                 + obj.coordinates(:, obj.faces(2, :)));

% Pull the new interface vertices onto Sigma, so the linear interpolation
% of the interface improves with the mesh instead of staying frozen at the
% base polyline. Only interface faces move: a box wall is straight and its
% midpoint is already exact.
if ~isempty(interface)
    on_interface = obj.faces(3, :) == obj.boundary_tags("Interface");
    if any(on_interface)
        midpoints(:, on_interface) = ...
            project_to_curve(interface, midpoints(:, on_interface));
    end
end
X = [obj.coordinates, midpoints];

V = obj.elements;                       % V(i, K) = vertex i of element K
M = n_vertices + obj.boundaries;        % M(i, K) = midpoint of local side i

% Local sides are e1 = x2->x3, e2 = x1->x3, e3 = x1->x2 (vertex_param in
% get_boundaries), so m_i is opposite vertex i: the two midpoints adjacent
% to a corner v_i are the ones not indexed by i. All four triples are
% listed in the parent's rotational sense, so orientation is preserved.
corner_1 = [V(1, :); M(3, :); M(2, :)];
corner_2 = [M(3, :); V(2, :); M(1, :)];
corner_3 = [M(2, :); M(1, :); V(3, :)];
middle   = [M(1, :); M(2, :); M(3, :)];

% Grouped by parent, so the child index is pure arithmetic: 4*(K-1) + c.
E = reshape(permute(cat(3, corner_1, corner_2, corner_3, middle), ...
    [1 3 2]), 3, []);

% "array" keeps coordinates and elements verbatim, which keeps the parent
% map above valid, but rebuilds the face table from scratch.
fine = SimplicialGrid(X, E, mesh_generator = "array");
fine.boundary_tags = obj.boundary_tags;   % apply_tag extends it per mesh

[fine, parent_face, alpha, beta] = inherit_face_data(obj, fine, ...
    n_vertices, n_faces, n_elements);

% A and b are stored per entity, never per type: after the first level a
% child's composed map depends on its whole ancestry. `child_type` groups
% the identical ones and is what the broadcast evaluation loops over.
child_type = repmat(1:4, 1, n_elements);
parent_element = repelem(1:n_elements, 4);
if isempty(interface)
    % Pure red refinement: every child is its parent's exact quarter, so
    % the four universal maps apply and no geometry needs consulting.
    A_types = child_maps_A();
    b_types = child_maps_b();
    A = A_types(:, :, child_type);
    b = b_types(:, child_type);
else
    % A projected child is not similar to its parent, so the maps are read
    % off the actual triangles. Exact for any pair.
    [A, b] = geometric_maps(fine, obj, parent_element);
    [~, ~, child_type] = unique([reshape(A, 4, []); b].', 'rows');
    child_type = child_type(:).';
end
info = struct( ...
    'parent_element', parent_element, ...
    'child_type',     child_type, ...
    'A', A, ...
    'b', b, ...
    'parent_face', parent_face, 'alpha', alpha, 'beta', beta, ...
    'n_types', max(child_type));
end

%=========================================================================%
function P = project_to_curve(interface, X)
% Closest point of Sigma, through whichever query the curve lineage
% provides. Both are global rather than local searches, which matters on a
% concave lobe where a local search takes the wrong branch.
if isa(interface, "LevelSetCurve")
    P = interface.project(X);
else
    P = zeros(size(X));
    for j = 1:size(X, 2)
        t = interface.find_closest_point(X(1, j), X(2, j));
        P(:, j) = interface.get_2d_points(t);
    end
end
end

%=========================================================================%
function [A, b] = geometric_maps(fine, coarse, parent_element)
% r_parent = A r_child + b, read off the actual triangles: with
% B = [v1 - v3, v2 - v3] the map is A = B_parent \ B_child and
% b = B_parent \ (x3_child - x3_parent).
Bc = element_basis(coarse);
Bf = element_basis(fine);
x3c = coarse.coordinates(:, coarse.elements(3, :));
x3f = fine.coordinates(:, fine.elements(3, :));
Bp = Bc(:, :, parent_element);
A = pagemldivide(Bp, Bf);
b = reshape(pagemldivide(Bp, reshape(x3f - x3c(:, parent_element), 2, 1, [])), 2, []);
end

function B = element_basis(mesh)
V = reshape(mesh.coordinates(:, mesh.elements), 2, 3, []);
B = cat(2, V(:, 1, :) - V(:, 3, :), V(:, 2, :) - V(:, 3, :));
end

%=========================================================================%
function [fine, parent_face, alpha, beta] = inherit_face_data( ...
    obj, fine, n_vertices, n_faces, n_elements)
% Restore the boundary markers and record the face parentage.
%
% SimplicialGrid.get_faces marks every boundary edge Dirichlet, so without
% this the interface faces (`intfaces`, used by the TPM transfer) would be
% relabelled. It also stores each
% edge sorted as (min, max), so neither the order nor the direction of the
% child faces is ours: they are matched on an orientation-free key, and the
% direction is recovered separately below.
midpoint_vertex = n_vertices + (1:n_faces);
sub_edges = [obj.faces(1, :), midpoint_vertex; ...
             midpoint_vertex,  obj.faces(2, :)];      % both halves, in
                                                      % the parent direction
[found, loc] = ismember(edge_key(sub_edges, size(fine.coordinates, 2)), ...
                        edge_key(fine.faces(1:2, :), size(fine.coordinates, 2)));
assert(all(found), "SimplicialGrid:uniform_refine:missingSubEdge", ...
    "A parent sub-edge is absent from the child face table.");

% Interior parents carry marker 0 and their children are interior anyway,
% so writing every parent's marker is a no-op there and overwrites the
% spurious Dirichlet on the boundary children.
fine.faces(3, loc) = repmat(obj.faces(3, :), 1, 2);

parent_face = zeros(1, size(fine.faces, 2));
parent_face(loc) = [1:n_faces, 1:n_faces];

% Face parameter map t_parent = alpha*t_child + beta. `half` selects which
% half of the parent the child covers; `flipped` records that get_faces
% stored it against the parent's direction (about half of them are).
flipped = fine.faces(1, loc) ~= sub_edges(1, :);
half    = [zeros(1, n_faces), ones(1, n_faces)];
alpha = zeros(1, size(fine.faces, 2));
beta  = zeros(1, size(fine.faces, 2));
alpha(loc) = 0.5 * (1 - 2 * flipped);
beta(loc)  = 0.5 * half + 0.5 * flipped;

assert(all(fine.faces(3, parent_face == 0) == 0), ...
    "SimplicialGrid:uniform_refine:spuriousMarker", ...
    "A face interior to a parent element carries a boundary marker.");
assert(size(fine.faces, 2) == 2*n_faces + 3*n_elements, ...
    "SimplicialGrid:uniform_refine:faceCount", ...
    "Face count %d does not match the expected %d.", ...
    size(fine.faces, 2), 2*n_faces + 3*n_elements);
end

%=========================================================================%
function key = edge_key(edges, n_vertices)
% Orientation-independent identifier, the same one get_faces sorts on.
key = min(edges, [], 1) + n_vertices * max(edges, [], 1);
end

%=========================================================================%
function A = child_maps_A()
% r_parent = A(:,:,c) * r_child + b(:,c), in the reference coordinates of
% reference_element(2): v1 = (1,0), v2 = (0,1), v3 = the origin.
%
% Every one is +-I/2, a homothety: the four maps are element-independent,
% and the enriched stress space stays nested under this refinement.
A = repmat(0.5 * eye(2), 1, 1, 4);
A(:, :, 4) = -0.5 * eye(2);
end

function b = child_maps_b()
b = [0.5 0 0 0.5; 0 0.5 0 0.5];
end

%=========================================================================%
function info = identity_info(obj)
% n_times = 0: every element is its own parent, the maps are the identity.
n_elements = size(obj.elements, 2);
n_faces    = size(obj.faces, 2);
info = struct( ...
    'parent_element', 1:n_elements, ...
    'child_type',     ones(1, n_elements), ...
    'A', repmat(eye(2), 1, 1, n_elements), ...
    'b', zeros(2, n_elements), ...
    'parent_face', 1:n_faces, ...
    'alpha', ones(1, n_faces), 'beta', zeros(1, n_faces), ...
    'n_types', 1);
end

%=========================================================================%
function info = compose(prev, step)
% Compose one refinement onto the accumulated maps. `prev` is indexed by
% the previous level's entities, `step` by the new level's.
%
% After the first level A and b are stored per entity rather than per type,
% because a child's composed map depends on the whole chain it came from.
% The distinct maps still number 4^n_times, and `child_type` keeps them
% addressable so the broadcast evaluation stays one call per distinct map.
up = step.parent_element;                        % new child -> prev entity
A_up = prev.A(:, :, up);
A = pagemtimes(A_up, step.A);
b = reshape(pagemtimes(A_up, reshape(step.b, 2, 1, [])), 2, []) ...
    + prev.b(:, up);

% A face interior to a previous-level element has no ancestor in the root
% skeleton, and that absence propagates.
parent_face = zeros(1, numel(step.parent_face));
alpha = zeros(1, numel(step.parent_face));
beta  = zeros(1, numel(step.parent_face));
inherited = step.parent_face > 0;
grandparent = prev.parent_face(step.parent_face(inherited));
keep = grandparent > 0;
idx = find(inherited);  idx = idx(keep);
parent_face(idx) = grandparent(keep);
alpha(idx) = prev.alpha(step.parent_face(idx)) .* step.alpha(idx);
beta(idx)  = prev.alpha(step.parent_face(idx)) .* step.beta(idx) ...
           + prev.beta(step.parent_face(idx));

% Distinct composed maps, so the caller can group by them.
[~, ~, child_type] = unique([reshape(A, 4, []); b].', 'rows');

info = struct( ...
    'parent_element', prev.parent_element(up), ...
    'child_type',     child_type(:).', ...
    'A', A, 'b', b, ...
    'parent_face', parent_face, 'alpha', alpha, 'beta', beta, ...
    'n_types', max(child_type));
end
