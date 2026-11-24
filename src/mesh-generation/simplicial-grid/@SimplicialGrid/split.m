function [mesh_1, mesh_2] = split(obj, interface, options)
%-------------------------------------------------------------------------
%Description: Carves this grid into two subdomains meeting at `interface`,
%             by selecting the elements on either side of it. The receiver
%             is the background grid and is not modified.
%
%             The two returned meshes are always ordered
%                 mesh_1  the outer / upper side (away from the enclosed
%                         region, or above the graph)
%                 mesh_2  the inner / lower side.
%
%             The elements Sigma crosses are discarded or handed to one
%             side, per `cut_elements`. The fitted constructions, where
%             each subdomain is meshed as its own polygon rather than
%             carved from a grid, live in fitted_split.
%
%Usage:       [outer, inner] = UnitSquareMesh(h).split(interface)
%             [outer, inner] = StructuredSimplicialGrid(1/h, 1/h, "crossed").split( ...
%                                  interface, cut_elements = "mesh_2")
%
%             Two background grids of different resolution: split each
%             grid and keep one half of each,
%                 [mesh_1, ~] = bg_1.split(interface, ...);
%                 [~, mesh_2] = bg_2.split(interface, ...);
%             Nothing checks that the two halves belong together, so keep
%             the two calls adjacent and identically optioned.
%
%See also: fitted_split, submesh, prune_infeasible_elements
%-------------------------------------------------------------------------
arguments
    obj SimplicialGrid
    interface AbstractCurve = FlatCurveHorizontal(0.5);

    % What happens to the elements Sigma crosses:
    %   "discard"  dropped, so both subdomains lie strictly inside their
    %              exact domains and the removed band is the gap;
    %   "mesh_1"   given to the outer/upper mesh;
    %   "mesh_2"   given to the inner/lower mesh.
    options.cut_elements (1,1) string ...
        {mustBeMember(options.cut_elements, ["discard", "mesh_1", "mesh_2"])} = "discard"

    % How those elements are found.
    %   "vertices"   an element is cut when its vertices straddle Sigma.
    %   "geometric"  additionally locates a dense sampling of Sigma in the
    %                background triangulation, which also catches a
    %                triangle crossed by a bulge with all vertices on one
    %                side.
    % Empty resolves to "geometric" for "discard" and "vertices" for the
    % assign families.
    options.cut_detection (1,1) string ...
        {mustBeMember(options.cut_detection, ["vertices", "geometric", ""])} = ""

    % Largest angle, in degrees, allowed between an interface face's
    % discrete normal and the exact normal at the point it maps to; a face
    % exceeding it makes its element infeasible and it is dropped. Empty
    % resolves to 45 for "discard" and 90 (off) for the assign families.
    % See prune_infeasible_elements.
    options.max_normal_deviation double = []
    % Same bound for mesh_1 (the "outside" side) only; empty (the default)
    % keeps max_normal_deviation on both sides.
    options.max_normal_deviation_mesh_1 double {mustBeScalarOrEmpty} = []

    % Number of samples of Sigma. Empty resolves from the grid's own face
    % size, see resolve_samples.
    options.n_interface_samples double = []
    % Drop elements with two or more carve-exposed faces. Empty resolves
    % to true for "discard" and false for the assign families, which keep
    % the cut elements. See prune_hanging_elements.
    options.prune_hanging logical {mustBeScalarOrEmpty} = []
    % The same switch for mesh_1 (the "outside" side) only; empty (the
    % default) keeps prune_hanging on both sides.
    options.prune_hanging_mesh_1 logical {mustBeScalarOrEmpty} = []
    % How close a face must lie to Sigma, in multiples of its own length,
    % to be judged by the feasibility filter. Forwarded to
    % prune_infeasible_elements; empty keeps its default (2).
    options.interface_band double {mustBeScalarOrEmpty} = []
    % Meshes on which every element touching Sigma through a diagonal face
    % (neither horizontal nor vertical) is dropped, after which hanging and
    % disconnected elements are re-pruned. Empty (the default) leaves the
    % carve untouched. A diagonal interface face can host an unstable mode
    % of the coupled operator.
    options.prune_diagonal_faces (1,:) string ...
        {mustBeMember(options.prune_diagonal_faces, ["mesh_1", "mesh_2"])} = string.empty
    % Repeat the diagonal prune (and its clean-up) up to this many times,
    % stopping when no diagonal carve face is left: the hanging/
    % disconnected clean-up of one pass can expose new diagonals.
    options.prune_diagonal_passes (1,1) double {mustBeInteger, mustBePositive} = 1
    % A face counts as diagonal when both of its components exceed this
    % fraction of its length (the default only spares exactly horizontal
    % and vertical faces; e.g. 0.1 also spares faces within ~6 degrees of
    % the axes, as on a slightly sheared grid)
    options.diagonal_tol (1,1) double {mustBeNonnegative} = 1e-8
end
detection = resolve_detection(options.cut_detection, options.cut_elements);
max_deviation = resolve_deviation(options.max_normal_deviation, options.cut_elements);
n_samples = resolve_samples(options.n_interface_samples, max(obj.area(:)), detection);
[curve, region] = enclosed_region(interface, n_samples, obj);

prune_hanging = options.prune_hanging;
if isempty(prune_hanging)
    prune_hanging = options.cut_elements == "discard";
end
extract_elements = @(side, dev, hang) get_active_elements(obj, interface, curve, region, ...
    detection, dev, side, options.cut_elements, hang, ...
    options.interface_band);

max_deviation_1 = options.max_normal_deviation_mesh_1;
if isempty(max_deviation_1), max_deviation_1 = max_deviation; end
prune_hanging_1 = options.prune_hanging_mesh_1;
if isempty(prune_hanging_1), prune_hanging_1 = prune_hanging; end
K_1 = extract_elements("outside", max_deviation_1, prune_hanging_1);
K_2 = extract_elements("inside", max_deviation, prune_hanging);
[K_1, K_2] = resolve_sealed_cavities(obj, K_1, K_2);
for pass = 1:options.prune_diagonal_passes
    n_before = [numel(K_1), numel(K_2)];
    if ismember("mesh_1", options.prune_diagonal_faces)
        K_1 = prune_diagonal_interface_elements(obj, K_1, prune_hanging_1, options.diagonal_tol);
    end
    if ismember("mesh_2", options.prune_diagonal_faces)
        K_2 = prune_diagonal_interface_elements(obj, K_2, prune_hanging, options.diagonal_tol);
    end
    if isequal(n_before, [numel(K_1), numel(K_2)]), break, end
end
% submesh retags every newly exposed face as "Interface", which is what
% makes the two halves usable without further tagging.
mesh_1 = submesh(obj, K_1(:).');
mesh_2 = submesh(obj, K_2(:).');
end

%=========================================================================%
function K_active = get_active_elements(bg_mesh, interface, curve, region, ...
    detection, max_deviation, side, cut_elements, prune_hanging, interface_band)
K_active = classify(bg_mesh, curve, region, detection, side, cut_elements);
% A discrete normal at 90 degrees to the exact one projects to zero
% interface measure, leaving trace DOFs unconstrained; such elements are
% dropped here, not repaired downstream.
if isempty(interface_band)
    K_active = prune_infeasible_elements(bg_mesh, K_active, interface, max_deviation);
else
    K_active = prune_infeasible_elements(bg_mesh, K_active, interface, ...
        max_deviation, interface_band = interface_band);
end
% Before the connectivity prune: shedding teeth can sever a neck and leave
% an island. After the feasibility filter, which passes teeth (their faces
% sit just inside its angle bound).
if prune_hanging
    K_active = prune_hanging_elements(bg_mesh, K_active);
end
% After the feasibility prune (which can itself sever a thin neck):
% islands share no faces with the kept component, so dropping them
% cannot re-expose infeasible faces and no re-iteration is needed.
K_active = prune_disconnected_elements(bg_mesh, K_active);
end

function K = prune_diagonal_interface_elements(bg_mesh, K, prune_hanging, diagonal_tol)
% Drops every element of K exposing a diagonal face to the carve: a face
% owned by exactly one element of K that is not on the background grid's
% own boundary (the faces submesh will tag "Interface"), and whose edge is
% neither horizontal nor vertical. One pass; the prunes that follow clean
% up any tooth or island it leaves.
K = K(:).';
B = bg_mesh.boundaries(:, K);                      % faces of the kept elements
counts_K = accumarray(B(:), 1, [size(bg_mesh.faces, 2), 1]);
counts_all = accumarray(bg_mesh.boundaries(:), 1, [size(bg_mesh.faces, 2), 1]);
exposed = (counts_K == 1) & (counts_all == 2);     % carve faces, not the domain boundary

V1 = bg_mesh.coordinates(:, bg_mesh.faces(1, :));
V2 = bg_mesh.coordinates(:, bg_mesh.faces(2, :));
d = abs(V2 - V1);
tol = diagonal_tol * vecnorm(d, 2, 1);
diagonal = (d(1, :) > tol & d(2, :) > tol).';

drop = any(exposed(B) & diagonal(B), 1);
if ~any(drop)
    return
end
K = K(~drop);
if prune_hanging
    K = prune_hanging_elements(bg_mesh, K);
end
K = prune_disconnected_elements(bg_mesh, K);
end

function [K_1, K_2] = resolve_sealed_cavities(bg_mesh, K_1, K_2)
% A connected component of the uncovered region (elements in neither
% mesh: the pruned band, dropped islands) whose element-neighbours all
% belong to one mesh is a sealed cavity inside it. Left alone, submesh
% would tag the cavity's ring of faces "Interface" and its projections
% would double-cover Sigma on top of that mesh's main interface loop
% (get_partitions:foldedProjection) -- so the cavity is handed to the
% enclosing mesh, where its exposed faces become interior faces. The
% gap band proper touches both meshes and is never touched.
uncovered = setdiff(1:size(bg_mesh.elements, 2), union(K_1, K_2));
if isempty(uncovered)
    return
end
B = bg_mesh.boundaries(:, uncovered);
[sorted_faces, order] = sort(B(:));
shared = find(diff(sorted_faces) == 0);
[~, e_first] = ind2sub(size(B), order(shared));
[~, e_second] = ind2sub(size(B), order(shared + 1));
components = conncomp(graph(e_first, e_second, [], numel(uncovered)));

for c = 1:max(components)
    cavity = uncovered(components == c);
    % elements outside the component sharing a face with it
    touching = bg_mesh.incidence(:, 1, bg_mesh.boundaries(:, cavity));
    neighbors = setdiff(touching(:).', cavity);
    if isempty(neighbors)
        continue
    end
    sides = [all(ismember(neighbors, K_1)), all(ismember(neighbors, K_2))];
    if ~any(sides)
        continue   % touches both meshes (or the domain boundary): a gap
    end
    warning("SimplicialGrid:split:sealedCavity", ...
        "Sealed cavity of %d uncovered element(s) handed to mesh %d, " + ...
        "where its faces become interior.", numel(cavity), find(sides, 1));
    if sides(1)
        K_1 = union(K_1, cavity);
    else
        K_2 = union(K_2, cavity);
    end
end
end

function elements = classify(bg_grid, curve, region, detection, side, cut_elements)
in_region = @(xv, yv) inpolygon(xv, yv, region(1, :), region(2, :));
vertices = {bg_grid.vertices(1, :, :), bg_grid.vertices(2, :, :)};
all_in = all(in_region(vertices{:}));
any_in = any(in_region(vertices{:}));

cut_T = find(any_in & ~all_in).';
if detection == "geometric"
    % Only the interface samples are located, not the segment closing an
    % open region, so elements along the far side of the grid are not
    % mistaken for cut ones.
    crossed = pointLocation(triangulation(bg_grid), curve.');
    cut_T = union(unique(crossed(~isnan(crossed))), cut_T);
end

if side == "inside"
    elements = setdiff(find(all_in).', cut_T);
    takes_cut = cut_elements == "mesh_2";
else
    elements = setdiff(find(~any_in).', cut_T);
    takes_cut = cut_elements == "mesh_1";
end

% union sorts, so the receiving mesh numbers its elements in
% background-grid order.
if takes_cut
    elements = union(elements, cut_T);
end
end

function [curve, region] = enclosed_region(interface, n_samples, bg_mesh)
% The polygon whose interior is the inner/lower subdomain, and the sampling
% of Sigma alone. A closed interface already encloses its region; an open
% one is closed off below the grid, using the grid's own bounding box
% rather than assuming the unit square, so this works for any background.
tt = linspace(interface.t_range(1), interface.t_range(2), n_samples);
curve = interface.get_2d_points(tt);
if interface.is_periodic()
    region = curve;
else
    box_min = min(bg_mesh.coordinates, [], 2);
    box_max = max(bg_mesh.coordinates, [], 2);
    % Drop below the grid so the closing segment cannot clip any element,
    % and span it in x so the region is closed whatever the curve's ends.
    floor_y = box_min(2) - (box_max(2) - box_min(2));
    region = [curve, [box_max(1); floor_y], [box_min(1); floor_y]];
end
end

function detection = resolve_detection(detection, cut_elements)
if detection == ""
    detection = "geometric";
    if cut_elements ~= "discard"
        detection = "vertices";
    end
end
end

function n_samples = resolve_samples(n_samples, h, detection)
if ~isempty(n_samples)
    return;
end
if detection == "geometric"
    % Spacing along Sigma well below h, so no crossed triangle is skipped.
    % h here is the largest face length of the background grid.
    n_samples = max(2000, ceil(200/h));
else
    n_samples = 250;
end
end

function max_deviation = resolve_deviation(max_deviation, cut_elements)
if isempty(max_deviation)
    max_deviation = 45 * (cut_elements == "discard") + 90 * (cut_elements ~= "discard");
end
if ~isscalar(max_deviation) || max_deviation < 0 || max_deviation > 90
    error("SimplicialGrid:split:badDeviation", ...
        "max_normal_deviation must be a scalar in [0, 90] degrees.");
end
end
