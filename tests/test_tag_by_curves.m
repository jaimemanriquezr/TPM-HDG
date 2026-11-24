%-------------------------------------------------------------------------
%Description: tag_by_curves: boundary faces assigned to curves by
%             closest-point projection, unresolved faces marked NaN, and
%             the remedies (bisect_boundary_face, force).
%
%             A fitted mesh must reproduce its own tags exactly. On a
%             carved flag the only conflicts are at the junctions of the
%             free sides and the clamped root, and after the remedies the
%             fluid and solid cover the same parameter range, so
%             get_partitions accepts them.
%-------------------------------------------------------------------------
n_failed = 0;
walls = {FlatCurveHorizontal(0, [0 1]), FlatCurveVertical(1, [0 1]), ...
         FlatCurveHorizontal(1, [0 1]), FlatCurveVertical(0, [0 1])};
wall_tags = repmat("Dirichlet", 1, 4);

%% Fitted circle: the induced tags are the existing ones ------------------
circle = CircularCurve([0.5; 0.5], 0.25);
[fluid, solid] = interface_meshes(circle, "fitted", 0.1);
[~, a] = tag_by_curves(fluid, [{circle}, walls], ["Interface", wall_tags]);
n_failed = n_failed + check("fitted fluid: every boundary face resolved", ...
    ~any(isnan(a.curve)) && all(a.status(a.curve > 0) == "ok"));
n_failed = n_failed + check("fitted fluid: no tag changed", isempty(a.changed));
n_failed = n_failed + check("fitted fluid: interior faces not judged", ...
    all(a.curve(fluid.faces(3, :) == 0) == 0));
[~, a] = tag_by_curves(solid, {circle}, "Interface");
n_failed = n_failed + check("fitted solid: resolved, no tag changed", ...
    ~any(isnan(a.curve)) && isempty(a.changed));

%% Carved flag: conflicts only at the junctions ---------------------------
x0 = 0.31;  x1 = 0.79;  y0 = 0.41;  y1 = 0.59;
corners = [x0 x1 x1 x0; y0 y0 y1 y1];
body = PolygonalCurve(corners);                       % closed, for carving
sides = PolygonalCurve(corners, closed = false);      % the free sides, t in [0, 3]
root = FlatCurveVertical(x0, [y0 y1]);                % the clamped side
% two non-nested grids, so the solid's range ends fall inside fluid faces
grid = StructuredSimplicialGrid(20, 20, 0, 1, 0, 1, style = "crossed");
[fluid, ~] = grid.split(body, cut_elements = "discard");
[~, solid] = StructuredSimplicialGrid(27, 27, 0, 1, 0, 1, style = "crossed").split(body, cut_elements = "discard");

[~, a] = tag_by_curves(solid, {sides, root}, ["Interface", "Clamp"]);
bad = find(isnan(a.curve));
n_failed = n_failed + check(sprintf("carved solid: %d NaN faces, all ambiguous", numel(bad)), ...
    ~isempty(bad) && all(a.status(bad) == "ambiguous"));
% remedy: split at the equidistant line of the two curves, then re-tag
g = @(X) sides.distance(X) - root.distance(X);
solid = bisect_boundary_face(solid, bad, repmat({g}, size(bad)));
[solid, a] = tag_by_curves(solid, {sides, root}, ["Interface", "Clamp"]);
n_failed = n_failed + check("carved solid: resolved after the equidistant split", ~any(isnan(a.curve)));
n_failed = n_failed + check("carved solid: clamp and interface both present", ...
    ~isempty(get_face_by_tag(solid, "Clamp")) && ~isempty(solid.intfaces));

curves = [{sides, root}, walls];
tags = ["Interface", "Clamp", wall_tags];
[~, a] = tag_by_curves(fluid, curves, tags);
bad = find(isnan(a.curve));
mid = 0.5*(fluid.coordinates(:, fluid.faces(1, bad)) + fluid.coordinates(:, fluid.faces(2, bad)));
near_junction = min(vecnorm(permute(mid, [1 2 3]) - permute(corners(:, [1 4]), [1 3 2])), [], 3) < 2*0.05;
n_failed = n_failed + check(sprintf("carved fluid: %d NaN faces, all at the root junctions", numel(bad)), ...
    ~isempty(bad) && all(near_junction) && all(a.status(bad) == "ambiguous"));
n_failed = n_failed + check("carved fluid: the rest resolved", ...
    all(a.curve(setdiff(find(a.curve ~= 0), bad)) > 0));

%% Forcing alone clears every NaN -----------------------------------------
[forced, b] = tag_by_curves(fluid, curves, tags, force = [bad; a.candidates(1, bad)]);
n_failed = n_failed + check("force: no NaN left, status forced", ...
    ~any(isnan(b.curve)) && all(b.status(bad) == "forced"));
n_failed = n_failed + check("force: the forced faces carry the candidate's tag", ...
    all(forced.faces(3, bad) == forced.boundary_tags(tags(a.candidates(1, bad)))));

%% Remedy: split at the solid's range ends, then the parameter rule -------
seg = map_mesh_to_interface(sides, solid);
t_solid = [min(seg.t(1, :)), max(seg.t(2, :))];
n_split = 0;
for t_end = t_solid
    g = @(X) sides.map_points_to_interface(X) - t_end;
    [~, a] = tag_by_curves(fluid, curves, tags);
    judged = find(a.curve == 1 | isnan(a.curve));     % faces facing the sides
    ga = g(fluid.coordinates(:, fluid.faces(1, judged)));
    gb = g(fluid.coordinates(:, fluid.faces(2, judged)));
    straddle = judged(ga.*gb < 0 & min(abs(ga), abs(gb)) > 1e-9);
    fluid = bisect_boundary_face(fluid, straddle, repmat({g}, size(straddle)));
    n_split = n_split + numel(straddle);
end
n_failed = n_failed + check(sprintf("remedy: %d fluid faces split at the solid's range ends", n_split), n_split > 0);
% faces whose image lies outside the solid's range, and the junction ties,
% belong to the root
[~, a] = tag_by_curves(fluid, curves, tags);
on_sides = find(a.curve == 1);
t_mid = sides.map_points_to_interface(0.5*(fluid.coordinates(:, fluid.faces(1, on_sides)) ...
    + fluid.coordinates(:, fluid.faces(2, on_sides))));
outside = on_sides(t_mid < t_solid(1) | t_mid > t_solid(2));
to_root = [find(isnan(a.curve)), outside];
[fluid, a] = tag_by_curves(fluid, curves, tags, force = [to_root; 2*ones(size(to_root))]);
n_failed = n_failed + check("remedy: no NaN left", ~any(isnan(a.curve)));
try
    get_partitions(sides, fluid, solid);
    ok = true;
catch err
    fprintf(2, "  get_partitions: %s\n", err.message);
    ok = false;
end
n_failed = n_failed + check("remedy: get_partitions accepts the pair", ok);

%% bisect_boundary_face: one face, one new element, tags kept -------------
mesh = StructuredSimplicialGrid(4, 4);
mesh = apply_tag(mesh, "Neumann", @(X) abs(X(1, :) - 1) < 1e-12);
f = find(mesh.faces(3, :) == mesh.boundary_tags("Neumann"), 1);
fine = bisect_boundary_face(mesh, f, {@(X) X(2, :) - mean(mesh.coordinates(2, mesh.faces(1:2, f))) - 0.01});
n_failed = n_failed + check("bisect: one more element, same area", ...
    size(fine.elements, 2) == size(mesh.elements, 2) + 1 && abs(sum(fine.volume) - sum(mesh.volume)) < 1e-14);
n_failed = n_failed + check("bisect: orientation kept", all(fine.volume > 0) == all(mesh.volume > 0));
n_failed = n_failed + check("bisect: Neumann faces + 1, Dirichlet unchanged", ...
    numel(get_face_by_tag(fine, "Neumann")) == numel(get_face_by_tag(mesh, "Neumann")) + 1 ...
    && numel(get_face_by_tag(fine, "Dirichlet")) == numel(get_face_by_tag(mesh, "Dirichlet")));
try
    bisect_boundary_face(mesh, f, {@(X) X(1, :) + 1});
    raised = false;
catch err
    raised = err.identifier == "SimplicialGrid:bisect_boundary_face:noCrossing";
end
n_failed = n_failed + check("bisect: no sign change is an error", raised);

%% Overlapping curves warn --------------------------------------------------
% a stretch of the bottom wall given as a curve of its own
lastwarn("");
tag_by_curves(grid, [walls, {FlatCurveHorizontal(0, [0.2 0.6])}], [wall_tags, "Neumann"]);
[~, id] = lastwarn;
n_failed = n_failed + check("overlap: warning raised", id == "SimplicialGrid:tag_by_curves:overlappingCurves");

%% Side check: a curve behind the face is not assigned ---------------------
mesh = StructuredSimplicialGrid(10, 10, style = "crossed");
inner_line = FlatCurveHorizontal(0.9, [0 1]);
[~, a] = tag_by_curves(mesh, [walls([1 2 4]), {inner_line}], ["Dirichlet", "Dirichlet", "Dirichlet", "Neumann"]);
% the two corner faces of the top are nearest to the side walls, legitimately
X = reshape(mesh.coordinates(:, mesh.faces(1:2, :)), 2, 2, []);
top = find(a.curve ~= 0 & reshape(all(X(2, :, :) == 1 & X(1, :, :) > 0.05 & X(1, :, :) < 0.95, 2), 1, []));
n_failed = n_failed + check(sprintf("side: the %d inner top faces are NaN, status side", numel(top)), ...
    numel(top) == 8 && all(isnan(a.curve(top))) && all(a.status(top) == "side"));

fprintf("\n%d check(s) failed.\n", n_failed);

%=========================================================================%
function failed = check(name, condition)
failed = ~condition;
if failed
    fprintf(2, "  FAIL  %s\n", name);
else
    fprintf("  ok    %s\n", name);
end
end
