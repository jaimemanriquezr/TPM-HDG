%-------------------------------------------------------------------------
%Description: Red refinement: conformity, geometry, marker inheritance and
%             the parent maps.
%
%             The counts are the conformity certificate. A midpoint that
%             failed to be shared between two neighbours would show up as a
%             vertex and face surplus long before it showed up as a wrong
%             answer, and the per-parent area check is what catches a
%             mistaken child triple -- the global area sum does not, since
%             any partition of the parent preserves it.
%-------------------------------------------------------------------------
n_failed = 0;

%% A mesh with more than one tag on its boundary --------------------------
mesh = StructuredSimplicialGrid(2, 2);
mesh = apply_tag(mesh, "Neumann", @(X) isapprox(X(1, :), 1));
mesh = add_stabilization(mesh, 3.0);

n_el_0 = size(mesh.elements, 2);
n_fc_0 = size(mesh.faces, 2);
n_vt_0 = size(mesh.coordinates, 2);

[fine, info] = uniform_refine(mesh);

%% Counts: conformity ------------------------------------------------------
n_failed = n_failed + check("elements x4", ...
    size(fine.elements, 2) == 4 * n_el_0);
n_failed = n_failed + check("vertices + one per face", ...
    size(fine.coordinates, 2) == n_vt_0 + n_fc_0);
n_failed = n_failed + check("faces 2F + 3E", ...
    size(fine.faces, 2) == 2 * n_fc_0 + 3 * n_el_0);

%% Geometry ---------------------------------------------------------------
per_parent = accumarray(info.parent_element(:), fine.volume(:));
n_failed = n_failed + check("per-parent area", ...
    max(abs(per_parent(:) - mesh.volume(:))) < 1e-12);
n_failed = n_failed + check("orientation preserved", ...
    all(sign(fine.volume) == sign(mesh.volume(info.parent_element))));
n_failed = n_failed + check("h halves", ...
    abs(max(fine.area) - 0.5 * max(mesh.area)) < 1e-12);

%% Markers ----------------------------------------------------------------
% The trap: the constructor stamps Dirichlet on every boundary edge, so an
% un-inherited interface would silently become a Dirichlet boundary.
for tag = keys(mesh.boundary_tags).'
    coarse_count = numel(get_face_by_tag(mesh, tag));
    fine_count   = numel(get_face_by_tag(fine, tag));
    n_failed = n_failed + check(sprintf("tag %s doubles (%d -> %d)", ...
        tag, coarse_count, fine_count), fine_count == 2 * coarse_count);
end
n_failed = n_failed + check("faces without a parent are interior", ...
    all(fine.faces(3, info.parent_face == 0) == 0));
n_failed = n_failed + check("boundary_tags carried over", ...
    isequal(keys(fine.boundary_tags), keys(mesh.boundary_tags)));

%% The reference-to-reference maps ----------------------------------------
% The map is what the error routine relies on, so test it the way it is
% used: the same physical point reached through the child and through the
% parent must coincide.
rs = [0.2 0.6 0.1; 0.3 0.1 0.7];
worst = 0;
for j = 1:size(fine.elements, 2)
    x_child  = get_2d_physical(fine, rs, j);
    rs_parent = info.A(:, :, j) * rs + info.b(:, j);
    x_parent = get_2d_physical(mesh, rs_parent, info.parent_element(j));
    worst = max(worst, max(abs(x_child(:) - x_parent(:))));
end
n_failed = n_failed + check(sprintf("child->parent map (%.1e)", worst), ...
    worst < 1e-13);

%% The face parameter map -------------------------------------------------
t = [0.15 0.5 0.85];
worst_face = 0;
inherited = find(info.parent_face > 0);
for f = inherited
    x_child  = get_1d_physical(fine, t, f);
    t_parent = info.alpha(f) * t + info.beta(f);
    x_parent = get_1d_physical(mesh, t_parent, info.parent_face(f));
    worst_face = max(worst_face, max(abs(x_child(:) - x_parent(:))));
end
n_failed = n_failed + check(sprintf("child->parent face map (%.1e)", worst_face), ...
    worst_face < 1e-13);
n_failed = n_failed + check("every parent face has two children", ...
    numel(inherited) == 2 * n_fc_0);

%% Repeated refinement ----------------------------------------------------
[fine3, info3] = uniform_refine(mesh, 3);
n_failed = n_failed + check("3 levels: element count", ...
    size(fine3.elements, 2) == 4^3 * n_el_0);
per_parent3 = accumarray(info3.parent_element(:), fine3.volume(:));
n_failed = n_failed + check("3 levels: per-parent area", ...
    max(abs(per_parent3(:) - mesh.volume(:))) < 1e-12);

worst3 = 0;
for j = 1:37:size(fine3.elements, 2)
    x_child  = get_2d_physical(fine3, rs, j);
    x_parent = get_2d_physical(mesh, ...
        info3.A(:, :, j) * rs + info3.b(:, j), info3.parent_element(j));
    worst3 = max(worst3, max(abs(x_child(:) - x_parent(:))));
end
n_failed = n_failed + check(sprintf("3 levels: composed map (%.1e)", worst3), ...
    worst3 < 1e-13);
n_failed = n_failed + check("3 levels: distinct maps = 4^3", ...
    info3.n_types == 64);

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
