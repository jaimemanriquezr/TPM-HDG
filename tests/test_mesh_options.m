% test_mesh_options.m
% The one-sided options of SimplicialGrid.split and interface_meshes act
% on their own side only, and size_function grades the background:
%   - split max_normal_deviation_mesh_1 / prune_hanging_mesh_1 change
%     mesh_1 and leave mesh_2 element-for-element unchanged;
%   - interface_meshes max_normal_deviation_fluid / prune_hanging_fluid
%     change the fluid mesh and leave the solid mesh unchanged;
%   - interface_meshes size_function: a constant h/2 gives about four
%     times the elements of h, and a graded one gives smaller elements
%     where it is small.
% Geometry: the sine graph 0.5 + 0.1 sin(2 pi t) on the unit square,
% fluid above (split's mesh_1), h = 0.1, where each option changes the
% carve.
clc; clear;

n_fail = 0;
n_el = @(m) size(m.elements, 2);
curve = SmoothCurveHorizontal("0.5 + .1*sin(2*pi*t)", [0 1]);
h = 0.1;

%% split: the mesh_1 options
bg = StructuredSimplicialGrid(10, 10, 0, 1, 0, 1, style = "crossed");
[m1, m2] = bg.split(curve, "cut_elements", "discard");
[d1, d2] = bg.split(curve, "cut_elements", "discard", max_normal_deviation_mesh_1 = 20);
n_fail = n_fail + one_sided("split max_normal_deviation_mesh_1", m1, m2, d1, d2);
[p1, p2] = bg.split(curve, "cut_elements", "discard", prune_hanging_mesh_1 = false);
n_fail = n_fail + one_sided("split prune_hanging_mesh_1 = false", m1, m2, p1, p2);
[q1, q2] = bg.split(curve, "cut_elements", "mesh_2");
[r1, r2] = bg.split(curve, "cut_elements", "mesh_2", prune_hanging_mesh_1 = true);
n_fail = n_fail + one_sided("split prune_hanging_mesh_1 = true", q1, q2, r1, r2);

%% interface_meshes: the fluid options
[fluid, solid] = interface_meshes(curve, "interior", h, fluid_side = "upper");
[fluid_d, solid_d] = interface_meshes(curve, "interior", h, fluid_side = "upper", ...
    max_normal_deviation_fluid = 20);
n_fail = n_fail + one_sided("interface_meshes max_normal_deviation_fluid", fluid, solid, fluid_d, solid_d);
[fluid_p, solid_p] = interface_meshes(curve, "interior", h, fluid_side = "upper", ...
    prune_hanging_fluid = false);
n_fail = n_fail + one_sided("interface_meshes prune_hanging_fluid", fluid, solid, fluid_p, solid_p);

%% interface_meshes: size_function
[fluid_h, solid_h] = interface_meshes(curve, "interior", h, fluid_side = "upper", ...
    size_function = @(P) h + 0*P(:, 1));
[fluid_h2, solid_h2] = interface_meshes(curve, "interior", h, fluid_side = "upper", ...
    size_function = @(P) h/2 + 0*P(:, 1));
ratio = (n_el(fluid_h2) + n_el(solid_h2)) / (n_el(fluid_h) + n_el(solid_h));
fprintf("size_function h/2 vs h: element ratio %.2f\n", ratio);
n_fail = n_fail + check(ratio > 3 && ratio < 6, sprintf("size_function h/2 vs h: ratio %.2f", ratio));
graded = interface_meshes(curve, "interior", h, fluid_side = "upper", ...
    size_function = @(P) 0.03 + 0.1*P(:, 1));
x_c = mean(reshape(graded.coordinates(1, graded.elements), 3, []));
area_ratio = mean(graded.volume(x_c > 0.7)) / mean(graded.volume(x_c < 0.3));
fprintf("graded size_function: mean area x > 0.7 / x < 0.3 = %.2f\n", area_ratio);
n_fail = n_fail + check(area_ratio > 3, sprintf("graded size_function: area ratio %.2f", area_ratio));

fprintf("%d check(s) failed.\n", n_fail);
if n_fail > 0
    error("test_mesh_options:failed", "%d checks failed", n_fail);
end

function n = one_sided(what, a1, a2, b1, b2)
% the option changed side 1 and left side 2 untouched
changed = ~isequal(a1.elements, b1.elements);
kept = isequal(a2.elements, b2.elements) && isequal(a2.coordinates, b2.coordinates);
fprintf("%-46s side 1: %d -> %d elements, side 2 unchanged: %d\n", what, ...
    size(a1.elements, 2), size(b1.elements, 2), kept);
n = check(changed && kept, what);
end
