% TPM geometry on a circular interface carved from one background grid:
% first colours the elements outside, inside and cut by the circle, then
% plots the two submeshes (cut elements to the inner side) with their
% induced partitions, transfer paths and quadrature points.
clc; clear; close all;
unit_square = [[0;0], [1;0], [1;1], [0;1]];

h = .25/4;
bg_grid = SimplicialGrid(unit_square, element_diameter = h);

xy0 = [.5;.65];
r0 = .2;
I = CircularCurve(xy0, r0);
tt = linspace(0, 1);
phi0 = I.get_2d_points(tt);

in_region = @(xv, yv) inpolygon(xv, yv, phi0(1,:), phi0(2,:));

vertices = {bg_grid.vertices(1, :, :), bg_grid.vertices(2, :, :)};
int_T = find(all(in_region(vertices{:})));
cut_T = setdiff(find(any(in_region(vertices{:}))), int_T);
out_T = setdiff((1:size(bg_grid.elements, 2)).', [int_T; cut_T]);

ax = axes();
plot(bg_grid, ax);
plot(bg_grid, ax, ElementList=out_T, FaceColor='r');
plot(bg_grid, ax, ElementList=int_T, FaceColor='b');
plot(bg_grid, ax, ElementList=cut_T, FaceColor='y');
plot(ax, phi0(1, :), phi0(2, :), 'k', 'LineWidth', 3.0);

figure()
ax = axes();
ax.NextPlot = 'add';
axis(ax, "equal")
outer_mesh = submesh(bg_grid, out_T);
inner_mesh = submesh(bg_grid, [cut_T; int_T]);
plot(outer_mesh, ax, 'FaceAlpha', 0.5, 'FaceColor', 'r');
plot(inner_mesh, ax, 'FaceAlpha', 0.5, 'FaceColor', 'b');   

plot(I, ax, 'LineWidth', 4, 'Color', 'k');
    
[outer_partition, inner_partition] = get_partitions(I, outer_mesh, inner_mesh);
plot_vertices(ax, outer_mesh, outer_partition)
plot_paths(ax, outer_mesh, outer_partition)
plot_patch(ax, I, outer_mesh, outer_partition)
plot_gauss_points(ax, outer_mesh, outer_partition)
plot_gauss_points_paths(ax, outer_mesh, outer_partition)

plot_vertices(ax, inner_mesh, inner_partition)
plot_paths(ax, inner_mesh, inner_partition)
plot_patch(ax, I, inner_mesh, inner_partition)
plot_gauss_points(ax, inner_mesh, inner_partition)
plot_gauss_points_paths(ax, inner_mesh, inner_partition)

xlim([0 1])
ylim([0 1])


function plot_vertices(ax, mesh, partition)
    v_x = mesh.get_1d_physical(partition.vertices, partition.faces);
    for e = 1:size(v_x, 3)
        plot(ax, v_x(1, :, e), v_x(2, :, e), 'ko', 'LineWidth', 1.2);
    end
end

function plot_patch(ax, interface, mesh, partition)
    v_x = mesh.get_1d_physical(partition.vertices, partition.faces);
    v_I = interface.get_2d_points(interface.map_points_to_interface(v_x));
    for e = 1:size(v_x, 3)
        for j = 1:size(v_x, 2);
            tp = [v_x(:, j, e), v_I(:, j, e)];
            plot(ax, tp(1, :), tp(2, :), 'y-', 'LineWidth', 1.2);
        end
    end
end

function plot_paths(ax, mesh, partition)
    v_1 = partition.quadrature.physical;
    v_2 = partition.quadrature.interface;
    for e = 1:size(v_1, 3)
        for j = 1:size(v_1, 2);
            tp = [v_1(:, j, e), v_2(:, j, e)];
            plot(ax, tp(1, 2), tp(2, 2), 'ko', 'MarkerSize', 6);
        end
    end
end

function plot_gauss_points(ax, mesh, partition)
    v_x = mesh.get_1d_physical(partition.quadrature.local, partition.faces);
    for e = 1:size(v_x, 3)
        plot(ax, v_x(1, :, e), v_x(2, :, e), 'bo', 'LineWidth', 1.2);
    end
end

function plot_gauss_points_paths(ax, mesh, partition)
    v_x = partition.quadrature.path;
    for e = 1:size(v_x, 3)
        x = v_x(1, :, e, :);
        y = v_x(2, :, e, :);
        plot(ax, x(:), y(:), 'rs', 'LineWidth', 1.2, 'MarkerSize', 4);
    end
end