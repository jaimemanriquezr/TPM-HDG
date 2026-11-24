% TPM geometry on a circular interface: fitted meshes of the unit square
% with a hole and of the disk. Plots both meshes, the induced partition
% vertices, transfer paths and quadrature points on each side.
clc; clear;
h   = 1/11;
d = 1e-1;
y0  = .5;
close all;
I = CircularCurve([.4; .5], .3);

[outer_mesh, inner_mesh] = fitted_split(I, h, 0);

ax = axes();
ax.NextPlot = 'add';
axis(ax, "equal")
plot(outer_mesh, ax, 'FaceAlpha', 0.5);
plot(inner_mesh, ax, 'FaceAlpha', 0.5);   

plot(I, ax, 'LineWidth', 4, 'Color', 'm');
    
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
            plot(ax, tp(1, :), tp(2, :), 'k-', 'LineWidth', 1.2);
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
