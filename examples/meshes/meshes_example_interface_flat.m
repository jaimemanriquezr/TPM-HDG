% TPM geometry on a flat interface y = y0: two rectangles separated by a
% gap of 2d, meshed with dissimilar sizes. Plots both meshes, the induced
% partition vertices, transfer paths and quadrature points on each side.
    h   = 1/16;
    d = 1e-1;
    y0  = .5;
close all;
I = FlatCurveHorizontal(y0);

lower_rectangle = [[0; 0], [1; 0], [1; y0-d], [0; y0-d]];
lower_conn = [[[1; 2]; 1], [[2; 3]; 1], [[3; 4]; 3], [[4; 1]; 1]];
lower_domain = {lower_rectangle; lower_conn};
lower_mesh = SimplicialGrid(lower_domain{:}, ...
        element_diameter = h, ...
        mesh_generator = "mesh2d");

upper_rectangle = [[0; y0+d], [1; y0+d], [1; 1], [0; 1]];
upper_conn = [[[1; 2]; 3], [[2; 3]; 1], [[3; 4]; 1], [[4; 1]; 1]];
upper_domain = {upper_rectangle; upper_conn};
upper_mesh = SimplicialGrid(upper_domain{:}, ...
        element_diameter = sqrt(2)*h, ...
        mesh_generator = "mesh2d");

ax = axes();
ax.NextPlot = 'add';
plot(upper_mesh, ax);
plot(lower_mesh, ax);   

plot(I, ax, 'LineWidth', 2);
    
[upper_partition, lower_partition] = get_partitions(I, upper_mesh, lower_mesh);
plot_vertices(ax, upper_mesh, upper_partition)
plot_paths(ax, upper_mesh, upper_partition)
plot_patch(ax, I, upper_mesh, upper_partition)
plot_gauss_points(ax, upper_mesh, upper_partition)
plot_gauss_points_paths(ax, upper_mesh, upper_partition)

plot_vertices(ax, lower_mesh, lower_partition)
plot_paths(ax, lower_mesh, lower_partition)
plot_patch(ax, I, lower_mesh, lower_partition)
plot_gauss_points(ax, lower_mesh, lower_partition)
plot_gauss_points_paths(ax, lower_mesh, lower_partition)



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
