function plot_transfer_paths(ax, interface, mesh, partition, options)
%-------------------------------------------------------------------------
%Description: Draws the TPM induced partition of one mesh on an existing
%             axes: patch endpoints on the mesh boundary, the transfer
%             paths joining them to the interface, the patch quadrature
%             nodes, and the path quadrature nodes.
%
%             Legend (defaults; each is overridable):
%               k+  patch endpoints on the mesh boundary  EndpointColor
%                                                           EndpointMarker
%               k-  transfer paths (endpoint -> interface)  PathColor
%               ko  quadrature landings on the interface    InterfaceColor
%               bo  patch quadrature nodes                  QuadratureColor
%                                                           QuadratureMarker
%               rs  path quadrature nodes                   PathQuadratureColor
%                                                           PathQuadratureMarker
%-------------------------------------------------------------------------
arguments (Input)
    ax                              matlab.graphics.axis.Axes
    interface   (1,1)               AbstractCurve
    mesh                            SimplicialGrid
    partition   (1,1)               struct
    options.PathColor                        = 'k'   % transfer paths
    options.InterfaceColor                   = 'k'   % landings on Sigma
    options.EndpointColor                    = 'k'   % patch ends on the mesh
    % Distinct from the interface landings by default
    options.EndpointMarker       (1,1) string = "+"
    options.QuadratureColor                  = 'b'   % patch quadrature nodes
    options.PathQuadratureColor              = 'r'   % path quadrature nodes
    options.QuadratureMarker     (1,1) string = "o"
    options.PathQuadratureMarker (1,1) string = "s"
    options.LineWidth            (1,1) double = 1.2
end
w = options.LineWidth;
% Accumulate onto the caller's axes; restore its hold state on exit.
held_before = ishold(ax);
hold(ax, "on")
restore_hold = onCleanup(@() hold_restore(ax, held_before));

% Patch endpoints on the mesh boundary, and the paths to the interface.
v_x = mesh.get_1d_physical(partition.vertices, partition.faces);
% Use the patch's own interface parameters, shared by both meshes and used
% by the solver, not a closest-point projection of the mesh endpoint.
v_I = interface.get_2d_points(partition.vertices_interface);
for e = 1:size(v_x, 3)
    plot(ax, v_x(1, :, e), v_x(2, :, e), options.EndpointMarker, ...
        'Color', options.EndpointColor, 'LineWidth', w);
    for j = 1:size(v_x, 2)
        tp = [v_x(:, j, e), v_I(:, j, e)];
        plot(ax, tp(1, :), tp(2, :), '-', 'Color', options.PathColor, 'LineWidth', w);
    end
end

% Where the quadrature nodes land on the interface.
x_interface = partition.quadrature.interface;
for e = 1:size(x_interface, 3)
    plot(ax, x_interface(1, :, e), x_interface(2, :, e), 'o', ...
        'Color', options.InterfaceColor, 'MarkerSize', 6);
end

% Quadrature nodes: on the mesh boundary, and along the transfer paths.
x_quad = mesh.get_1d_physical(partition.quadrature.local, partition.faces);
x_path = partition.quadrature.path;
for e = 1:size(x_quad, 3)
    plot(ax, x_quad(1, :, e), x_quad(2, :, e), options.QuadratureMarker, ...
        'Color', options.QuadratureColor, 'LineWidth', w);
    x = x_path(1, :, e, :);
    y = x_path(2, :, e, :);
    plot(ax, x(:), y(:), options.PathQuadratureMarker, ...
        'Color', options.PathQuadratureColor, 'LineWidth', w, 'MarkerSize', 4);
end
end

function hold_restore(ax, was_held)
if ~was_held && isvalid(ax)
    hold(ax, "off")
end
end
