function draw_mesh_couplings(h, interface, splitting, ax)
% Plots a mesh pair produced by `splitting`, the interface and the
% transfer paths of both induced partitions on the axes `ax`.
arguments
    h double = 0.1;
    interface ParametrizedCurve = SmoothCurveHorizontal()
    splitting function_handle = @(I, h) UnitSquareMesh(h).split(I, cut_elements = "mesh_1")
    ax = axes();
end
[mesh_1, mesh_2] = splitting(interface, h);

ax.NextPlot = 'add';
axis(ax, "square");

plot(mesh_1, ax, FaceColor='g', FaceAlpha=0.5);
plot(mesh_2, ax, FaceColor='y', FaceAlpha=0.5);
plot(interface, ax, LineWidth=4, Color='k');

[partition_1, partition_2] = get_partitions(interface, mesh_1, mesh_2);
plot_transfer_paths(ax, interface, mesh_1, partition_1);
plot_transfer_paths(ax, interface, mesh_2, partition_2);
end
