function plot(obj, mesh, options)
%PLOT  Plot a Datum over the bounding box of a mesh.
%
%   plot(datum, mesh) samples datum on an n_grid_points x n_grid_points
%   grid over the bounding box of mesh and draws each component as a
%   surface seen from above, one tile per component of a tiledlayout.
%
%   Name-value options (defaults):
%     parent          tiledlayout to draw in (the one in the current
%                     figure, or a new 'flow' layout when there is none)
%     n_grid_points   grid points per direction (50)
%     colormap        colormap of each tile ("jet")
%     name            variable name for the titles (the datum's name);
%                     empty gives no title
%
%See also: FEField/plot, FEFunction/plot, SimplicialGrid/plot
arguments
    obj Datum
    mesh AbstractGrid
    options.parent = findobj(gcf().Children, 'Type', 'tiledlayout');
    options.n_grid_points = 50;
    options.colormap = "jet";
    options.name = obj.name;
end
if isempty(options.parent)
    parent = tiledlayout('flow');
else
    parent = options.parent;
end

x1 = min(mesh.coordinates(1, :));
x2 = max(mesh.coordinates(1, :));

y1 = min(mesh.coordinates(2, :));
y2 = max(mesh.coordinates(2, :));

[X, Y] = meshgrid(linspace(x1, x2, options.n_grid_points), ...
                  linspace(y1, y2, options.n_grid_points));
Z = obj.evaluate([X(:), Y(:)].', coordinate_system = "physical");

out_sz = size(Z, 3:4);
for i = 1:out_sz(1)
    for j = 1:out_sz(2)
        if isa(parent, 'matlab.graphics.layout.TiledChartLayout')
            ax = nexttile(parent);
        else
            ax = parent;
        end

        surf(ax, X, Y, reshape(Z(:, :, i, j), size(X)));
        if ~isempty(options.name)
            if out_sz(2) > 1
                title(ax, sprintf("$%s^{%i%i}$", options.name, i, j), ...
                'Interpreter', 'latex');
            elseif out_sz(1)*out_sz(2) ~= 1
                title(ax, sprintf("$%s^{%i}$", options.name, i), ...
                'Interpreter', 'latex');
            else
                title(ax, sprintf("$%s$", options.name), ...
                'Interpreter', 'latex');
            end
        end
        view(ax, 2)

        shading(ax, 'interp')
        colormap(ax, options.colormap)
    end
end
end