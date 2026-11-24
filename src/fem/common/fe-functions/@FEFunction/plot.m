function parent = plot(obj, mesh, options)
%PLOT  Plot an FEFunction on a mesh.
%
%   parent = plot(f, mesh) evaluates f at the vertices of the elements of
%   mesh and draws each component as a coloured surface, one tile per
%   component of a tiledlayout.
%
%   Name-value options (defaults):
%     parent  tiledlayout to draw in (the one in the current figure, or a
%             new 'flow' layout when there is none)
%     name    variable name for the titles (empty: no title)
%
%   parent is the layout drawn in.
%
%See also: FEField/plot, Datum/plot, SimplicialGrid/plot
    arguments
        obj FEFunction
        mesh AbstractGrid

        options.parent = findobj(gcf().Children, 'Type', 'tiledlayout');
        options.name = string.empty;
    end

    if isempty(options.parent)
        parent = tiledlayout('flow');
    else
        parent = options.parent;
    end
    [~, geo_dim] = read_signature(obj);
    vertices_reference = reference_element(geo_dim, mesh.ElementType);
    P = obj.evaluate(vertices_reference, mesh);
    switch geo_dim
        case 1
            plot_scalar = @plot_scalar_1d;
            x = reshape(mesh.coordinates(1, mesh.faces(1:2,:)), 2, []);
            y = reshape(mesh.coordinates(2, mesh.faces(1:2,:)), 2, []);
        case 2
            plot_scalar = @plot_scalar_2d;
            n_v = size(mesh.elements, 1);   % vertices per element
            x = reshape(mesh.coordinates(1, mesh.elements), n_v, []);
            y = reshape(mesh.coordinates(2, mesh.elements), n_v, []);
    end

    if isequal(obj.GeometricOutputDimension, [1 1])
        if isempty(options.name)
            name_str = string.empty;
        else
            name_str = "$" + options.name + "$";
        end
        plot_scalar(P, x, y, name_str, name_str, parent)
    else
        n = obj.GeometricOutputDimension(1);
        m = obj.GeometricOutputDimension(2);
        c_str = {string(1:n), string(1:m)};
        [c_str{cellfun(@(str) isequal(str, "1"), c_str)}] = "";
        component_str = c_str{1}.' + c_str{2};

        for i = 1:n
            for j = 1:m
                k = sub2ind([n m], i, j);
                if isempty(options.name)
                    name_str = string.empty;
                    tag = "component_" + component_str(k);
                else
                    name_str = "$" + options.name + "^{" + component_str(k) + "}$";
                    tag = name_str;
                end
                plot_scalar(P(:, :, i, j, :), x, y, name_str, tag, parent)
            end
        end
    end
end

function plot_scalar_1d(p, x, y, var_name, tag, parent)
    assert(size(x,2) == size(y,2))

    if isa(parent,'matlab.graphics.layout.TiledChartLayout')
        if isempty(findobj(parent, 'Type', 'axes', 'Tag', tag))
            parent = nexttile(parent);
            parent.Tag = tag;
        end
    end
    
    parent.NextPlot = 'add';
    for e = 1:size(x, 2)
        x_p = [x(:, e) x(:,e)];
        y_p = [y(:, e) y(:,e)];
        z_p = [p(:, e) p(:,e)];
        c_p = z_p;
        surface(x_p, y_p, z_p, c_p, ...
            'EdgeColor', 'interp', 'FaceColor', 'none', ...
            'Parent', parent);
    end
    c = colorbar(parent);
    axis tight; 
    colormap(parent, 'jet');
    title(parent, var_name, 'interpreter','latex','fontsize',14)
    set([parent c], "TickLabelInterpreter", 'latex');
    set(parent, "color", "none");
end

function plot_scalar_2d(p, x, y, var_name, tag, parent)
    assert(size(x,2) == size(y,2))

    existing_axes = findobj(parent, 'Type', 'axes', 'Tag', tag);
    if ~isempty(existing_axes)
        parent = existing_axes;
    else
        parent = nexttile(parent);
        parent.Tag = tag;
    end

    parent.NextPlot = 'add';
    for K = 1:size(x, 2)
        x_p = x(:, K);
        y_p = y(:, K);
        z_p = p(:, :, K);
        % the patch trisurf(1:3, ...) draws, for any number of vertices
        patch('Faces', 1:numel(x_p), 'Vertices', [x_p y_p z_p(:)], ...
            'FaceVertexCData', z_p(:), ...
            'FaceColor', get(parent, 'DefaultSurfaceFaceColor'), ...
            'EdgeColor', get(parent, 'DefaultSurfaceEdgeColor'), ...
            'Parent', parent);
    end
    c = colorbar(parent);
    shading(parent,'flat'); 
    axis tight; 
    colormap(parent, 'jet');
    title(parent, var_name, 'interpreter','latex','fontsize',14);
    set([parent c], "TickLabelInterpreter", 'latex');
    set(parent, "color", "none");
end
