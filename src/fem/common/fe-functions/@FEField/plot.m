function parent = plot(obj, mesh, options)
%PLOT  Plot a discrete field on a mesh.
%
%   parent = plot(field, mesh) draws each component of field as a
%   coloured surface over mesh, one tile per component of a tiledlayout.
%   mesh defaults to the field's own mesh.
%
%   Name-value options (defaults):
%     name        variable name for the titles (the field's name); empty
%                 gives no title
%     parent      tiledlayout to draw in (the one in the current figure,
%                 or a new 'flow' layout when there is none)
%     components  which components of a vector- or tensor-valued field are
%                 drawn ("all"):
%                   "all"    every component, one tile each, in parent
%                   "tiled"  every component, in a nested layout of its
%                            own, so the variable stays one block of
%                            parent
%                   "first"  the (1,1) component alone
%
%   parent is the layout drawn in.
%
%See also: Datum/plot, FEFunction/plot, SimplicialGrid/plot
    arguments
        obj FEField
        mesh AbstractGrid = obj.mesh;

        options.name = obj.name;
        options.parent = findobj(gcf().Children, 'Type', 'tiledlayout');

        % Which components of a vector- or tensor-valued field are drawn,
        % and how they are arranged:
        %   "all"    every component, one tile each, in the parent layout
        %   "tiled"  every component, in a nested layout of its own, so
        %            one variable stays one block of the parent layout
        %   "first"  the (1,1) component alone
        options.components (1,1) string ...
            {mustBeMember(options.components, ["all", "tiled", "first"])} = "all"
    end
    if isempty(options.parent)
        parent = tiledlayout('flow');
    else
        parent = options.parent;
    end
    geo_dim = obj.basis.GeometricInputDimension;
    vertices_reference = reference_element(geo_dim, mesh.ElementType);
    P = obj.evaluate(vertices_reference, coordinate_system = "vertex");
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
    
    is_scalar = isa(obj.basis, 'ScalarFEBasis');
    if is_scalar
        [n, m] = deal(1, 1);
    else
        % component shape from the basis (not size(obj), the array size)
        shape = obj.basis.GeometricOutputDimension;
        if isscalar(shape); shape = [shape, 1]; end
        [n, m] = deal(shape(1), shape(2));
        if options.components == "first"
            [n, m] = deal(1, 1);
        end
    end

    % One nested layout per variable, so its components stay together as a
    % single block of the parent layout. A nested layout does not flow (it
    % defaults to tile 1), so it is placed after whatever the parent holds.
    target = parent;
    if options.components == "tiled" ...
            && isa(parent, 'matlab.graphics.layout.TiledChartLayout')
        occupied = numel(parent.Children);
        target = tiledlayout(parent, n, m, ...
            TileSpacing = "compact", Padding = "compact");
        target.Layout.Tile = occupied + 1;
    end

    if is_scalar
        if isempty(options.name)
            name_str = string.empty;
        else
            name_str = "$" + options.name + "$";
        end
        plot_scalar(P, x, y, name_str, name_str, target)
    else
        c_str = {string(1:n), string(1:m)};
        for c_ix = find(cellfun(@(str) isequal(str, "1"), c_str))
            c_str{c_ix} = "";
        end
        component_str = c_str{1}.' + c_str{2};

        for i = 1:n
            for j = 1:m
                k = sub2ind([n m], i, j);
                if isempty(options.name)
                    % No title, but the tag must still separate the
                    % components: axes are reused by tag.
                    name_str = string.empty;
                    tag = "component_" + component_str(k);
                else
                    name_str = "$" + options.name + "^{" + component_str(k) + "}$";
                    tag = name_str;
                end
                plot_scalar(P(:, :, i, j, :), x, y, name_str, tag, target)
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
    elseif isa(parent, 'matlab.graphics.layout.TiledChartLayout')
        parent = nexttile(parent);
        parent.Tag = tag;
    else
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
