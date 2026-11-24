function p = plot(obj, y, options)
%PLOT  Draw a SimplicialGrid, or an FE function on it.
%
%   p = plot(mesh) draws the triangles of mesh as patches on the current
%   axes, with a 2-D view.
%   p = plot(mesh, ax) draws on the axes ax.
%   plot(mesh, y) with y an FEFunction plots y on mesh through its own
%   plot method, instead of drawing the grid.
%
%   Name-value options (defaults), applied to the patch:
%     Parent       target axes (gca())
%     ElementList  indices of the elements drawn (all)
%     LineWidth    (0.8)
%     LineStyle    ('-')
%     EdgeColor    ([0 0 0])
%     FaceColor    ('w')
%     FaceAlpha    (1.0)
%
%   p is the patch object (not set when y is given).
%
%See also: FEField/plot, plot_transfer_paths
arguments
    obj SimplicialGrid
    y = FEFunction.empty;
    options.Parent = gca()
    options.ElementList = 1:size(obj.elements, 2);
    options.LineWidth = .8;
    options.LineStyle = '-';
    options.EdgeColor = [0 0 0];
    options.FaceColor = 'w';
    options.FaceAlpha = 1.0;
end
% plot(mesh, ax, ...) : axes passed positionally
if isa(y, 'matlab.graphics.axis.Axes')
    options.Parent = y;
    y = FEFunction.empty;
end
ax = options.Parent;

ax.NextPlot = 'add';
axes(ax);
plot_options = namedargs2cell(rmfield(options, "ElementList"));

connectivity = obj.elements(:, options.ElementList).';

if isempty(y)
    p = patch(...
        "faces", connectivity, ...
        "vertices", obj.coordinates.', ...
        plot_options{:} ...
        );
else
    plot(y, obj);
    ax = gca();
end
view(ax, 2)
end
