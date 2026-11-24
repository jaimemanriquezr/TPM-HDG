function plot(obj, ax, options)
%PLOT  Draw a level-set curve.
%
%   plot(curve) draws the zero level set phi = 0 of the curve on the
%   current axes with fimplicit, over curve.box; plot(curve, ax) uses ax.
%
%   Name-value options (defaults), forwarded to the plotting call:
%     LineWidth  (2.0)
%     Color      ('k')
%
%See also: plot_transfer_paths, SimplicialGrid/plot
arguments
    obj LevelSetCurve
    ax = gca();
    options.LineWidth = 2.0;
    options.Color = 'k';
end
args = namedargs2cell(options);
fimplicit(ax, @(x, y) reshape(obj.evaluate([x(:).'; y(:).']), size(x)), ...
    obj.box, args{:})
end
