function plot(obj, ax, options)
%PLOT  Draw a parametrized curve.
%
%   plot(curve) draws the curve on the current axes; plot(curve, ax) on ax.
%
%   Name-value options (defaults), forwarded to the plotting call:
%     LineWidth  (2.0)
%     Color      ('k')
%
%See also: plot_transfer_paths, SimplicialGrid/plot
arguments
    obj ParametrizedCurve
    ax = gca();
    options.LineWidth = 2.0;
    options.Color = 'k';
end
t = linspace(obj.t_range(1), obj.t_range(2));
args = namedargs2cell(options);
plot(ax, obj.x_coord(t), obj.y_coord(t), args{:})
end
