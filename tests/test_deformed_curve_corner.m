% test_deformed_curve_corner.m
% DeformedCurve.find_closest_point at a corner of the curve: a point in
% the normal cone of a convex corner (beyond both one-sided normals) has
% the corner as its closest point. The Gauss-Newton polish after the
% global search used to step along a one-sided tangent past the kink and
% walk away (the carved FSI1 fluid: the tip-corner vertex (0.60333,
% 0.18667) of the flag's sides mapped to t = 0.8333, x = 0.54, instead of
% the corner t = 1, and get_partitions reported "Mesh 1 folds over the
% interface").
%   1. the flag's free sides (open polygon (0.24899, 0.19), (0.6, 0.19),
%      (0.6, 0.21), (0.24899, 0.21)) as a DeformedCurve with d = 0: points
%      off both corners, at offsets 0.001 ... 0.005 along directions
%      inside the normal cones, map to t = 1 and t = 2 (1e-7: the global
%      search's fminbnd stops at about sqrt(eps) |t| on the kink, where
%      the polish cannot improve it; before the fix the error was 0.25);
%      points
%      over the sides map to their feet (1e-10).
%   2. the same with a small smooth d (the corners move): points off the
%      moved corners map to t = 1 and 2.
n_fail = 0;
xc = 0.2 + sqrt(0.05^2 - 0.01^2);
sides = PolygonalCurve([xc 0.6 0.6 xc; 0.19 0.19 0.21 0.21], closed = false);
for variant = 1:2
    if variant == 1
        d = @(t) zeros(2, numel(t));
        dp = @(t) zeros(2, numel(t));
    else
        d = @(t) 1e-3*[sin(t); t.^2/4];
        dp = @(t) 1e-3*[cos(t); t/2];
    end
    curve = DeformedCurve(sides, d, parametric = true, derivative = dp);
    corners = curve.get_2d_points([1 2]);
    err = 0;
    for r = [0.001 0.0033 0.005]
        for angle = [10 45 80]
            % bottom-right corner: cone between (1, 0) and (0, -1)
            p = corners(:, 1) + r*[cosd(angle); -sind(angle)];
            err = max(err, abs(curve.map_points_to_interface(p) - 1));
            % top-right corner: cone between (1, 0) and (0, 1)
            p = corners(:, 2) + r*[cosd(angle); sind(angle)];
            err = max(err, abs(curve.map_points_to_interface(p) - 2));
        end
    end
    tq = [0.3 0.7 1.4 1.6 2.3 2.8];
    N = reshape(curve.unit_normal(tq), 2, []);
    feet = curve.map_points_to_interface(curve.get_2d_points(tq) + 0.002*N);
    err_sides = max(abs(feet - tq));
    fprintf("variant %d: corner points max |t - t_corner| %.1e; side points max |t - t_foot| %.1e\n", ...
        variant, err, err_sides);
    n_fail = n_fail + check(err <= 1e-7, sprintf("variant %d: corners", variant));
    n_fail = n_fail + check(err_sides <= 1e-10, sprintf("variant %d: sides", variant));
end
% the FSI1 vertex itself
curve = DeformedCurve(sides, @(t) zeros(2, numel(t)), parametric = true, derivative = @(t) zeros(2, numel(t)));
t = curve.map_points_to_interface([0.60333; 0.18667]);
fprintf("FSI1 tip vertex (0.60333, 0.18667): t = %.6f (corner 1)\n", t);
n_fail = n_fail + check(abs(t - 1) <= 1e-7, "FSI1 tip vertex");
if n_fail == 0
    fprintf("test_deformed_curve_corner: all checks passed\n");
else
    fprintf("test_deformed_curve_corner: %d check(s) failed\n", n_fail);
end
