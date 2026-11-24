% test_deformed_curve_parametric.m
% DeformedCurve with parametric = true and a derivative handle, as
% lagrangian_carved_march builds Sigma_h (x(t) = X(t) + d(t), tangent
% X'(t) + d'(t)), against analytic deformed curves:
%   1. closed: the circle of radius 0.4 about (0.5, 0.5) under a linear
%      map (rotation and stretch), a translation and a radial ripple of
%      three waves: position, tangent, outward unit normal, Nanson speed
%      ratio |x'|/|X'| to 1e-10, arc length to the integral's 1e-6, and
%      the closest-point map (points on the curve and 0.005 off it on
%      either side, also next to the seam t = 0 = 1) to 1e-10 in t;
%   2. open: y = 0.5 under a sinusoidal bump with a tangential component
%      vanishing at the walls, as a WallEndedCurve: the same quantities
%      (upward normal), and points on the walls x = 0, 1 map to the ends;
%   3. the march's construction: an affine solid displacement projected
%      on a structured solid mesh below y = 0.5 (exact at k = 1) and
%      evaluated on Sigma_0 by solid_interface_displacement (piecewise
%      d, finite-difference d'): position, tangent and normal of the
%      affine image to 1e-10.
clc; clear;

n_fail = 0;
t = linspace(0, 1, 201);

%% 1. closed curve
c = [0.5; 0.5];  R = 0.4;
M = [1.1, -0.2; 0.15, 0.9];  b = [0.01; -0.02];  ripple = 0.01;
X = @(t) c + R*[cos(2*pi*t); sin(2*pi*t)];
X_prime = @(t) 2*pi*R*[-sin(2*pi*t); cos(2*pi*t)];
x = @(t) c + b + M*(X(t) - c) + ripple*sin(6*pi*t).*(X(t) - c)/R;
x_prime = @(t) M*X_prime(t) + ripple*(6*pi*cos(6*pi*t).*(X(t) - c) + sin(6*pi*t).*X_prime(t))/R;
turn = @(v, sense) sense*[-v(2, :); v(1, :)]./vecnorm(v);   % the tangent turned by sense*90 degrees
outward = @(t) turn(x_prime(t), -1);                        % counterclockwise circle: by -90 degrees
curve = DeformedCurve(CircularCurve(c, R), @(t) x(t) - X(t), parametric = true, ...
    derivative = @(t) x_prime(t) - X_prime(t));
n_fail = n_fail + compare(curve, t, x, x_prime, outward, 2*pi*R, "closed", [], []);
tq = [linspace(0.01, 0.99, 25), 0.9995, 3e-4];
for delta = [0, 0.005, -0.005]
    tq_found = curve.map_points_to_interface(x(tq) + delta*outward(tq));
    e = max(abs(mod(tq_found - tq + 0.5, 1) - 0.5));
    fprintf("closed: closest point, offset %+.3f: max |t - t_exact| %.1e\n", delta, e);
    n_fail = n_fail + check(e < 1e-10, sprintf("closed: closest-point map off by %.1e at offset %g", e, delta));
end

%% 2. open curve with ends on the walls
X = @(t) [t; 0.5 + 0*t];
x = @(t) X(t) + [0.02*sin(2*pi*t); 0.05*sin(pi*t) + 0.01];
x_prime = @(t) [1 + 0.04*pi*cos(2*pi*t); 0.05*pi*cos(pi*t)];
upward = @(t) turn(x_prime(t), 1);
curve = WallEndedCurve(FlatCurveHorizontal(0.5), @(t) x(t) - X(t), [0 1], "parametric", true, ...
    "derivative", @(t) x_prime(t) - [1; 0]);
ends = x([0, 1]);
walls = [0, 1; ends(2, 1) + 0.003, ends(2, 2) - 0.004];
n_fail = n_fail + compare(curve, t, x, x_prime, upward, 1, "open", walls, [0, 1]);
tq = linspace(0.05, 0.95, 19);
for delta = [0, 0.005, -0.005]
    e = max(abs(curve.map_points_to_interface(x(tq) + delta*upward(tq)) - tq));
    fprintf("open: closest point, offset %+.3f: max |t - t_exact| %.1e\n", delta, e);
    n_fail = n_fail + check(e < 1e-10, sprintf("open: closest-point map off by %.1e at offset %g", e, delta));
end

%% 3. the march's construction from solid dofs
k = 1;  N = 8;
solid_mesh = apply_tag(StructuredSimplicialGrid(1.5*N, 0.75*N, 0, 1, 0, 0.5), ...
    "Interface", @(P) abs(P(2, :) - 0.5) < 1e-14);
M = [1.05, 0.1; -0.08, 0.97];  b = [0; 0.01];
u = {@(x, y) (M(1, 1) - 1)*x + M(1, 2)*y + b(1) + 0*x; @(x, y) M(2, 1)*x + (M(2, 2) - 1)*y + b(2) + 0*x};
dofs = L2proj2D(solid_mesh, u, k);
reference = FlatCurveHorizontal(0.5);
[d, d_prime] = solid_interface_displacement(reference, solid_mesh, dofs, k);
curve = DeformedCurve(reference, d, parametric = true, derivative = d_prime);
x = @(t) M*[t; 0.5 + 0*t] + b;
x_prime = @(t) M(:, 1) + 0*t;
upward = @(t) turn(x_prime(t), 1);
n_fail = n_fail + compare(curve, t, x, x_prime, upward, 1, "from solid dofs", [], []);

fprintf("%d check(s) failed.\n", n_fail);
if n_fail > 0
    error("test_deformed_curve_parametric:failed", "%d check(s) failed", n_fail);
end

function n_fail = compare(curve, t, x, x_prime, normal, reference_speed, name, wall_points, wall_t)
% position, tangent, unit normal, Nanson ratio and arc length of curve
% against the analytic ones; points on the walls mapped to wall_t
e = [max(abs(curve.get_2d_points(t) - x(t)), [], "all"), ...
     max(abs(curve.tangent(t) - x_prime(t)), [], "all"), ...
     max(abs(reshape(curve.unit_normal(t), 2, []) - normal(t)), [], "all"), ...
     max(abs(vecnorm(curve.tangent(t))./vecnorm(curve.reference.tangent(t) + 0*t) ...
         - vecnorm(x_prime(t))/reference_speed))];
L = integral(@(s) vecnorm(x_prime(s)), 0, 1, AbsTol = 1e-14, RelTol = 1e-13);
e_length = abs(curve.arc_length(0, 1) - L)/L;
fprintf("%s: position %.1e, tangent %.1e, normal %.1e, Nanson %.1e; arc length %.6f (relative error %.1e)\n", ...
    name, e, L, e_length);
n_fail = check(max(e) < 1e-10, name + ": position / tangent / normal / Nanson ratio off the analytic ones");
n_fail = n_fail + check(e_length < 1e-6, name + ": arc length off the analytic one");
if ~isempty(wall_points)
    found = curve.map_points_to_interface(wall_points);
    fprintf("%s: wall points map to t = %s\n", name, mat2str(found));
    n_fail = n_fail + check(isequal(found, wall_t), name + ": wall points not mapped to the ends");
end
end
