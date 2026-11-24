%% PiecewiseSplineCurve: an open three-piece curve (a rectangle's bottom,
% right and top sides, bent by a cubic displacement). Not-a-knot cubic
% splines reproduce cubics exactly, so position, tangent and arc length
% are checked against the closed forms; the normal must point out of the
% rectangle, and the closest-point mapping must recover the parameter.
clc; clear;
n_fail = 0;
check = @(name, err, tol) fprintf("%-40s %.3e  %s\n", name, err, ...
    string(ifelse(err <= tol, "ok", "FAILED")));

a = 0.25;  b = 0.6;  y0 = 0.19;  y1 = 0.21;  c = 2e-2;
bend = @(x) c*(x - a).^2.*(x - a + 0.1);          % cubic in x, zero at x = a
dbend = @(x) c*(2*(x - a).*(x - a + 0.1) + (x - a).^2);
sb = linspace(a, b, 5);  st = linspace(y0, y1, 4);  su = linspace(b, a, 6);
pieces = struct( ...
    "s", {sb, st, su}, ...
    "P", {[sb; y0 + bend(sb)], [b + 0*st; st + bend(b)], [su; y1 + bend(su)]});
curve = PiecewiseSplineCurve(pieces);

% exact parametrization: piece i, t in [i-1, i] -> s affine
exact = {@(s) [s; y0 + bend(s)], @(s) [b + 0*s; s + bend(b)], @(s) [s; y1 + bend(s)]};
dexact = {@(s) [1 + 0*s; dbend(s)], @(s) [0*s; 1 + 0*s], @(s) [1 + 0*s; dbend(s)]};
ends = [a b; y0 y1; b a];

t = linspace(0, 3, 601);
[e_pos, e_tan] = deal(0);
for i = 1:3
    ti = t(t >= i - 1 & t <= i);
    ti = ti(ti < i | i == 3);
    s = ends(i, 1) + (ti - (i - 1))*(ends(i, 2) - ends(i, 1));
    e_pos = max(e_pos, max(abs(curve.get_2d_points(ti) - exact{i}(s)), [], "all"));
    e_tan = max(e_tan, max(abs(curve.tangent(ti) - dexact{i}(s)*(ends(i, 2) - ends(i, 1))), [], "all"));
end
check("position vs exact cubic", e_pos, 1e-13);   n_fail = n_fail + (e_pos > 1e-13);
check("tangent vs exact derivative", e_tan, 1e-12); n_fail = n_fail + (e_tan > 1e-12);

% arc length: closed-form speed integrated per piece
L_exact = integral(@(x) sqrt(1 + dbend(x).^2), a, b, AbsTol = 1e-15)*2 + (y1 - y0);
e_len = abs(curve.arc_length(0, 3) - L_exact);
check("arc length of the whole curve", e_len, 1e-11); n_fail = n_fail + (e_len > 1e-11);
e_sub = abs(curve.arc_length(0.5, 2.5) - curve.arc_length(0.5, 1) ...
    - curve.arc_length(1, 2) - curve.arc_length(2, 2.5));
check("arc length additivity across joins", e_sub, 1e-12); n_fail = n_fail + (e_sub > 1e-12);

% normals out of the rectangle: down on the bottom, right at the tip, up on top
n = curve.unit_normal([0.5, 1.5, 2.5]);
e_n = max(abs(n - [0 1 0; -1 0 1]), [], "all");
check("outward normals (bend-tilted)", e_n, 0.05);   n_fail = n_fail + (e_n > 0.05);
e_unit = max(abs(vecnorm(curve.unit_normal(t)) - 1));
check("unit normal length", e_unit, 1e-14);          n_fail = n_fail + (e_unit > 1e-14);

% closest-point mapping recovers the parameter of points on the curve
tq = [0.13 0.77 1.4 2.05 2.9];
tq_back = curve.map_points_to_interface(curve.get_2d_points(tq));
e_map = max(abs(tq_back(:).' - tq));
check("closest-point mapping", e_map, 1e-6);         n_fail = n_fail + (e_map > 1e-6);

% the corners: continuous position, jumping tangent
jump = norm(curve.get_2d_points(1 - 1e-12) - curve.get_2d_points(1));
check("continuity at the first corner", jump, 1e-11); n_fail = n_fail + (jump > 1e-11);

% closed: the same three pieces closed by the left side x = a (bend(a) = 0)
sl = linspace(y1, y0, 3);
loop = PiecewiseSplineCurve([pieces, struct("s", sl, "P", [a + 0*sl; sl])], closed = true);
e_per = max([norm(loop.get_2d_points(4) - loop.get_2d_points(0)), ...
             norm(loop.get_2d_points(-0.5) - loop.get_2d_points(3.5)), ...
             norm(loop.get_2d_points(4.3) - loop.get_2d_points(0.3))]);
check("closed: parameter wraps", e_per, 1e-14);      n_fail = n_fail + (e_per > 1e-14);
fail_periodic = ~loop.is_periodic() || curve.is_periodic();
check("closed: is_periodic (open: false)", double(fail_periodic), 0);  n_fail = n_fail + fail_periodic;
e_open = max(abs(curve.get_2d_points(t) - PiecewiseSplineCurve(pieces).get_2d_points(t)), [], "all");
check("open default unchanged by the option", e_open, 0);  n_fail = n_fail + (e_open > 0);
e_len = abs(loop.arc_length(0, 4) - (L_exact + (y1 - y0)));
check("closed: perimeter", e_len, 1e-11);            n_fail = n_fail + (e_len > 1e-11);
e_n = max(abs(loop.unit_normal(3.5) - [-1; 0]));
check("closed: outward normal on the closing side", e_n, 1e-12);  n_fail = n_fail + (e_n > 1e-12);
tq = [0.02 1.5 3.4 3.97];
tq_back = loop.map_points_to_interface(loop.get_2d_points(tq) + 1e-4*loop.unit_normal(tq));
e_map = max(abs(tq_back(:).' - tq));
check("closed: closest point across the seam", e_map, 1e-6);  n_fail = n_fail + (e_map > 1e-6);

fprintf("\n%d check(s) failed\n", n_fail);

function v = ifelse(c, a, b)
if c, v = a; else, v = b; end
end
