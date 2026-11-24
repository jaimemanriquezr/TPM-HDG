% test_eulerian_interface_reconstruction.m
% eulerian_interface_displacement and interface_curve_metrics (tests/helpers):
% the interface of an Eulerian solid on Sigma_0, D(s) = x_h(s) - X(s) with
% x_h - rigid - u_K(x_h) = X(s), on the solid side of a carved background
% grid ("left" grid, Sigma = {y = 0.45 + d}; "crossed" grid, the circle of
% radius 0.2 centred at (0.5 + d, 0.5)), d = 0.037, rigid = d e_y / d e_x.
%   a. polynomial reproduction: a global polynomial u of degree p = k + 1
%      (k = 1, 2), projected elementwise onto Dubiner degree p (exact),
%      sliding the material tangentially by 0.08 (more than h = 1/16) with
%      a gradient of size 0.3-0.4: x_h equals the exact solution of
%      x - rigid - u(x) = X(s) (fixed point with the exact u) to 1e-10 at
%      about 200 parameters (breaks +- 1e-9 included), the jumps are below
%      1e-11 and D_prime equals (I - grad u)^-1 X' - X' to 1e-9; flat and
%      circle, N = 16, unsnapped and snapped.
%   b. continuity correction (p = 2, N = 16, snapped): with the Dubiner
%      dofs perturbed by 1e-5 noise (the basis is O(10) and its gradient
%      O(10/h) on an element, so 1e-3 already perturbs grad u by O(1) and
%      folds the pieces) the raw curve jumps (max |jump| > 1e-5);
%      continuous = true makes the one-sided values agree at every
%      interior break (and across the wrap of the circle) to 1e-12 and
%      changes D by at most the largest jump; the tangent of the
%      DeformedCurve built from (D, D_prime) matches a central difference
%      of its points at the piece midpoints to 1e-6.
%   c. convergence without a solve: u = [0.1 sin(pi x) cos(y);
%      0.2 (y - 0.45 - d) cos(x - y)] (flat, zero normal component on
%      Sigma), and u = xi - Q(phi) xi, xi = (x - 0.5 - d, y - 0.5), Q the
%      rotation by phi = 0.3 (1 + 0.5 x y) (circle; it maps the translated
%      circle onto Sigma_0, so the exact set is the translated circle while
%      the material parameter slides), L2-projected onto degree p = k + 1:
%      the distance to the exact set and the error against the exact
%      material points converge at order >= p + 0.6 on the last pair
%      (expected p + 1), on the flat curve the measure ratio at the same
%      order (or within 1e-12 of 1; printed only for the circle, whose two
%      levels give a noisy rate), and the continuous reconstruction is
%      monotone. The raw (discontinuous) curve folds back at the breaks
%      where the tangential part of the jump points backwards, by a
%      fraction of the jump: its monotonicity is printed, slightly
%      negative. Flat N = 8, 16, 32; circle N = 16, 32. Snapped.
%   d. the metrics of the exact curve: distance and point error below
%      1e-12, measure ratio 1 to 1e-8 (the Clenshaw-Curtis rule of
%      interface_curve_metrics gives about 1e-14 here), monotone.
clc; clear;

n_fail = 0;
d = 0.037;
setups = {flat_setup(d), circle_setup(d)};

%% a. polynomial reproduction
for g = 1:2
    G = setups{g};
    for k = 1:2
        p = k + 1;
        [u, grad_u] = polynomial_field(G, p);
        for snap = [false, true]
            sm = solid_carve(16, G, snap);
            dofs = L2proj2D(sm, u, p);
            [D, D_prime, info] = eulerian_interface_displacement(G.reference, G.curve, sm, ...
                dofs, p, G.rigid);
            s = test_parameters(G, info.breaks);
            X = G.reference.get_2d_points(s);
            x = exact_points(X, G.rigid, u);
            e_x = max(vecnorm(X + D(s) - x, 2, 1));
            T = G.reference.tangent(s) + 0*s;
            e_d = max(vecnorm(D_prime(s) - (solve_jacobian(grad_u(x), T) - T), 2, 1));
            jump = max(vecnorm(info.jumps, 2, 1));
            slide = max(abs(D(s) - G.rigid), [], "all");
            fprintf("a. %-6s p = %d snap %d: |x_h - x| %.2e  |D' - exact| %.2e  max jump %.2e  " + ...
                "(%d pieces, slide %.3f)\n", G.name, p, snap, e_x, e_d, jump, numel(info.elements), slide);
            n_fail = n_fail + check(e_x <= 1e-10 && e_d <= 1e-9 && jump <= 1e-11, ...
                sprintf("polynomial reproduction, %s p = %d snap %d", G.name, p, snap));
        end
    end
end

%% b. continuity correction
rng(1);
for g = 1:2
    G = setups{g};
    p = 2;
    u = polynomial_field(G, p);
    sm = solid_carve(16, G, true);
    dofs = L2proj2D(sm, u, p);
    dofs = dofs + 1e-5*randn(size(dofs));
    [D, ~, raw] = eulerian_interface_displacement(G.reference, G.curve, sm, dofs, p, G.rigid);
    [Dc, Dc_prime, info] = eulerian_interface_displacement(G.reference, G.curve, sm, dofs, p, ...
        G.rigid, continuous = true);
    max_jump = max(vecnorm(raw.jumps, 2, 1));
    b = info.breaks;
    n_pieces = numel(info.elements);
    if G.closed
        j = 1:n_pieces;
        before = [n_pieces, 1:n_pieces - 1];
    else
        j = 2:n_pieces;
        before = j - 1;
    end
    one_sided = max(vecnorm(info.on_piece(b(j), before) - info.on_piece(b(j), j), 2, 1));
    s = test_parameters(G, b);
    change = max(vecnorm(Dc(s) - D(s), 2, 1));
    curve = DeformedCurve(G.reference, Dc, parametric = true, derivative = Dc_prime);
    if G.closed, ends = [b(2:end), b(1) + 1]; else, ends = b(2:end); end
    mid = (b(1:n_pieces) + ends)/2;
    h = 1e-6;
    fd = (curve.get_2d_points(mid + h) - curve.get_2d_points(mid - h))/(2*h);
    e_t = max(vecnorm(curve.tangent(mid) - fd, 2, 1));
    fprintf("b. %-6s raw max jump %.2e; corrected one-sided mismatch %.2e; max |D_c - D| %.2e; " + ...
        "|tangent - FD| %.2e\n", G.name, max_jump, one_sided, change, e_t);
    n_fail = n_fail + check(max_jump > 1e-5, sprintf("noise makes no jump, %s", G.name));
    n_fail = n_fail + check(one_sided <= 1e-12, sprintf("corrected curve discontinuous, %s", G.name));
    n_fail = n_fail + check(change <= max_jump + 1e-12, sprintf("correction above the jump, %s", G.name));
    n_fail = n_fail + check(e_t <= 1e-6, sprintf("DeformedCurve tangent, %s", G.name));
end

%% c. convergence of the reconstruction
levels = {[8 16 32], [16 32]};
for g = 1:2
    G = setups{g};
    u = smooth_field(G);
    for k = 1:2
        p = k + 1;
        N = levels{g};
        E = zeros(numel(N), 5);
        for i = 1:numel(N)
            sm = solid_carve(N(i), G, true);
            dofs = L2proj2D(sm, u, p);
            [D, ~, info] = eulerian_interface_displacement(G.reference, G.curve, sm, dofs, p, G.rigid);
            exact = @(s) exact_points(G.reference.get_2d_points(s), G.rigid, u);
            m = interface_curve_metrics(@(s) G.reference.get_2d_points(s) + D(s), info.breaks, ...
                exact, G.geometry);
            Dc = eulerian_interface_displacement(G.reference, G.curve, sm, dofs, p, G.rigid, ...
                continuous = true);
            mc = interface_curve_metrics(@(s) G.reference.get_2d_points(s) + Dc(s), info.breaks, ...
                [], G.geometry);
            E(i, :) = [m.distance, m.point_error, abs(m.measure_ratio - 1), m.monotonicity, ...
                mc.monotonicity];
        end
        rate = [nan(1, 3); log2(E(1:end - 1, 1:3)./E(2:end, 1:3))];
        fprintf("c. %s, p = %d (expected order %d)\n", G.name, p, p + 1);
        fprintf("    N   distance  rate   point err  rate   |ratio-1|  rate   monotone (raw, continuous)\n");
        for i = 1:numel(N)
            fprintf("   %2d   %.2e  %5.2f  %.2e  %5.2f  %.2e  %5.2f  %9.2e  %.2e\n", N(i), ...
                E(i, 1), rate(i, 1), E(i, 2), rate(i, 2), E(i, 3), rate(i, 3), E(i, 4), E(i, 5));
        end
        n_fail = n_fail + check(all(rate(end, 1:2) >= p + 0.6), ...
            sprintf("%s p = %d: distance / point error order below p + 0.6", G.name, p));
        if ~G.closed
            n_fail = n_fail + check(rate(end, 3) >= p + 0.6 || E(end, 3) <= 1e-12, ...
                sprintf("%s p = %d: measure order below p + 0.6", G.name, p));
        end
        n_fail = n_fail + check(all(E(:, 5) > 0), sprintf("%s p = %d: not monotone", G.name, p));
    end
end

%% d. the metrics of the exact curve
for g = 1:2
    G = setups{g};
    sm = solid_carve(16, G, true);
    [~, ~, info] = eulerian_interface_displacement(G.reference, G.curve, sm, ...
        L2proj2D(sm, smooth_field(G), 2), 2, G.rigid);
    exact = @(s) G.curve.get_2d_points(s);
    m = interface_curve_metrics(exact, info.breaks, exact, G.geometry);
    fprintf("d. %-6s exact curve: distance %.2e  point error %.2e  |ratio - 1| %.2e  monotone %.2e\n", ...
        G.name, m.distance, m.point_error, abs(m.measure_ratio - 1), m.monotonicity);
    n_fail = n_fail + check(m.distance < 1e-12 && m.point_error < 1e-12 && ...
        abs(m.measure_ratio - 1) <= 1e-8 && m.monotonicity > 0, sprintf("exact metrics, %s", G.name));
end

fprintf("%d check(s) failed.\n", n_fail);
if n_fail > 0
    error("test_eulerian_interface_reconstruction:failed", "%d check(s) failed", n_fail);
end

function G = flat_setup(d)
G.name = "flat";
G.closed = false;
G.curve = FlatCurveHorizontal(0.45 + d);
G.reference = FlatCurveHorizontal(0.45);
G.rigid = [0; d];
G.style = "left";
G.geometry = struct("type", "flat", "height", 0.45 + d);
G.d = d;
end

function G = circle_setup(d)
G.name = "circle";
G.closed = true;
G.curve = CircularCurve([0.5 + d; 0.5], 0.2);
G.reference = CircularCurve([0.5; 0.5], 0.2);
G.rigid = [d; 0];
G.style = "crossed";
G.geometry = struct("type", "circle", "center", [0.5 + d; 0.5], "radius", 0.2);
G.d = d;
end

function sm = solid_carve(N, G, snap)
% the solid side of the background grid, as eulerian_solid_march carves it
grid = StructuredSimplicialGrid(N, N, G.style);
if snap
    [~, sm] = grid.split(G.curve, max_normal_deviation = 90);
    sm = snap_to_curve(sm, G.curve, normal_sign = -1);
else
    [~, sm] = grid.split(G.curve);
end
end

function [u, grad_u] = polynomial_field(G, p)
% a global polynomial of degree p sliding the material by about 0.08;
% grad_u(x) returns the rows u1_x, u2_x, u1_y, u2_y
if ~G.closed
    yc = 0.45 + G.d;
    c3 = double(p == 3);
    u = {@(x, y) 0.08 + 0.2*x.*(y - yc) - 0.1*(x - 0.5).^2 + c3*0.1*(x - 0.5).^3; ...
         @(x, y) 0.3*(y - yc).*(1 - x) + 0.01*x.^2 + c3*0.05*x.^2.*(y - yc)};
    grad_u = @(x) [0.2*(x(2, :) - yc) - 0.2*(x(1, :) - 0.5) + c3*0.3*(x(1, :) - 0.5).^2; ...
        -0.3*(x(2, :) - yc) + 0.02*x(1, :) + c3*0.1*x(1, :).*(x(2, :) - yc); ...
        0.2*x(1, :); ...
        0.3*(1 - x(1, :)) + c3*0.05*x(1, :).^2];
else
    a = 0.5 + G.d;  b = 0.5;
    c3 = double(p == 3);
    u = {@(x, y) -0.4*(y - b) + 0.3*(x - a).*(y - b) + c3*0.5*(x - a).^3; ...
         @(x, y) 0.4*(x - a) + 0.2*((x - a).^2 - (y - b).^2) - c3*0.5*(y - b).^3};
    grad_u = @(x) [0.3*(x(2, :) - b) + c3*1.5*(x(1, :) - a).^2; ...
        0.4 + 0.4*(x(1, :) - a); ...
        -0.4 + 0.3*(x(1, :) - a); ...
        -0.4*(x(2, :) - b) - c3*1.5*(x(2, :) - b).^2];
end
end

function u = smooth_field(G)
d = G.d;
if ~G.closed
    u = {@(x, y) 0.1*sin(pi*x).*cos(y); @(x, y) 0.2*(y - 0.45 - d).*cos(x - y)};
else
    phi = @(x, y) 0.3*(1 + 0.5*x.*y);
    u = {@(x, y) (x - 0.5 - d) - (cos(phi(x, y)).*(x - 0.5 - d) - sin(phi(x, y)).*(y - 0.5)); ...
         @(x, y) (y - 0.5) - (sin(phi(x, y)).*(x - 0.5 - d) + cos(phi(x, y)).*(y - 0.5))};
end
end

function x = exact_points(X, rigid, u)
% the solution of x - rigid - u(x) = X by the fixed-point iteration
% (u is a contraction, |grad u| < 1/2)
x = X + rigid;
for it = 1:200
    x_new = X + rigid + [u{1}(x(1, :), x(2, :)); u{2}(x(1, :), x(2, :))];
    change = max(abs(x_new - x), [], "all");
    x = x_new;
    if change < 1e-16, break; end
end
end

function z = solve_jacobian(grad_u, T)
% (I - grad u) \ T, grad u as the rows u1_x, u2_x, u1_y, u2_y
J = [1 - grad_u(1, :); -grad_u(2, :); -grad_u(3, :); 1 - grad_u(4, :)];
det = J(1, :).*J(4, :) - J(3, :).*J(2, :);
z = [J(4, :).*T(1, :) - J(3, :).*T(2, :); J(1, :).*T(2, :) - J(2, :).*T(1, :)]./det;
end

function s = test_parameters(G, breaks)
% 150 parameters over the range and the interior breaks +- 1e-9 (at most 50)
if G.closed
    s = linspace(0, 1, 151);
    s = s(1:end - 1);
    b = breaks;
else
    s = linspace(G.reference.t_range(1), G.reference.t_range(2), 150);
    b = breaks(2:end - 1);
end
b = b(round(linspace(1, numel(b), min(numel(b), 25))));
s = [s, b - 1e-9, b + 1e-9];
end
