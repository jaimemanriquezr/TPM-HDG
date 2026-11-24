function out = cylinder_benchmark(k, n_y, options)
% Schäfer–Turek 2D-1 laminar flow around a cylinder (Re = 20) on an unfitted
% mesh (Solano & Vargas 2022, Example 5): channel [0, 2.2] x [0, 0.41],
% cylinder of radius 0.05 centred at (0.2, 0.2), nu = 1e-3, parabolic
% velocity 4 U_m y (H - y)/H^2 with U_m = 0.3 at both x = 0 and x = 2.2 (as
% in the paper), no-slip elsewhere. The background grid of the channel is
% carved against the cylinder (immerse_mesh), so the hole boundary Gamma_h
% is a staircase at distance O(h) from the circle; the no-slip datum is
% transferred to it along closest-point paths. Navier–Stokes by Picard
% (solve_navier_stokes, beta = u*), velocity-gradient form.
%
% Drag and lift are integrated on the true circle C from the fields
% extrapolated along the paths (the same polynomial extension the transfer
% uses), in the benchmark's form
%   F_d = int_C (nu d(u.t)/dn n_2 - p n_1) dC,  F_l = -int_C (nu d(u.t)/dn n_1 + p n_2) dC,
% n the normal pointing into the fluid, t = (n_2, -n_1);
% c_d = 2 F_d / (ubar^2 D), ubar = 0.2, D = 0.1. Reference (Schäfer–Turek):
% c_d = 5.5795, c_l = 0.0106; the paper reports 5.5791 and 0.0110 (k = 3).
arguments
    k (1,1) double = 2
    n_y (1,1) double = 41           % background cells across the channel (h = 0.41/n_y)
    options.tolerance (1,1) double = 1e-8
    options.verbose (1,1) logical = true
    options.convection (1,1) logical = true   % false: the Stokes flow (one solve, no wake)
    % "dirichlet": the parabola at x = 2.2 too (Solano & Vargas); "neumann":
    % the DFG benchmark's do-nothing outlet nu du/dn - p n = 0
    options.outlet (1,1) string {mustBeMember(options.outlet, ["dirichlet", "neumann"])} = "dirichlet"
end
L = 2.2;  H = 0.41;  c = [0.2; 0.2];  R = 0.05;  nu = 1e-3;  U_m = 0.3;
n_x = round(n_y*L/H);
h = H/n_y;

background = StructuredSimplicialGrid(n_x, n_y, 0, L, 0, H);
phi = @(X) R - sqrt((X(1, :, :) - c(1)).^2 + (X(2, :, :) - c(2)).^2);   % < 0 outside the cylinder
mesh = immerse_mesh(background, phi);
% closest point on the circle for the hole faces, identity on the fitted walls
mapping = @(X) closest_point(X, c, R, 3*h);

if options.outlet == "neumann"
    mesh = apply_tag(mesh, "Neumann", @(X) isapprox(X(1, :), L));
    pressure_mean = "pinned";
else
    pressure_mean = "constrained";
end
inflow = {@(x, y) 4*U_m*y.*(H - y)/H^2 .* (abs(x) < 1e-12 | abs(x - L) < 1e-12);
          @(x, y) 0*x};
problem = OseenProblem(k, viscosity = nu, datum_dirichlet = inflow, transfer_mapping = mapping, pressure_mean = pressure_mean);
if options.verbose
    fprintf("cylinder benchmark: k = %d, n_y = %d (h = %.4f, %d elements)\n", k, n_y, h, size(mesh.elements, 2));
end
timer = tic;
if options.convection
    [solution, history, u_star] = solve_navier_stokes(problem, mesh, k, tolerance = options.tolerance, verbose = options.verbose);
else
    solution = solve(problem, mesh, Order = k);  history = [];
    u_star = divergence_free_velocity(solution);
end
[c_d, c_l] = drag_lift(solution, mesh, k, c, R, nu, mapping);
% DFG pressure difference p(0.15, 0.2) - p(0.25, 0.2): both points sit on
% the circle, inside the carved band, so p_h is extrapolated from the
% element whose transfer path reaches each point (as the drag is)
p_diff = point_pressure(solution, mesh, k, [0.15; 0.2], mapping) - point_pressure(solution, mesh, k, [0.25; 0.2], mapping);
out = struct("c_d", c_d, "c_l", c_l, "p_diff", p_diff, "picard", numel(history), "h", h, ...
    "elements", size(mesh.elements, 2), "seconds", toc(timer), "solution", solution, "u_star", {u_star}, ...
    "L_dof", solution.Unknowns.field(1).dof, "u_dof", solution.Unknowns.field(2).dof, "p_dof", solution.Unknowns.field(3).dof, ...
    "outlet", options.outlet, "convection", options.convection);
if options.verbose
    fprintf("c_d = %.4f (ref 5.5795)   c_l = %.4f (ref 0.0106)   p_diff = %.4f (ref 0.1175)   Picard %d it, %.0f s\n", ...
        c_d, c_l, p_diff, numel(history), out.seconds);
end
end

function Xb = closest_point(X, c, R, band)
r = vecnorm(X - c);
near = r < R + band;
Xb = X;
Xb(:, near) = c + R*(X(:, near) - c)./r(near);
end

function [c_d, c_l] = drag_lift(solution, mesh, k, c, R, nu, mapping)
% integrate on the circle: each hole face projects to an arc [th1, th2];
% Gauss points on the arc, fields extrapolated from the adjacent element
d2 = nchoosek(k + 2, 2);
Ld = solution.Unknowns.field(1).dof;  pd = solution.Unknowns.field(3).dof;
L_h = {Ld(1:d2, :), Ld(2*d2+1:3*d2, :); Ld(d2+1:2*d2, :), Ld(3*d2+1:4*d2, :)};
faces = mesh.dirfaces(:).';
P = mesh.coordinates(:, mesh.faces(1:2, faces));           % 2 x 2*n
P = reshape(P, 2, 2, []);
mid = squeeze(mean(P, 2));
hole = vecnorm(mid - c) < 2*R;
faces = faces(hole);  P = P(:, :, hole);
[K, ~] = get_element_by_face(mesh, faces);
if iscell(K), K = cellfun(@(C) C(1), K); end
[s, w] = get_1d_reference_quadrature(k + 1);
F_d = 0;  F_l = 0;
for b = 1:numel(faces)
    e1 = mapping(P(:, 1, b)) - c;  e2 = mapping(P(:, 2, b)) - c;
    th1 = atan2(e1(2), e1(1));  th2 = atan2(e2(2), e2(1));
    dth = mod(th2 - th1 + pi, 2*pi) - pi;                   % shortest arc
    th = th1 + dth*s(:).';                                  % 1 x nq
    xb = c + R*[cos(th); sin(th)];
    n = [cos(th); sin(th)];  t = [n(2, :); -n(1, :)];       % into the fluid; tangent
    rs = mesh.get_2d_local(reshape(xb, 2, [], 1), K(b));
    Pk = dubiner_basis(k, rs(1, :).', rs(2, :).');          % nq x d2 (extrapolated)
    Lq = cellfun(@(d) (Pk*d(:, K(b))).', L_h, 'UniformOutput', false);
    pq = (Pk*pd(:, K(b))).';
    % d(u.t)/dn = t^T L n
    dut_dn = t(1, :).*(Lq{1,1}.*n(1, :) + Lq{1,2}.*n(2, :)) + t(2, :).*(Lq{2,1}.*n(1, :) + Lq{2,2}.*n(2, :));
    dC = R*abs(dth)*w(:).';
    F_d = F_d + sum((nu*dut_dn.*n(2, :) - pq.*n(1, :)).*dC);
    F_l = F_l - sum((nu*dut_dn.*n(1, :) + pq.*n(2, :)).*dC);
end
ubar = 0.2;  D = 2*R;
c_d = 2*F_d/(ubar^2*D);
c_l = 2*F_l/(ubar^2*D);
end

function p = point_pressure(solution, mesh, k, x, mapping)
% p_h at a point on the circle, extrapolated from the element adjacent to
% the hole face whose closest-point image is nearest to x
d2 = nchoosek(k + 2, 2);
pd = solution.Unknowns.field(3).dof;
faces = mesh.dirfaces(:).';
P = reshape(mesh.coordinates(:, mesh.faces(1:2, faces)), 2, 2, []);
mid = squeeze(mean(P, 2));
xb = mapping(mid);
[~, b] = min(vecnorm(xb - x));
[K, ~] = get_element_by_face(mesh, faces(b));
if iscell(K), K = K{1}(1); end
rs = mesh.get_2d_local(reshape(x, 2, 1, 1), K(1));
p = dubiner_basis(k, rs(1), rs(2)) * pd(:, K(1));
end
