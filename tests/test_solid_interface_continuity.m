% test_solid_interface_continuity.m
% solid_interface_displacement with continuous = true: the solid field on
% Sigma_0, evaluated piecewise from the element polynomials, jumps at the
% solid's interface vertices; the continuous variant corrects each piece
% linearly to the average of the one-sided values. Disc of radius 0.2
% (mesh2d, boundary vertices on the circle), k = 1 and 2, element-wise
% random fields (discontinuous across elements).
%   1. continuous = true: the jumps at the breakpoints vanish (round-off);
%   2. the change is at most the largest jump;
%   3. d_prime is the derivative of d inside the pieces (central
%      differences);
%   4. a smooth field (a quadratic, exact in every element at k = 2) has
%      no jumps and is unchanged.
clc; clear;

n_fail = 0;
center = [0.6; 0.5];  R = 0.2;  h = 0.05;
n = ceil(2*pi*R/h);
theta = 2*pi*(0:n - 1)/n;
mesh = SimplicialGrid(center + R*[cos(theta); sin(theta)], [], element_diameter = h, ...
    refine_boundaries = false);
mesh = apply_tag(mesh, "Interface", @(P) abs(vecnorm(P - center) - R) < 1e-10);
reference = CircularCurve(center, R);
segments = reference.map_mesh_to_interface(mesh);
breaks = unique(segments.t(:)).';
s = linspace(0, 1, 1001);
inside = min(abs(s.' - breaks), [], 2).' > 1e-4;
rng(1);
for k = [1 2]
    d2 = nchoosek(k + 2, 2);
    dofs = 1e-2*randn(2*d2, size(mesh.elements, 2));
    [d0, p0] = solid_interface_displacement(reference, mesh, dofs, k);
    [d1, p1] = solid_interface_displacement(reference, mesh, dofs, k, continuous = true);
    jump = @(d) vecnorm(d(breaks - 1e-10) - d(breaks + 1e-10));
    J0 = max(jump(d0));  J1 = max(jump(d1));
    change = max(vecnorm(d1(s) - d0(s)));
    fd = (d1(s + 1e-6) - d1(s - 1e-6))/2e-6;
    e_derivative = max(vecnorm(p1(s(inside)) - fd(:, inside)));
    fprintf("k = %d: jumps %.2e -> %.2e, change %.2e, |d'' - FD| %.2e\n", k, J0, J1, change, e_derivative);
    n_fail = n_fail + check(J1 < 1e-8, sprintf("k = %d: jumps remain", k));
    n_fail = n_fail + check(change <= J0 + 1e-14, sprintf("k = %d: change above the largest jump", k));
    n_fail = n_fail + check(e_derivative < 1e-6, sprintf("k = %d: d_prime off the derivative", k));
end

%% 4. a quadratic field, exact in every element at k = 2
k = 2;  d2 = nchoosek(k + 2, 2);
u = {@(x, y) 0.1*x.^2 - 0.05*x.*y, @(x, y) 0.02*y.^2 + 0.03*x};
dofs = zeros(2*d2, size(mesh.elements, 2));
lambda = [1 0 0 .5 .5 0 1/3 .2 .6 .2; 0 1 0 .5 0 .5 1/3 .2 .2 .6];   % barycentric samples
for K = 1:size(mesh.elements, 2)
    % least-squares fit through the mesh's own local coordinates (exact:
    % the field is in P_k)
    V = mesh.coordinates(:, mesh.elements(:, K));
    x = V(:, 3) + (V(:, 1) - V(:, 3))*lambda(1, :) + (V(:, 2) - V(:, 3))*lambda(2, :);
    rs = mesh.get_2d_local(reshape(x, 2, 1, []), K*ones(size(x, 2), 1));
    B = dubiner_basis(k, reshape(rs(1, 1, :), [], 1), reshape(rs(2, 1, :), [], 1));
    for c = 1:2
        dofs((c - 1)*d2 + (1:d2), K) = B \ u{c}(x(1, :), x(2, :)).';
    end
end
[d0, ~] = solid_interface_displacement(reference, mesh, dofs, k);
[d1, ~] = solid_interface_displacement(reference, mesh, dofs, k, continuous = true);
change = max(vecnorm(d1(s) - d0(s)));
fprintf("quadratic field, k = 2: change %.2e\n", change);
n_fail = n_fail + check(change < 1e-12, "a continuous field is changed");

fprintf("%d check(s) failed.\n", n_fail);
if n_fail > 0
    error("test_solid_interface_continuity:failed", "%d check(s) failed", n_fail);
end
