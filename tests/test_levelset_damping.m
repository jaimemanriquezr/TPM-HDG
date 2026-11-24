%-------------------------------------------------------------------------
%Description: LevelSetCurve.project with the damping option.
%
%             damping = 1 (the default) must reproduce plain Newton
%             x <- x - phi grad/|grad|^2 exactly. For damping < 1 the
%             iterates must still land on {phi = 0}, and, on the ellipse
%             x^2/a2 + y^2 = 1, approach the end of the field line of
%             grad phi through the start point (y = C x^a2) as the
%             damping goes to zero. Invalid damping values are rejected.
%-------------------------------------------------------------------------
n_fail = 0;
a2 = 6.25;
phi  = @(X) X(1, :).^2/a2 + X(2, :).^2 - 1;
gphi = @(X) [2*X(1, :)/a2; 2*X(2, :)];
curve = LevelSetCurve(phi, gradient = gphi, box = [-3 3 -1.5 1.5]);

% start points off the ellipse (the radially scaled ellipse, as in Fig. B)
u = [20 38 56];
X0 = 1.28*[sqrt(a2)*cosd(u); sind(u)];

% 1. default = plain Newton, bit for bit
ref = X0;
for it = 1:50
    z = phi(ref);
    if all(abs(z) <= 1e-12), break; end
    G = gphi(ref);
    ref = ref - z .* G ./ sum(G.^2, 1);
end
n_fail = n_fail + check(isequal(curve.project(X0), ref), "default project equals plain Newton");
n_fail = n_fail + check(isequal(curve.project(X0, damping = 1), ref), "damping = 1 equals plain Newton");

% 2. damped iterates land on the curve
for lam = [0.5 0.25 0.05]
    P = curve.project(X0, damping = lam);
    n_fail = n_fail + check(all(abs(phi(P)) <= 1e-12), sprintf("damping = %.2f lands on phi = 0", lam));
end

% 3. damped landing -> end of the field line of grad phi, at rate O(damping)
Lfl = zeros(2, numel(u));
for j = 1:numel(u)
    C = X0(2, j)/X0(1, j)^a2;
    xf = fzero(@(x) phi([x; C*x^a2]), [1e-3, X0(1, j)]);
    Lfl(:, j) = [xf; C*xf^a2];
end
e1 = max(vecnorm(curve.project(X0) - Lfl));
e2 = max(vecnorm(curve.project(X0, damping = 0.05) - Lfl));
e3 = max(vecnorm(curve.project(X0, damping = 0.01) - Lfl));
fprintf("distance to the field-line end: damping 1: %.2e, 0.05: %.2e, 0.01: %.2e\n", e1, e2, e3);
n_fail = n_fail + check(e2 < 0.1*e1, "damping = 0.05 is much closer to the field line than plain Newton");
n_fail = n_fail + check(e3 < e2 && e3/e2 < 0.3, "the gap shrinks roughly like the damping");

% 4. invalid damping
for bad = [0 -0.5 1.5]
    try
        curve.project(X0, damping = bad);
        n_fail = n_fail + check(false, sprintf("damping = %g is rejected", bad));
    catch
    end
end

% 5. the projection_damping property reaches the TPM query surface
n_fail = n_fail + check(curve.projection_damping == 1, "projection_damping defaults to 1");
curve = discretize(curve);
damped = curve;
damped.projection_damping = 0.05;
n_fail = n_fail + check(isequal(damped.project(X0), curve.project(X0, damping = 0.05)), ...
    "projection_damping is the default damping of project");
% vertex parameters (map_points_to_interface) move toward the field-line end
t_newton = curve.map_points_to_interface(X0);
t_damped = damped.map_points_to_interface(X0);
t_line = curve.map_points_to_interface(Lfl);
fprintf("parameter distance to the field-line end: damping 1: %.2e, 0.05: %.2e\n", ...
    max(abs(t_newton - t_line)), max(abs(t_damped - t_line)));
n_fail = n_fail + check(max(abs(t_damped - t_line)) < 0.1*max(abs(t_newton - t_line)), ...
    "map_points_to_interface with damping 0.05 lands near the field-line end");
n_fail = n_fail + check(isequal(t_damped, curve.map_points_to_interface(curve.project(X0, damping = 0.05))), ...
    "map_points_to_interface uses the damped projection");
% canonical interface points stay on the curve
damped.projection_damping = 0.25;
P = damped.get_2d_points(linspace(0, damped.t_range(2), 50));
n_fail = n_fail + check(all(abs(phi(P)) <= 1e-12), "get_2d_points with damping 0.25 lands on phi = 0");
% construction option and validation
n_fail = n_fail + check(LevelSetCurve(phi, projection_damping = 0.5).projection_damping == 0.5, ...
    "projection_damping is a constructor option");
for bad = [0 -0.5 1.5]
    try
        damped.projection_damping = bad;
        n_fail = n_fail + check(false, sprintf("projection_damping = %g is rejected", bad));
    catch
    end
end

fprintf("%d check(s) failed.\n", n_fail);
