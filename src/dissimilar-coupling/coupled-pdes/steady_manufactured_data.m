function data = steady_manufactured_data(curve, box, uu_fluid, p_fluid, u_solid, ...
    mu_F, E, nu, dt, rho_f, rho_s)
%-------------------------------------------------------------------------
%Description: Manufactured data for the steady coupled problem on a closed
%             interface, with the stationary-reaction convention of
%             get_manufactured_solution:
%
%               f_F = rho_f/dt   u_F - div sigma_F
%               f_S = rho_s/dt^2 u_S - div sigma_S
%
%             matching the mass terms the problems assemble from an at-rest
%             history. Construct the solid problem with time = Inf; the
%             default time = 0 gives 2*rho/dt^2 instead of rho/dt^2.
%
%             The traction jump T_1 uses the analytic normal
%             grad(phi)/|grad(phi)| on a level-set curve, and a bracketed
%             closest-point normal on a closed parametrized curve (the
%             unconstrained query can land on the opposite surface where
%             two surfaces pass within a mesh size of each other).
%
%Usage:       data = steady_manufactured_data(kidney, [-1 1 -1 1], ...
%                 "sin(pi*x)*sin(pi*y)", "10*(x^2 + y^2)", ...
%                 ["sin(pi*x)*cos(pi*y)"; "cos(pi*x)*sin(pi*y)"], ...
%                 0.4, 1.0, 0.2, 1e-2, 1.0, 1.0)
%-------------------------------------------------------------------------
arguments
    curve AbstractCurve
    box (1,4) double
    uu_fluid, p_fluid, u_solid
    mu_F, E, nu, dt, rho_f, rho_s
end
assert(isa(curve, "LevelSetCurve") || is_periodic(curve), ...
    "steady_manufactured_data:unsupportedCurve", ...
    "Needs a level set or a closed parametrized curve.");

syms x y
strain = @(v) .5*(jacobian(v, [x y]) + jacobian(v, [x y]).');
shear  = @(v) -.5*(-diff(v(1), y) + diff(v(2), x));
Div    = @(v) [divergence(v(1, :), [x y]); divergence(v(2, :), [x y])];
mu_S = E / (2*(1 + nu));
lambda_S = E*nu / ((1 + nu)*(1 - 2*nu));
zeta = @(v) 2*mu_S*v + lambda_S*trace(v)*eye(2);

u_F = [-diff(str2sym(uu_fluid), y); diff(str2sym(uu_fluid), x)];   % curl of the
p0  = str2sym(p_fluid);                                            % streamfunction
u_S = str2sym(u_solid);
sigma_F = 2*mu_F*strain(u_F) - p0*eye(2);
sigma_S = zeta(strain(u_S));

to_h = @(s) matlabFunction(s, Vars = [x y]);
wrap = @(f) @(xx, yy) f(xx, yy) + 0*xx;
cellify = @(S) arrayfun(@(s) wrap(to_h(s)), S, UniformOutput = false);

data.u_F = cellify(u_F);
data.f_F = cellify(rho_f/dt*u_F - Div(sigma_F));
data.sigma_F = cellify(sigma_F);  data.gamma_F = wrap(to_h(shear(u_F)));
data.u_S = cellify(u_S);
data.f_S = cellify(rho_s/dt^2*u_S - Div(sigma_S));
data.sigma_S = cellify(sigma_S);  data.gamma_S = wrap(to_h(shear(u_S)));

if isa(curve, "LevelSetCurve")
    % Analytic outward normal, straight from the level set.
    phi = curve.evaluate([x; y]);
    grad_phi = [diff(phi, x); diff(phi, y)];
    n_S = grad_phi / sqrt(grad_phi(1)^2 + grad_phi(2)^2);
    data.T_1 = cellify((sigma_S - sigma_F) * n_S);
else
    % Closed parametrized curve: normal from a bracketed closest-point
    % search (see header).
    sig_diff = cellify(sigma_S - sigma_F);
    t_ref = linspace(curve.t_range(1), curve.t_range(2), 4001);
    P_ref = curve.get_2d_points(t_ref);
    data.T_1 = {traction_handle(curve, sig_diff, 1, t_ref, P_ref);
        traction_handle(curve, sig_diff, 2, t_ref, P_ref)};
end
data.T_2 = cellify(u_F - u_S/dt);

tt = linspace(curve.t_range(1), curve.t_range(2), 4001);
P = curve.get_2d_points(tt);
data.measure_F = (box(2) - box(1))*(box(4) - box(3)) - polyarea(P(1, :), P(2, :));
end

%=========================================================================%
function f = traction_handle(curve, sig_diff, comp, t_ref, P_ref)
f = @(xx, yy) traction_component(curve, sig_diff, comp, t_ref, P_ref, xx, yy);
end

function v = traction_component(curve, sig_diff, comp, t_ref, P_ref, xx, yy)
sz = size(xx);
X = [xx(:).'; yy(:).'];
n = zeros(2, size(X, 2));
opts = optimset('TolX', 1e-12);
dt_ref = t_ref(2) - t_ref(1);
for j = 1:size(X, 2)
    [~, i] = min((P_ref(1, :) - X(1, j)).^2 + (P_ref(2, :) - X(2, j)).^2);
    distance = @(t) sum((X(:, j) - curve.get_2d_points(t)).^2);
    t = fminbnd(distance, t_ref(i) - dt_ref, t_ref(i) + dt_ref, opts);
    n(:, j) = curve.unit_normal(t);
end
v = sig_diff{comp, 1}(X(1, :), X(2, :)) .* n(1, :) ...
    + sig_diff{comp, 2}(X(1, :), X(2, :)) .* n(2, :);
v = reshape(v, sz);
end
