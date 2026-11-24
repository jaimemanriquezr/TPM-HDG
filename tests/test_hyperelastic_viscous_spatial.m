%-------------------------------------------------------------------------
%Description: The spatial Kelvin–Voigt viscosity of HyperelasticProblem
%             (Data.viscous_model = "spatial"): σ = P(F) + 2μ_v J D F⁻ᵀ,
%             D = sym(Ḟ F⁻¹), Ḟ the Newmark rate of the discrete gradient,
%             through HyperelasticProblem.viscous_stress (the function the
%             Newton step calls).
%             Checks:
%               (a) the tangent dP of viscous_stress matches central
%                   differences of P along Ġ = rate G + c (20 random G
%                   with det F > 0): relative error < 1e-6, for "spatial"
%                   and "material";
%               (b) objectivity, F = Q(θ) U, U constant, θ̇ = ω: with the
%                   exact rate Ḟ = ω Q'(θ) U both models vanish (< 1e-12);
%                   with the Newmark rate (β = 1/4, γ = 1/2) from the exact
%                   history at θ₀, Ḟ = (2/Δt)(F₁ − F₀) − Ḟ₀, the rate of
%                   deformation is exactly a I, a = (2(1 − cos φ) − φ sin φ)/Δt,
%                   φ = ωΔt, so P = 2μ_v J₁ a F₁⁻ᵀ = O(Δt³): match to 1e-6
%                   and ratio under halving Δt in [7.5, 8.5];
%               (c) small strain: with G = s G₀, Ġ = s Ġ₀ the two models
%                   differ by O(s²) (ratio under halving s in [3.6, 4.4]),
%                   and their tangents coincide at G = Ġ = 0;
%               (d) μ_v = 0 gives P = 0 and dP = 0 exactly; unknown
%                   models are rejected;
%               (e) Newmark march, manufactured u = t⁴ a (sin(πx) y,
%                   x cos(πy)), a = 0.05, neo-Hookean, E = 1, ν = 0.3,
%                   ρ = 1, μ_v = 0.1, T = 1, data with the spatial viscous
%                   stress (as in test_hyperelastic_viscous): k = 1,
%                   τ = 4/h, N = 4, 8, 16, Δt = 1/N: EOC of u_h, û_h ≥ 1.8
%                   and of σ_h, γ_h ≥ 0.8; Newton converged at every step
%                   and quadratic (log r_{n+1}/log r_n ≥ 1.5, as in
%                   test_hyperelastic_viscous); the option reaches the
%                   solver (one step with "material" differs).
%-------------------------------------------------------------------------
n_failed = 0;
rng(7);
mu_v = 0.1;
row = @(M) M(:).';                                    % column-major [11 21 12 22]
stress = @(G, G_dot, viscosity, rate, model) ...
    HyperelasticProblem.viscous_stress(G, G_dot, viscosity, rate, model);

%% (a) tangent against central differences
n = 20; rate = 37; h = 1e-6;
G = 0.3*(rand(n, 4) - 0.5);
c = rand(n, 4) - 0.5;
for model = ["spatial", "material"]
    [~, dP] = stress(G, rate*G + c, mu_v, rate, model);
    dP_fd = zeros(n, 4, 4);
    for b = 1:4
        H = zeros(1, 4); H(b) = h;
        P_plus = stress(G + H, rate*(G + H) + c, mu_v, rate, model);
        P_minus = stress(G - H, rate*(G - H) + c, mu_v, rate, model);
        dP_fd(:, :, b) = (P_plus - P_minus)/(2*h);
    end
    err = max(abs(dP - dP_fd), [], "all")/max(abs(dP), [], "all");
    fprintf("(a) %s: tangent vs central differences, relative error %.2e\n", model, err);
    n_failed = n_failed + check(err < 1e-6, sprintf("(a) %s tangent", model));
end

%% (b) objectivity
U = [1.2 0.3; 0.3 0.9];
Q = @(th) [cos(th) -sin(th); sin(th) cos(th)];
dQ = @(th) [-sin(th) -cos(th); cos(th) -sin(th)];
theta0 = 0.4; omega = 2;
F0 = Q(theta0)*U;
F0_dot = omega*dQ(theta0)*U;
for model = ["spatial", "material"]
    P = stress(row(F0 - eye(2)), row(F0_dot), mu_v, 1, model);
    fprintf("(b) %s: |P_v| at a rigid rotation rate %.2e\n", model, norm(P));
    n_failed = n_failed + check(norm(P) < 1e-12, sprintf("(b) %s vanishes on rigid rotation", model));
end
dts = [0.1 0.05 0.025];
P_norm = zeros(size(dts));
for i = 1:numel(dts)
    dt = dts(i); phi = omega*dt; rate = 2/dt;
    F1 = Q(theta0 + phi)*U;
    P = stress(row(F1 - eye(2)), row(rate*(F1 - F0) - F0_dot), mu_v, rate, "spatial");
    a = (2*(1 - cos(phi)) - phi*sin(phi))/dt;
    P_exact = row(2*mu_v*det(F1)*a*inv(F1).');
    err = norm(P - P_exact)/norm(P_exact);
    P_norm(i) = norm(P);
    fprintf("(b) Newmark rate, dt = %.3f: |P_v| = %.3e, vs 2 mu_v J a F^-T: %.2e\n", dt, P_norm(i), err);
    n_failed = n_failed + check(err < 1e-6, sprintf("(b) Newmark rotation, dt = %g", dt));
end
ratio = P_norm(1:end-1)./P_norm(2:end);
fprintf("(b) ratios under halving dt: %s\n", sprintf("%.3f ", ratio));
n_failed = n_failed + check(all(ratio > 7.5 & ratio < 8.5), "(b) Newmark rotation stress O(dt^3)");

%% (c) small strain
G0 = [0.3 -0.2 0.1 0.25];
G0_dot = [0.5 0.4 -0.3 0.2];
rate = 37;
s = [0.04 0.02 0.01];
d = zeros(size(s));
for i = 1:numel(s)
    d(i) = norm(stress(s(i)*G0, s(i)*G0_dot, mu_v, rate, "spatial") ...
        - stress(s(i)*G0, s(i)*G0_dot, mu_v, rate, "material"));
end
ratio = d(1:end-1)./d(2:end);
fprintf("(c) |P_spatial - P_material| at s = %s: %s, ratios %s\n", sprintf("%g ", s), ...
    sprintf("%.3e ", d), sprintf("%.3f ", ratio));
n_failed = n_failed + check(all(ratio > 3.6 & ratio < 4.4), "(c) models agree to O(s^2)");
[~, dP_s] = stress(zeros(1, 4), zeros(1, 4), mu_v, rate, "spatial");
[~, dP_m] = stress(zeros(1, 4), zeros(1, 4), mu_v, rate, "material");
err = max(abs(dP_s - dP_m), [], "all")/max(abs(dP_m), [], "all");
n_failed = n_failed + check(err < 1e-14, "(c) tangents coincide at G = 0");

%% (d) zero viscosity, option validation
[P, dP] = stress(G, 37*G + c, 0, 37, "spatial");
n_failed = n_failed + check(all(P == 0, "all") && all(dP == 0, "all"), "(d) mu_v = 0 gives zero");
n_failed = n_failed + check(throws(@() stress(G, G, mu_v, 1, "deviatoric")), ...
    "(d) viscous_stress rejects an unknown model");
n_failed = n_failed + check(throws(@() HyperelasticProblem(1, viscous_model = "deviatoric")), ...
    "(d) the constructor rejects an unknown model");

%% (e) Newmark march on spatial manufactured data
E = 1; nu = 0.3; rho = 1; T = 1; k = 1;
u_str = ["t^4*0.05*sin(pi*x)*y"; "t^4*0.05*x*cos(pi*y)"];
[f, u, exact] = time_data(u_str, E, nu, rho, mu_v, T);
make = @(dt, model) HyperelasticProblem(k, material = "neo_hookean", young_modulus = E, ...
    poisson_ratio = nu, density = rho, time_step = dt, solid_viscosity = mu_v, ...
    viscous_model = model, source_term = at(f, 0), datum_dirichlet = at(u, 0));
Ns = [4 8 16];
e = zeros(numel(Ns), 4);
converged = true;
orders = [];
for i = 1:numel(Ns)
    mesh = StructuredSimplicialGrid(Ns(i), Ns(i)).add_stabilization(4*Ns(i));
    [solution, ok, r] = march(make(T/Ns(i), "spatial"), mesh, k, f, u, T/Ns(i), Ns(i));
    converged = converged && ok;
    orders = [orders, r]; %#ok<AGROW>
    [e_vol, e_hat] = solution.compute_l2_error({exact.sigma, exact.u, exact.gamma}, {exact.u});
    e(i, :) = [e_vol(:).', e_hat];
end
eoc = [nan(1, 4); log2(e(1:end-1, :)./e(2:end, :))];
fprintf("\n(e) k = 1, tau = 4/h, dt = 1/N, spatial\n   N |  e(σ_h)    EOC |  e(u_h)    EOC |  e(γ_h)    EOC |  e(û_h)    EOC\n");
for i = 1:numel(Ns)
    fprintf(" %3d | %.3e %5.2f | %.3e %5.2f | %.3e %5.2f | %.3e %5.2f\n", Ns(i), ...
        e(i, 1), eoc(i, 1), e(i, 2), eoc(i, 2), e(i, 3), eoc(i, 3), e(i, 4), eoc(i, 4));
end
target = [k - 0.2, k + 0.8, k - 0.2, k + 0.8];
names = ["sigma", "u", "gamma", "u_hat"];
for j = 1:4
    n_failed = n_failed + check(eoc(end, j) >= target(j), sprintf("(e) EOC %s >= %.1f", names(j), target(j)));
end
n_failed = n_failed + check(converged, "(e) Newton converged at every step");
fprintf("(e) observed Newton orders log r_{n+1}/log r_n: min %.2f over %d pairs\n", min(orders), numel(orders));
n_failed = n_failed + check(~isempty(orders) && all(orders >= 1.5), "(e) quadratic Newton convergence");
mesh = StructuredSimplicialGrid(4, 4).add_stabilization(16);
solution = march(make(T/4, "spatial"), mesh, k, f, u, T/4, 4);
U_s = solution.Unknowns.field(2).dof;
solution = march(make(T/4, "material"), mesh, k, f, u, T/4, 4);
U_m = solution.Unknowns.field(2).dof;
gap = max(abs(U_s - U_m), [], "all")/max(abs(U_s), [], "all");
fprintf("(e) N = 4: max|u_spatial - u_material|/max|u_spatial| = %.2e\n", gap);
n_failed = n_failed + check(gap > 1e-8, "(e) viscous_model reaches the solver");
fprintf("\n%d check(s) failed.\n", n_failed);
if n_failed > 0
    error("test_hyperelastic_viscous_spatial:failed", "%d check(s) failed", n_failed);
end

function [solution, converged, orders] = march(pde, mesh, k, f, u, dt, n_steps)
% Newmark march from initial_state; the solution at T, Newton convergence
% at every step and the observed orders of its consecutive residuals
state = initial_state(pde, mesh, k);
converged = true;
orders = [];
for n = 1:n_steps
    t = n*dt;
    pde.Data.source_term = at(f, t);
    pde.Data.datum_dirichlet = at(u, t);
    pde.Data.state = state;
    [solution, ~, newton] = solve(pde, mesh, k);
    converged = converged && newton.converged;
    r = newton.residual;
    pairs = find(r(1:end-1) < 1e-3 & r(2:end) > 100*r(end));
    orders = [orders, log(r(pairs + 1))./log(r(pairs))]; %#ok<AGROW>
    state = get_kinematic_state(pde, solution, state);
end
end

function h = at(g, t)
% the spatial handles of the space-time handles g at time t
h = {@(x, y) g{1}(x, y, t) + 0*x; @(x, y) g{2}(x, y, t) + 0*x};
end

function [f, u, exact] = time_data(u_str, E, nu, rho, mu_v, T)
% space-time handles of f = ρ ü − ∇·(P + 2μ_v J D F⁻ᵀ) and u
% (neo-Hookean), and the exact σ, u, γ at t = T
[lambda, mu] = LinearElasticityProblem.get_lame_parameters(E, nu);
syms x y t real
us = str2sym(u_str);
G = jacobian(us, [x y]);
F = eye(2) + G;
F_inv = inv(F);
L = diff(G, t)*F_inv;
P = HyperelasticProblem.first_piola_symbolic(G, lambda, mu, "neo_hookean") ...
    + mu_v*det(F)*(L + L.')*F_inv.';
fs = rho*diff(us, t, 2) - [divergence(P(1, :), [x y]); divergence(P(2, :), [x y])];
f = {matlabFunction(fs(1), Vars = [x y t]); matlabFunction(fs(2), Vars = [x y t])};
u = {matlabFunction(us(1), Vars = [x y t]); matlabFunction(us(2), Vars = [x y t])};
gamma = (diff(us(1), y) - diff(us(2), x))/2;
exact = struct("sigma", {handles(subs(P, t, T))}, "u", {handles(subs(us, t, T))}, ...
    "gamma", handles(subs(gamma, t, T)));
end

function h = handles(e)
% vectorized @(x, y) handles of a symbolic array (a single handle for a scalar)
syms x y real
h = cell(size(e));
for i = 1:numel(e)
    g = matlabFunction(e(i), Vars = [x, y]);
    h{i} = @(x, y) g(x, y) + 0*x;
end
if isscalar(h), h = h{1}; end
end

function tf = throws(fn)
% true when fn() raises an error
try
    fn();
    tf = false;
catch
    tf = true;
end
end
