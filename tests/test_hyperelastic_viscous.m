%-------------------------------------------------------------------------
%Description: The Kelvin–Voigt solid viscosity of HyperelasticProblem
%             (Data.solid_viscosity = μ_v): σ = P(F) + F S_v,
%             S_v = μ_v(ḞᵀF + FᵀḞ), Ḟ the Newmark rate of the discrete
%             gradient.
%
%             Manufactured displacement u = t⁴ a (sin(πx) y, x cos(πy)),
%             a = 0.05, on the unit square, neo-Hookean, E = 1, ν = 0.3,
%             ρ = 1, μ_v = 0.1 (the viscous stress is of the size of the
%             elastic one at T = 1), Dirichlet data u(t) on the whole
%             boundary, average acceleration (β = 1/4, γ = 1/2), marched to
%             T = 1 from initial_state (u, u̇, ü vanish at t = 0), f(t_{n+1})
%             and u(t_{n+1}) set before each step.
%             Checks:
%               (1) space, τ = 4/h, k = 1, 2, N = 4, 8, 16, Δt = 1/N^((k+1)/2)
%                   (rounded; the time error below the spatial one): at T,
%                   EOC of u_h, û_h ≥ k + 0.8 and of σ_h, γ_h ≥ k − 0.2 against
%                   the exact (σ = P + F S_v, u, γ, u);
%               (2) time, k = 2, N = 4, τ = 1, Δt = 1/8 ... 1/64:
%                   self-convergence EOC of the volume displacement ≥ 1.8;
%               (3) Newton converged at every step of (1) and (2), and is
%                   quadratic: consecutive residuals with r_n < 1e-3 and
%                   r_{n+1} > 100 r_end (above the round-off floor of the
%                   step, which is high relative to the small data of the
%                   first steps) have log r_{n+1} / log r_n ≥ 1.5;
%               (4) the term is active: the scheme with μ_v = 0 on the same
%                   data (k = 2, N = 8) has a u_h error above 10x that of
%                   the viscous scheme.
%-------------------------------------------------------------------------
n_failed = 0;
E = 1; nu = 0.3; rho = 1; mu_v = 0.1; T = 1;
u_str = ["t^4*0.05*sin(pi*x)*y"; "t^4*0.05*x*cos(pi*y)"];
[f, u] = time_data(u_str, E, nu, rho, mu_v);
exact = HyperelasticProblem.get_manufactured_solution(u_str, E, nu, 0, rho, T, "neo_hookean", mu_v);
make = @(k, dt, mu_v) HyperelasticProblem(k, material = "neo_hookean", young_modulus = E, ...
    poisson_ratio = nu, density = rho, time_step = dt, solid_viscosity = mu_v, ...
    source_term = at(f, 0), datum_dirichlet = at(u, 0));
orders = [];

%% (1) space
Ns = [4 8 16];
for k = 1:2
    e = zeros(numel(Ns), 4);
    converged = true;
    for i = 1:numel(Ns)
        n_steps = round(Ns(i)^((k + 1)/2));
        mesh = StructuredSimplicialGrid(Ns(i), Ns(i)).add_stabilization(4*Ns(i));
        [solution, ok, r] = march(make(k, T/n_steps, mu_v), mesh, k, f, u, T/n_steps, n_steps);
        converged = converged && ok;
        orders = [orders, r]; %#ok<AGROW>
        [e_vol, e_hat] = solution.compute_l2_error({exact.sigma, exact.u, exact.gamma}, {exact.u});
        e(i, :) = [e_vol(:).', e_hat];
    end
    eoc = [nan(1, 4); log2(e(1:end-1, :)./e(2:end, :))];
    fprintf("\nk = %d, tau = 4/h, dt = 1/N^((k+1)/2)\n   N |  e(σ_h)    EOC |  e(u_h)    EOC |  e(γ_h)    EOC |  e(û_h)    EOC\n", k);
    for i = 1:numel(Ns)
        fprintf(" %3d | %.3e %5.2f | %.3e %5.2f | %.3e %5.2f | %.3e %5.2f\n", Ns(i), ...
            e(i, 1), eoc(i, 1), e(i, 2), eoc(i, 2), e(i, 3), eoc(i, 3), e(i, 4), eoc(i, 4));
    end
    target = [k - 0.2, k + 0.8, k - 0.2, k + 0.8];
    names = ["sigma", "u", "gamma", "u_hat"];
    for j = 1:4
        n_failed = n_failed + check(eoc(end, j) >= target(j), ...
            sprintf("k = %d: EOC %s >= %.1f", k, names(j), target(j)));
    end
    n_failed = n_failed + check(converged, sprintf("k = %d: Newton converged at every step", k));
end

%% (2) time
k = 2;
mesh = StructuredSimplicialGrid(4, 4);
steps = [8 16 32 64];
U = cell(size(steps));
converged = true;
for i = 1:numel(steps)
    [solution, ok, r] = march(make(k, T/steps(i), mu_v), mesh, k, f, u, T/steps(i), steps(i));
    U{i} = solution.Unknowns.field(2).dof;
    converged = converged && ok;
    orders = [orders, r]; %#ok<AGROW>
end
e = arrayfun(@(i) max(abs(U{i} - U{i + 1}), [], "all"), 1:numel(steps) - 1);
eoc = [nan, log2(e(1:end-1)./e(2:end))];
fprintf("\nk = 2, N = 4, tau = 1\n    dt    | max|u_dt - u_dt/2|  EOC\n");
for i = 1:numel(e)
    fprintf("  1/%-4d |     %.3e     %5.2f\n", steps(i), e(i), eoc(i));
end
n_failed = n_failed + check(eoc(end) >= 1.8, "temporal EOC >= 1.8");
n_failed = n_failed + check(converged, "time: Newton converged at every step");

%% (3) Newton
fprintf("\nobserved Newton orders log r_{n+1}/log r_n: min %.2f over %d pairs\n", min(orders), numel(orders));
n_failed = n_failed + check(~isempty(orders) && all(orders >= 1.5), "quadratic Newton convergence");

%% (4) the viscous term is active
mesh = StructuredSimplicialGrid(8, 8).add_stabilization(32);
e_on = march(make(2, T/23, mu_v), mesh, 2, f, u, T/23, 23).compute_l2_error({exact.sigma, exact.u, exact.gamma});
e_off = march(make(2, T/23, 0), mesh, 2, f, u, T/23, 23).compute_l2_error({exact.sigma, exact.u, exact.gamma});
fprintf("\nk = 2, N = 8: e(u_h) with mu_v = %g: %.3e, with mu_v = 0: %.3e\n", mu_v, e_on(2), e_off(2));
n_failed = n_failed + check(e_off(2) > 10*e_on(2), "the viscous term is active");
fprintf("\n%d check(s) failed.\n", n_failed);

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

function [f, u] = time_data(u_str, E, nu, rho, mu_v)
% space-time handles of f = ρ ü − ∇·(P + F S_v) and u (neo-Hookean)
[lambda, mu] = LinearElasticityProblem.get_lame_parameters(E, nu);
syms x y t real
us = str2sym(u_str);
G = jacobian(us, [x y]);
F = eye(2) + G;
F_dot = diff(G, t);
P = HyperelasticProblem.first_piola_symbolic(G, lambda, mu, "neo_hookean") ...
    + F*(mu_v*(F_dot.'*F + F.'*F_dot));
fs = rho*diff(us, t, 2) - [divergence(P(1, :), [x y]); divergence(P(2, :), [x y])];
f = {matlabFunction(fs(1), Vars = [x y t]); matlabFunction(fs(2), Vars = [x y t])};
u = {matlabFunction(us(1), Vars = [x y t]); matlabFunction(us(2), Vars = [x y t])};
end
