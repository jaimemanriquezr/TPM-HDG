%-------------------------------------------------------------------------
%Description: Temporal order of the Newmark march of
%             HyperelasticGradientProblem.
%
%             u = t⁴ a (sin(πx) y, x cos(πy)), a = 0.05, Saint
%             Venant–Kirchhoff and neo-Hookean (E = 1, ν = 0.3, ρ = 1),
%             Dirichlet data u(t), k = 2, N = 4, default stabilization,
%             β = 1/4, γ = 1/2, T = 1, Δt = 1/8 ... 1/64, from
%             initial_state (the data start from rest). The loop of the
%             class help: f(t_{n+1}), u(t_{n+1}) and the state set before
%             each step, get_kinematic_state after it.
%             Checks:
%               (1) Newton converged at every step;
%               (2) self-convergence of u_h(T): EOC ≥ 1.8;
%               (3) the march is nonlinear: the linear law marched on the
%                   same data (Δt = 1/64) differs from it at T by more
%                   than 10 e(1/32).
%-------------------------------------------------------------------------
n_failed = 0;
E = 1; nu = 0.3; rho = 1; k = 2; N = 4; T = 1;
u_str = ["t^4*0.05*sin(pi*x)*y"; "t^4*0.05*x*cos(pi*y)"];
mesh = StructuredSimplicialGrid(N, N);
steps = [8 16 32 64];
for material = ["saint_venant_kirchhoff", "neo_hookean"]
    [f, u] = time_data(u_str, E, nu, rho, material);
    U = cell(size(steps));
    converged = true;
    for i = 1:numel(steps)
        [U{i}, ok] = march(material, k, E, nu, rho, mesh, f, u, T/steps(i), steps(i));
        converged = converged && ok;
    end
    e = arrayfun(@(i) max(abs(U{i} - U{i + 1}), [], 'all'), 1:numel(steps) - 1);
    r = [nan, log(e(1:end-1)./e(2:end))/log(2)];
    U_lin = march("linear", k, E, nu, rho, mesh, f, u, T/steps(end), steps(end));
    d_lin = max(abs(U_lin - U{end}), [], 'all');
    fprintf("\n%s\n    dt    | max|u_dt - u_dt/2|  EOC\n", material);
    for i = 1:numel(e)
        fprintf("  1/%-4d |     %.3e     %5.2f\n", steps(i), e(i), r(i));
    end
    fprintf("  linear law at dt = 1/%d: max|u_lin - u| = %.3e\n", steps(end), d_lin);
    n_failed = n_failed + check(converged, material + ": Newton converged at every step");
    n_failed = n_failed + check(r(end) >= 1.8, material + ": temporal EOC >= 1.8");
    n_failed = n_failed + check(d_lin > 10*e(end), material + ": the march is nonlinear");
end
fprintf("\n%d check(s) failed.\n", n_failed);

function [u_T, converged] = march(material, k, E, nu, rho, mesh, f, u, dt, n_steps)
pde = HyperelasticGradientProblem(k, material = material, young_modulus = E, poisson_ratio = nu, ...
    density = rho, time_step = dt, source_term = at(f, 0), datum_dirichlet = at(u, 0));
state = initial_state(pde, mesh, k);
converged = true;
for n = 1:n_steps
    pde.Data.source_term = at(f, n*dt);
    pde.Data.datum_dirichlet = at(u, n*dt);
    pde.Data.state = state;
    [solution, ~, newton] = solve(pde, mesh, k);
    converged = converged && newton.converged;
    state = get_kinematic_state(pde, solution, state);
end
u_T = state.solid.displacement;
end

function h = at(g, t)
h = {@(x, y) g{1}(x, y, t) + 0*x; @(x, y) g{2}(x, y, t) + 0*x};
end

function [f, u] = time_data(u_str, E, nu, rho, material)
% space-time handles of f = ρ ü − ∇·P(I + ∇u) and u
[lambda, mu] = LinearElasticityProblem.get_lame_parameters(E, nu);
syms x y t real
us = str2sym(u_str);
P = HyperelasticProblem.first_piola_symbolic(jacobian(us, [x y]), lambda, mu, material);
fs = rho*diff(us, t, 2) - [divergence(P(1, :), [x y]); divergence(P(2, :), [x y])];
f = {matlabFunction(fs(1), Vars = [x y t]); matlabFunction(fs(2), Vars = [x y t])};
u = {matlabFunction(us(1), Vars = [x y t]); matlabFunction(us(2), Vars = [x y t])};
end
