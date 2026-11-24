% test_biot_dynamic_ladder.m
% Verification ladder, rungs L1 and L2, for the fully dynamic Biot element
% (BiotProblem with density = rho_P > 0), k = 0, 1, 2, 3, with a patch
% test on every rung.
%
%   L1: unit square, every face "Dirichlet", strong Dirichlet data.
%   L2: the strip [0, 1] x [0, 1 - delta] (delta = h^2 and delta = h, h =
%       1/N), every face "Dirichlet", the data transferred along vertical
%       paths from y = 1 to the top face y = 1 - delta (BiotProblem's
%       transfer_mapping branch, biot_path_gradient); the paths of the
%       other faces have length zero (weak Dirichlet).
%
% Every run starts from the exact state at t0 = 0.5 (initial_state) and
% marches n_steps steps of dt = 0.25, with data linear in t: the backward
% differences of the inertia (v = (eta - eta*)/cdt, a = (v - v*)/cdt) and
% of the mass equation are then exact, so the errors are spatial only. The
% history (eta*, v*, sigma*, p*) is volume data of the strip mesh at t_n,
% the transferred datum is boundary data on y = 1 at t_{n+1}.
%
%   (a) patch test: eta, p polynomials of degree k in space, linear in t;
%       sigma, gamma, q, v then lie in the discrete spaces and every
%       field is reproduced to round-off: max over fields and steps of
%       ||e||/max(||exact||, 1) <= 1e-10, N = 2, 4, schemes euler and
%       bdf2 (its Euler start included), 3 steps; L1 and L2 with both gaps.
%   (b) convergence: the smooth solution of test_biot_dynamic (c), euler;
%       errors of sigma, eta, q, p, v after the first step and after the
%       last step; rate coarsest-finest >= k + 1 - 0.3 (eta, p, v) and
%       >= k + 1 - 0.5 (sigma, q) on L1 and on L2 with delta = h^2; the
%       gap delta = h is reported, not asserted.
% E = 1, nu = 0.3, alpha = 1, s0 = 1, kappa = I, rho_P = 1.
clc; clear;
n_fail = 0;

E = 1; nu = 0.3; alpha = 1; s0 = 1; kappa = eye(2); rho = 1;
params = {"young_modulus", E, "poisson_ratio", nu, "biot_alpha", alpha, ...
    "storage", s0, "diffusion", kappa};
t0 = 0.5; dt = 0.25;
geometries = ["L1", "L2 gap h^2", "L2 gap h"];
names = ["sigma", "eta", "q", "p", "v"];

%% (a) patch test
n_steps = 3;
patch_tol = 1e-10;
for k = 0:3
    patch = BiotProblem.get_manufactured_solution( ...
        [sprintf("(1 + t)*(1 + x + 2*y)^%d", k); sprintf("(2 - t)*(x - y + 1/2)^%d + t", k)], ...
        sprintf("(1 + 2*t)*(3/10 + x + y/2)^%d", k), params{:}, density = rho);
    for geometry = geometries
        for scheme = ["euler", "bdf2"]
            worst = 0;
            for N = [2 4]
                [mesh, mapping] = strip_mesh(geometry, N);
                rel = march(patch, mesh, k, params, rho, scheme, t0, dt, n_steps, mapping, true);
                worst = max(worst, max(rel, [], "all"));
            end
            fprintf("(a) patch k = %d, %-11s %-5s: max relative error %.3e\n", ...
                k, geometry, scheme, worst);
            n_fail = n_fail + check(worst <= patch_tol, ...
                sprintf("(a) patch k = %d, %s, %s", k, geometry, scheme));
        end
    end
end

%% (b) convergence
n_steps = 4;
smooth = BiotProblem.get_manufactured_solution( ...
    ["(1 + t)*sin(pi*x)*cos(pi*y)"; "(2 - t)*x*exp(y)"], ...
    "(1 + 2*t)*cos(pi*x)*sin(pi*y) + t*x*y", params{:}, density = rho);
for k = 0:3
    if k == 0
        levels = [8 16 32 64];
    else
        levels = [4 8 16 32];
    end
    for geometry = geometries
        [err_first, err_last] = deal(zeros(numel(levels), 5));
        for i = 1:numel(levels)
            [mesh, mapping] = strip_mesh(geometry, levels(i));
            [err_last(i, :), err_first(i, :)] = march(smooth, mesh, k, params, rho, ...
                "euler", t0, dt, n_steps, mapping, false);
        end
        rate_first = print_table(sprintf("(b) k = %d, %s, one step", k, geometry), ...
            levels, err_first, names);
        rate_last = print_table(sprintf("(b) k = %d, %s, %d steps", k, geometry, n_steps), ...
            levels, err_last, names);
        if geometry ~= "L2 gap h"
            bound = k + 1 - [0.5 0.3 0.5 0.3 0.3];
            ok = all(rate_first >= bound) && all(rate_last >= bound);
            n_fail = n_fail + check(ok, sprintf("(b) rates k = %d, %s", k, geometry));
        end
    end
end

fprintf("\n%d check(s) failed.\n", n_fail);

%% local functions
function n = check(ok, label)
if ok
    fprintf("  ok    %s\n", label);
    n = 0;
else
    fprintf("  FAIL  %s\n", label);
    n = 1;
end
end

function [mesh, mapping] = strip_mesh(geometry, N)
% N x N mesh of the unit square (L1) or of the strip [0, 1] x [0, 1 -
% delta] with the vertical transfer from y = 1 - delta to y = 1 (L2)
h = 1/N;
switch geometry
    case "L1"
        mesh = StructuredSimplicialGrid(N, N);
        mapping = function_handle.empty;
        return
    case "L2 gap h^2"
        delta = h^2;
    case "L2 gap h"
        delta = h;
end
mesh = StructuredSimplicialGrid(N, N, 0, 1, 0, 1 - delta);
mapping = @(X) [X(1, :); X(2, :) + delta*(X(2, :) > 1 - delta - 1e-12)];
end

function [err, err_first] = march(mms, mesh, k, params, rho, scheme, t0, dt, n_steps, mapping, relative)
% n_steps steps from the exact state at t0; err (1 x 5) the errors of
% sigma, eta, q, p, v at the last step (relative: the maximum over the
% steps of ||e||/max(||exact||, 1)), err_first those of the first step
options = {"density", rho, "time_step", dt, "time_scheme", scheme};
if ~isempty(mapping)
    options = [options, {"transfer_mapping", mapping}];
end
problem = BiotProblem(k, params{:}, options{:});
d = mms.at(t0);
problem.Data.time = t0;
problem.Data.initial_displacement = d.u_P;
problem.Data.initial_velocity = d.v_P;
problem.Data.initial_pressure = d.p_P;
state = initial_state(problem, mesh, k);
err = zeros(1, 5);
for n = 1:n_steps
    d = mms.at(state.time + dt);
    problem.Data.source_term = d.f_P;
    problem.Data.source_mass = d.g_P;
    problem.Data.datum_dirichlet = [d.u_P; {d.p_P}];
    problem.Data.state = state;
    solution = solve(problem, mesh, k);
    state = get_kinematic_state(problem, solution, state);
    velocity = solution.Unknowns.field(2);
    velocity.dof = state.biot.velocity;
    fields = [solution.Unknowns.field([1 2 4 5]), velocity];
    exact = {d.sigma_P, d.u_P, d.q_P, d.p_P, d.v_P};
    e = zeros(1, 5);
    for j = 1:5
        e(j) = compute_l2_error(fields(j), mesh, Datum(exact{j}));
        if relative
            zero = fields(j);
            zero.dof = 0*zero.dof;
            e(j) = e(j)/max(compute_l2_error(zero, mesh, Datum(exact{j})), 1);
        end
    end
    if n == 1, err_first = e; end
    if relative
        err = max(err, e);
    else
        err = e;
    end
end
end

function rate = print_table(title, levels, err, names)
% errors and EOC per level, and the rate between the coarsest and finest
eoc = [nan(1, numel(names)); log2(err(1:end-1, :)./err(2:end, :))];
fprintf("\n%s\n  %4s  %s\n", title, "N", ...
    strjoin(compose("%10s  %5s", names.', repmat("eoc", numel(names), 1)).', "  "));
for i = 1:numel(levels)
    fprintf("  %4d  %s\n", levels(i), ...
        strjoin(compose("%10.3e  %5.2f", err(i, :).', eoc(i, :).').', "  "));
end
rate = log(err(1, :)./err(end, :))/log(levels(end)/levels(1));
fprintf("  rate coarsest-finest: %s\n", strjoin(names + " " + compose("%.2f", rate), ", "));
end
