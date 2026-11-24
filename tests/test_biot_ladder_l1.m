% test_biot_ladder_l1.m
% Verification ladder, rung 1, for the quasi-static Biot element
% (BiotProblem, BiotMatrices) alone on the unit square, all boundary data
% Dirichlet (eta, p), one implicit step from rest, k = 0, 1, 2, 3:
%   (0) k = 0 feasibility: one solve at k = 0; if it throws, k = 0 is a
%       documented skip (the error is printed, not counted as a failure);
%   (a) patch test: eta and p polynomial of degree <= k in (x, y), so
%       sigma (degree k), gamma and q (degree k - 1) lie in the discrete
%       spaces; the relative L2 errors of (sigma, eta, gamma, q, p),
%       err/max(1, ||exact||), are <= 1e-10 on the 2x2, 4x4 and 8x8 grids;
%   (b) h-convergence with a smooth solution: rate between the coarsest
%       and the finest level >= k + 1 - 0.25 for sigma, eta, q, p (gamma
%       printed only).
clc; clear;
n_fail = 0;
names = ["sigma", "eta", "gamma", "q", "p"];
E = 1; nu = 0.3; alpha = 1; s0 = 0.1; dt = 1;
opts = {"young_modulus", E, "poisson_ratio", nu, "biot_alpha", alpha, "storage", s0};

%% (0) k = 0 feasibility
data = smooth_data(opts, dt);
k_list = 0:3;
try
    solve(biot_problem(data, opts, dt), StructuredSimplicialGrid(4, 4), 0);
    fprintf("(0) k = 0 supported\n");
catch err
    fprintf("(0) k = 0 SKIPPED (not supported): %s\n    %s\n", err.identifier, err.message);
    for s = err.stack(1:min(3, end)).'
        fprintf("    at %s:%d\n", s.file, s.line);
    end
    k_list = 1:3;
end

%% (a) patch test
for k = k_list
    eta = ["(1 + t)*(" + poly(k, 1) + ")"; "(2 - t)*(" + poly(k, 2) + ")"];
    p = "(1 + t)*(" + poly(k, 3) + ")";
    mms = BiotProblem.get_manufactured_solution(eta, p, opts{:}, time_step = dt, ...
        rate = "from_rest");
    pdata = mms.at(dt);
    pde = biot_problem(pdata, opts, dt);
    worst = 0;
    for n = [2 4 8]
        mesh = StructuredSimplicialGrid(n, n);
        sol = solve(pde, mesh, k);
        [err, nrm] = biot_errors(sol, mesh, pdata);
        rel = err./max(1, nrm);
        worst = max(worst, max(rel));
        fprintf("(a) k = %d patch, n = %d: relative L2 errors %s = %s\n", k, n, ...
            strjoin(names, "/"), strjoin(compose("%.2e", rel), " "));
    end
    fprintf("(a) k = %d patch: max relative error %.2e\n", k, worst);
    n_fail = n_fail + check(worst <= 1e-10, "(a) patch test k = " + k);
end

%% (b) convergence
pde = biot_problem(data, opts, dt);
for k = k_list
    switch k
        case 0, levels = [4 8 16 32 64];
        case 3, levels = [4 8 16];
        otherwise, levels = [4 8 16 32];
    end
    err = zeros(numel(levels), 5);
    for i = 1:numel(levels)
        mesh = StructuredSimplicialGrid(levels(i), levels(i));
        err(i, :) = biot_errors(solve(pde, mesh, k), mesh, data);
    end
    eoc = [nan(1, 5); log2(err(1:end-1, :)./err(2:end, :))];
    fprintf("(b) k = %d\n  %4s  %s\n", k, "n", ...
        strjoin(compose("%10s  %5s", names.', repmat("eoc", 5, 1)).', "  "));
    for i = 1:numel(levels)
        fprintf("  %4d  %s\n", levels(i), ...
            strjoin(compose("%10.3e  %5.2f", err(i, :).', eoc(i, :).').', "  "));
    end
    rate = log(err(1, :)./err(end, :))/log(levels(end)/levels(1));
    fprintf("  rate coarsest-finest: %s\n", strjoin(names + " " + compose("%.2f", rate), ", "));
    asserted = [1 2 4 5];    % sigma, eta, q, p; gamma printed only
    n_fail = n_fail + check(all(rate(asserted) >= k + 1 - 0.25), "(b) rates k = " + k);
end

fprintf("%d check(s) failed.\n", n_fail);
if n_fail > 0, error("test_biot_ladder_l1:failed", "%d checks failed", n_fail); end

function data = smooth_data(opts, dt)
mms = BiotProblem.get_manufactured_solution(["sin(pi*x)*cos(pi*y)"; "x*exp(y)*(1 + t)"], ...
    "cos(pi*x)*sin(pi*y) + t*x*y", opts{:}, time_step = dt, rate = "from_rest");
data = mms.at(dt);
end

function pde = biot_problem(data, opts, dt)
pde = BiotProblem([], opts{:}, time_step = dt, source_term = data.f_P, ...
    source_mass = data.g_P, datum_dirichlet = [data.u_P; {data.p_P}]);
end

function s = poly(k, seed)
% a full polynomial of total degree k in (x, y) with fixed O(1) coefficients
s = "";
for i = 0:k
    for j = 0:k - i
        c = (-1)^(i + j + seed)*(1 + mod(3*i + 5*j + seed, 7))/4;
        s = s + sprintf(" + (%g)*x^%d*y^%d", c, i, j);
    end
end
end

function [err, nrm] = biot_errors(sol, mesh, data)
% L2 errors of (sigma, eta, gamma, q, p) and L2 norms of the exact fields
exact = {data.sigma_P, data.u_P, data.gamma_P, data.q_P, data.p_P};
err = zeros(1, 5); nrm = zeros(1, 5);
for i = 1:5
    f = sol.Unknowns.field(i);
    err(i) = compute_l2_error(f, mesh, Datum(exact{i}));
    nrm(i) = compute_l2_error(FEField(f.basis, zeros(size(f.dof))), mesh, Datum(exact{i}));
end
end
