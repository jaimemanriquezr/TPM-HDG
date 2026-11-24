% test_transfer_bc_incompressible_elasticity.m
% TPM transfer of the Dirichlet velocity for the Stokes fluid in stress
% form (IncompressibleElasticityProblem, PEERS; transfer_mapping ->
% apply_transfer_conditions, path gradient grad u = A sigma + gamma with
% A sigma = dev(sigma)/(2 mu)). The domain is the unit square; the mesh
% covers [0,1] x [0,1-delta] and the datum on y = 1 is transferred along
% vertical paths of length delta. mu = 1, tau = 1, pressure with zero mean
% (the manufactured pressures have zero mean in x, hence on the strip).
%
% Stage 1 (consistency): fitted square (delta = 0), identity mapping ->
%   the transfer reduces to the projected Dirichlet datum; dof for dof
%   equal to the plain Dirichlet solve, k = 0..3, polynomial data.
% Stage 2 (patch test): u, sigma, gamma polynomial of degree <= k (stream
%   function of degree k + 1, pressure of degree <= k), delta = h^2 and
%   delta = h, n = 4, 8, 16, k = 0..3: the exact solution is in the
%   discrete space and the Taylor transfer is exact, so the relative L2
%   errors (normalised by the largest field norm) must be <= 1e-10.
% Stage 3 (rates): U = sin(pi x) sin(pi y), p = cos(pi x) sin(pi y),
%   n = 4..32 (k = 0, 1 also 64), k = 0..3. delta = h^2: the EOC of the
%   last pair must be >= k + 1 - 0.3 for u and sigma. delta = h: reported
%   only (not expected to be optimal).
%   gamma is printed, not checked: at k = 0 PEERSBasis(0) has no bubble
%   modes and gamma_h does not converge, already on the fitted square
%   without transfer (an element property, not a transfer one). With
%   delta = h and k = 3 the sigma and gamma errors grow from n = 16 on.
clc; clear;
n_fail = 0;
mu = 1.0;
tol_patch = 1e-10;

U_patch = ["x + 2*y", "x^2 + x*y + y^2", "x^3 + x^2*y + y^3", ...
           "x^4 + x^2*y^2 + x*y^3 + y^4"];
P_patch = ["0", "x - 1/2", "(x - 1/2)*(1 + y)", "(x - 1/2)*y^2 + (x^2 - 1/3)*y"];

%% Stage 1: fitted square, identity mapping == plain Dirichlet
mesh = StructuredSimplicialGrid(8, 8);
for k = 0:3
    d = IncompressibleElasticityProblem.get_manufactured_solution(U_patch(k+1), P_patch(k+1), mu, 0, 0);
    pb = IncompressibleElasticityProblem(k, viscosity = mu, source_term = d.f, datum_dirichlet = d.u);
    s0 = solve(pb, mesh, Order = k);
    pb.Data.transfer_mapping = @(X) X;
    s1 = solve(pb, mesh, Order = k);
    e = 0;
    for f = 1:3
        e = max(e, max(abs(s0.Unknowns.field(f).dof - s1.Unknowns.field(f).dof), [], 'all'));
    end
    fprintf("Stage 1 k = %d: identity transfer vs Dirichlet, max dof diff = %.2e\n", k, e);
    if ~(e < 1e-9)
        fprintf("Stage 1 k = %d FAILED\n", k); n_fail = n_fail + 1;
    end
end

%% Stage 2: patch test on the strip, both gaps
for k = 0:3
    d = IncompressibleElasticityProblem.get_manufactured_solution(U_patch(k+1), P_patch(k+1), mu, 0, 0);
    for gap = ["h2", "h"]
        for n = [4 8 16]
            [sol, mesh] = solve_strip(k, n, gap, mu, d);
            [e, r] = field_errors(sol, mesh, d);
            fprintf("Stage 2 k = %d gap %-2s n = %2d: rel err sigma %.1e  u %.1e  gamma %.1e\n", ...
                k, gap, n, r);
            if ~(max(r) <= tol_patch)
                fprintf("Stage 2 k = %d gap %s n = %d FAILED (abs %.1e %.1e %.1e)\n", k, gap, n, e);
                n_fail = n_fail + 1;
            end
        end
    end
end

%% Stage 3: smooth solution, rates
d = IncompressibleElasticityProblem.get_manufactured_solution( ...
    "sin(pi*x)*sin(pi*y)", "cos(pi*x)*sin(pi*y)", mu, 0, 0);
for k = 0:3
    ns = [4 8 16 32];
    if k <= 1, ns = [ns 64]; end %#ok<AGROW>
    for gap = ["h2", "h"]
        err = zeros(numel(ns), 3);
        for i = 1:numel(ns)
            [sol, mesh] = solve_strip(k, ns(i), gap, mu, d);
            err(i, :) = field_errors(sol, mesh, d);
        end
        eoc = [nan(1, 3); compute_eoc(err, 1./ns(:))];
        fprintf("\n=== k = %d, gap %s ===\n", k, gap);
        fprintf("   n   delta      e(sigma)  eoc    e(u)      eoc    e(gamma)  eoc\n");
        for i = 1:numel(ns)
            fprintf("%4d  %.2e  %.3e  %5.2f  %.3e  %5.2f  %.3e  %5.2f\n", ns(i), ...
                gap_size(gap, 1/ns(i)), err(i, 1), eoc(i, 1), err(i, 2), eoc(i, 2), err(i, 3), eoc(i, 3));
        end
        if gap == "h2" && min(eoc(end, 1:2)) < k + 1 - 0.3
            fprintf("Stage 3 k = %d gap h2 rate FAILED\n", k); n_fail = n_fail + 1;
        end
    end
end

fprintf("%d check(s) failed.\n", n_fail);

%% ------------------------------------------------------------------ %%
function delta = gap_size(gap, h)
if gap == "h2", delta = h^2; else, delta = h; end
end

function [sol, mesh] = solve_strip(k, n, gap, mu, d)
% mesh of [0,1] x [0,1-delta]; points on y = 1 - delta go vertically to y = 1
delta = gap_size(gap, 1/n);
y_top = 1 - delta;
mesh = StructuredSimplicialGrid(n, n, 0, 1, 0, y_top);
mapping = @(X) [X(1, :); X(2, :) + (1 - X(2, :)) .* (X(2, :) > y_top - 1e-12)];
pb = IncompressibleElasticityProblem(k, viscosity = mu, source_term = d.f, ...
    datum_dirichlet = d.u, domain_measure = y_top, transfer_mapping = mapping);
sol = solve(pb, mesh, Order = k);
end

function [e, r] = field_errors(sol, mesh, d)
% absolute L2 errors of (sigma, u, gamma) and the same relative to the
% largest discrete field norm (sigma or gamma may vanish identically)
e = sol.compute_l2_error({d.sigma, d.u, d.gamma});
e = e(1:3).';
nrm = arrayfun(@(f) compute_l2_error(sol.Unknowns.field(f), mesh), 1:3);
r = e / max(nrm);
end
