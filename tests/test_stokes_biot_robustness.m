% test_stokes_biot_robustness.m
% Stokes-Biot (HDG + TPM): robustness in the small and large parameters.
%
% Flat interface, k = 2, h = 1/8 and 1/16, one Euler step from rest
% (Delta t = 1), friction_bjs = 1, gap h^2, divergence-free manufactured
% displacement eta = (sin(pi x) cos(pi y), -cos(pi x) sin(pi y)) (so
% lambda div eta stays bounded as lambda grows), other fields as in
% StokesBiotCoupling.get_manufactured_solution. Parameter sets:
%   baseline     kappa = I, s0 = 1, nu = 0.3 (E = 1, alpha = 1)
%   s0 = 1e-4, s0 = 0
%   kappa0 = 1e-2, 1e-4, each with stabilization_darcy = 1 (tau_p =
%       kappa0 tau) and = 1/kappa0 (tau_p = tau)
%   nu = 0.4999
%   (kappa0, s0, lambda, mu) = (1e-4, 1e-4, 1e6, 1)
% Asserted per set: finite solution and residue <= 1e-6 ||rhs|| on both
% meshes; errors of u_F, u_P, p_P at h = 1/16 <= 20 x baseline; error ratio
% e(1/8) / e(1/16) >= 2^k for the same fields. For the two small-kappa
% sets the assertions apply to the default stabilization_darcy = 1
% (tau_p = kappa0 tau); the 1/kappa0 variant is printed for comparison.
% condest of the free block is printed for every run.
%
% Skips (prints SKIPPED, 0 checks) while the Stokes-Biot code is not on
% the path.
clc; clear;

k = 2;
hs = [1/8, 1/16];
fields = ["u_F", "u_P", "p_P"];
common = {'friction_bjs', 1.0, 'displacement', ...
    ["sin(pi*x)*cos(pi*y)"; "-cos(pi*x)*sin(pi*y)"]};
lam = 1e6; mu = 1;          % E, nu of (lambda, mu) = (1e6, 1)
E_lm = mu*(3*lam + 2*mu)/(lam + mu);
nu_lm = lam/(2*(lam + mu));

n_fail = 0;
need = ["BiotProblem", "StokesBiotCoupling.get_manufactured_solution", ...
        "StokesBiotCoupling.solve"];
missing = need(arrayfun(@(f) isempty(which(f)), need));
if ~isempty(missing)
    fprintf("SKIPPED test_stokes_biot_robustness: missing %s\n", strjoin(missing, ", "));
    fprintf("%d check(s) failed.\n", n_fail);
    return
end

% name, kappa0, s0, E, nu, stabilization_darcy variants
sets = struct( ...
    "name",  {"baseline", "s0 = 1e-4", "s0 = 0", "kappa0 = 1e-2", "kappa0 = 1e-4", ...
              "nu = 0.4999", "kappa0 = s0 = 1e-4, lambda = 1e6"}, ...
    "kappa", {1, 1, 1, 1e-2, 1e-4, 1, 1e-4}, ...
    "s0",    {1, 1e-4, 0, 1, 1, 1, 1e-4}, ...
    "E",     {1, 1, 1, 1, 1, 1, E_lm}, ...
    "nu",    {0.3, 0.3, 0.3, 0.3, 0.3, 0.4999, nu_lm}, ...
    "stab",  {1, 1, 1, [1, 1e2], [1, 1e4], 1, 1});

e_base = [];
for S = sets
    n_var = numel(S.stab);
    E = zeros(numel(hs), numel(fields), n_var);
    ok_finite = true(1, n_var);  res = zeros(numel(hs), n_var);
    for v = 1:n_var
        for i = 1:numel(hs)
            out = stokes_biot_case("flat", k, hs(i), common{:}, ...
                diffusion = S.kappa*eye(2), storage = S.s0, ...
                young_modulus = S.E, poisson_ratio = S.nu, ...
                stabilization_darcy = S.stab(v));
            E(i, :, v) = arrayfun(@(f) out.errors.(f), fields);
            ok_finite(v) = ok_finite(v) && solution_is_finite(out);
            res(i, v) = out.residue;
            gs = out.solutions.global_system;
            c = condest(gs.matrix(gs.dof, gs.dof));
            fprintf("%-32s stab %-6g h = 1/%-2d  residue %.1e  condest %.2e  %s\n", ...
                S.name, S.stab(v), round(1/hs(i)), res(i, v), c, ...
                strjoin(compose("%s %.3e", fields.', E(i, :, v).'), "  "));
        end
    end
    if isempty(e_base), e_base = E(end, :, 1); end    % the first set is the baseline
    % assert the default (the first variant); print the others' errors
    b = 1;
    for v = 2:n_var
        fprintf("%-32s stabilization_darcy = %g: e/e_base %s (printed only)\n", S.name, ...
            S.stab(v), strjoin(compose("%.2f", E(end, :, v) ./ e_base), " "));
    end
    e_fine = E(end, :, b);  ratio = E(1, :, b) ./ E(end, :, b);
    fprintf("%-32s e/e_base %s   ratio %s\n", S.name, ...
        strjoin(compose("%.2f", e_fine ./ e_base), " "), strjoin(compose("%.2f", ratio), " "));
    n_fail = n_fail + check(ok_finite(b), sprintf("%s: non-finite solution", S.name));
    n_fail = n_fail + check(all(res(:, b) <= 1e-6), ...
        sprintf("%s: residue %.1e > 1e-6 ||rhs||", S.name, max(res(:, b))));
    for j = 1:numel(fields)
        n_fail = n_fail + check(e_fine(j) <= 20*e_base(j), ...
            sprintf("%s: e(%s) at h = 1/16 is %.1f x baseline (> 20)", ...
                S.name, fields(j), e_fine(j)/e_base(j)));
        n_fail = n_fail + check(ratio(j) >= 2^k, ...
            sprintf("%s: e(%s) ratio %.2f < 2^k = %d", S.name, fields(j), ratio(j), 2^k));
    end
end

fprintf("%d check(s) failed.\n", n_fail);
if n_fail > 0
    error("test_stokes_biot_robustness:failed", "%d checks failed", n_fail);
end

function ok = solution_is_finite(out)
% the coupled vector and the six recovered fields
ok = all(isfinite(out.solutions.global_system.solution));
for i = [1, 2]
    ok = ok && all(isfinite(out.solutions.fluid.Unknowns.field(i).dof), "all");
end
for i = [1, 2, 4, 5]
    ok = ok && all(isfinite(out.solutions.biot.Unknowns.field(i).dof), "all");
end
end
