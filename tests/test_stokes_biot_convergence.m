% test_stokes_biot_convergence.m
% Stokes-Biot (HDG + TPM, InterfaceProblem path): h-convergence of the
% coupled manufactured solution, k = 2, one Euler step from rest (Delta t
% = 1), defaults of StokesBiotCoupling.get_manufactured_solution with
% s0 = 1 and friction_bjs = 1, gap h^2, all-Dirichlet outer boundaries.
%
%   flat   (Sigma: y = 1/2), h = 1/4, 1/8, 1/16: rate over all levels
%          (coarsest to finest) >= k + 1 - 0.3 for u_F, u_P, p_P and
%          >= k + 1 - 0.5 for sigma_F, sigma_P, q_P;
%   circle (R = 1/4, Biot inside), h = 1/10, 1/20, 1/40 (Sigma refined
%          at every level: 16, 31, 63 Biot interface faces): the same with
%          0.4 / 0.6 (the interface normal enters the rows only here, so
%          this block is the one that sees a flipped nu_p);
%   circle s0 = 0: the circle ladder with storage s0 = 0 (alpha = 1, no
%          Biot boundary: the pore-pressure level is free and fixed by
%          pore_pressure_mean = the mean of the exact p over the Biot
%          mesh), the same bounds; on every level of the three ladders
%          the level scalar must be kept exactly when s0 = 0;
%   bdf2   flat, h = 1/8, exact time derivative, Delta t = 0.1: the error
%          of every field after 3 steps is within 2x of the error after
%          the first step (an Euler step from initial_state): the
%          history and its transport do not amplify the error.
%
% Skips (prints SKIPPED, 0 checks) while the Stokes-Biot code is not on
% the path.
clc; clear;

k = 2;
fields = ["u_F", "u_P", "p_P", "sigma_F", "sigma_P", "q_P"];
n_primal = 3;                            % the first three get the tight bound
base = {'storage', 1.0, 'friction_bjs', 1.0};

n_fail = 0;
need = ["BiotProblem", "StokesBiotCoupling.get_manufactured_solution", ...
        "StokesBiotCoupling.solve"];
missing = need(arrayfun(@(f) isempty(which(f)), need));
if ~isempty(missing)
    fprintf("SKIPPED test_stokes_biot_convergence: missing %s\n", strjoin(missing, ", "));
    fprintf("%d check(s) failed.\n", n_fail);
    return
end

%% h-convergence, flat and circle
ladders = struct("geometry", {"flat", "circle", "circle"}, ...
    "name", {"flat", "circle", "circle s0 = 0"}, ...
    "h", {[1/4, 1/8, 1/16], [1/10, 1/20, 1/40], [1/10, 1/20, 1/40]}, ...
    "slack", {[0.3, 0.5], [0.4, 0.6], [0.4, 0.6]}, ...
    "extra", {{}, {}, {'storage', 0.0, 'pore_pressure_mean', "exact"}});
for L = ladders
    n_h = numel(L.h);
    E = zeros(n_h, numel(fields));
    for i = 1:n_h
        out = stokes_biot_case(L.geometry, k, L.h(i), base{:}, L.extra{:});
        E(i, :) = arrayfun(@(f) out.errors.(f), fields);
        fprintf("%-6s h = 1/%-3d  residue %.1e  %s\n", L.name, round(1/L.h(i)), ...
            out.residue, strjoin(compose("%s %.3e", fields.', E(i, :).'), "  "));
        % the pore-pressure level scalar is kept only for s0 = 0 (the last ladder)
        gs = out.solutions.global_system;
        n_fail = n_fail + check(ismember(numel(gs.solution), gs.dof) == ~isempty(L.extra), ...
            sprintf("%s: pore-pressure level scalar kept/eliminated wrongly", L.name));
    end
    rate = log(E(1, :) ./ E(end, :)) / log(L.h(1) / L.h(end));
    for j = 1:numel(fields)
        bound = k + 1 - L.slack(1 + (j > n_primal));
        fprintf("%-6s rate(%-7s) = %.2f (pass >= %.2f)\n", L.name, fields(j), rate(j), bound);
        n_fail = n_fail + check(isfinite(rate(j)) && rate(j) >= bound, ...
            sprintf("%s: rate of %s %.2f < %.2f", L.name, fields(j), rate(j), bound));
    end
end

%% bdf2: 3 steps against the first step
need = ["biot_step_history", "BiotProblem/initial_state"];
has = [~isempty(which("biot_step_history")), ...
       ismember("initial_state", methods("BiotProblem"))];
if ~all(has)
    fprintf("SKIPPED bdf2 block: missing %s\n", strjoin(need(~has), ", "));
else
    marching = {'time_scheme', "bdf2", 'time_step', 0.1, 'rate', "exact"};
    out_1 = stokes_biot_case("flat", k, 1/8, base{:}, marching{:}, n_steps = 1);
    out_3 = stokes_biot_case("flat", k, 1/8, base{:}, marching{:}, n_steps = 3);
    for f = fields
        e1 = out_1.errors.(f);  e3 = out_3.errors.(f);
        fprintf("bdf2   %-7s step 1 %.3e  step 3 %.3e  ratio %.2f (pass <= 2)\n", ...
            f, e1, e3, e3/e1);
        n_fail = n_fail + check(isfinite(e3) && e3 <= 2*e1, ...
            sprintf("bdf2: error of %s at step 3 (%.3e) > 2x step 1 (%.3e)", f, e3, e1));
    end
end

fprintf("%d check(s) failed.\n", n_fail);
if n_fail > 0
    error("test_stokes_biot_convergence:failed", "%d checks failed", n_fail);
end
