% test_fsi_lagrangian_pair.m
% The per-side interface curve pair {curve_fluid, curve_solid} of the fsi
% coupling (InterfaceProblem.fsi, get_partitions option interface_2):
%   1. the pair {C, C} reproduces the one-curve coupling bit for bit;
%   2. DeformedCurve's finite-difference tangent, normal and Nanson factor
%      match the symbolic ones of lagrangian_pair_problems;
%   3. steady manufactured FSI with the fluid on the deformed flat
%      interface: k + 1 in u_F and u_S (k = 1), and no convergence when
%      the data omit the Nanson factor (negative control);
%   4. a pair rejects the slip laws.
% The full convergence tables are experiments/moving-interface/
% lagrangian_pair_convergence.m.
clc; clear;

n_fail = 0;

%% 1. {C, C} == C, steady circle with a gap, k = 2
k = 2;
circle = CircularCurve([0.5; 0.5], 0.4);
[fm, sm] = fitted_split(circle, 1/8, 1/64);
[fp, sp, src] = FluidStructureCoupling.load_problems(k, circle, 0, 1e-3, ...
    "sin(pi*x)*sin(pi*y)", "10*x^2", ["sin(pi*x)*cos(pi*y)"; "cos(pi*x)*sin(pi*y)"], ...
    0, 0.4, 0, 1.0, 0.2, 0);
[~, ~, ~, g1] = FluidStructureCoupling.solve(fp, sp, fm, sm, circle, src{:});
[~, ~, ~, g2] = FluidStructureCoupling.solve(fp, sp, fm, sm, {circle, circle}, src{:});
n_fail = n_fail + check(isequal(g1.solution, g2.solution), "pair {C, C} differs from C");
fprintf("{C, C} vs C: max |dx| = %g\n", max(abs(g1.solution - g2.solution)));

%% 2. DeformedCurve against the symbolic curve
t = linspace(0, 1, 101);
for geometry = ["flat", "circle"]
    [~, ~, ~, pd, curves] = lagrangian_pair_problems(1, geometry);
    c = curves{1};
    e_tangent = max(abs(c.tangent(t) - pd.tangent(t)), [], "all");
    e_normal = max(abs(-reshape(c.unit_normal(t), 2, []) - pd.normal_fluid(t)), [], "all");
    e_nanson = max(abs(vecnorm(c.tangent(t))./vecnorm(curves{2}.tangent(t)) - pd.nanson(t)));
    fprintf("%s: tangent %.1e, normal %.1e, J %.1e\n", geometry, e_tangent, e_normal, e_nanson);
    n_fail = n_fail + check(max([e_tangent, e_normal, e_nanson]) < 1e-10, ...
        geometry + ": DeformedCurve tangent / normal / J off the symbolic ones");
end

%% 3. flat steady convergence, k = 1, and the no-Nanson control
for nanson = [true, false]
    [fp, sp, src, pd, curves, mesher] = lagrangian_pair_problems(1, "flat", nanson = nanson);
    e = zeros(2, 2);
    for i = 1:2
        [fm, sm] = mesher(1/(8*i));
        [fs, ss] = FluidStructureCoupling.solve(fp, sp, fm, sm, curves, src{:});
        eF = fs.compute_l2_error({pd.sigma_F, pd.u_F, pd.gamma_F});
        eS = ss.compute_l2_error({pd.sigma_S, pd.u_S, pd.gamma_S});
        e(i, :) = [eF(2), eS(2)];
    end
    eoc = log2(e(1, :) ./ e(2, :));
    fprintf("flat, Nanson in the data %d: |u_F - u| %.2e -> %.2e, |u_S - u| %.2e -> %.2e, EOC %.2f %.2f\n", ...
        nanson, e(:, 1), e(:, 2), eoc);
    if nanson
        n_fail = n_fail + check(all(eoc > 1.9), "flat pair: EOC of u below k + 1");
    else
        n_fail = n_fail + check(all(eoc < 0.5), "flat pair: the data without J still converge");
    end
end

%% 4. the slip laws are rejected with a pair
[fm, sm] = mesher(1/8);
n_fail = n_fail + check(throws(@() FluidStructureCoupling.solve(fp, sp, fm, sm, curves, src{:}, ...
    interface_law = "full"), "InterfaceProblem:curvePairSlip"), "a pair accepted interface_law = ""full""");

fprintf("%d check(s) failed.\n", n_fail);
if n_fail > 0
    error("test_fsi_lagrangian_pair:failed", "%d check(s) failed", n_fail);
end
