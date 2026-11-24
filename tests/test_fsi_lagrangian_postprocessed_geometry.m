% test_fsi_lagrangian_postprocessed_geometry.m
% Sigma_h built from the postprocessed solid displacement u* (degree k+1,
% postprocessed_displacement of LinearElasticityProblem, and of
% HyperelasticProblem from its discrete gradient L_h) instead of u_h
% (lagrangian_carved_march, interface_geometry = "postprocessed"). Flat
% Sigma_0, linear solid, k = 1, coupling data on the exact curve
% (lagrangian_march_problems, exact_geometry = true), N = 8 and 16 with
% dt = 0.025, 0.0125 (dt ~ h), T = 0.5, snapped fluid, continuous_interface
% = false (the discontinuous Sigma_h, whose limit item 3 shows).
%   1. interface_geometry = "postprocessed" given explicitly reproduces
%      the default march bit for bit;
%   2. the geometry error |Sigma_h - (X + u_S)| with u* is at most 1/10 of
%      that with u_h on both meshes;
%   3. with u*, the fluid stress and velocity errors converge at order
%      > 2 (k + 1) between N = 8 and 16, while with u_h the fluid stress
%      order is below 1.9 (the geometric limit O(h^k) of a degree-k
%      Sigma_h, its normal and Nanson factor);
%   4. hyperelastic (SVK) solid, N = 8: the march with u* converges
%      (Newton) and its geometry error is at most that with u_h;
%   5. for a linear solution, the hyperelastic postprocessing (from L_h,
%      the S_h projection of the linear reconstructed gradient) equals the
%      linear one (from sigma and gamma) to round-off: grad P_k+1 lies in
%      S_h, so the projection does not change the right-hand side.
clc; clear;
warning("off", "MATLAB:nearlySingularMatrix");   % trace extension on sliver patches

n_fail = 0;
k = 1;  T = 0.5;  A = 0.03;
Ns = [8 16];  dts = [0.025 0.0125];
runs = cell(2, 2);           % (geometry: displacement / postprocessed, level)
geometries = ["displacement", "postprocessed"];
for level = 1:2
    N = Ns(level);  dt = dts(level);
    background = StructuredSimplicialGrid(N, N, "left");
    solid_mesh = apply_tag(StructuredSimplicialGrid(1.5*N, 0.75*N, 0, 1, 0, 0.5), ...
        "Interface", @(P) abs(P(2, :) - 0.5) < 1e-14);
    problems = @(t) lagrangian_march_problems(k, "flat", t, dt, A, exact_geometry = true);
    for g = 1:2
        runs{g, level} = lagrangian_carved_march(k, "flat", background, solid_mesh, A, dt, ...
            round(T/dt), "snap", true, "problems", problems, "interface_geometry", geometries(g), ...
            "continuous_interface", false);
    end
    if level == 1
        %% 1. the default, given explicitly
        default = lagrangian_carved_march(k, "flat", background, solid_mesh, A, dt, round(T/dt), ...
            "snap", true, "problems", problems, "continuous_interface", false);
        a = runs{2, 1};
        same = isequal(a.eF, default.eF) && isequal(a.eS, default.eS) && isequal(a.history, default.history);
        fprintf("interface_geometry = ""postprocessed"" bit-identical to the default: %d\n", same);
        n_fail = n_fail + check(same, "explicit interface_geometry = ""postprocessed"" changes the march");
        %% 4. hyperelastic solid
        hyper = cell(1, 2);
        for g = 1:2
            hyper{g} = lagrangian_carved_march(k, "flat", background, solid_mesh, A, dt, round(T/dt), ...
                "snap", true, "solid_model", "hyperelastic", "interface_geometry", geometries(g), ...
                "continuous_interface", false, ...
                "problems", @(t) lagrangian_march_problems(k, "flat", t, dt, A, ...
                solid_model = "hyperelastic", exact_geometry = true));
        end
        g_h = cellfun(@(r) max(r.history(:, 2)), hyper);
        fprintf("hyperelastic, N = %d: geometry error u_h %.3e, u* %.3e; Newton converged %d\n", ...
            N, g_h, all([hyper{2}.newton.converged]));
        n_fail = n_fail + check(all([hyper{2}.newton.converged]), "hyperelastic march with u*: Newton failed");
        n_fail = n_fail + check(g_h(2) <= g_h(1), "hyperelastic: u* geometry error above u_h's");
        %% 5. hyperelastic vs linear postprocessing of one linear solution
        u = default.ss.Unknowns;
        s = struct("displacement", u.field(2).dof, "stress", u.field(1).dof, ...
            "rotation", {u.field(3).dof}, "displacement_skeleton", u.hybrid(1).dof);
        U_lin = LinearElasticityProblem(k, young_modulus = 10, poisson_ratio = 0.3) ...
            .postprocessed_displacement(solid_mesh, k, s);
        U_hyp = HyperelasticProblem(k, young_modulus = 10, poisson_ratio = 0.3) ...
            .postprocessed_displacement(solid_mesh, k, s);
        d2 = nchoosek(k + 2, 2);  d2p = nchoosek(k + 3, 2);
        U_k = zeros(size(U_lin));  U_k([1:d2, d2p + (1:d2)], :) = s.displacement;
        ratio = max(abs(U_hyp - U_lin), [], "all")/max(abs(U_lin - U_k), [], "all");
        fprintf("|u*(L_h) - u*(sigma, gamma)| / |u* - u_h| = %.3e\n", ratio);
        n_fail = n_fail + check(ratio < 1e-10, "hyperelastic and linear postprocessings disagree");
    end
end

%% 2.-3. geometry error and orders
fprintf("\n%-14s %-3s %-10s %-10s %-10s %-10s %-10s\n", "Sigma_h from", "N", "|sF-s|", "|uF-u|", ...
    "|sS-s|", "|uS-u|", "geometry");
order = zeros(2, 4);
for g = 1:2
    for level = 1:2
        r = runs{g, level};
        fprintf("%-14s %-3d %-10.3e %-10.3e %-10.3e %-10.3e %-10.3e\n", geometries(g), Ns(level), ...
            r.eF(1), r.eF(2), r.eS(1), r.eS(2), max(r.history(:, 2)));
    end
    e = @(level) [runs{g, level}.eF(1:2), runs{g, level}.eS(1:2)];
    order(g, :) = log2(e(1) ./ e(2));
    fprintf("  order (sF uF sS uS): %s\n", sprintf("%.2f ", order(g, :)));
end
geometry = cellfun(@(r) max(r.history(:, 2)), runs);
n_fail = n_fail + check(all(geometry(2, :) <= geometry(1, :)/10), ...
    "u* does not reduce the geometry error tenfold");
n_fail = n_fail + check(all(order(2, 1:2) > 2), "fluid order with u* not above 2");
n_fail = n_fail + check(order(1, 1) < 1.9, "fluid stress order with u_h not below 1.9 (negative control)");

fprintf("%d check(s) failed.\n", n_fail);
if n_fail > 0
    error("test_fsi_lagrangian_postprocessed_geometry:failed", "%d check(s) failed", n_fail);
end
