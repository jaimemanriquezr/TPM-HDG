% test_fsi_lagrangian_navier_stokes.m
% The FSI march with a Navier-Stokes fluid (lagrangian_carved_march,
% fluid_model = "navier_stokes": beta the divergence-free u* of earlier
% solves, remapped to each new carve): flat clamped Sigma_0, fluid carved
% from a "left" background grid and snapped, k = 1, data from
% lagrangian_march_problems with the convective source rho_f (u_F.grad) u_F.
%   1. a fluid of density 0 convects nothing: the march is identical
%      (isequal) to the Stokes march;
%   2. Picard: with picard_iterations = 20 every step stops at the
%      tolerance, each increment at least 10 times below the previous one,
%      with the linear solid and merged into the coupled Newton loop of the
%      hyperelastic solid (every step converged);
%   3. time order, N = 8, T = 1, dt = 0.1, 0.05, 0.025: the semi-implicit
%      march (one Oseen solve per step, beta = 2 u*^n - u*^{n-1}) and the
%      Picard march self-converge at order 2 in u_F and u_S, and their
%      difference is O(dt^2); the constant extrapolation beta = u*^n drops
%      to order 1 (the check is sensitive to the extrapolation);
%   4. space, semi-implicit, T = 1, dt = 0.0125, N = 8, 16: EOC of
%      (sigma_F, u_F, sigma_S, u_S) above 1.8 (k + 1 = 2).
clc; clear;
warning("off", "MATLAB:nearlySingularMatrix");   % trace extension on sliver patches

n_fail = 0;
k = 1;  A = 0.03;
meshes = @(N) deal(StructuredSimplicialGrid(N, N, "left"), ...
    apply_tag(StructuredSimplicialGrid(1.5*N, 0.75*N, 0, 1, 0, 0.5), ...
    "Interface", @(P) abs(P(2, :) - 0.5) < 1e-14));
[background, solid_mesh] = meshes(8);
march = @(dt, T, varargin) lagrangian_carved_march(k, "flat", background, solid_mesh, A, dt, ...
    round(T/dt), "snap", true, varargin{:});
ns = {"fluid_model", "navier_stokes"};
u_dof = @(r, side) r.(side).Unknowns.field(2).dof;

%% 1. density 0: identical to Stokes
steady = @(t) lagrangian_march_problems(k, "flat", t, 0.1, A, rho_f = 0, convection = true);
stokes = march(0.1, 0.3, "problems", steady);
navier = march(0.1, 0.3, "problems", steady, ns{:}, "picard_iterations", 5);
same = isequal(stokes.fs.Unknowns.field(2).dof, navier.fs.Unknowns.field(2).dof) ...
    && isequal(stokes.ss.Unknowns.field(2).dof, navier.ss.Unknowns.field(2).dof) ...
    && isequal(stokes.history, navier.history);
fprintf("rho_f = 0: Navier-Stokes march identical to Stokes: %d\n", same);
n_fail = n_fail + check(same, "rho_f = 0 does not reproduce the Stokes march");

%% 2. Picard convergence (linear solid, and merged with the hyperelastic Newton loop)
tolerance = 1e-8;
for solid = ["linear", "hyperelastic"]
    r = march(0.1, 1, ns{:}, "picard_iterations", 20, "picard_tolerance", tolerance, ...
        "solid_model", solid);
    contracts = true;
    for n = 1:numel(r.picard)
        h = r.picard(n).increment;
        contracts = contracts && h(end) <= tolerance && all(h(2:end) <= h(1:end-1)/10);
    end
    fprintf("%s solid, Picard updates per step: %s; last step increments %s\n", solid, ...
        sprintf("%d ", [r.picard.updates]), sprintf("%.1e ", r.picard(end).increment));
    n_fail = n_fail + check(contracts, solid + ": Picard did not contract to the tolerance");
    if solid == "hyperelastic"
        fprintf("  Newton solves per step: %s\n", sprintf("%d ", [r.newton.iterations]));
        n_fail = n_fail + check(all([r.newton.converged]), "a merged Newton-Picard step did not converge");
    end
end

%% 3. time order of the semi-implicit and Picard marches
dts = [0.1, 0.05, 0.025];
schemes = {{}, {"picard_iterations", 20, "picard_tolerance", 1e-10}, ...
    {"convective_extrapolation", "constant"}};
names = ["semi-implicit", "Picard", "constant beta"];
runs = cell(3, 3);
for i = 1:3
    for j = 1:3
        runs{i, j} = march(dts(i), 1, ns{:}, schemes{j}{:});
    end
end
order = zeros(3, 2);
for j = 1:3
    D = zeros(2, 2);
    for i = 1:2
        D(i, :) = [max(abs(u_dof(runs{i, j}, "fs") - u_dof(runs{i + 1, j}, "fs")), [], "all"), ...
                   max(abs(u_dof(runs{i, j}, "ss") - u_dof(runs{i + 1, j}, "ss")), [], "all")];
    end
    order(j, :) = log2(D(1, :)./D(2, :));
    fprintf("%-13s: self-convergence order (u_F, u_S) %.2f %.2f\n", names(j), order(j, :));
end
gap = arrayfun(@(i) max(abs(u_dof(runs{i, 1}, "fs") - u_dof(runs{i, 2}, "fs")), [], "all"), 1:3);
gap_order = log2(gap(1:2)./gap(2:3));
fprintf("|semi-implicit - Picard| (u_F): %s, orders %s\n", sprintf("%.2e ", gap), sprintf("%.2f ", gap_order));
n_fail = n_fail + check(all(order(1:2, :) > 1.8, "all"), "semi-implicit or Picard march below order 2");
n_fail = n_fail + check(all(gap_order > 1.8), "semi-implicit - Picard is not O(dt^2)");
n_fail = n_fail + check(all(order(3, :) < 1.3), "constant extrapolation not first order");

%% 4. spatial convergence, semi-implicit
E = zeros(2, 4);
for level = 1:2
    [background, solid_mesh] = meshes(8*level);
    r = lagrangian_carved_march(k, "flat", background, solid_mesh, A, 0.0125, 80, "snap", true, ns{:});
    E(level, :) = [r.eF(1:2), r.eS(1:2)];
end
eoc = log2(E(1, :)./E(2, :));
fprintf("N = 8, 16: (sigma_F, u_F, sigma_S, u_S) %s / %s; EOC %s\n", sprintf("%.2e ", E(1, :)), ...
    sprintf("%.2e ", E(2, :)), sprintf("%.2f ", eoc));
n_fail = n_fail + check(all(eoc > 1.8), "spatial EOC below k + 1");

fprintf("%d check(s) failed.\n", n_fail);
if n_fail > 0
    error("test_fsi_lagrangian_navier_stokes:failed", "%d check(s) failed", n_fail);
end
