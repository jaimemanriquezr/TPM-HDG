% test_interface_stokes_stokes_lup.m
% Stokes-Stokes interface problem in the velocity gradient-velocity-
% pressure (L-u-p) formulation, transmission rows spelled entirely as
% symbolic ExtendedFunction forms. Checks EOC k+1, and that
% InterfaceProblem.stokes_stokes gives the same global solution
% (dof for dof) as these explicit forms, its oracle:
%   T1: extrapolated traction balance,
%       sum_i < nu E(L_i) n_i - E(p_i) n_i - nu tau (u_i - u_hat_i), mu_1 > = 0
%       (tau pairs on the face);
%   T2: transferred velocity continuity,
%       < (u_hat_1 + int L_1 . m) - (u_hat_2 + int L_2 . m), mu_2 > = 0.
% Both sides run the constrained pressure closure; the manufactured
% pressure has zero x-average so the closure is exact on any strip.
clc; clear;

nu = 0.7;
y0 = 0.5;
interface = FlatCurveHorizontal(y0);

syms x y
psi = sin(pi*x) .* sin(pi*y);
u_sym = [-diff(psi, y); diff(psi, x)];
p_sym = cos(pi*x) .* sin(pi*y);
f_sym = -nu*[laplacian(u_sym(1), [x y]); laplacian(u_sym(2), [x y])] ...
      + gradient(p_sym, [x y]);
u_ = Datum(u_sym);
f_han = {matlabFunction(f_sym(1), "Vars", [x y]); ...
         matlabFunction(f_sym(2), "Vars", [x y])};

n_fail = 0;
for k = [1, 2]
    Pk = DubinerBasis(k);
    n_ref = 3;
    [err_u1, err_u2, h_list] = deal(zeros(n_ref, 1));
    for i = 1:n_ref
        h = 1/4 * (1/2)^(i-1);
        h_list(i) = h;
        [upper_mesh, lower_mesh] = fitted_split(interface, h, h^2);

        ms = {upper_mesh, lower_mesh};
        pb = cell(1, 2);
        for j = 1:2
            pb{j} = FEProblem( ...
                Spaces  = [FESpace([Pk.to_tensor(), Pk.to_vector(), Pk, DubinerBasis(0)], Condense=true), FESpace([LegendreBasis(k).to_vector(2), DubinerBasis(0), ConstantBasis()])], ...
                Order = k);
            pb{j}.Mesh = ms{j};
            [lhs, flux, rhs] = StokesMatrices(k, ms{j}, nu, f_han);
            pb{j}.Operator = AffineOperator(lhs + rhs);
            pb{j}.Operator = pb{j}.Operator + (flux);
            pb{j}.BoundaryConditions = WeakEssentialBoundaryConditions(u_, "Dirichlet");
        end

        [L1, u1, p1, ~, uhat1, ~, ~] = get_trial_functions(pb{1});
        [L2, u2, p2, ~, uhat2, ~, ~] = get_trial_functions(pb{2});
        [~, ~, ~, ~, mu1, ~, ~] = get_test_functions(pb{1});
        [~, ~, ~, ~, mu2, ~, ~] = get_test_functions(pb{2});
        n1_ = InterfaceNormal(upper_mesh);
        n2_ = InterfaceNormal(lower_mesh);
        tau_ = Stabilization();
        dI1 = Measure(upper_mesh, "dI");
        dI2 = Measure(lower_mesh, "dI");

        T1 = ( nu*(Extrapolate(L1)*n1_) - Extrapolate(p1)*n1_ ...
                 - (nu*tau_)*(u1 - uhat1) ...
             + nu*(Extrapolate(L2)*n2_) - Extrapolate(p2)*n2_ ...
                 - (nu*tau_)*(u2 - uhat2) ) * mu1 * dI1;
        T2 = ( uhat1 + PathIntegral(L1) ...
             - uhat2 - PathIntegral(L2) ) * mu2 * dI2;

        coupled = InterfaceProblem([pb{1}, pb{2}], interface);
        coupled.transmission = [T1, T2];
        [sol1, sol2, gs] = solve(coupled, k);

        % the constructor against the explicit forms
        [~, ~, gs_c] = solve(InterfaceProblem.stokes_stokes(pb{1}, pb{2}, interface, nu, nu), k);
        diff_c = max(abs(gs_c.solution - gs.solution));
        fprintf("k = %d, h = %.4f: stokes_stokes vs explicit forms, max abs diff = %.3e (bitwise equal: %d)\n", ...
            k, h, diff_c, isequal(gs_c.solution, gs.solution));
        if ~(diff_c <= 1e-12 * max(abs(gs.solution)))
            fprintf("k = %d, h = %.4f constructor MISMATCH\n", k, h); n_fail = n_fail + 1;
        end

        err_u1(i) = compute_l2_error(sol1.Unknowns.field(2), upper_mesh, u_);
        err_u2(i) = compute_l2_error(sol2.Unknowns.field(2), lower_mesh, u_);
        fprintf("k = %d, h = %.4f: err u = (%.3e, %.3e), residue = %.1e\n", ...
            k, h, err_u1(i), err_u2(i), gs.residue);
    end
    err_u = max([err_u1, err_u2], [], 2);
    eoc_u = compute_eoc(err_u, h_list);
    fprintf("k = %d: EOC(u) = %s (expected ~ %d)\n\n", k, num2str(eoc_u.', "%.2f "), k+1);
    if mean(eoc_u) < k + 0.5
        fprintf("k = %d EOC FAILED\n", k); n_fail = n_fail + 1;
    end
end
fprintf("%d check(s) failed.\n", n_fail);
if n_fail > 0
    error("test_interface_stokes_stokes_lup:failed", "%d checks failed", n_fail);
end
