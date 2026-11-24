% test_interface_darcy_darcy.m
% Darcy-Darcy TPM interface problem on InterfaceProblem: flat interface, fitted meshes retreating h^2 from Sigma,
% one smooth manufactured p over both subdomains (so the transmission data
% vanish), non-symmetric kappa to catch transpose bugs. Expected EOC k+1.
clc; clear;

kappa = [2, 1; 0, 3];
y0 = 0.5;
interface = FlatCurveHorizontal(y0);

syms x y
p_sym = sin(pi*x) .* sin(pi*y);
q_sym = -kappa * gradient(p_sym, [x y]);
f_sym = divergence(q_sym, [x y]);
p_ = Datum(p_sym);
q_ = Datum(q_sym);
f_ = Datum(f_sym);

dx = Measure(preset = "dx");

n_fail = 0;
for k = [1, 2]
    Pk = DubinerBasis(k); VPk = Pk.to_vector(); Mk = LegendreBasis(k);
    n_ref = 3;
    [err_p1, err_p2, err_q1, err_q2, h_list] = deal(zeros(n_ref, 1));
    for i = 1:n_ref
        h = 1/4 * (1/2)^(i-1);
        h_list(i) = h;
        [upper_mesh, lower_mesh] = fitted_split(interface, h, h^2);

        make = @(mesh) configure(mesh, k, kappa, f_, p_, VPk, Pk, Mk, dx);
        pb1 = make(upper_mesh);
        pb2 = make(lower_mesh);

        coupled = InterfaceProblem.darcy_darcy(pb1, pb2, interface, kappa, kappa);
        [sol1, sol2, gs] = solve(coupled, k);

        %% symbolic transmission rows: same coupling spelled as forms,
        %% checked dof-for-dof against darcy_darcy_rows
        [q1, p1, phat1] = get_trial_functions(pb1);
        [q2, p2, phat2] = get_trial_functions(pb2);
        [~, ~, mu1] = get_test_functions(pb1);
        [~, ~, mu2] = get_test_functions(pb2);
        n1_ = InterfaceNormal(upper_mesh);
        n2_ = InterfaceNormal(lower_mesh);
        tau_ = Stabilization();
        dI1 = Measure(upper_mesh, "dI");
        dI2 = Measure(lower_mesh, "dI");
        T1 = ( dot(Extrapolate(q1), n1_) + tau_*(p1 - phat1) ...
             + dot(Extrapolate(q2), n2_) + tau_*(p2 - phat2) ...
             ) * mu1 * dI1;
        T2 = ( (phat1 - PathIntegral(kappa \ q1)) ...
             - (phat2 - PathIntegral(kappa \ q2)) ) * mu2 * dI2;
        coupled_sym = coupled;
        coupled_sym.transmission = [T1, T2];
        [sol1s, sol2s] = solve(coupled_sym, k);
        d_sym = 0;
        for pair = {{sol1, sol1s}, {sol2, sol2s}}
            a = pair{1}{1}; b = pair{1}{2};
            for fld = 1:2
                d_sym = max(d_sym, max(abs(a.Unknowns.field(fld).dof ...
                                         - b.Unknowns.field(fld).dof), [], 'all'));
            end
            d_sym = max(d_sym, max(abs(a.Unknowns.hybrid(1).dof ...
                                     - b.Unknowns.hybrid(1).dof), [], 'all'));
        end
        fprintf("k = %d, h = %.4f: symbolic vs darcy_darcy_rows max dof diff = %.2e\n", ...
            k, h, d_sym);
        if d_sym > 1e-10
            fprintf("k = %d symbolic transmission FAILED\n", k); n_fail = n_fail + 1;
        end

        err_p1(i) = compute_l2_error(sol1.Unknowns.field(2), upper_mesh, p_);
        err_p2(i) = compute_l2_error(sol2.Unknowns.field(2), lower_mesh, p_);
        err_q1(i) = compute_l2_error(sol1.Unknowns.field(1), upper_mesh, q_);
        err_q2(i) = compute_l2_error(sol2.Unknowns.field(1), lower_mesh, q_);
        fprintf("k = %d, h = %.4f: err p = (%.3e, %.3e), err q = (%.3e, %.3e), residue = %.1e\n", ...
            k, h, err_p1(i), err_p2(i), err_q1(i), err_q2(i), gs.residue);
    end
    err_p = max([err_p1, err_p2], [], 2);
    err_q = max([err_q1, err_q2], [], 2);
    eoc_p = compute_eoc(err_p, h_list);
    eoc_q = compute_eoc(err_q, h_list);
    fprintf("k = %d: EOC(p) = %s, EOC(q) = %s (expected ~ %d)\n\n", ...
        k, num2str(eoc_p.', "%.2f "), num2str(eoc_q.', "%.2f "), k+1);
    if mean(eoc_p) < k + 0.5 || mean(eoc_q) < k + 0.5
        fprintf("k = %d EOC FAILED\n", k); n_fail = n_fail + 1;
    end
end
fprintf("%d check(s) failed.\n", n_fail);
if n_fail > 0
    error("test_interface_darcy_darcy:failed", "%d checks failed", n_fail);
end

function pb = configure(mesh, k, kappa, f_, p_, VPk, Pk, Mk, dx)
pb = FEProblem(Spaces= [FESpace([VPk, Pk], Condense=true), FESpace(Mk)], ...
    Order = k);
pb.Mesh = mesh;
[lhs_op, flux_op] = DarcyMatrices(k, mesh, kappa);
pb.Operator = AffineOperator(lhs_op, flux_op);
[~, w, ~] = get_test_functions(pb.Space);
F = assemble(f_ * w * dx, mesh);
pb.Operator = pb.Operator + F;
pb.BoundaryConditions = WeakEssentialBoundaryConditions(p_, "Dirichlet");
end
