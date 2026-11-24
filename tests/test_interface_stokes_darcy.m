% test_interface_stokes_darcy.m
% Stokes-Darcy (velocity-gradient Stokes over Darcy, BJS friction) on
% InterfaceProblem, compared dof-for-dof against the low-level
% +StokesDarcyCoupling solver. Same meshes, same manufactured data (interface
% sources included), same quadrature orders (k+1/k+1, the low-level
% defaults), strong Dirichlet on both solvers. Also checks the symbolic
% transmission rows against stokes_darcy_rows.
clc; clear;

y0 = 0.5;
nu = 1.0;
kappa = eye(2);
wbjs = 10.0;
interface = FlatCurveHorizontal(y0);

data = StokesDarcyCoupling.get_manufactured_solution( ...
    "x^2*y^2", "x^2 + y^2", "x^3 + y^3", y0, nu, kappa, wbjs);

n_fail = 0;
for k = [1, 2]
    Pk = DubinerBasis(k); VPk = Pk.to_vector(); Mk = LegendreBasis(k);
    for h = [1/8, 1/16]
        d_gap = h^2;
        [upper_mesh, lower_mesh] = fitted_split(interface, [h, h], d_gap);

        %% low-level reference
        sp2 = StokesProblem(viscosity = data.nu, source_term = data.f_S, ...
            datum_dirichlet = data.u_S, domain_measure = data.domain_measure);
        dp2 = DarcyProblem(diffusion = data.kappa, source = data.f_D, ...
            datum_dirichlet = data.p_D);
        [ss2, ds2] = StokesDarcyCoupling.solve(sp2, dp2, ...
            upper_mesh, lower_mesh, interface, wbjs, data.T_1, data.T_2, ...
            Order = k, quiet = true);

        %% InterfaceProblem
        sp1 = FEProblem( ...
            Spaces  = [FESpace([Pk.to_tensor(), VPk, Pk, DubinerBasis(0)], Condense=true), FESpace([LegendreBasis(k).to_vector(2), DubinerBasis(0), ConstantBasis()])], ...
            Order = k);
        sp1.Mesh = upper_mesh;
        [lhs, flux, rhs] = StokesMatrices(k, upper_mesh, nu, data.f_S);
        sp1.Operator = AffineOperator(lhs + rhs);
        sp1.Operator = sp1.Operator + (flux);
        sp1.BoundaryConditions = StrongEssentialBoundaryConditions(Datum(data.u_S), "Dirichlet", projection_order = k);

        dp1 = FEProblem(Spaces= [FESpace([VPk, Pk], Condense=true), FESpace(Mk)], ...
            Order = k);
        dp1.Mesh = lower_mesh;
        [lhs_op, flux_op] = DarcyMatrices(k, lower_mesh, kappa);
dp1.Operator = AffineOperator(lhs_op, flux_op);
        [~, w, ~] = get_test_functions(dp1.Space);
        dp1.Operator = dp1.Operator ...
            + assemble(Datum(data.f_D) * w * Measure(preset = "dx"), lower_mesh);
        dp1.BoundaryConditions = StrongEssentialBoundaryConditions(Datum(data.p_D), "Dirichlet", projection_order = k);

        coupled = InterfaceProblem.stokes_darcy(sp1, dp1, interface, ...
            viscosity = nu, diffusion = kappa, friction_bjs = wbjs, ...
            domain_measure = data.domain_measure, ...
            source_neumann = data.T_1, source_dirichlet = data.T_2, ...
            quadrature_patch = k + 1, quadrature_path = k + 1);
        [ss1, ds1, gs] = solve(coupled, k);

        %% symbolic transmission: same coupling as forms,
        %% alpha = a 4th side-1 global space (ConstantBasis)
        sp1s = FEProblem( ...
            Spaces  = [FESpace([Pk.to_tensor(), VPk, Pk, DubinerBasis(0)], Condense=true), FESpace([LegendreBasis(k).to_vector(2), DubinerBasis(0), ...
                             ConstantBasis(), ConstantBasis()])], ...
            Order = k);
        sp1s.Mesh = upper_mesh;
        sp1s.Operator = sp1.Operator;
        sp1s.BoundaryConditions = sp1.BoundaryConditions;
        [L1, u1, p1, ~, uhat1, ~, ~, alpha1] = get_trial_functions(sp1s);
        [~, ~, ~, ~, mu1, ~, ~, xa] = get_test_functions(sp1s);
        [q2, p2, phat2] = get_trial_functions(dp1);
        [~, ~, mu2] = get_test_functions(dp1);
        n1_ = InterfaceNormal(upper_mesh);
        n2_ = InterfaceNormal(lower_mesh);
        tau_ = Stabilization();
        dI1 = Measure(upper_mesh, "dI");
        dI2 = Measure(lower_mesh, "dI");
        cA = dubiner_basis(0, 0, 0);
        nu_ = data.nu;
        T1_ = Datum(data.T_1); T2_ = Datum(data.T_2);
        T1s = ( nu_*(Extrapolate(L1)*n1_) - Extrapolate(p1)*n1_ ...
              - (nu_*tau_)*(u1 - uhat1) ...
              + wbjs*TangentialProjection(uhat1 + PathIntegral(L1)) ...
              - cA*(alpha1*n1_) ...
              + PathIntegral(kappa\q2)*n2_ - phat2*n2_ ...
              - T1_ ) * mu1 * dI1;
        T2s = ( dot(uhat1, n1_) + dot(PathIntegral(L1), n1_) ...
              + dot(Extrapolate(q2), n2_) + tau_*p2 - tau_*phat2 ...
              - T2_ ) * mu2 * dI2;
        Ta  = ( cA*PathIntegral(p1) ...
              + (data.domain_measure*cA*cA)*alpha1 ) * xa * dI1;
        coupled_sym = InterfaceProblem([sp1s, dp1], interface, ...
            struct("kind", "none"), ...
            quadrature_patch = k + 1, quadrature_path = k + 1);
        coupled_sym.transmission = [T1s, T2s, Ta];
        [ss1s, ds1s, gss] = solve(coupled_sym, k);
        nF1 = size(upper_mesh.faces, 2); nel1 = size(upper_mesh.elements, 2);
        alpha_h = gss.solution(2*(k+1)*nF1 + nel1 + 2);
        ss1s.Unknowns.field(3).dof(1, :) = ss1s.Unknowns.field(3).dof(1, :) + alpha_h;
        d_sym = 0;
        for fld = 1:3
            d_sym = max(d_sym, max(abs(ss1s.Unknowns.field(fld).dof ...
                                     - ss1.Unknowns.field(fld).dof), [], 'all'));
            if fld < 3
                d_sym = max(d_sym, max(abs(ds1s.Unknowns.field(fld).dof ...
                                         - ds1.Unknowns.field(fld).dof), [], 'all'));
            end
        end
        d_sym = max(d_sym, max(abs(ss1s.Unknowns.hybrid(1).dof ...
                                 - ss1.Unknowns.hybrid(1).dof), [], 'all'));
        d_sym = max(d_sym, max(abs(ds1s.Unknowns.hybrid(1).dof ...
                                 - ds1.Unknowns.hybrid(1).dof), [], 'all'));
        fprintf("k = %d, h = %.4f: symbolic vs stokes_darcy_rows max dof diff = %.2e\n", ...
            k, h, d_sym);
        if d_sym > 1e-10
            fprintf("k = %d symbolic transmission FAILED\n", k); n_fail = n_fail + 1;
        end

        %% dof-for-dof comparison
        dL = max(abs(ss1.Unknowns.field(1).dof - ss2.Unknowns.field(1).dof), [], "all");
        du = max(abs(ss1.Unknowns.field(2).dof - ss2.Unknowns.field(2).dof), [], "all");
        dp = max(abs(ss1.Unknowns.field(3).dof - ss2.Unknowns.field(3).dof), [], 'all');
        dq = max(abs(ds1.Unknowns.field(1).dof - ds2.Unknowns.field(1).dof), [], "all");
        dpd = max(abs(ds1.Unknowns.field(2).dof - ds2.Unknowns.field(2).dof), [], 'all');
        d1n = k + 1;
        uh1 = reshape(ss1.Unknowns.hybrid(1).dof, 2*d1n, []);
        duh = max(abs(uh1 - ss2.Unknowns.hybrid(1).dof), [], "all");
        drho = max(abs(reshape(ss1.Unknowns.hybrid(2).dof, 1, []) ...
                     - reshape(ss2.Unknowns.hybrid(2).dof, 1, [])), [], 'all');
        dph = max(abs(reshape(ds1.Unknowns.hybrid(1).dof, d1n, []) ...
                    - ds2.Unknowns.hybrid(1).dof), [], 'all');

        worst = max([dL, du, dp, dq, dpd, duh, drho, dph]);
        fprintf("k = %d, h = %.4f: worst dof diff = %.3e (L %.1e, u %.1e, p %.1e, q %.1e, pD %.1e, uhat %.1e, rho %.1e, phat %.1e)\n", ...
            k, h, worst, dL, du, dp, dq, dpd, duh, drho, dph);
        if worst > 1e-11
            fprintf("k = %d, h = %.4f ORACLE COMPARISON FAILED\n", k, h); n_fail = n_fail + 1;
        end
    end
end
fprintf("%d check(s) failed.\n", n_fail);
if n_fail > 0
    error("test_interface_stokes_darcy:failed", "%d checks failed", n_fail);
end
