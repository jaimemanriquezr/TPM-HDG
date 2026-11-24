% test_dsl_stokes.m
% HDG Stokes on the heterogeneous global layout (u_hat skeleton + rho
% element means + mu scalar zero-mean multiplier), assembled from
% StokesMatrices on the declarative path, must match the low-level
% @StokesProblem dof-for-dof. Also checks strong vs weak zero Dirichlet and
% the fully symbolic spelling of the same system.
clc; clear;

nu = 0.7;
tau = 1;
f = {@(x, y) sin(pi*x).*cos(pi*y); @(x, y) exp(x).*y - exp(1)/2 + 1/2};
u_zero = Datum({@(x, y) 0*x; @(x, y) 0*x});

n_fail = 0;
for k = [1, 2]
    d1 = k + 1;
    mesh = StructuredSimplicialGrid(8, 8).add_stabilization(tau);

    %% low-level reference
    pb2 = StokesProblem(k, viscosity = nu, source_term = f);
    sol2 = solve(pb2, mesh, Order= k);

    %% declarative: heterogeneous global layout
    Pk  = DubinerBasis(k);
    problem = FEProblem( ...
        Spaces  = [FESpace([Pk.to_tensor(), Pk.to_vector(), Pk, DubinerBasis(0)], Condense=true), FESpace([LegendreBasis(k).to_vector(2), DubinerBasis(0), ConstantBasis()])]);
    problem.Mesh = mesh;
    [lhs, flux, rhs] = StokesMatrices(k, mesh, nu, f);
    problem.Operator = AffineOperator(lhs + rhs);
    problem.Operator = problem.Operator + (flux);
    problem.BoundaryConditions = WeakEssentialBoundaryConditions(u_zero, "Dirichlet");

    sol1 = solve(problem, mesh, k);

    %% Strong essential BCs on a vector trace: with the same zero datum,
    %% strong and weak treatments must agree
    problem.BoundaryConditions = StrongEssentialBoundaryConditions(u_zero, "Dirichlet");
    sol_strong = solve(problem, mesh, k);
    d_strong = max(abs(sol_strong.Unknowns.field(2).dof ...
                     - sol1.Unknowns.field(2).dof), [], 'all');
    fprintf("k = %d: strong-vs-weak zero-datum |du| = %.2e\n", k, d_strong);
    if d_strong > 1e-10
        fprintf("k = %d strong BC comparison FAILED\n", k); n_fail = n_fail + 1;
    end

    %% fully symbolic spelling: the same system written as variational
    %% forms, block-for-block from
    %% StokesMatrices (no integration by parts), closures as integrals,
    %% tau via Stabilization(). Must match the factory path dof-for-dof.
    [L, u, p, xi, u_hat, rho, mu_c] = get_trial_functions(problem.Space);
    [T, v, q, e, v_hat, s_rho, s_mu] = get_test_functions(problem.Space);
    dx = Measure(preset = "dx");
    dS = Measure(preset = "dS");
    % the factory RHS (test_elements) integrates at order k; pin the
    % source moments to the same rule for the dof-for-dof identity
    dxk = Measure(preset = "dx");
    dxk.order = k;
    n_ = ExteriorNormal();
    tau_ = Stabilization();
    f_ = Datum(f);

    sym_problem = problem;
    sym_problem.BoundaryConditions = WeakEssentialBoundaryConditions(u_zero, "Dirichlet");
    sym_problem.Operator = AffineOperator(...
          dot(L, T)*dx + dot(u, divergence(T))*dx - dot(u_hat*n_.', T)*dS ...
        - dot(divergence(L), nu*v)*dx + dot((nu*tau_)*u, v)*dS ...
            + dot(gradient(p), v)*dx - dot((nu*tau_)*u_hat, v)*dS + dot(f_, v)*dxk ...
        + dot(u, gradient(q))*dx + xi*q*dx - dot(u_hat, n_)*q*dS ...
        + p*e*dx - rho*e*dx);
    sym_problem.Operator = sym_problem.Operator + (...
          dot(L*n_, nu*v_hat)*dS - dot((nu*tau_)*u, v_hat)*dS ...
        - p*dot(v_hat, n_)*dS + dot((nu*tau_)*u_hat, v_hat)*dS ...
        + dot(u_hat, n_)*s_rho*dS + mu_c*s_rho*dx + rho*s_mu*dx);
    sol_sym = solve(sym_problem, mesh, k);
    d_sym = 0;
    for fld = 1:3
        d_sym = max(d_sym, max(abs(sol_sym.Unknowns.field(fld).dof ...
                                 - sol1.Unknowns.field(fld).dof), [], 'all'));
    end
    for hyb = 1:2
        d_sym = max(d_sym, max(abs(sol_sym.Unknowns.hybrid(hyb).dof ...
                                 - sol1.Unknowns.hybrid(hyb).dof), [], 'all'));
    end
    fprintf("k = %d: symbolic vs StokesMatrices max dof diff = %.2e\n", k, d_sym);
    if d_sym > 1e-10
        fprintf("k = %d symbolic spelling FAILED\n", k); n_fail = n_fail + 1;
    end

    %% dof-for-dof comparison
    L2 = sol2.Unknowns.field(1).dof;     % [L11; L21; L12; L22] entity layout
    u2 = sol2.Unknowns.field(2).dof;
    p2 = sol2.Unknowns.field(3).dof;
    uhat2 = sol2.Unknowns.hybrid(1).dof; % 2*d1 x nf entity layout
    rho2  = sol2.Unknowns.hybrid(2).dof;

    L1 = sol1.Unknowns.field(1).dof;
    u1 = sol1.Unknowns.field(2).dof;
    p1 = sol1.Unknowns.field(3).dof;
    uhat1 = reshape(sol1.Unknowns.hybrid(1).dof, 2*d1, []);
    rho1  = reshape(sol1.Unknowns.hybrid(2).dof, 1, []);

    dL = max(abs(L1 - L2), [], 'all');
    du = max(abs(u1 - u2), [], 'all');
    dp = max(abs(p1 - p2), [], 'all');
    duh = max(abs(uhat1 - uhat2), [], "all");
    drho = max(abs(rho1 - reshape(rho2, 1, [])), [], 'all');

    scale = max(abs(u2), [], 'all');
    fprintf("k = %d: |dL| = %.2e  |du| = %.2e  |dp| = %.2e  |duhat| = %.2e  |drho| = %.2e  (|u| ~ %.2e)\n", ...
        k, dL, du, dp, duh, drho, scale);
    if max([dL, du, dp, duh, drho]) > 1e-10 * max(1, scale)
        fprintf("k = %d dof comparison FAILED\n", k); n_fail = n_fail + 1;
    end
end
fprintf("%d check(s) failed.\n", n_fail);
if n_fail > 0
    error("test_dsl_stokes:failed", "%d checks failed", n_fail);
end
