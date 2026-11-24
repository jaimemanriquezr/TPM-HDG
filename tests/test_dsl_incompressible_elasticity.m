% test_dsl_incompressible_elasticity.m
% Stress-velocity (PEERS) incompressible elasticity on a heterogeneous
% global layout: IncompressibleElasticityMatrices on the declarative path vs
% the low-level @IncompressibleElasticityProblem, dof-for-dof, under both
% pressure closures ("constrained" carries the scalar multiplier space,
% "pinned" omits it). Then the same system as variational forms vs the factory.
clc; clear;

mu = 0.05;
tau = 1;
f = {@(x, y) sin(pi*x).*cos(pi*y); @(x, y) exp(x).*y - exp(1)/2 + 1/2};
u_zero = Datum({@(x, y) 0*x; @(x, y) 0*x});

n_fail = 0;
for closure = ["constrained", "pinned"]
    for k = [1, 2]
        d1 = k + 1;
        mesh = StructuredSimplicialGrid(8, 8).add_stabilization(tau);
        if closure == "pinned"
            % pinned needs the boundary to fix the pressure level:
            % traction (Neumann) faces on x = 1, zero datum on both paths
            mesh = apply_tag(mesh, "Neumann", @(X) isapprox(X(1, :), 1));
        end

        %% low-level reference
        pb2 = IncompressibleElasticityProblem(k, viscosity = mu, ...
            source_term = f, pressure_mean = closure);
        sol2 = solve(pb2, mesh, Order= k);

        %% declarative
        Pk = DubinerBasis(k);
        gspaces = [LegendreBasis(k).to_vector(2), DubinerBasis(0)];
        if closure == "constrained"
            gspaces = [gspaces, ConstantBasis()]; %#ok<AGROW>
        end
        problem = FEProblem( ...
            Spaces  = [FESpace([PEERSBasis(k), Pk.to_vector(), Pk, DubinerBasis(0)], Condense=true), FESpace(gspaces)]);
        problem.Mesh = mesh;
        [lhs, flux, rhs] = IncompressibleElasticityMatrices(k, mesh, mu, f, ...
            pressure_mean = closure);
        problem.Operator = AffineOperator(lhs + rhs);
        problem.Operator = problem.Operator + (flux);
        problem.BoundaryConditions = WeakEssentialBoundaryConditions(u_zero, "Dirichlet");

        sol1 = solve(problem, mesh, k);

        %% dof-for-dof comparison
        s2 = sol2.Unknowns.field(1).dof;
        u2 = sol2.Unknowns.field(2).dof;
        g2 = sol2.Unknowns.field(3).dof;
        uhat2 = sol2.Unknowns.hybrid(1).dof; % 2*d1 x nf entity layout
        rho2  = sol2.Unknowns.hybrid(2).dof;

        s1 = sol1.Unknowns.field(1).dof;
        u1 = sol1.Unknowns.field(2).dof;
        g1 = sol1.Unknowns.field(3).dof;
        uhat1 = reshape(sol1.Unknowns.hybrid(1).dof, 2*d1, []);
        rho1  = reshape(sol1.Unknowns.hybrid(2).dof, 1, []);

        ds  = max(abs(s1 - s2), [], 'all');
        du  = max(abs(u1 - u2), [], 'all');
        dg  = max(abs(g1 - g2), [], 'all');
        duh = max(abs(uhat1 - uhat2), [], "all");
        drho = max(abs(rho1 - reshape(rho2, 1, [])), [], 'all');

        scale = max(abs(u2), [], 'all');
        fprintf("%s k = %d: |ds| = %.2e |du| = %.2e |dg| = %.2e |duhat| = %.2e |drho| = %.2e (|u| ~ %.2e)\n", ...
            closure, k, ds, du, dg, duh, drho, scale);
        if max([ds, du, dg, duh, drho]) > 1e-9 * max(1, scale)
            fprintf("%s k = %d dof comparison FAILED\n", closure, k); n_fail = n_fail + 1;
        end
    end
end
%% The same system written entirely as variational forms (both closures,
%% plus the theta-scheme momentum mass/history terms, the histories through
%% the factory's state option) vs the factory.
J = [0, 1; -1, 0];
dt = 1e-2;
for cfg = ["constrained-steady", "pinned-steady", "euler", "cranknicolson"]
    closure = "constrained";
    if cfg == "pinned-steady", closure = "pinned"; end
    for k = [1, 2]
        d1 = k + 1; d2 = nchoosek(k + 2, 2);
        c1 = 1/(2*mu); c2 = 1/(4*mu);
        mesh = StructuredSimplicialGrid(8, 8).add_stabilization(tau);
        if closure == "pinned"
            mesh = apply_tag(mesh, "Neumann", @(X) isapprox(X(1, :), 1));
        end
        nel = size(mesh.elements, 2);
        wave = @(c) reshape(sin(c + (1:2*d2*nel)), 2*d2, nel);
        u0 = wave(0.1); a0 = wave(0.2);
        switch cfg
            case {"constrained-steady", "pinned-steady"}
                fl_opts = {}; cm = 0; hist = [];
            case "euler"
                fl_opts = {'density', 1.1, 'time_step', dt, 'theta', 1.0, ...
                    'state', struct('velocity', u0, 'acceleration', [])};
                cm = 1.1/dt; hist = 1.1*u0/dt;
            case "cranknicolson"
                fl_opts = {'density', 1.1, 'time_step', dt, 'theta', 0.5, ...
                    'state', struct('velocity', u0, 'acceleration', a0)};
                cm = 1.1/(0.5*dt); hist = 1.1*(u0/(0.5*dt) + a0);
        end

        Pk = DubinerBasis(k);
        gspaces = [LegendreBasis(k).to_vector(2), DubinerBasis(0)];
        if closure == "constrained"
            gspaces = [gspaces, ConstantBasis()]; %#ok<AGROW>
        end
        mk = @() FEProblem( ...
            Spaces  = [FESpace([PEERSBasis(k), Pk.to_vector(), Pk, DubinerBasis(0)], Condense=true), FESpace(gspaces)], Order= k);

        pbA = mk(); pbA.Mesh = mesh;
        [lhs, flux, frhs] = IncompressibleElasticityMatrices(k, mesh, mu, f, ...
            "pressure_mean", closure, fl_opts{:});
        pbA.Operator = AffineOperator(lhs + frhs);
        pbA.Operator = pbA.Operator + (flux);
        pbA.BoundaryConditions = WeakEssentialBoundaryConditions(u_zero, "Dirichlet");
        solA = solve(pbA, mesh, k);

        pbB = mk(); pbB.Mesh = mesh;
        if closure == "constrained"
            [sg, u, gm, xi, u_hat, rho, mu_c] = get_trial_functions(pbB.Space);
            [T, v, e, z, v_hat, s_rho, s_mu] = get_test_functions(pbB.Space);
        else
            [sg, u, gm, xi, u_hat, rho] = get_trial_functions(pbB.Space);
            [T, v, e, z, v_hat, s_rho] = get_test_functions(pbB.Space);
        end
        dx = Measure(preset = "dx"); dS = Measure(preset = "dS");
        dxk = Measure(preset = "dx"); dxk.order = k;   % factory rhs rule
        n_ = ExteriorNormal(); tau_ = Stabilization();
        f_ = Datum(f);
        local = ...
              dot(sg, c1*T)*dx - trace(sg)*(c2*trace(T))*dx ...
            + dot(u, divergence(T))*dx + dot(gm*J, T)*dx + xi*trace(T)*dx ...
            - dot(u_hat*n_.', T)*dS ...
            - dot(divergence(sg), v)*dx + dot((mu*tau_)*u, v)*dS ...
            - dot((mu*tau_)*u_hat, v)*dS + dot(f_, v)*dxk ...
            + dot(sg, e*J)*dx ...
            + trace(sg)*z*dx - rho*z*dx;
        if cm > 0
            hist_f = FEField(FESpace(Pk.to_vector(), mesh), hist);
            local = local + dot(cm*u, v)*dx + dot(hist_f, v)*dx;
        end
        pbB.Operator = AffineOperator(local);
        glob = ...
              dot(sg*n_, v_hat)*dS - dot((mu*tau_)*u, v_hat)*dS ...
            + dot((mu*tau_)*u_hat, v_hat)*dS ...
            + dot(u_hat, n_)*s_rho*dS;
        if closure == "constrained"
            glob = glob + mu_c*s_rho*dx + rho*s_mu*dx;
        end
        pbB.Operator = pbB.Operator + (glob);
        pbB.BoundaryConditions = WeakEssentialBoundaryConditions(u_zero, "Dirichlet");
        solB = solve(pbB, mesh, k);

        d_sym = 0;
        for fld = 1:4
            d_sym = max(d_sym, max(abs(solA.Unknowns.field(fld).dof ...
                                     - solB.Unknowns.field(fld).dof), [], 'all'));
        end
        for hyb = 1:2
            d_sym = max(d_sym, max(abs(solA.Unknowns.hybrid(hyb).dof ...
                                     - solB.Unknowns.hybrid(hyb).dof), [], 'all'));
        end
        fprintf("%-19s k = %d: symbolic vs factory max dof diff = %.2e\n", cfg, k, d_sym);
        if d_sym > 1e-9
            fprintf("%s k = %d symbolic spelling FAILED\n", cfg, k); n_fail = n_fail + 1;
        end
    end
end

fprintf("%d check(s) failed.\n", n_fail);
if n_fail > 0
    error("test_dsl_incompressible_elasticity:failed", "%d checks failed", n_fail);
end
