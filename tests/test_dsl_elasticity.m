% test_dsl_elasticity.m
% The PEERS linear elasticity system written entirely as variational forms,
% including the density time terms (newmark, onestep start, marching;
% histories as FEField moments), must match LinearElasticityMatrices
% dof-for-dof.
% Polynomial source so every quadrature rule in play is exact.
clc; clear;

E = 1.0; nu = 0.25;
lambda = E*nu/((1 + nu)*(1 - 2*nu));
mu = E/2/(1 + nu);
c1 = 1/(2*mu);
c2 = lambda/(4*mu*(lambda + mu));
rho_s = 1.3; dt = 1e-2; beta_nm = 1/4;
f = {@(x, y) x.^2 - y; @(x, y) x.*y + 1};
u_zero = Datum({@(x, y) 0*x; @(x, y) 0*x});
J = [0, 1; -1, 0];

n_fail = 0;
for scheme = ["steady", "onestep", "marching", "newmark"]
    for k = [1, 2]
        d1 = k + 1; d2 = nchoosek(k + 2, 2);
        mesh = StructuredSimplicialGrid(8, 8).add_stabilization(1);
        nel = size(mesh.elements, 2);
        wave = @(c) reshape(sin(c + (1:2*d2*nel)), 2*d2, nel);
        u0 = wave(0.1); v0 = wave(0.2); a0 = wave(0.3); u00 = wave(0.4);

        % factory options and the matching (coefficient, history dofs)
        switch scheme
            case "steady"
                solid_opts = {};
                cm = 0; hist = [];
            case "onestep"
                solid_opts = {'density', rho_s, 'time_step', dt, ...
                    'state', struct('displacement', u0, 'velocity', v0)};
                cm = 2*rho_s/dt^2;
                hist = rho_s*(2*u0/dt^2 + 2*v0/dt);
            case "marching"
                solid_opts = {'density', rho_s, 'time_step', dt, ...
                    'state', struct('displacement', u0, 'displacement_previous', u00)};
                cm = rho_s/dt^2;
                hist = rho_s*(2*u0 - u00)/dt^2;
            case "newmark"
                solid_opts = {'density', rho_s, 'time_step', dt, ...
                    'time_scheme', "newmark", 'newmark_beta', beta_nm, ...
                    'state', struct('displacement', u0, 'velocity', v0, ...
                    'acceleration', a0)};
                cm = rho_s/(beta_nm*dt^2);
                hist = rho_s*((u0 + dt*v0)/(beta_nm*dt^2) + (1/2 - beta_nm)/beta_nm*a0);
        end

        Pk = DubinerBasis(k);
        mk = @() FEProblem( ...
            Spaces  = [FESpace([PEERSBasis(k), Pk.to_vector(), Pk], Condense=true), FESpace(LegendreBasis(k).to_vector(2))], ...
            Order = k);

        %% factory path
        pbA = mk(); pbA.Mesh = mesh;
        [lhs, flux, frhs] = LinearElasticityMatrices(k, mesh, E, nu, solid_opts{:});
        pbA.Operator = AffineOperator(lhs);
        pbA.Operator = pbA.Operator + (flux);
        if ~isempty(frhs), pbA.Operator = pbA.Operator + frhs; end
        w = TestFunction(Pk);
        F1 = double(Datum(f{1}) * w * Measure(mesh, 'dx'));
        F2 = double(Datum(f{2}) * w * Measure(mesh, 'dx'));
        pbA.Operator = pbA.Operator ...
            + IndexedMatrix(IndexedArray(cat(1, F1, F2), [-1, 2]));
        pbA.BoundaryConditions = WeakEssentialBoundaryConditions(u_zero, "Dirichlet");
        solA = solve(pbA, mesh, k);

        %% symbolic path
        pbB = mk(); pbB.Mesh = mesh;
        [sg, u, gm, u_hat] = get_trial_functions(pbB.Space);
        [T, v, e, v_hat] = get_test_functions(pbB.Space);
        dx = Measure(preset = "dx"); dS = Measure(preset = "dS");
        n_ = ExteriorNormal(); tau_ = Stabilization();
        f_ = Datum(f);
        local = ...
              dot(sg, c1*T)*dx - trace(sg)*(c2*trace(T))*dx ...
            + dot(u, divergence(T))*dx + dot(gm*J, T)*dx - dot(u_hat*n_.', T)*dS ...
            - dot(divergence(sg), v)*dx + dot(tau_*u, v)*dS ...
            - dot(tau_*u_hat, v)*dS + dot(f_, v)*dx ...
            + dot(sg, e*J)*dx;
        if cm > 0
            hist_f = FEField(FESpace(Pk.to_vector(), mesh), hist);
            local = local + dot(cm*u, v)*dx + dot(hist_f, v)*dx;
        end
        pbB.Operator = AffineOperator(local);
        pbB.Operator = pbB.Operator + (...
              dot(sg*n_, v_hat)*dS - dot(tau_*u, v_hat)*dS + dot(tau_*u_hat, v_hat)*dS);
        pbB.BoundaryConditions = WeakEssentialBoundaryConditions(u_zero, "Dirichlet");
        solB = solve(pbB, mesh, k);

        %% dof-for-dof
        d_sym = 0;
        for fld = 1:3
            d_sym = max(d_sym, max(abs(solA.Unknowns.field(fld).dof ...
                                     - solB.Unknowns.field(fld).dof), [], 'all'));
        end
        d_sym = max(d_sym, max(abs(solA.Unknowns.hybrid(1).dof ...
                                 - solB.Unknowns.hybrid(1).dof), [], 'all'));
        fprintf("%-8s k = %d: symbolic vs LinearElasticityMatrices max dof diff = %.2e\n", ...
            scheme, k, d_sym);
        if d_sym > 1e-9
            fprintf("%s k = %d symbolic spelling FAILED\n", scheme, k); n_fail = n_fail + 1;
        end
    end
end
fprintf("%d check(s) failed.\n", n_fail);
if n_fail > 0
    error("test_dsl_elasticity:failed", "%d checks failed", n_fail);
end
