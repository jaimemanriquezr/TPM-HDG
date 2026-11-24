function FSI_dsl_flat_convergence(k, s_list, max_iterations, h0)
% Convergence study of steady FSI on a flat interface with a gap closing
% as h^s (fitted meshes retreating from y0; TPM bridges the band), with
% both subdomains written in the form language as in FSI_dsl_example.m.
% For each s prints the L2 errors (fluid and solid combined) of sigma_h,
% u_h, gamma_h, u_hat and their EOC, plus the h-weighted skeleton error
% ||P_M u - u_hat||_h (P_M the L2 projection onto the trace space), which
% superconverges at k+2.
arguments
    k               = 2
    s_list          = [1, 2]
    max_iterations  = 4
    h0              = 1/4
end
y0 = 0.5;
viscosity = 1.0;  young_modulus = 1.0;  poisson_ratio = 0.33;  dt = 1.0;

interface = FlatCurveHorizontal(y0);
[~, ~, interface_sources, pd] = FluidStructureCoupling.load_problems( ...
    k, interface, 0.0, dt, "sin(pi*x)*sin(pi*y)", "0", ["cos(pi*x)"; "sin(pi*y)"], ...
    0.0, viscosity, 0.0, young_modulus, poisson_ratio, 0.0, ...
    pressure_mean = "pinned");
material = struct(viscosity = viscosity, young_modulus = young_modulus, ...
    poisson_ratio = poisson_ratio, time_step = dt);

for s = s_list
    [error_sigma, error_u, error_gamma, error_hat, error_proj, h_list] = deal(zeros(max_iterations, 1));
    for i = 1:max_iterations
        disp("s = " + s + ": iteration " + i + "/" + max_iterations);
        h = h0*(1/2)^(i - 1);
        [fluid_mesh, solid_mesh] = fitted_split(interface, h, h^s);

        coupled = fsi_dsl_problem(k, fluid_mesh, solid_mesh, interface, ...
            material, pd, interface_sources);
        [fluid_solution, solid_solution] = solve(coupled, k);

        e_F = field_errors(fluid_solution, fluid_mesh, pd.sigma_F, pd.u_F, pd.gamma_F);
        e_S = field_errors(solid_solution, solid_mesh, pd.sigma_S, pd.u_S, pd.gamma_S);
        e = sqrt(e_F.^2 + e_S.^2);

        n_elements = size(fluid_mesh.elements, 2) + size(solid_mesh.elements, 2);
        h_list(i)      = 1/sqrt(n_elements);
        error_sigma(i) = e(1);
        error_u(i)     = e(2);
        error_gamma(i) = e(3);
        error_hat(i)   = e(4);
        error_proj(i)  = e(5);
    end

    r_sigma = [NaN; compute_eoc(error_sigma, h_list)];
    r_u     = [NaN; compute_eoc(error_u,     h_list)];
    r_gamma = [NaN; compute_eoc(error_gamma, h_list)];
    r_hat   = [NaN; compute_eoc(error_hat,   h_list)];
    r_proj  = [NaN; compute_eoc(error_proj,  h_list)];

    fprintf("\nk = %d, gap = h^%g\n", k, s);
    hdr = "     h    |  e(σ_h)   | r(σ_h) |   e(u_h)  | r(u_h) |  e(γ_h)   | r(γ_h) |  e(û_h)   | r(û_h) | e(Pu-û_h) | r(Pu-û_h)";
    fprintf("%s\n%s\n", hdr, repmat('-', 1, strlength(hdr)));
    for i = 1:max_iterations
        fprintf(" %8.2e | %9.3e | %6.2f | %9.3e | %6.2f | %9.3e | %6.2f | %9.3e | %6.2f | %9.3e | %6.2f\n", ...
            h_list(i), error_sigma(i), r_sigma(i), error_u(i), r_u(i), ...
            error_gamma(i), r_gamma(i), error_hat(i), r_hat(i), error_proj(i), r_proj(i));
    end
end
end

%=========================================================================%
function E = field_errors(solution, mesh, sigma, u, gamma)
% L2 errors of sigma_h, u_h, gamma_h and u_hat, and ||P_M u - u_hat||_h;
% the multiplier spaces have no exact solution. On face fields
% compute_l2_error is h-weighted (sum_e h_e ||.||_e^2).
fields = solution.Unknowns.field;
u_hat = solution.Unknowns.hybrid(1);
P_u = Datum(u).project(mesh, u_hat.basis);
E = [compute_l2_error(fields(1), mesh, Datum(sigma));
     compute_l2_error(fields(2), mesh, Datum(u));
     compute_l2_error(fields(3), mesh, Datum(gamma));
     compute_l2_error(u_hat, mesh, Datum(u));
     compute_l2_error(FEField(FESpace(u_hat.basis, mesh), u_hat.dof - P_u.dof), mesh)];
end

%=========================================================================%
function coupled = fsi_dsl_problem(k, fluid_mesh, solid_mesh, interface, material, pd, interface_sources)
% The no-slip steady FSI problem of FSI_dsl_example.m on the given meshes.
J  = [0, 1; -1, 0];
Pk = DubinerBasis(k);
dx = Measure(preset = "dx");
dS = Measure(preset = "dS");
n_ = ExteriorNormal();
tau_ = Stabilization();

% --- fluid: incompressible elasticity (PEERS stress / velocity) ---------
mu = material.viscosity;
c1 = 1/(2*mu);  c2 = 1/(4*mu);
fluid = FEProblem(Order = k, ...
    Spaces = [FESpace([PEERSBasis(k), Pk.to_vector(), Pk, DubinerBasis(0)], Condense = true), ...
              FESpace([LegendreBasis(k).to_vector(2), DubinerBasis(0)])]);
fluid.Mesh = fluid_mesh;
[sg, u, gm, xi, u_hat, rho] = get_trial_functions(fluid.Space);
[T, v, e, z, v_hat, s_rho]  = get_test_functions(fluid.Space);
fluid.Operator = AffineOperator( ...
      dot(sg, c1*T)*dx - trace(sg)*(c2*trace(T))*dx ...
    + dot(u, divergence(T))*dx + dot(gm*J, T)*dx + xi*trace(T)*dx ...
    - dot(u_hat*n_.', T)*dS ...
    - dot(divergence(sg), v)*dx + dot((mu*tau_)*u, v)*dS ...
    - dot((mu*tau_)*u_hat, v)*dS + dot(Datum(pd.f_F), v)*dx ...
    + dot(sg, e*J)*dx ...
    + trace(sg)*z*dx - rho*z*dx);
fluid.Operator = fluid.Operator + ( ...
      dot(sg*n_, v_hat)*dS - dot((mu*tau_)*u, v_hat)*dS ...
    + dot((mu*tau_)*u_hat, v_hat)*dS ...
    + dot(u_hat, n_)*s_rho*dS);
fluid.BoundaryConditions = StrongEssentialBoundaryConditions(Datum(pd.u_F), ...
    "Dirichlet", projection_order = k);

% --- solid: linear elasticity (PEERS stress / displacement) -------------
E_S = material.young_modulus;  nu_S = material.poisson_ratio;
lambda = E_S*nu_S/((1 + nu_S)*(1 - 2*nu_S));
mu_S   = E_S/(2*(1 + nu_S));
c1_S = 1/(2*mu_S);  c2_S = lambda/(4*mu_S*(lambda + mu_S));
solid = FEProblem(Order = k, ...
    Spaces = [FESpace([PEERSBasis(k), Pk.to_vector(), Pk], Condense = true), ...
              FESpace(LegendreBasis(k).to_vector(2))]);
solid.Mesh = solid_mesh;
[sg, u, gm, u_hat] = get_trial_functions(solid.Space);
[T, v, e, v_hat]   = get_test_functions(solid.Space);
solid.Operator = AffineOperator( ...
      dot(sg, c1_S*T)*dx - trace(sg)*(c2_S*trace(T))*dx ...
    + dot(u, divergence(T))*dx + dot(gm*J, T)*dx - dot(u_hat*n_.', T)*dS ...
    - dot(divergence(sg), v)*dx + dot(tau_*u, v)*dS ...
    - dot(tau_*u_hat, v)*dS + dot(Datum(pd.f_S), v)*dx ...
    + dot(sg, e*J)*dx);
solid.Operator = solid.Operator + ( ...
      dot(sg*n_, v_hat)*dS - dot(tau_*u, v_hat)*dS + dot(tau_*u_hat, v_hat)*dS);
solid.BoundaryConditions = StrongEssentialBoundaryConditions(Datum(pd.u_S), ...
    "Dirichlet", projection_order = k);

% --- coupling -----------------------------------------------------------
coupled = InterfaceProblem.fsi(fluid, solid, interface, ...
    viscosity = mu, young_modulus = E_S, poisson_ratio = nu_S, ...
    domain_measure = pd.measure_F, ...
    scheme = "steady", time_step = material.time_step, pressure_mean = "pinned", ...
    source_neumann = interface_sources{1}, ...
    source_dirichlet = interface_sources{2});
end
